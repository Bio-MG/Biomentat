# =============================================================================
# tools/pkg_skeleton/R/bulk_meta_analysis.R — Horizontal multi-study DE
#                                             meta-analysis engine (BK-07)
# =============================================================================
# Pure domain logic (zero Shiny symbols, no top-level library()/require()).
#
# Combines differential expression summary tables across >= 2 independent
# studies without pooling raw count matrices (horizontal meta-analysis).
# Harmonizes genes observed in at least `min_studies` studies, clamps
# p-values to [1e-300, 1], and combines evidence using one of four standard
# meta-analysis combiners:
#   - "stouffer": weighted/unweighted Stouffer Z-score (direction-aware or
#                 unsigned)
#   - "fisher"  : Fisher's chi-squared combination of log p-values
#   - "maxP"    : Wilkinson conjunctive test (requires significance in all
#                 studies; p_max^K)
#   - "minP"    : Tippett disjunctive test (1 - (1 - p_min)^K)
#
# Deterministic: no random sampling.
# =============================================================================

#' Horizontal multi-study differential expression meta-analysis
#'
#' Harmonizes per-study DE result tables and combines p-values across studies
#' using Stouffer, Fisher, maxP (Wilkinson), or minP (Tippett).
#'
#' @param study_tables List of >= 2 data.frames (named or unnamed; unnamed
#'   entries are auto-named `study_1`, `study_2`, ...). Each table must provide
#'   gene identifiers (via `gene`, `gene_id`, `symbol` column or explicit
#'   `rownames`), a log2 fold-change column (`log2FC`, `log2FoldChange`, or
#'   `logFC`), and a p-value column (`pvalue`, `PValue`, `P.Value`, or `pval`).
#' @param method Character scalar: one of `"stouffer"`, `"fisher"`, `"maxP"`,
#'   `"minP"`.
#' @param min_studies Integer scalar >= 1: minimum number of studies in which
#'   a gene must be present to be included in the meta-analysis (default `2L`).
#' @param fdr_method Character scalar passed to [stats::p.adjust()] (default
#'   `"BH"`).
#' @param weights Optional positive numeric vector of study weights (length
#'   equal to `length(study_tables)`), used by `"stouffer"`. If named with
#'   study names, weights are aligned by name. Default `NULL` (equal weights).
#' @param direction_aware Logical scalar (default `TRUE`). For `"stouffer"`,
#'   when `TRUE`, study p-values are converted to signed two-sided Z-scores
#'   using `sign(log2FC)`; when `FALSE`, unsigned one-sided Z-scores are used.
#' @return A list with elements:
#'   - `table`: `data.frame` sorted by `FDR`, `meta_pval`, with columns
#'     `gene`, `meta_pval`, `FDR`, `mean_log2FC`, `n_studies`,
#'     `direction_concordance`.
#'   - `method`: combiner method used.
#'   - `direction_aware`: logical flag used.
#'   - `n_studies`: total number of input studies.
#'   - `n_genes`: number of harmonized genes in `table`.
#'   - `study_names`: character vector of study names.
#'   - `provenance`: `list(engine = "biomentat_horizontal_meta_de", method = method)`.
run_horizontal_meta_de <- function(study_tables,
                                   method = c("stouffer", "fisher", "maxP", "minP"),
                                   min_studies = 2L,
                                   fdr_method = "BH",
                                   weights = NULL,
                                   direction_aware = TRUE) {
  .stop <- function(msg) {
    stop(errorCondition(msg, class = "bulk_meta_analysis_error"))
  }

  # --- 1. Validate inputs ---------------------------------------------------
  if (is.null(study_tables) || !is.list(study_tables) ||
      is.data.frame(study_tables) || length(study_tables) < 2L) {
    .stop("Horizontal meta-analysis requires a list of at least 2 study data.frames.")
  }

  n_studies_total <- length(study_tables)
  study_names <- names(study_tables)
  if (is.null(study_names)) {
    study_names <- paste0("study_", seq_len(n_studies_total))
  } else {
    missing_nm <- is.na(study_names) | !nzchar(trimws(study_names))
    if (any(missing_nm)) {
      study_names[missing_nm] <- paste0("study_", which(missing_nm))
    }
    if (anyDuplicated(study_names)) {
      study_names <- make.unique(study_names, sep = "_")
    }
  }
  names(study_tables) <- study_names

  if (!is.character(method) || length(method) == 0L || anyNA(method)) {
    .stop("`method` must be one of: 'stouffer', 'fisher', 'maxP', 'minP'.")
  }
  valid_methods <- c("stouffer", "fisher", "maxP", "minP")
  if (length(method) == 1L && !(method %in% valid_methods)) {
    .stop(sprintf(
      "Unsupported meta-analysis method '%s'. Expected one of: %s.",
      method, paste(valid_methods, collapse = ", ")
    ))
  }
  method <- match.arg(method, choices = valid_methods)

  if (!is.numeric(min_studies) || length(min_studies) != 1L ||
      is.na(min_studies) || !is.finite(min_studies) || min_studies < 1) {
    .stop("`min_studies` must be a single positive integer >= 1.")
  }
  min_studies <- as.integer(min_studies)

  if (!is.character(fdr_method) || length(fdr_method) != 1L ||
      is.na(fdr_method) || !(fdr_method %in% stats::p.adjust.methods)) {
    .stop(sprintf(
      "Invalid `fdr_method` '%s'. Expected one of stats::p.adjust.methods.",
      as.character(fdr_method)[1L]
    ))
  }

  if (!is.logical(direction_aware) || length(direction_aware) != 1L ||
      is.na(direction_aware)) {
    .stop("`direction_aware` must be a single non-NA logical value.")
  }

  if (is.null(weights)) {
    w_vec <- rep(1, n_studies_total)
    names(w_vec) <- study_names
  } else {
    if (!is.numeric(weights) || length(weights) != n_studies_total ||
        anyNA(weights) || any(!is.finite(weights)) || any(weights <= 0)) {
      .stop(sprintf(
        "`weights` must be a positive finite numeric vector of length %d.",
        n_studies_total
      ))
    }
    if (!is.null(names(weights)) && all(study_names %in% names(weights))) {
      w_vec <- as.numeric(weights[study_names])
      names(w_vec) <- study_names
    } else {
      w_vec <- as.numeric(weights)
      names(w_vec) <- study_names
    }
  }

  # --- 2. Extract and harmonize per-study tables ----------------------------
  gene_col_candidates <- c("gene", "gene_id", "symbol")
  lfc_col_candidates  <- c("log2FC", "log2FoldChange", "logFC")
  pval_col_candidates <- c("pvalue", "PValue", "P.Value", "pval")

  extracted <- vector("list", n_studies_total)
  names(extracted) <- study_names

  for (i in seq_len(n_studies_total)) {
    sname <- study_names[i]
    df <- study_tables[[i]]
    if (is.null(df) || !is.data.frame(df)) {
      .stop(sprintf("Study '%s' must be a data.frame.", sname))
    }

    cols <- colnames(df)
    lfc_col <- lfc_col_candidates[lfc_col_candidates %in% cols][1L]
    if (is.na(lfc_col)) {
      .stop(sprintf(
        "Study '%s' is missing a log2 fold-change column (expected one of: %s).",
        sname, paste(lfc_col_candidates, collapse = ", ")
      ))
    }

    pval_col <- pval_col_candidates[pval_col_candidates %in% cols][1L]
    if (is.na(pval_col)) {
      .stop(sprintf(
        "Study '%s' is missing a p-value column (expected one of: %s).",
        sname, paste(pval_col_candidates, collapse = ", ")
      ))
    }

    gene_col <- gene_col_candidates[gene_col_candidates %in% cols][1L]
    if (!is.na(gene_col)) {
      genes <- as.character(df[[gene_col]])
    } else if (.row_names_info(df) > 0L) {
      genes <- as.character(rownames(df))
    } else {
      .stop(sprintf(
        "Study '%s' is missing gene identifiers (expected column %s or explicit rownames).",
        sname, paste(gene_col_candidates, collapse = "/")
      ))
    }

    lfc <- suppressWarnings(as.numeric(df[[lfc_col]]))
    pval <- suppressWarnings(as.numeric(df[[pval_col]]))

    keep <- !is.na(genes) & nzchar(trimws(genes)) &
      is.finite(lfc) & is.finite(pval)
    genes <- trimws(genes[keep])
    lfc <- lfc[keep]
    pval <- pmin(pmax(pval[keep], 1e-300), 1)

    if (anyDuplicated(genes)) {
      # Keep most significant entry per gene deterministically.
      ord_dup <- order(pval, seq_along(pval))
      genes <- genes[ord_dup]
      lfc <- lfc[ord_dup]
      pval <- pval[ord_dup]
      first_idx <- !duplicated(genes)
      genes <- genes[first_idx]
      lfc <- lfc[first_idx]
      pval <- pval[first_idx]
    }

    extracted[[i]] <- list(gene = genes, log2FC = lfc, pvalue = pval)
  }

  all_genes <- unique(unlist(lapply(extracted, `[[`, "gene"), use.names = FALSE))

  empty_result <- function() {
    empty_df <- data.frame(
      gene = character(0),
      meta_pval = numeric(0),
      FDR = numeric(0),
      mean_log2FC = numeric(0),
      n_studies = integer(0),
      direction_concordance = numeric(0),
      stringsAsFactors = FALSE
    )
    list(
      table = empty_df,
      method = method,
      direction_aware = direction_aware,
      n_studies = n_studies_total,
      n_genes = 0L,
      study_names = study_names,
      provenance = list(
        engine = "biomentat_horizontal_meta_de",
        method = method
      )
    )
  }

  if (length(all_genes) == 0L || min_studies > n_studies_total) {
    return(empty_result())
  }

  n_all <- length(all_genes)
  lfc_mat <- matrix(NA_real_, nrow = n_all, ncol = n_studies_total)
  pval_mat <- matrix(NA_real_, nrow = n_all, ncol = n_studies_total)

  for (j in seq_len(n_studies_total)) {
    ex <- extracted[[j]]
    if (length(ex$gene) == 0L) next
    idx <- match(ex$gene, all_genes)
    lfc_mat[idx, j] <- ex$log2FC
    pval_mat[idx, j] <- ex$pvalue
  }

  present <- !is.na(pval_mat)
  k_vec <- as.integer(rowSums(present))
  keep_genes <- k_vec >= min_studies

  if (!any(keep_genes)) {
    return(empty_result())
  }

  genes_kept <- all_genes[keep_genes]
  lfc_mat <- lfc_mat[keep_genes, , drop = FALSE]
  pval_mat <- pval_mat[keep_genes, , drop = FALSE]
  present <- present[keep_genes, , drop = FALSE]
  k_vec <- k_vec[keep_genes]
  n_g <- length(genes_kept)

  # --- 3. Combine p-values across studies -----------------------------------
  if (method == "stouffer") {
    w_mat <- matrix(w_vec, nrow = n_g, ncol = n_studies_total, byrow = TRUE)
    w_mat[!present] <- 0
    denom <- sqrt(rowSums(w_mat^2))

    if (direction_aware) {
      z_mat <- sign(lfc_mat) * stats::qnorm(pval_mat / 2, lower.tail = FALSE)
      z_mat[!present] <- 0
      z_comb <- rowSums(w_mat * z_mat) / denom
      meta_pval <- pmin(pmax(2 * stats::pnorm(-abs(z_comb)), 0), 1)
    } else {
      p_eff <- pmin(pval_mat, 1 - 1e-15)
      z_mat <- stats::qnorm(p_eff, lower.tail = FALSE)
      z_mat[!present] <- 0
      z_comb <- rowSums(w_mat * z_mat) / denom
      meta_pval <- pmin(pmax(stats::pnorm(z_comb, lower.tail = FALSE), 0), 1)
    }
  } else if (method == "fisher") {
    log_p <- log(pval_mat)
    log_p[!present] <- 0
    chi2 <- -2 * rowSums(log_p)
    meta_pval <- stats::pchisq(chi2, df = 2 * k_vec, lower.tail = FALSE)
  } else if (method == "maxP") {
    p_max <- rep(-Inf, n_g)
    for (j in seq_len(n_studies_total)) {
      idx <- present[, j]
      if (any(idx)) {
        p_max[idx] <- pmax(p_max[idx], pval_mat[idx, j])
      }
    }
    meta_pval <- p_max^k_vec
  } else {
    # method == "minP" (Tippett: 1 - (1 - p_min)^K)
    p_min <- rep(Inf, n_g)
    for (j in seq_len(n_studies_total)) {
      idx <- present[, j]
      if (any(idx)) {
        p_min[idx] <- pmin(p_min[idx], pval_mat[idx, j])
      }
    }
    meta_pval <- -expm1(k_vec * log1p(-p_min))
  }

  # --- 4. Compute summary metrics & sort ------------------------------------
  mean_lfc <- rowMeans(lfc_mat, na.rm = TRUE)
  all_zero <- rowSums(lfc_mat != 0, na.rm = TRUE) == 0L
  same_sign <- rowSums(sign(lfc_mat) == sign(mean_lfc), na.rm = TRUE)
  dir_concordance <- ifelse(all_zero, 1, same_sign / k_vec)

  fdr <- stats::p.adjust(meta_pval, method = fdr_method)

  res_df <- data.frame(
    gene = genes_kept,
    meta_pval = as.numeric(meta_pval),
    FDR = as.numeric(fdr),
    mean_log2FC = as.numeric(mean_lfc),
    n_studies = as.integer(k_vec),
    direction_concordance = as.numeric(dir_concordance),
    stringsAsFactors = FALSE
  )

  ord <- order(res_df$FDR, res_df$meta_pval, res_df$gene)
  res_df <- res_df[ord, , drop = FALSE]
  rownames(res_df) <- NULL

  list(
    table = res_df,
    method = method,
    direction_aware = direction_aware,
    n_studies = n_studies_total,
    n_genes = nrow(res_df),
    study_names = study_names,
    provenance = list(
      engine = "biomentat_horizontal_meta_de",
      method = method
    )
  )
}
