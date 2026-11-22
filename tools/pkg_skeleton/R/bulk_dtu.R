# =============================================================================
# tools/pkg_skeleton/R/bulk_dtu.R — Differential Transcript Usage (DTU) engine
# =============================================================================
# Pure domain logic (zero Shiny symbols). Implements transcript-level isoform
# switching / Differential Transcript Usage (DTU) analysis via
# edgeR::diffSpliceDGE() and edgeR::topSpliceDGE().
#
# Dependencies: edgeR (checked via requireNamespace), stats, utils.
# =============================================================================

#' Run Differential Transcript Usage (DTU) analysis via edgeR::diffSpliceDGE
#'
#' @param counts Numeric/integer matrix (or coercible data.frame) of transcript
#'   counts with rownames = transcript IDs and columns = samples.
#' @param tx2gene data.frame mapping transcripts to genes, containing
#'   `transcript_id` and `gene_id` (also accepts `tx_id` / `gene_name`).
#' @param design_or_groups Either a numeric design matrix (`nrow == ncol(counts)`)
#'   or a factor/character vector of sample group labels (`length == ncol(counts)`,
#'   at least 2 levels).
#' @param coef Integer or character coefficient index/name in the design matrix
#'   to test (default `2L`).
#' @param min_count Minimum count per transcript passed to `edgeR::filterByExpr`.
#' @param min_total_count Minimum total count per transcript passed to
#'   `edgeR::filterByExpr`.
#' @param robust Logical; whether to use robust empirical Bayes squeezing in
#'   `edgeR::glmQLFit`.
#' @param test_method Character; gene-level test method: `"Simes"`, `"F"` (maps
#'   to `"gene"` in `edgeR::topSpliceDGE`), or `"t"` (maps to `"exon"`).
#' @return A list with `gene_table`, `transcript_table`, `n_genes_tested`,
#'   `n_transcripts_tested`, `test_method`, and `provenance`.
run_bulk_dtu <- function(counts,
                         tx2gene,
                         design_or_groups,
                         coef = 2L,
                         min_count = 10L,
                         min_total_count = 15L,
                         robust = TRUE,
                         test_method = c("Simes", "F", "t")) {
  .stop <- function(msg) {
    errorCondition(msg, class = "bulk_dtu_error")
  }

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop(.stop("Package 'edgeR' is required for DTU analysis but is not installed."))
  }

  test_method <- match.arg(test_method)

  # --- Validate counts ---
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  if (is.null(counts) || !is.matrix(counts) || !is.numeric(counts)) {
    stop(.stop(sprintf(
      "DTU counts must be a numeric/integer matrix or data.frame (got: %s).",
      if (is.null(counts)) "NULL" else paste(class(counts), collapse = "/")
    )))
  }
  if (nrow(counts) < 2L || ncol(counts) < 2L) {
    stop(.stop("DTU counts matrix must have at least 2 transcripts and 2 samples."))
  }
  tx_names <- rownames(counts)
  if (is.null(tx_names) || anyNA(tx_names) || any(!nzchar(tx_names))) {
    stop(.stop("DTU counts matrix must have non-empty rownames corresponding to transcript IDs."))
  }
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop(.stop("DTU counts matrix must contain finite non-negative values without NAs."))
  }

  # --- Validate tx2gene ---
  if (is.null(tx2gene) || !is.data.frame(tx2gene)) {
    stop(.stop("tx2gene must be a data.frame mapping transcripts to genes."))
  }
  tx_col <- if ("transcript_id" %in% colnames(tx2gene)) {
    "transcript_id"
  } else if ("tx_id" %in% colnames(tx2gene)) {
    "tx_id"
  } else {
    NULL
  }
  gene_col <- if ("gene_id" %in% colnames(tx2gene)) {
    "gene_id"
  } else if ("gene_name" %in% colnames(tx2gene)) {
    "gene_name"
  } else {
    NULL
  }
  if (is.null(tx_col) || is.null(gene_col)) {
    stop(.stop("tx2gene must contain 'transcript_id' (or 'tx_id') and 'gene_id' (or 'gene_name') columns."))
  }

  map_df <- data.frame(
    transcript_id = as.character(tx2gene[[tx_col]]),
    gene_id = as.character(tx2gene[[gene_col]]),
    stringsAsFactors = FALSE
  )
  valid_map <- !is.na(map_df$transcript_id) & nzchar(map_df$transcript_id) &
    !is.na(map_df$gene_id) & nzchar(map_df$gene_id)
  map_df <- map_df[valid_map, , drop = FALSE]
  map_df <- map_df[!duplicated(map_df$transcript_id), , drop = FALSE]

  # --- Validate design_or_groups ---
  if (is.matrix(design_or_groups) && is.numeric(design_or_groups)) {
    design <- design_or_groups
    if (nrow(design) != ncol(counts)) {
      stop(.stop(sprintf(
        "Design matrix row count (%d) does not match sample count (%d).",
        nrow(design), ncol(counts)
      )))
    }
    if (ncol(design) < 2L || anyNA(design)) {
      stop(.stop("Design matrix must have at least 2 columns and no missing values."))
    }
  } else if (is.factor(design_or_groups) || is.character(design_or_groups)) {
    if (length(design_or_groups) != ncol(counts)) {
      stop(.stop(sprintf(
        "Group vector length (%d) does not match sample count (%d).",
        length(design_or_groups), ncol(counts)
      )))
    }
    if (anyNA(design_or_groups) || any(!nzchar(as.character(design_or_groups)))) {
      stop(.stop("Group vector must not contain NA or empty values."))
    }
    grp <- if (is.factor(design_or_groups)) {
      droplevels(design_or_groups)
    } else {
      factor(design_or_groups)
    }
    if (nlevels(grp) < 2L) {
      stop(.stop("Group vector must have at least 2 distinct levels for DTU testing."))
    }
    design <- stats::model.matrix(~ grp)
  } else {
    stop(.stop("design_or_groups must be a numeric design matrix or a factor/character group vector."))
  }

  # --- Intersect transcripts and enforce >= 2 isoforms per gene ---
  common_tx <- intersect(tx_names, map_df$transcript_id)
  if (length(common_tx) < 2L) {
    stop(.stop("Fewer than 2 transcripts overlap between counts rownames and tx2gene."))
  }

  tx_to_gene <- stats::setNames(map_df$gene_id, map_df$transcript_id)
  gene_vec <- unname(tx_to_gene[common_tx])

  gene_counts <- table(gene_vec)
  multi_genes <- names(gene_counts)[gene_counts >= 2L]
  if (length(multi_genes) == 0L) {
    stop(.stop("No multi-isoform genes (>= 2 transcripts per gene) found after intersecting counts with tx2gene."))
  }

  keep_multi <- gene_vec %in% multi_genes
  common_tx <- common_tx[keep_multi]
  gene_vec <- gene_vec[keep_multi]
  counts_sub <- counts[common_tx, , drop = FALSE]

  # --- Build DGEList and filter ---
  dge <- edgeR::DGEList(counts = counts_sub)
  keep_expr <- edgeR::filterByExpr(
    dge,
    design = design,
    min.count = min_count,
    min.total.count = min_total_count
  )
  if (sum(keep_expr) < 2L) {
    stop(.stop("Fewer than 2 transcripts remain after edgeR::filterByExpr filtering."))
  }

  dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
  common_tx <- common_tx[keep_expr]
  gene_vec <- gene_vec[keep_expr]

  # Re-enforce >= 2 transcripts per gene after expression filtering
  gene_counts_post <- table(gene_vec)
  multi_genes_post <- names(gene_counts_post)[gene_counts_post >= 2L]
  if (length(multi_genes_post) == 0L) {
    stop(.stop("No multi-isoform genes (>= 2 transcripts per gene) remain after expression filtering."))
  }
  keep_post <- gene_vec %in% multi_genes_post
  if (!all(keep_post)) {
    dge <- dge[keep_post, , keep.lib.sizes = FALSE]
    common_tx <- common_tx[keep_post]
    gene_vec <- gene_vec[keep_post]
  }

  # --- Fit edgeR QLF + diffSpliceDGE ---
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = isTRUE(robust))
  sp <- edgeR::diffSpliceDGE(
    fit,
    coef = coef,
    geneid = gene_vec,
    exonid = common_tx,
    verbose = FALSE
  )

  top_test <- switch(
    test_method,
    Simes = "Simes",
    F = "gene",
    t = "exon"
  )

  gene_table <- edgeR::topSpliceDGE(sp, test = top_test, number = Inf)
  transcript_table <- edgeR::topSpliceDGE(sp, test = "exon", number = Inf)

  list(
    gene_table = gene_table,
    transcript_table = transcript_table,
    n_genes_tested = length(unique(gene_vec)),
    n_transcripts_tested = length(common_tx),
    test_method = test_method,
    provenance = list(
      engine = "edgeR::diffSpliceDGE",
      edgeR_version = as.character(utils::packageVersion("edgeR"))
    )
  )
}

#' Summarize transcript-level DTU results for a single gene
#'
#' @param dtu_res Result list returned by `run_bulk_dtu()`.
#' @param gene_id Character scalar specifying the gene ID to extract.
#' @return A `data.frame` of transcript-level DTU statistics for `gene_id`,
#'   sorted by `FDR` and `P.Value`.
summarize_dtu_gene <- function(dtu_res, gene_id) {
  .stop <- function(msg) {
    errorCondition(msg, class = "bulk_dtu_error")
  }

  if (!is.list(dtu_res) || is.null(dtu_res$transcript_table) ||
      !is.data.frame(dtu_res$transcript_table)) {
    stop(.stop("dtu_res must be a valid result list from run_bulk_dtu() containing transcript_table."))
  }
  if (!is.character(gene_id) || length(gene_id) != 1L || is.na(gene_id) || !nzchar(gene_id)) {
    stop(.stop("gene_id must be a single non-empty character string."))
  }

  tt <- dtu_res$transcript_table
  gene_col <- if ("GeneID" %in% colnames(tt)) {
    "GeneID"
  } else if ("gene_id" %in% colnames(tt)) {
    "gene_id"
  } else {
    stop(.stop("transcript_table does not contain a 'GeneID' or 'gene_id' column."))
  }

  sub <- tt[tt[[gene_col]] == gene_id, , drop = FALSE]
  if (nrow(sub) > 1L) {
    fdr_vals <- if ("FDR" %in% colnames(sub)) sub$FDR else rep(NA_real_, nrow(sub))
    pval_vals <- if ("P.Value" %in% colnames(sub)) sub$P.Value else rep(NA_real_, nrow(sub))
    ord <- order(fdr_vals, pval_vals, na.last = TRUE)
    sub <- sub[ord, , drop = FALSE]
  }
  rownames(sub) <- NULL
  sub
}
