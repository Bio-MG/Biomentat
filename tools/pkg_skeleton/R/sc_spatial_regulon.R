# =============================================================================
# tools/pkg_skeleton/R/sc_spatial_regulon.R — P3.3: decoupleR SC/Spatial regulons
# =============================================================================
# Pure-R TF/pathway regulon activity inference via decoupleR (ULM / MLM)
# supporting both Single-Cell and Spatial normalized expression matrices with
# optional cluster/niche group aggregation. Zero network egress.
# =============================================================================

.regulon_stop <- function(msg) {
  stop(errorCondition(msg, class = "sc_spatial_regulon_error"))
}

#' Infer TF/pathway regulon activities on SC or Spatial expression matrices
#'
#' @param expr_mat Gene-by-cell/spot expression matrix (dense or sparse Matrix)
#'   with rownames (gene symbols) and colnames (cell/spot IDs).
#' @param network Data.frame with columns `source`, `target`, and `mor`.
#' @param method Character, `"ulm"` (default) or `"mlm"`.
#' @param minsize Integer minimum target regulon size (default 5).
#' @param group_by Optional character/factor vector of cluster/niche labels
#'   aligned to `colnames(expr_mat)` (named or positional).
#' @return List with `scores`, `score_matrix`, `group_summary`, `method`, `n_regulons`.
run_sc_spatial_regulons <- function(expr_mat,
                                    network,
                                    method = "ulm",
                                    minsize = 5L,
                                    group_by = NULL) {
  if (is.null(expr_mat) || length(dim(expr_mat)) != 2L ||
      is.null(rownames(expr_mat)) || is.null(colnames(expr_mat))) {
    .regulon_stop("run_sc_spatial_regulons() : expr_mat doit etre une matrice 2D avec rownames et colnames.")
  }
  if (!is.data.frame(network) || !all(c("source", "target", "mor") %in% names(network))) {
    .regulon_stop("run_sc_spatial_regulons() : network doit contenir les colonnes source, target et mor.")
  }
  if (!is.character(method) || length(method) != 1L || !(method %in% c("ulm", "mlm"))) {
    .regulon_stop("run_sc_spatial_regulons() : method doit etre 'ulm' ou 'mlm'.")
  }
  ms <- as.integer(minsize)
  if (is.na(ms) || ms < 1L) {
    .regulon_stop("run_sc_spatial_regulons() : minsize doit etre >= 1.")
  }
  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    .regulon_stop("run_sc_spatial_regulons() : le package 'decoupleR' est requis.")
  }

  mat_in <- if (inherits(expr_mat, "Matrix")) as.matrix(expr_mat) else expr_mat
  fn <- if (identical(method, "ulm")) decoupleR::run_ulm else decoupleR::run_mlm
  raw_res <- fn(
    mat = mat_in,
    network = network,
    .source = "source",
    .target = "target",
    .mor = "mor",
    minsize = ms
  )
  scores_df <- as.data.frame(raw_res, stringsAsFactors = FALSE)
  if (nrow(scores_df) == 0L) {
    .regulon_stop("run_sc_spatial_regulons() : aucun regulon actif apres filtrage minsize.")
  }

  sources <- sort(unique(as.character(scores_df$source)))
  cells <- colnames(expr_mat)
  score_mat <- matrix(NA_real_, nrow = length(sources), ncol = length(cells),
                      dimnames = list(sources, cells))
  r_idx <- match(as.character(scores_df$source), sources)
  c_idx <- match(as.character(scores_df$condition), cells)
  ok <- !is.na(r_idx) & !is.na(c_idx)
  score_mat[cbind(r_idx[ok], c_idx[ok])] <- as.numeric(scores_df$score[ok])

  group_summary <- NULL
  if (!is.null(group_by)) {
    grp <- if (!is.null(names(group_by))) group_by[cells] else group_by
    if (length(grp) != length(cells)) {
      .regulon_stop("run_sc_spatial_regulons() : group_by doit correspondre aux colonnes de expr_mat.")
    }
    grp_chr <- as.character(grp)
    uniq_grp <- sort(unique(grp_chr[!is.na(grp_chr)]))
    rows <- list()
    for (g in uniq_grp) {
      g_idx <- which(grp_chr == g)
      rows[[length(rows) + 1L]] <- data.frame(
        source = sources,
        group = g,
        mean_score = rowMeans(score_mat[, g_idx, drop = FALSE], na.rm = TRUE),
        n_cells = length(g_idx),
        stringsAsFactors = FALSE
      )
    }
    group_summary <- do.call(rbind, rows)
    rownames(group_summary) <- NULL
  }

  list(
    scores = scores_df,
    score_matrix = score_mat,
    group_summary = group_summary,
    method = method,
    n_regulons = length(sources)
  )
}
