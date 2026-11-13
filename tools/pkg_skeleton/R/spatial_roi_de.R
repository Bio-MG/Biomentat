# =============================================================================
# tools/pkg_skeleton/R/spatial_roi_de.R — P2.4 (i): Reactive BPCells ROI DE
# =============================================================================
# Pure-R differential expression for spatial ROIs backed by BPCells on-disk
# IterableMatrix or in-memory Matrix/matrix. Bounds RAM by deterministic
# subsampling BEFORE slicing the on-disk matrix.
# =============================================================================

.roi_de_stop <- function(msg) {
  stop(errorCondition(msg, class = "spatial_roi_de_error"))
}

.roi_de_subsample <- function(ids, max_n, seed) {
  if (length(ids) <= max_n) return(ids)
  set.seed(seed)
  sort(sample(ids, size = max_n, replace = FALSE))
}

.roi_de_slice_sparse <- function(counts_mat, cell_ids) {
  sub_mat <- counts_mat[, cell_ids, drop = FALSE]
  if (inherits(sub_mat, "IterableMatrix")) {
    sub_mat <- methods::as(sub_mat, "dgCMatrix")
  } else if (!inherits(sub_mat, "dgCMatrix")) {
    sub_mat <- Matrix::Matrix(sub_mat, sparse = TRUE)
  }
  sub_mat
}

.roi_de_wilcox_row <- function(x1, x2) {
  if (length(unique(c(x1, x2))) <= 1L) return(1.0)
  suppressWarnings(stats::wilcox.test(x1, x2)$p.value)
}

#' Differential expression for a spatial ROI against a control/background set
#'
#' @param counts_mat Gene-by-cell count matrix (BPCells IterableMatrix, dgCMatrix,
#'   or dense matrix) with row and column names.
#' @param roi_ids Character vector of cell/spot IDs in the ROI.
#' @param control_ids Optional character vector of background/control IDs.
#'   Defaults to all cells in `counts_mat` not in `roi_ids`.
#' @param max_cells_per_group Integer ceiling per group before materialization.
#' @param min_pct Minimum fraction of expressing cells (> 0) in either group.
#' @param logfc_threshold Minimum absolute log2 fold-change to retain.
#' @param seed Integer seed for deterministic subsampling.
#' @return List with `markers` data.frame and group/backend metadata.
spatial_roi_de_bpcells <- function(counts_mat,
                                   roi_ids,
                                   control_ids = NULL,
                                   max_cells_per_group = 5000L,
                                   min_pct = 0.1,
                                   logfc_threshold = 0.25,
                                   seed = 42L) {
  if (is.null(counts_mat) || length(dim(counts_mat)) != 2L ||
      is.null(rownames(counts_mat)) || is.null(colnames(counts_mat))) {
    .roi_de_stop("spatial_roi_de_bpcells() : counts_mat doit etre une matrice 2D avec rownames et colnames.")
  }
  if (is.null(roi_ids) || !is.character(roi_ids) || length(roi_ids) == 0L) {
    .roi_de_stop("spatial_roi_de_bpcells() : roi_ids doit etre un vecteur de caracteres non vide.")
  }

  all_cells <- colnames(counts_mat)
  roi_clean <- intersect(unique(roi_ids), all_cells)
  if (length(roi_clean) < 3L) {
    .roi_de_stop("spatial_roi_de_bpcells() : au moins 3 cellules ROI valides sont requises.")
  }

  ctrl_clean <- if (is.null(control_ids)) {
    setdiff(all_cells, roi_clean)
  } else {
    setdiff(intersect(unique(control_ids), all_cells), roi_clean)
  }
  if (length(ctrl_clean) < 3L) {
    .roi_de_stop("spatial_roi_de_bpcells() : au moins 3 cellules de controle valides sont requises.")
  }

  max_n <- as.integer(max_cells_per_group)
  if (is.na(max_n) || max_n < 3L) {
    .roi_de_stop("spatial_roi_de_bpcells() : max_cells_per_group doit etre >= 3.")
  }

  backend <- if (inherits(counts_mat, "IterableMatrix")) "BPCells" else "Matrix"
  roi_used <- .roi_de_subsample(roi_clean, max_n, seed)
  ctrl_used <- .roi_de_subsample(ctrl_clean, max_n, seed + 1L)

  sub_mat <- .roi_de_slice_sparse(counts_mat, c(roi_used, ctrl_used))
  n1 <- length(roi_used)
  n2 <- length(ctrl_used)
  idx1 <- seq_len(n1)
  idx2 <- n1 + seq_len(n2)

  # Library-size normalize to CPM / 1e4 (Seurat-compatible log1p scale)
  lib_sizes <- Matrix::colSums(sub_mat)
  lib_sizes[lib_sizes <= 0] <- 1
  norm_mat <- Matrix::t(Matrix::t(sub_mat) / lib_sizes) * 1e4
  log_mat <- log1p(norm_mat)

  m1 <- sub_mat[, idx1, drop = FALSE]
  m2 <- sub_mat[, idx2, drop = FALSE]
  pct1 <- Matrix::rowSums(m1 > 0) / n1
  pct2 <- Matrix::rowSums(m2 > 0) / n2

  mean1 <- Matrix::rowMeans(expm1(log_mat[, idx1, drop = FALSE]))
  mean2 <- Matrix::rowMeans(expm1(log_mat[, idx2, drop = FALSE]))
  avg_log2fc <- log2((mean1 + 1) / (mean2 + 1))

  keep <- (pmax(pct1, pct2) >= min_pct) & (abs(avg_log2fc) >= logfc_threshold)
  genes <- rownames(sub_mat)

  if (!any(keep)) {
    empty_df <- data.frame(
      gene = character(0), avg_log2FC = numeric(0),
      pct.1 = numeric(0), pct.2 = numeric(0),
      p_val = numeric(0), p_val_adj = numeric(0),
      stringsAsFactors = FALSE
    )
    return(list(
      markers = empty_df,
      n_roi_total = length(roi_clean), n_roi_used = n1,
      n_ctrl_total = length(ctrl_clean), n_ctrl_used = n2,
      backend = backend
    ))
  }

  keep_idx <- which(keep)
  dense_log <- as.matrix(log_mat[keep_idx, , drop = FALSE])
  p_vals <- vapply(seq_len(nrow(dense_log)), function(i) {
    .roi_de_wilcox_row(dense_log[i, idx1], dense_log[i, idx2])
  }, numeric(1))
  p_vals[!is.finite(p_vals)] <- 1.0
  p_adj <- stats::p.adjust(p_vals, method = "BH")

  markers <- data.frame(
    gene = genes[keep_idx],
    avg_log2FC = unname(avg_log2fc[keep_idx]),
    pct.1 = unname(pct1[keep_idx]),
    pct.2 = unname(pct2[keep_idx]),
    p_val = p_vals,
    p_val_adj = p_adj,
    stringsAsFactors = FALSE
  )
  ord <- order(markers$p_val_adj, -abs(markers$avg_log2FC))
  markers <- markers[ord, , drop = FALSE]
  rownames(markers) <- NULL

  list(
    markers = markers,
    n_roi_total = length(roi_clean),
    n_roi_used = n1,
    n_ctrl_total = length(ctrl_clean),
    n_ctrl_used = n2,
    backend = backend
  )
}
