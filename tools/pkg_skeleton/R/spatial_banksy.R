# =============================================================================
# tools/pkg_skeleton/R/spatial_banksy.R — P2.4 (ii): BANKSY AGF + Multi-sample
# =============================================================================
# Pure-R spatial feature augmentation extending BANKSY-lite with:
#   1. Optional 1st-order Azimuthal Gabor Filter (AGF) gradient magnitude block
#      capturing local spatial expression anisotropy (pure base R + RANN + Matrix).
#   2. Optional per-sample neighborhood isolation (`sample_ids`) so multi-section
#      datasets never connect spots across distinct tissue sections.
# =============================================================================

.banksy_stop <- function(msg) {
  stop(errorCondition(msg, class = "spatial_banksy_error"))
}

.banksy_scale_cols <- function(mat) {
  out <- scale(mat)
  out[!is.finite(out)] <- 0
  out
}

.banksy_sample_blocks <- function(own_sub, xy_sub, k_geom, use_agf) {
  n <- nrow(xy_sub)
  kk <- min(as.integer(k_geom), n - 1L)
  if (kk < 1L) {
    zero_mat <- matrix(0, nrow = n, ncol = ncol(own_sub), dimnames = dimnames(own_sub))
    return(list(nbr = zero_mat, agf = if (use_agf) zero_mat else NULL))
  }

  nn <- RANN::nn2(xy_sub, k = kk + 1L)
  nbr_idx <- nn$nn.idx[, -1L, drop = FALSE]
  nbr_dist <- nn$nn.dists[, -1L, drop = FALSE]

  # Gaussian kernel weights (sigma = median 1-NN distance, floor 1e-6)
  sigma <- stats::median(nbr_dist[, 1L], na.rm = TRUE)
  if (!is.finite(sigma) || sigma <= 1e-6) sigma <- 1.0
  raw_w <- exp(-(nbr_dist^2) / (2 * sigma^2))
  row_s <- rowSums(raw_w)
  row_s[row_s <= 0] <- 1
  w_norm <- raw_w / row_s

  W <- Matrix::sparseMatrix(
    i = rep(seq_len(n), each = kk),
    j = as.vector(t(nbr_idx)),
    x = as.vector(t(w_norm)),
    dims = c(n, n)
  )
  nbr_mat <- as.matrix(W %*% own_sub)
  dimnames(nbr_mat) <- dimnames(own_sub)

  agf_mat <- NULL
  if (isTRUE(use_agf)) {
    # 1st-order azimuthal Gabor gradient: |sum_j w_ij * exp(i * theta_ij) * x_j|
    i_vec <- rep(seq_len(n), each = kk)
    j_vec <- as.vector(t(nbr_idx))
    dx <- xy_sub[j_vec, 1L] - xy_sub[i_vec, 1L]
    dy <- xy_sub[j_vec, 2L] - xy_sub[i_vec, 2L]
    theta <- atan2(dy, dx)
    w_vec <- as.vector(t(w_norm))

    W_cos <- Matrix::sparseMatrix(i = i_vec, j = j_vec, x = w_vec * cos(theta), dims = c(n, n))
    W_sin <- Matrix::sparseMatrix(i = i_vec, j = j_vec, x = w_vec * sin(theta), dims = c(n, n))
    re_part <- as.matrix(W_cos %*% own_sub)
    im_part <- as.matrix(W_sin %*% own_sub)
    agf_mat <- sqrt(re_part^2 + im_part^2)
    dimnames(agf_mat) <- dimnames(own_sub)
  }

  list(nbr = nbr_mat, agf = agf_mat)
}

#' Build BANKSY augmented spatial expression matrix (with optional AGF & multi-sample)
#'
#' @param expr_mat Normalized gene-by-spot matrix with column names.
#' @param coords Data.frame containing `x`, `y` (and `id` or rownames matching `colnames(expr_mat)`).
#' @param k_geom Integer number of spatial nearest neighbors (default 18).
#' @param lambda Numeric spatial mixing weight in `[0, 1]` (default 0.2).
#' @param use_agf Logical; if TRUE, append 1st-order Azimuthal Gabor Filter block.
#' @param agf_weight Numeric fraction of `lambda` allocated to AGF vs neighborhood mean (default 0.5).
#' @param sample_ids Optional character/factor vector of sample labels per spot (isolates kNN per sample).
#' @return List with `augmented` (spots x augmented features), `n_spots`, `n_features`, `use_agf`, `n_samples`.
spatial_banksy_augment <- function(expr_mat,
                                   coords,
                                   k_geom = 18L,
                                   lambda = 0.2,
                                   use_agf = TRUE,
                                   agf_weight = 0.5,
                                   sample_ids = NULL) {
  if (is.null(expr_mat) || length(dim(expr_mat)) != 2L || is.null(colnames(expr_mat))) {
    .banksy_stop("spatial_banksy_augment() : expr_mat doit etre une matrice 2D avec colnames.")
  }
  if (!is.data.frame(coords) || !all(c("x", "y") %in% names(coords))) {
    .banksy_stop("spatial_banksy_augment() : coords doit etre un data.frame avec colonnes x et y.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) || lambda < 0 || lambda > 1) {
    .banksy_stop("spatial_banksy_augment() : lambda doit etre dans [0, 1].")
  }
  k_int <- as.integer(k_geom)
  if (is.na(k_int) || k_int < 1L) {
    .banksy_stop("spatial_banksy_augment() : k_geom doit etre >= 1.")
  }
  if (!requireNamespace("RANN", quietly = TRUE)) {
    .banksy_stop("spatial_banksy_augment() : le package 'RANN' est requis.")
  }

  spots <- colnames(expr_mat)
  coord_ids <- if ("id" %in% names(coords)) coords$id else rownames(coords)
  idx <- match(spots, coord_ids)
  if (anyNA(idx)) {
    .banksy_stop("spatial_banksy_augment() : certains spots de expr_mat sont absents de coords.")
  }
  xy <- as.matrix(coords[idx, c("x", "y"), drop = FALSE])
  if (any(!is.finite(xy))) {
    .banksy_stop("spatial_banksy_augment() : coordonnees x/y non finies detectees.")
  }

  own_mat <- t(as.matrix(expr_mat))
  if (is.null(sample_ids)) {
    grp <- rep("sample_1", nrow(own_mat))
  } else {
    if (length(sample_ids) != nrow(own_mat)) {
      .banksy_stop("spatial_banksy_augment() : length(sample_ids) doit egaler ncol(expr_mat).")
    }
    grp <- as.character(sample_ids)
  }

  nbr_mat <- matrix(0, nrow = nrow(own_mat), ncol = ncol(own_mat), dimnames = dimnames(own_mat))
  agf_mat <- if (isTRUE(use_agf)) {
    matrix(0, nrow = nrow(own_mat), ncol = ncol(own_mat), dimnames = dimnames(own_mat))
  } else NULL

  uniq_samples <- unique(grp)
  for (s in uniq_samples) {
    sel <- which(grp == s)
    blk <- .banksy_sample_blocks(own_mat[sel, , drop = FALSE], xy[sel, , drop = FALSE], k_int, isTRUE(use_agf))
    nbr_mat[sel, ] <- blk$nbr
    if (isTRUE(use_agf)) agf_mat[sel, ] <- blk$agf
  }

  own_sc <- .banksy_scale_cols(own_mat)
  nbr_sc <- .banksy_scale_cols(nbr_mat)

  if (isTRUE(use_agf)) {
    agf_sc <- .banksy_scale_cols(agf_mat)
    w_own <- sqrt(1 - lambda)
    w_nbr <- sqrt(lambda * (1 - agf_weight))
    w_agf <- sqrt(lambda * agf_weight)
    augmented <- cbind(w_own * own_sc, w_nbr * nbr_sc, w_agf * agf_sc)
  } else {
    augmented <- cbind(sqrt(1 - lambda) * own_sc, sqrt(lambda) * nbr_sc)
  }

  list(
    augmented = augmented,
    nbr = nbr_mat,
    agf = agf_mat,
    n_spots = nrow(augmented),
    n_features = ncol(augmented),
    use_agf = isTRUE(use_agf),
    n_samples = length(uniq_samples)
  )
}
