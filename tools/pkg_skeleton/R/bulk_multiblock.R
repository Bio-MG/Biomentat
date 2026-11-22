# =============================================================================
# tools/pkg_skeleton/R/bulk_multiblock.R — Multi-Block DIABLO & NEMO Engine (BK-06)
# =============================================================================
# Pure domain logic (zero Shiny symbols). Implements:
#   1. Supervised multi-block integration (DIABLO-style coupled sGCCA / PLS-DA
#      with L1 feature selection, block correlation diagnostics, and stratified
#      cross-validation; optional mixOmics::block.splsda dispatch when installed).
#   2. Partial-block unsupervised neighborhood fusion (NEMO-style adaptive
#      Gaussian affinity averaging across partially overlapping omics blocks +
#      spectral clustering).
#
# Dependencies: base R + stats only (mixOmics optional via requireNamespace).
# =============================================================================

#' Internal helper to construct a typed bulk_multiblock_error condition
#'
#' @param msg Character scalar error message.
#' @return An error condition of class `"bulk_multiblock_error"`.
#' @keywords internal
.bulk_multiblock_stop <- function(msg) {
  errorCondition(msg, class = "bulk_multiblock_error")
}

#' Validate multi-block omics input list
#'
#' @param blocks Named list of numeric matrices (features in rows, samples in
#'   columns).
#' @param min_blocks Minimum number of blocks required (default 2L).
#' @param min_samples Minimum number of aligned/union samples required.
#' @param partial Logical; when TRUE, allows partial sample overlap (union of
#'   samples). When FALSE, requires >= `min_samples` common samples across all
#'   blocks.
#' @return A list with `blocks`, `common_samples`, and `union_samples`.
#' @keywords internal
.validate_multiblock_blocks <- function(blocks,
                                        min_blocks = 2L,
                                        min_samples = 3L,
                                        partial = FALSE) {
  if (is.null(blocks) || !is.list(blocks) || length(blocks) < min_blocks) {
    stop(.bulk_multiblock_stop(sprintf(
      "Input 'blocks' must be a named list of at least %d numeric matrices.",
      min_blocks
    )))
  }
  bnames <- names(blocks)
  if (is.null(bnames) || anyNA(bnames) || any(!nzchar(bnames)) || anyDuplicated(bnames)) {
    stop(.bulk_multiblock_stop("Input 'blocks' must have unique, non-empty block names."))
  }

  cleaned <- vector("list", length(blocks))
  names(cleaned) <- bnames

  for (i in seq_along(blocks)) {
    mat <- blocks[[i]]
    nm <- bnames[i]
    if (is.data.frame(mat)) {
      if (!all(vapply(mat, is.numeric, logical(1L)))) {
        stop(.bulk_multiblock_stop(sprintf(
          "Block '%s' data.frame must contain only numeric columns.", nm
        )))
      }
      mat <- as.matrix(mat)
    }
    if (is.null(mat) || !is.matrix(mat) || !is.numeric(mat)) {
      stop(.bulk_multiblock_stop(sprintf(
        "Block '%s' must be a numeric matrix or numeric data.frame.", nm
      )))
    }
    if (nrow(mat) < 1L || ncol(mat) < 2L) {
      stop(.bulk_multiblock_stop(sprintf(
        "Block '%s' must have at least 1 feature (row) and 2 samples (columns).", nm
      )))
    }
    rn <- rownames(mat)
    cn <- colnames(mat)
    if (is.null(rn) || anyNA(rn) || any(!nzchar(rn))) {
      stop(.bulk_multiblock_stop(sprintf(
        "Block '%s' must have non-empty feature rownames.", nm
      )))
    }
    if (is.null(cn) || anyNA(cn) || any(!nzchar(cn)) || anyDuplicated(cn)) {
      stop(.bulk_multiblock_stop(sprintf(
        "Block '%s' must have unique, non-empty sample colnames.", nm
      )))
    }
    if (any(is.infinite(mat))) {
      stop(.bulk_multiblock_stop(sprintf(
        "Block '%s' must not contain Inf/-Inf values.", nm
      )))
    }
    cleaned[[i]] <- mat
  }

  sample_lists <- lapply(cleaned, colnames)
  common_samples <- Reduce(intersect, sample_lists)
  union_samples <- unique(unlist(sample_lists, use.names = FALSE))

  if (!isTRUE(partial) && length(common_samples) < min_samples) {
    stop(.bulk_multiblock_stop(sprintf(
      "Need at least %d common samples across all blocks (found %d).",
      min_samples, length(common_samples)
    )))
  }
  if (isTRUE(partial) && length(union_samples) < min_samples) {
    stop(.bulk_multiblock_stop(sprintf(
      "Need at least %d total union samples across blocks (found %d).",
      min_samples, length(union_samples)
    )))
  }

  list(
    blocks = cleaned,
    common_samples = common_samples,
    union_samples = union_samples
  )
}

#' Filter zero-variance features, impute isolated NAs, and standardize block
#'
#' @param mat Numeric matrix (features x samples).
#' @param block_name Character label used in error messages.
#' @return Standardized sample-by-feature matrix (`samples x features`).
#' @keywords internal
.filter_and_scale_block <- function(mat, block_name = "block") {
  x <- t(mat)
  # Impute isolated NAs with feature mean before variance check
  if (anyNA(x)) {
    for (j in seq_len(ncol(x))) {
      na_idx <- is.na(x[, j])
      if (all(na_idx)) {
        x[, j] <- 0
      } else if (any(na_idx)) {
        x[na_idx, j] <- mean(x[!na_idx, j])
      }
    }
  }
  sds <- apply(x, 2L, stats::sd)
  keep <- is.finite(sds) & (sds > 1e-12)
  if (!any(keep)) {
    stop(.bulk_multiblock_stop(sprintf(
      "Block '%s' has zero non-constant features on the selected samples.",
      block_name
    )))
  }
  x_sub <- x[, keep, drop = FALSE]
  means <- colMeans(x_sub)
  sds_sub <- apply(x_sub, 2L, stats::sd)
  scaled <- sweep(x_sub, 2L, means, "-")
  scaled <- sweep(scaled, 2L, sds_sub, "/")
  scaled
}

#' Sparse unit-norm projection keeping at most `keep_x` top loadings
#'
#' @param w Numeric loading vector.
#' @param keep_x Optional positive integer number of non-zero entries to retain.
#' @return L2-normalized numeric vector.
#' @keywords internal
.soft_threshold_keepx <- function(w, keep_x = NULL) {
  if (!is.null(keep_x) && is.finite(keep_x) && keep_x >= 1L && keep_x < length(w)) {
    k <- as.integer(keep_x)
    ord <- order(abs(w), decreasing = TRUE)
    w_sparse <- numeric(length(w))
    names(w_sparse) <- names(w)
    keep_idx <- ord[seq_len(k)]
    w_sparse[keep_idx] <- w[keep_idx]
    w <- w_sparse
  }
  nrm <- sqrt(sum(w * w))
  if (!is.finite(nrm) || nrm <= 1e-12) {
    out <- rep(0, length(w))
    names(out) <- names(w)
    if (length(out) >= 1L) out[1L] <- 1
    return(out)
  }
  w / nrm
}

#' Align outcome vector y to common sample IDs
#'
#' @param y Factor or character vector of class labels.
#' @param common_samples Character vector of common sample IDs.
#' @return Factor aligned to `common_samples` with >= 2 levels.
#' @keywords internal
.align_diablo_outcome <- function(y, common_samples) {
  if (is.null(y) || (!is.factor(y) && !is.character(y) && !is.numeric(y))) {
    stop(.bulk_multiblock_stop("Input 'y' must be a factor or character vector of class labels."))
  }
  if (!is.null(names(y))) {
    missing_s <- setdiff(common_samples, names(y))
    if (length(missing_s) > 0L) {
      stop(.bulk_multiblock_stop(sprintf(
        "Outcome 'y' is missing %d common sample(s) (e.g. '%s').",
        length(missing_s), missing_s[1L]
      )))
    }
    y_aligned <- y[common_samples]
  } else {
    if (length(y) != length(common_samples)) {
      stop(.bulk_multiblock_stop(sprintf(
        "Unnamed outcome 'y' length (%d) must equal number of common samples (%d).",
        length(y), length(common_samples)
      )))
    }
    y_aligned <- y
  }
  y_chr <- as.character(y_aligned)
  if (anyNA(y_chr) || any(!nzchar(y_chr))) {
    stop(.bulk_multiblock_stop("Outcome 'y' must not contain NA or empty labels."))
  }
  y_fac <- factor(y_chr)
  if (nlevels(y_fac) < 2L) {
    stop(.bulk_multiblock_stop("Outcome 'y' must have at least 2 distinct classes."))
  }
  names(y_fac) <- common_samples
  y_fac
}

#' Fit pure-R coupled multi-block sGCCA / DIABLO model
#'
#' @keywords internal
.fit_pure_r_diablo <- function(x_list,
                               y_fac,
                               ncomp,
                               keep_vec,
                               design_weight,
                               max_iter,
                               tol) {
  k_blocks <- length(x_list)
  bnames <- names(x_list)
  n_samples <- nrow(x_list[[1L]])
  snames <- rownames(x_list[[1L]])

  # Centered one-hot outcome indicator matrix
  y_ind <- stats::model.matrix(~ y_fac - 1)
  y_ind <- scale(y_ind, center = TRUE, scale = FALSE)

  x_res <- x_list
  y_res <- y_ind

  variates <- lapply(x_list, function(x) {
    m <- matrix(0, nrow = n_samples, ncol = ncomp)
    rownames(m) <- snames
    colnames(m) <- paste0("comp", seq_len(ncomp))
    m
  })
  names(variates) <- bnames

  loadings <- lapply(x_list, function(x) {
    m <- matrix(0, nrow = ncol(x), ncol = ncomp)
    rownames(m) <- colnames(x)
    colnames(m) <- paste0("comp", seq_len(ncomp))
    m
  })
  names(loadings) <- bnames

  for (h in seq_len(ncomp)) {
    w_list <- vector("list", k_blocks)
    t_list <- vector("list", k_blocks)

    # Initialize weights from leading block-outcome covariance
    for (k in seq_len(k_blocks)) {
      cov_xy <- crossprod(x_res[[k]], y_res)
      sv <- svd(cov_xy, nu = 1L, nv = 1L)
      w0 <- as.numeric(sv$u[, 1L])
      names(w0) <- colnames(x_res[[k]])
      w_list[[k]] <- .soft_threshold_keepx(w0, keep_vec[k])
      t_list[[k]] <- as.numeric(x_res[[k]] %*% w_list[[k]])
    }

    # Iterative RGCCA-style coupling between blocks and outcome
    for (it in seq_len(max_iter)) {
      max_delta <- 0
      for (k in seq_len(k_blocks)) {
        # Outcome direction weighted against y_res
        y_proj <- y_res %*% crossprod(y_res, t_list[[k]])
        y_norm <- sqrt(sum(y_proj * y_proj))
        if (y_norm > 1e-12) y_proj <- y_proj / y_norm

        # Inter-block coupling sum
        block_couple <- numeric(n_samples)
        if (design_weight > 0 && k_blocks > 1L) {
          for (j in seq_len(k_blocks)) {
            if (j != k) {
              block_couple <- block_couple + design_weight * t_list[[j]]
            }
          }
        }
        target_vec <- as.numeric(y_proj) + block_couple
        w_raw <- as.numeric(crossprod(x_res[[k]], target_vec))
        names(w_raw) <- colnames(x_res[[k]])
        w_new <- .soft_threshold_keepx(w_raw, keep_vec[k])
        delta <- max(abs(w_new - w_list[[k]]))
        if (delta > max_delta) max_delta <- delta
        w_list[[k]] <- w_new
        t_list[[k]] <- as.numeric(x_res[[k]] %*% w_new)
      }
      if (max_delta < tol) break
    }

    # Store variates/loadings and deflate residuals
    for (k in seq_len(k_blocks)) {
      tk <- t_list[[k]]
      variates[[k]][, h] <- tk
      loadings[[k]][, h] <- w_list[[k]]
      denom <- sum(tk * tk)
      if (denom > 1e-12) {
        pk <- as.numeric(crossprod(x_res[[k]], tk)) / denom
        x_res[[k]] <- x_res[[k]] - tcrossprod(tk, pk)
      }
    }
    t_mean <- Reduce("+", t_list) / k_blocks
    denom_y <- sum(t_mean * t_mean)
    if (denom_y > 1e-12) {
      qy <- as.numeric(crossprod(y_res, t_mean)) / denom_y
      y_res <- y_res - tcrossprod(t_mean, qy)
    }
  }

  consensus <- Reduce("+", variates) / k_blocks
  rownames(consensus) <- snames
  colnames(consensus) <- paste0("comp", seq_len(ncomp))

  list(
    variates = variates,
    loadings = loadings,
    consensus_variates = consensus
  )
}

#' Predict classes on held-out samples using fitted DIABLO block loadings
#'
#' @keywords internal
.predict_diablo_nearest_centroid <- function(fit, x_train_list, y_train, x_test_list) {
  k_blocks <- length(x_train_list)
  n_test <- nrow(x_test_list[[1L]])
  ncomp <- ncol(fit$consensus_variates)

  test_variates <- lapply(seq_len(k_blocks), function(k) {
    x_test_list[[k]] %*% fit$loadings[[k]]
  })
  test_consensus <- Reduce("+", test_variates) / k_blocks

  lvls <- levels(y_train)
  centroids <- matrix(0, nrow = length(lvls), ncol = ncomp)
  for (c_idx in seq_along(lvls)) {
    idx <- which(y_train == lvls[c_idx])
    centroids[c_idx, ] <- colMeans(fit$consensus_variates[idx, , drop = FALSE])
  }

  preds <- character(n_test)
  for (i in seq_len(n_test)) {
    dists <- apply(centroids, 1L, function(ct) sum((test_consensus[i, ] - ct)^2))
    preds[i] <- lvls[which.min(dists)]
  }
  factor(preds, levels = lvls)
}

#' Run supervised multi-block integration (DIABLO-style, BK-06)
#'
#' Integrates two or more bulk/omics matrices measured on overlapping samples
#' against a categorical outcome `y`, extracting correlated latent variates,
#' sparse discriminative feature loadings, inter-block correlations, and
#' stratified cross-validation accuracy.
#'
#' @param blocks Named list of numeric matrices (features in rows, samples in
#'   columns).
#' @param y Factor or character vector of class labels (named by sample IDs or
#'   matching common sample count).
#' @param ncomp Positive integer number of latent components (default `2L`).
#' @param keepX Optional named integer vector (or scalar) specifying maximum
#'   features retained per block per component.
#' @param design_weight Numeric coupling weight in `[0, 1]` between omics blocks
#'   (default `0.1`).
#' @param cv_folds Positive integer number of stratified CV folds (default `5L`).
#' @param max_iter Maximum inner coupling iterations (default `50L`).
#' @param tol Convergence tolerance on loading vectors (default `1e-6`).
#' @param use_mixomics Logical; when `TRUE` and `mixOmics` is installed, calls
#'   `mixOmics::block.splsda`. Defaults to `FALSE` (pure-R engine).
#' @param seed Optional integer seed for stratified CV fold assignment.
#' @return A list of class `"biomentat_diablo"` containing `variates`,
#'   `loadings`, `consensus_variates`, `selected_features`,
#'   `block_correlations`, `cv_accuracy`, `cv_balanced_accuracy`, `predictions`,
#'   `common_samples`, `ncomp`, and `method_used`.
#' @export
run_multiblock_diablo <- function(blocks,
                                  y,
                                  ncomp = 2L,
                                  keepX = NULL,
                                  design_weight = 0.1,
                                  cv_folds = 5L,
                                  max_iter = 50L,
                                  tol = 1e-6,
                                  use_mixomics = FALSE,
                                  seed = 42L) {
  val <- .validate_multiblock_blocks(blocks, min_blocks = 2L, min_samples = 3L, partial = FALSE)
  common_s <- val$common_samples
  y_fac <- .align_diablo_outcome(y, common_s)

  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) || ncomp < 1L) {
    stop(.bulk_multiblock_stop("'ncomp' must be a positive integer."))
  }
  if (!is.numeric(design_weight) || length(design_weight) != 1L ||
      is.na(design_weight) || design_weight < 0 || design_weight > 1) {
    stop(.bulk_multiblock_stop("'design_weight' must be a scalar in [0, 1]."))
  }

  # Subset blocks to common samples and standardize features
  bnames <- names(val$blocks)
  x_list <- lapply(seq_along(val$blocks), function(k) {
    sub_mat <- val$blocks[[k]][, common_s, drop = FALSE]
    .filter_and_scale_block(sub_mat, block_name = bnames[k])
  })
  names(x_list) <- bnames

  min_p <- min(vapply(x_list, ncol, integer(1L)))
  ncomp_eff <- as.integer(min(as.integer(ncomp), length(common_s) - 1L, min_p))
  if (ncomp_eff < 1L) {
    stop(.bulk_multiblock_stop("Insufficient samples or features to extract 1 component."))
  }

  keep_vec <- rep(NA_integer_, length(x_list))
  names(keep_vec) <- bnames
  if (!is.null(keepX)) {
    if (length(keepX) == 1L && is.null(names(keepX))) {
      keep_vec[] <- as.integer(keepX)
    } else {
      for (bn in bnames) {
        if (!is.null(keepX[bn]) && !is.na(keepX[bn])) {
          keep_vec[bn] <- as.integer(keepX[bn])
        }
      }
    }
  }

  method_used <- "pure_r_sgcca_diablo"
  fit <- NULL

  if (isTRUE(use_mixomics) && requireNamespace("mixOmics", quietly = TRUE)) {
    design_mat <- matrix(design_weight, nrow = length(x_list), ncol = length(x_list))
    diag(design_mat) <- 0
    keep_arg <- lapply(seq_along(x_list), function(k) {
      kv <- if (is.na(keep_vec[k])) ncol(x_list[[k]]) else min(keep_vec[k], ncol(x_list[[k]]))
      rep(kv, ncomp_eff)
    })
    names(keep_arg) <- bnames
    mx <- mixOmics::block.splsda(
      X = x_list, Y = y_fac, ncomp = ncomp_eff,
      keepX = keep_arg, design = design_mat
    )
    fit <- list(
      variates = mx$variates[bnames],
      loadings = mx$loadings[bnames],
      consensus_variates = Reduce("+", mx$variates[bnames]) / length(bnames)
    )
    method_used <- "mixOmics_block_splsda"
  } else {
    fit <- .fit_pure_r_diablo(
      x_list = x_list,
      y_fac = y_fac,
      ncomp = ncomp_eff,
      keep_vec = keep_vec,
      design_weight = design_weight,
      max_iter = max_iter,
      tol = tol
    )
  }

  selected_features <- lapply(fit$loadings, function(ld) {
    nonzero <- apply(abs(ld) > 1e-10, 1L, any)
    rownames(ld)[nonzero]
  })

  # Pairwise Pearson correlations between block variates on component 1
  k_blocks <- length(bnames)
  block_cor <- matrix(1, nrow = k_blocks, ncol = k_blocks,
                      dimnames = list(bnames, bnames))
  for (i in seq_len(k_blocks)) {
    for (j in seq_len(k_blocks)) {
      if (i < j) {
        v1 <- fit$variates[[i]][, 1L]
        v2 <- fit$variates[[j]][, 1L]
        r_val <- if (stats::sd(v1) > 1e-12 && stats::sd(v2) > 1e-12) {
          stats::cor(v1, v2)
        } else {
          0
        }
        block_cor[i, j] <- r_val
        block_cor[j, i] <- r_val
      }
    }
  }

  # Stratified cross-validation on consensus centroid classifier
  min_class_n <- min(as.integer(table(y_fac)))
  k_folds <- min(as.integer(cv_folds), min_class_n)
  cv_acc <- NA_real_
  cv_bal_acc <- NA_real_
  cv_preds <- factor(rep(NA_character_, length(y_fac)), levels = levels(y_fac))
  names(cv_preds) <- common_s

  if (k_folds >= 2L) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    if (!is.null(seed)) set.seed(as.integer(seed))

    fold_id <- integer(length(y_fac))
    for (lvl in levels(y_fac)) {
      idx <- which(y_fac == lvl)
      perm <- sample(idx)
      fold_id[perm] <- rep(seq_len(k_folds), length.out = length(perm))
    }

    for (f in seq_len(k_folds)) {
      tr_idx <- which(fold_id != f)
      te_idx <- which(fold_id == f)
      x_tr <- lapply(seq_along(val$blocks), function(k) {
        .filter_and_scale_block(val$blocks[[k]][, common_s[tr_idx], drop = FALSE], bnames[k])
      })
      names(x_tr) <- bnames
      x_te <- lapply(seq_along(val$blocks), function(k) {
        feats <- colnames(x_tr[[k]])
        sub_te <- t(val$blocks[[k]][feats, common_s[te_idx], drop = FALSE])
        tr_raw <- t(val$blocks[[k]][feats, common_s[tr_idx], drop = FALSE])
        tr_means <- colMeans(tr_raw, na.rm = TRUE)
        tr_means[!is.finite(tr_means)] <- 0
        if (anyNA(sub_te)) {
          for (j in seq_len(ncol(sub_te))) {
            na_j <- is.na(sub_te[, j])
            if (any(na_j)) sub_te[na_j, j] <- tr_means[j]
          }
        }
        tr_sds <- apply(tr_raw, 2L, stats::sd, na.rm = TRUE)
        tr_sds[tr_sds <= 1e-12 | !is.finite(tr_sds)] <- 1
        sweep(sweep(sub_te, 2L, tr_means, "-"), 2L, tr_sds, "/")
      })
      names(x_te) <- bnames

      fold_ncomp <- min(ncomp_eff, length(tr_idx) - 1L, min(vapply(x_tr, ncol, integer(1L))))
      fold_fit <- .fit_pure_r_diablo(
        x_list = x_tr,
        y_fac = droplevels(y_fac[tr_idx]),
        ncomp = max(1L, fold_ncomp),
        keep_vec = keep_vec,
        design_weight = design_weight,
        max_iter = max_iter,
        tol = tol
      )
      cv_preds[te_idx] <- .predict_diablo_nearest_centroid(
        fold_fit, x_tr, y_fac[tr_idx], x_te
      )
    }

    valid <- !is.na(cv_preds)
    if (any(valid)) {
      cv_acc <- mean(cv_preds[valid] == y_fac[valid])
      recalls <- vapply(levels(y_fac), function(lvl) {
        idx <- which(y_fac == lvl & valid)
        if (!length(idx)) return(NA_real_)
        mean(cv_preds[idx] == lvl)
      }, numeric(1L))
      cv_bal_acc <- mean(recalls, na.rm = TRUE)
    }
  }

  structure(
    list(
      variates = fit$variates,
      loadings = fit$loadings,
      consensus_variates = fit$consensus_variates,
      selected_features = selected_features,
      block_correlations = block_cor,
      cv_accuracy = cv_acc,
      cv_balanced_accuracy = cv_bal_acc,
      predictions = cv_preds,
      y = y_fac,
      common_samples = common_s,
      ncomp = ncomp_eff,
      method_used = method_used
    ),
    class = "biomentat_diablo"
  )
}

#' Compute per-block NEMO local Gaussian relative similarity matrix
#'
#' @param x_scaled Standardized sample-by-feature matrix (`samples x features`).
#' @param num_neighbors Neighborhood size for adaptive local bandwidth.
#' @return Row-normalized relative similarity matrix (`samples x samples`).
#' @keywords internal
.nemo_block_relative_similarity <- function(x_scaled, num_neighbors = NULL) {
  n <- nrow(x_scaled)
  snames <- rownames(x_scaled)
  d_mat <- as.matrix(stats::dist(x_scaled, method = "euclidean"))
  k_nn <- if (is.null(num_neighbors)) {
    min(10L, max(2L, floor(n / 6)))
  } else {
    min(as.integer(num_neighbors), n - 1L)
  }
  k_nn <- max(1L, k_nn)

  sigma <- numeric(n)
  nn_mask <- matrix(FALSE, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    other_idx <- setdiff(seq_len(n), i)
    ord_i <- other_idx[order(d_mat[i, other_idx])]
    nbrs <- ord_i[seq_len(min(k_nn, length(ord_i)))]
    nn_mask[i, nbrs] <- TRUE
    sigma[i] <- mean(d_mat[i, nbrs])
  }
  sigma[!is.finite(sigma) | sigma <= 1e-12] <- stats::median(d_mat[d_mat > 0])
  sigma[!is.finite(sigma) | sigma <= 1e-12] <- 1

  sig_prod <- outer(sigma, sigma, "*")
  w <- exp(-(d_mat^2) / (sig_prod + 1e-12))
  diag(w) <- 0
  # Restrict each row to its k-nearest neighbors (NEMO Eq. 1-2 relative similarity)
  w[!nn_mask] <- 0

  row_sums <- rowSums(w)
  row_sums[row_sums <= 1e-12] <- 1
  p <- sweep(w, 1L, row_sums, "/")
  p <- 0.5 * (p + t(p))
  dimnames(p) <- list(snames, snames)
  p
}

#' Compute mean silhouette score from dissimilarity matrix and cluster vector
#'
#' @keywords internal
.mean_silhouette_from_dist <- function(dist_mat, clusters) {
  u_cl <- sort(unique(clusters))
  if (length(u_cl) < 2L || length(clusters) < 3L) return(NA_real_)
  n <- length(clusters)
  s_vals <- numeric(n)
  for (i in seq_len(n)) {
    ci <- clusters[i]
    same <- which(clusters == ci & seq_len(n) != i)
    a_i <- if (length(same) > 0L) mean(dist_mat[i, same]) else 0
    b_i <- Inf
    for (cj in u_cl) {
      if (cj != ci) {
        other <- which(clusters == cj)
        d_other <- mean(dist_mat[i, other])
        if (d_other < b_i) b_i <- d_other
      }
    }
    denom <- max(a_i, b_i)
    s_vals[i] <- if (denom > 1e-12) (b_i - a_i) / denom else 0
  }
  mean(s_vals)
}

#' Run NEMO-style partial-block multi-omics integration & spectral clustering (BK-06)
#'
#' Fuses multi-block omics matrices with partial sample coverage (some samples
#' may be missing from one or more blocks) by computing per-block adaptive
#' relative neighborhood affinities and averaging each pair `(i, j)` only over
#' blocks where both samples were measured.
#'
#' @param blocks Named list of numeric matrices (features in rows, samples in
#'   columns). Samples do not need to be identical across blocks.
#' @param k Optional integer number of clusters. When `NULL`, selected via
#'   eigengap heuristic on `[2, max_k]`.
#' @param max_k Maximum number of clusters considered by eigengap heuristic
#'   (default `6L`).
#' @param num_neighbors Optional neighborhood size for local Gaussian scaling.
#' @param seed Optional integer seed for deterministic spectral k-means.
#' @return A list of class `"biomentat_nemo"` containing `clusters`,
#'   `n_clusters`, `fused_similarity`, `n_shared_blocks`, `block_coverage`,
#'   `eigengap_scores`, `embedding`, `silhouette_mean`, `union_samples`, and
#'   `method_used`.
#' @export
run_multiblock_nemo <- function(blocks,
                                k = NULL,
                                max_k = 6L,
                                num_neighbors = NULL,
                                seed = 42L) {
  val <- .validate_multiblock_blocks(blocks, min_blocks = 2L, min_samples = 4L, partial = TRUE)
  union_s <- val$union_samples
  n_total <- length(union_s)
  bnames <- names(val$blocks)
  k_blocks <- length(bnames)

  if (!is.null(k)) {
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2L || k >= n_total) {
      stop(.bulk_multiblock_stop(sprintf(
        "'k' must be an integer in [2, %d].", n_total - 1L
      )))
    }
  }

  block_coverage <- matrix(FALSE, nrow = n_total, ncol = k_blocks,
                           dimnames = list(union_s, bnames))
  sim_sum <- matrix(0, nrow = n_total, ncol = n_total,
                    dimnames = list(union_s, union_s))
  n_shared <- matrix(0L, nrow = n_total, ncol = n_total,
                     dimnames = list(union_s, union_s))

  for (b in seq_len(k_blocks)) {
    mat_b <- val$blocks[[b]]
    s_b <- colnames(mat_b)
    block_coverage[s_b, b] <- TRUE
    x_b <- .filter_and_scale_block(mat_b, block_name = bnames[b])
    p_b <- .nemo_block_relative_similarity(x_b, num_neighbors = num_neighbors)
    sim_sum[s_b, s_b] <- sim_sum[s_b, s_b] + p_b
    n_shared[s_b, s_b] <- n_shared[s_b, s_b] + 1L
  }

  fused <- matrix(0, nrow = n_total, ncol = n_total,
                  dimnames = list(union_s, union_s))
  pos <- n_shared > 0L
  fused[pos] <- sim_sum[pos] / n_shared[pos]
  fused <- 0.5 * (fused + t(fused))
  diag(fused) <- 1

  # Normalized symmetric Laplacian L_sym = I - D^{-1/2} W D^{-1/2}
  w_no_diag <- fused
  diag(w_no_diag) <- 0
  deg <- rowSums(w_no_diag)
  d_inv_sqrt <- ifelse(deg > 1e-12, 1 / sqrt(deg), 0)
  l_norm <- diag(n_total) - (w_no_diag * outer(d_inv_sqrt, d_inv_sqrt, "*"))
  l_norm <- 0.5 * (l_norm + t(l_norm))

  eig <- eigen(l_norm, symmetric = TRUE)
  # Eigenvalues ascending order
  ord <- order(eig$values, decreasing = FALSE)
  evals <- eig$values[ord]
  evecs <- eig$vectors[, ord, drop = FALSE]

  max_k_eff <- min(as.integer(max_k), n_total - 1L)
  eigengaps <- numeric(0)
  if (max_k_eff >= 2L) {
    upper <- min(max_k_eff + 1L, length(evals))
    if (upper >= 3L) {
      eigengaps <- diff(evals[2L:upper])
      names(eigengaps) <- paste0("k=", seq(2L, length.out = length(eigengaps)))
    }
  }

  n_clusters <- if (!is.null(k)) {
    as.integer(k)
  } else if (length(eigengaps) >= 1L) {
    as.integer(which.max(eigengaps) + 1L)
  } else {
    2L
  }
  n_clusters <- max(2L, min(n_clusters, n_total - 1L))

  emb <- evecs[, seq_len(n_clusters), drop = FALSE]
  row_norms <- sqrt(rowSums(emb^2))
  row_norms[row_norms <= 1e-12] <- 1
  emb_norm <- sweep(emb, 1L, row_norms, "/")
  rownames(emb_norm) <- union_s
  colnames(emb_norm) <- paste0("spec", seq_len(n_clusters))

  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (!is.null(seed)) set.seed(as.integer(seed))

  km <- stats::kmeans(emb_norm, centers = n_clusters, nstart = 20L, iter.max = 100L)
  clusters <- as.integer(km$cluster)
  names(clusters) <- union_s

  emb_dist <- as.matrix(stats::dist(emb_norm, method = "euclidean"))
  sil_mean <- .mean_silhouette_from_dist(emb_dist, clusters)

  structure(
    list(
      clusters = clusters,
      n_clusters = n_clusters,
      fused_similarity = fused,
      n_shared_blocks = n_shared,
      block_coverage = block_coverage,
      eigengap_scores = eigengaps,
      embedding = emb_norm,
      silhouette_mean = sil_mean,
      union_samples = union_s,
      method_used = "pure_r_nemo_partial"
    ),
    class = "biomentat_nemo"
  )
}
