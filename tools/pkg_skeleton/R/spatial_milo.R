# =============================================================================
# tools/pkg_skeleton/R/spatial_milo.R — P2.4 (iii): Milo Spatial (DA sur voisinages)
# =============================================================================
# Pure-R spatial neighbourhood differential abundance engine:
#   1. Enforces biological-replicate design validation (blocks pseudoreplication:
#      >= 2 distinct samples per condition, unique sample -> condition mapping).
#   2. Builds k-NN neighbourhoods on a joint spatial/niche embedding (RANN::nn2)
#      with deterministic anchor sampling.
#   3. Aggregates spot counts per biological sample and tests differential
#      abundance (target - ref) via edgeR quasi-likelihood NB GLM.
# =============================================================================

.spatial_milo_stop <- function(msg) {
  stop(errorCondition(msg, class = "spatial_milo_error"))
}

.spatial_milo_validate_design <- function(metadata, sample_col, condition_col, ref_level, target_level) {
  if (!is.data.frame(metadata) || nrow(metadata) == 0L || is.null(rownames(metadata))) {
    .spatial_milo_stop("spatial_milo_da() : metadata doit etre un data.frame non vide avec rownames.")
  }
  if (!(sample_col %in% names(metadata)) || !(condition_col %in% names(metadata))) {
    .spatial_milo_stop("spatial_milo_da() : sample_col ou condition_col absent de metadata.")
  }
  s_vec <- as.character(metadata[[sample_col]])
  c_vec <- as.character(metadata[[condition_col]])
  keep <- !is.na(s_vec) & nzchar(s_vec) & !is.na(c_vec) & nzchar(c_vec)
  if (sum(keep) < 6L) {
    .spatial_milo_stop("spatial_milo_da() : trop peu de spots avec sample_id et condition valides.")
  }
  sc_tab <- unique(data.frame(sample = s_vec[keep], cond = c_vec[keep], stringsAsFactors = FALSE))
  if (anyDuplicated(sc_tab$sample) > 0L) {
    .spatial_milo_stop("spatial_milo_da() : un echantillon est associe a plusieurs conditions (design confondu).")
  }
  sc_sub <- sc_tab[sc_tab$cond %in% c(ref_level, target_level), , drop = FALSE]
  n_ref <- sum(sc_sub$cond == ref_level)
  n_tgt <- sum(sc_sub$cond == target_level)
  if (n_ref < 2L || n_tgt < 2L) {
    .spatial_milo_stop(paste0(
      "spatial_milo_da() : au moins 2 echantillons biologiques par condition sont requis ",
      "(ref=", n_ref, ", target=", n_tgt, ") — pseudoreplication bloquee."
    ))
  }
  list(sample_table = sc_sub, keep_ids = rownames(metadata)[keep & (c_vec %in% c(ref_level, target_level))])
}

#' Spatial neighbourhood differential abundance (Milo spatial)
#'
#' @param embedding Numeric matrix (spots x dims) of joint spatial/niche coordinates
#'   with rownames matching `rownames(metadata)`.
#' @param metadata Data.frame of spot metadata with `sample_col` and `condition_col`.
#' @param sample_col Name of the biological sample column.
#' @param condition_col Name of the experimental condition column.
#' @param ref_level Reference condition level.
#' @param target_level Target condition level.
#' @param k Integer number of nearest neighbours per neighbourhood (default 25).
#' @param prop Numeric fraction of spots sampled as neighbourhood anchors in `(0, 1]`.
#' @param niche_col Optional metadata column to summarize dominant niche per neighbourhood.
#' @param seed Integer seed for deterministic anchor sampling.
#' @return Canonical list with `type = "spatial_milo_da"`, `DA_table`, `sample_table`, and parameters.
spatial_milo_da <- function(embedding,
                            metadata,
                            sample_col,
                            condition_col,
                            ref_level,
                            target_level,
                            k = 25L,
                            prop = 0.2,
                            niche_col = "niche",
                            seed = 42L) {
  if (is.null(embedding) || length(dim(embedding)) != 2L || is.null(rownames(embedding))) {
    .spatial_milo_stop("spatial_milo_da() : embedding doit etre une matrice 2D avec rownames.")
  }
  if (!requireNamespace("RANN", quietly = TRUE)) {
    .spatial_milo_stop("spatial_milo_da() : le package 'RANN' est requis.")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    .spatial_milo_stop("spatial_milo_da() : le package 'edgeR' est requis.")
  }

  design <- .spatial_milo_validate_design(metadata, sample_col, condition_col, ref_level, target_level)
  common_ids <- intersect(design$keep_ids, rownames(embedding))
  if (length(common_ids) < 10L) {
    .spatial_milo_stop("spatial_milo_da() : moins de 10 spots communs entre embedding et metadata.")
  }

  emb <- as.matrix(embedding[common_ids, , drop = FALSE])
  meta <- metadata[common_ids, , drop = FALSE]
  n_spots <- nrow(emb)
  k_int <- min(as.integer(k), n_spots)
  if (is.na(k_int) || k_int < 3L) {
    .spatial_milo_stop("spatial_milo_da() : k doit etre >= 3.")
  }
  n_nhoods <- max(5L, min(n_spots, as.integer(ceiling(prop * n_spots))))

  set.seed(seed)
  anchors <- sort(sample.int(n_spots, size = n_nhoods, replace = FALSE))
  nn <- RANN::nn2(data = emb, query = emb[anchors, , drop = FALSE], k = k_int)

  # Sparse neighbourhood x spot incidence matrix
  nhood_spot <- Matrix::sparseMatrix(
    i = rep(seq_len(n_nhoods), each = k_int),
    j = as.vector(t(nn$nn.idx)),
    x = 1,
    dims = c(n_nhoods, n_spots)
  )

  # Spot x sample indicator matrix
  samp_vec <- factor(as.character(meta[[sample_col]]), levels = design$sample_table$sample)
  spot_samp <- Matrix::sparseMatrix(
    i = seq_len(n_spots),
    j = as.integer(samp_vec),
    x = 1,
    dims = c(n_spots, nlevels(samp_vec)),
    dimnames = list(NULL, levels(samp_vec))
  )

  nhood_counts <- as.matrix(nhood_spot %*% spot_samp)
  rownames(nhood_counts) <- paste0("Nhood_", seq_len(n_nhoods))

  cond_f <- factor(design$sample_table$cond, levels = c(ref_level, target_level))
  mm <- stats::model.matrix(~ cond_f)
  dge <- edgeR::DGEList(counts = nhood_counts)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  dge <- edgeR::estimateDisp(dge, mm)
  fit <- edgeR::glmQLFit(dge, mm, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = 2L)
  tab <- qlf$table
  tab$FDR <- stats::p.adjust(tab$PValue, method = "BH")

  niche_major <- rep(NA_character_, n_nhoods)
  if (!is.null(niche_col) && niche_col %in% names(meta)) {
    niche_vals <- as.character(meta[[niche_col]])
    niche_major <- vapply(seq_len(n_nhoods), function(r) {
      vals <- niche_vals[nn$nn.idx[r, ]]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0L) return(NA_character_)
      names(sort(table(vals), decreasing = TRUE))[1L]
    }, character(1))
  }

  da_df <- data.frame(
    Nhood = rownames(nhood_counts),
    n_spots = as.integer(rowSums(nhood_counts)),
    logFC = tab$logFC,
    logCPM = tab$logCPM,
    F = tab$F,
    PValue = tab$PValue,
    FDR = tab$FDR,
    niche_major = niche_major,
    stringsAsFactors = FALSE
  )
  da_df <- da_df[order(da_df$FDR, da_df$PValue), , drop = FALSE]
  rownames(da_df) <- NULL

  list(
    type = "spatial_milo_da",
    DA_table = da_df,
    sample_table = design$sample_table,
    parameters = list(
      sample_col = sample_col,
      condition_col = condition_col,
      ref_level = ref_level,
      target_level = target_level,
      k = k_int,
      prop = prop,
      seed = seed
    )
  )
}
