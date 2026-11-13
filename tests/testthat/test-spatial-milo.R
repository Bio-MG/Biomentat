# =============================================================================
# test-spatial-milo.R — P2.4 (iii): Milo Spatial (neighbourhood DA)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/spatial_milo.R")

.mk_spatial_milo_fixture <- function(n_samples_per_cond = 3L, spots_per_sample = 60L, seed = 42L) {
  set.seed(seed)
  samples <- c(paste0("CTL_", seq_len(n_samples_per_cond)),
               paste0("TRT_", seq_len(n_samples_per_cond)))
  conds <- c(rep("CTL", n_samples_per_cond), rep("TRT", n_samples_per_cond))

  emb_list <- list()
  meta_list <- list()
  for (i in seq_along(samples)) {
    s <- samples[i]
    cnd <- conds[i]
    # Two spatial niches: Niche1 centered at (0,0), Niche2 centered at (10,10).
    # In TRT, Niche2 is enriched (45 vs 15 spots) and Niche1 depleted (15 vs 45).
    n1 <- if (cnd == "CTL") as.integer(round(0.75 * spots_per_sample)) else as.integer(round(0.25 * spots_per_sample))
    n2 <- spots_per_sample - n1
    xy1 <- cbind(stats::rnorm(n1, 0, 0.8), stats::rnorm(n1, 0, 0.8))
    xy2 <- cbind(stats::rnorm(n2, 10, 0.8), stats::rnorm(n2, 10, 0.8))
    xy <- rbind(xy1, xy2)
    ids <- paste0(s, "_spot_", seq_len(spots_per_sample))
    rownames(xy) <- ids
    colnames(xy) <- c("dim1", "dim2")
    emb_list[[i]] <- xy
    meta_list[[i]] <- data.frame(
      sample_id = s,
      condition = cnd,
      niche = c(rep("Niche1", n1), rep("Niche2", n2)),
      row.names = ids,
      stringsAsFactors = FALSE
    )
  }

  list(
    embedding = do.call(rbind, emb_list),
    metadata = do.call(rbind, meta_list)
  )
}

test_that("spatial_milo_da blocks pseudoreplication and invalid designs", {
  fx <- .mk_spatial_milo_fixture(n_samples_per_cond = 1L)
  err1 <- tryCatch(
    spatial_milo_da(fx$embedding, fx$metadata, sample_col = "sample_id",
                    condition_col = "condition", ref_level = "CTL", target_level = "TRT"),
    error = function(e) e
  )
  expect_s3_class(err1, "spatial_milo_error")
  expect_match(conditionMessage(err1), "2 echantillons", fixed = TRUE)

  fx3 <- .mk_spatial_milo_fixture(n_samples_per_cond = 3L)
  err2 <- tryCatch(
    spatial_milo_da(fx3$embedding, fx3$metadata, sample_col = "missing_col",
                    condition_col = "condition", ref_level = "CTL", target_level = "TRT"),
    error = function(e) e
  )
  expect_s3_class(err2, "spatial_milo_error")
})

test_that("spatial_milo_da detects enriched and depleted spatial neighbourhoods via edgeR", {
  skip_if_not_installed("edgeR")
  fx <- .mk_spatial_milo_fixture(n_samples_per_cond = 3L, spots_per_sample = 60L)
  res <- spatial_milo_da(
    fx$embedding,
    fx$metadata,
    sample_col = "sample_id",
    condition_col = "condition",
    ref_level = "CTL",
    target_level = "TRT",
    k = 20L,
    prop = 0.25,
    seed = 42L
  )

  expect_identical(res$type, "spatial_milo_da")
  expect_true(is.data.frame(res$DA_table))
  expect_true(all(c("Nhood", "n_spots", "logFC", "PValue", "FDR", "niche_major") %in% names(res$DA_table)))
  expect_true(nrow(res$DA_table) >= 10L)

  n1_rows <- res$DA_table[res$DA_table$niche_major == "Niche1", ]
  n2_rows <- res$DA_table[res$DA_table$niche_major == "Niche2", ]
  expect_true(nrow(n1_rows) > 0L && nrow(n2_rows) > 0L)
  # Niche2 is enriched in TRT -> positive logFC; Niche1 is depleted -> negative logFC
  expect_gt(stats::median(n2_rows$logFC), 0.5)
  expect_lt(stats::median(n1_rows$logFC), -0.5)
  expect_true(any(res$DA_table$FDR < 0.05))
})

test_that("spatial_milo_da is deterministic under seed", {
  skip_if_not_installed("edgeR")
  fx <- .mk_spatial_milo_fixture(n_samples_per_cond = 3L, spots_per_sample = 40L)
  r1 <- spatial_milo_da(fx$embedding, fx$metadata, "sample_id", "condition", "CTL", "TRT", k = 15L, prop = 0.2, seed = 99L)
  r2 <- spatial_milo_da(fx$embedding, fx$metadata, "sample_id", "condition", "CTL", "TRT", k = 15L, prop = 0.2, seed = 99L)
  expect_identical(r1$DA_table, r2$DA_table)
})
