# =============================================================================
# test-spatial-banksy.R — P2.4 (ii): BANKSY AGF + Multi-sample
# =============================================================================

source_project_file("tools/pkg_skeleton/R/spatial_banksy.R")

.mk_banksy_fixture <- function(n_genes = 8L, n_spots = 30L, seed = 42L) {
  set.seed(seed)
  coords <- data.frame(
    id = paste0("S", seq_len(n_spots)),
    x = stats::runif(n_spots, 0, 100),
    y = stats::runif(n_spots, 0, 100),
    stringsAsFactors = FALSE
  )
  expr <- matrix(stats::rnorm(n_genes * n_spots), nrow = n_genes, ncol = n_spots)
  # Spatial gradient along x for gene 1 to induce directional AGF signal
  expr[1L, ] <- coords$x / 10 + stats::rnorm(n_spots, sd = 0.05)
  rownames(expr) <- paste0("G", seq_len(n_genes))
  colnames(expr) <- coords$id
  list(expr = expr, coords = coords)
}

test_that("spatial_banksy_augment validates coords, lambda, and k_geom", {
  fx <- .mk_banksy_fixture()

  err1 <- tryCatch(spatial_banksy_augment(fx$expr, fx$coords, lambda = 1.5), error = function(e) e)
  expect_s3_class(err1, "spatial_banksy_error")

  err2 <- tryCatch(spatial_banksy_augment(fx$expr, fx$coords[, c("id", "x")]), error = function(e) e)
  expect_s3_class(err2, "spatial_banksy_error")

  err3 <- tryCatch(spatial_banksy_augment(fx$expr, fx$coords, k_geom = 0L), error = function(e) e)
  expect_s3_class(err3, "spatial_banksy_error")
})

test_that("spatial_banksy_augment supports 2-block (BANKSY-lite) and 3-block (AGF) modes", {
  fx <- .mk_banksy_fixture()
  res_lite <- spatial_banksy_augment(fx$expr, fx$coords, k_geom = 6L, lambda = 0.2, use_agf = FALSE)
  res_agf  <- spatial_banksy_augment(fx$expr, fx$coords, k_geom = 6L, lambda = 0.2, use_agf = TRUE)

  expect_identical(nrow(res_lite$augmented), 30L)
  expect_identical(ncol(res_lite$augmented), 16L)
  expect_false(res_lite$use_agf)

  expect_identical(nrow(res_agf$augmented), 30L)
  expect_identical(ncol(res_agf$augmented), 24L)
  expect_true(res_agf$use_agf)
  expect_true(all(is.finite(res_agf$augmented)))
  # AGF block (cols 17..24) must be non-degenerate for spatially varying genes
  expect_gt(stats::sd(res_agf$augmented[, 17L]), 0)
})

test_that("spatial_banksy_augment isolates kNN neighborhoods per sample_ids", {
  fx_a <- .mk_banksy_fixture(n_genes = 6L, n_spots = 20L, seed = 10L)
  fx_b <- .mk_banksy_fixture(n_genes = 6L, n_spots = 20L, seed = 20L)
  # Deliberately overlap coordinates of sample B with sample A
  fx_b$coords$x <- fx_a$coords$x
  fx_b$coords$y <- fx_a$coords$y
  fx_b$coords$id <- paste0("B_", seq_len(20L))
  colnames(fx_b$expr) <- fx_b$coords$id

  expr_joint <- cbind(fx_a$expr, fx_b$expr)
  coords_joint <- rbind(fx_a$coords, fx_b$coords)
  samples <- c(rep("A", 20L), rep("B", 20L))

  res_joint <- spatial_banksy_augment(expr_joint, coords_joint, k_geom = 5L, use_agf = TRUE, sample_ids = samples)
  res_a_only <- spatial_banksy_augment(fx_a$expr, fx_a$coords, k_geom = 5L, use_agf = TRUE)

  # Sample A's unscaled nbr and agf blocks in joint run must match A-only run (zero cross-sample leakage)
  expect_equal(unname(res_joint$nbr[1:20, ]), unname(res_a_only$nbr), tolerance = 1e-12)
  expect_equal(unname(res_joint$agf[1:20, ]), unname(res_a_only$agf), tolerance = 1e-12)
  expect_identical(res_joint$n_samples, 2L)
})
