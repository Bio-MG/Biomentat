# =============================================================================
# test-spatial-roi-de.R — P2.4 (i): Reactive BPCells ROI DE
# =============================================================================

source_project_file("tools/pkg_skeleton/R/spatial_roi_de.R")

.mk_roi_fixture <- function(n_genes = 50L, n_roi = 40L, n_ctrl = 80L, seed = 42L) {
  set.seed(seed)
  n_cells <- n_roi + n_ctrl
  mat <- matrix(stats::rpois(n_genes * n_cells, lambda = 2), nrow = n_genes, ncol = n_cells)
  # Boost first 5 genes in ROI cells
  mat[1:5, seq_len(n_roi)] <- matrix(
    stats::rpois(5L * n_roi, lambda = 25), nrow = 5L, ncol = n_roi
  )
  rownames(mat) <- paste0("Gene", seq_len(n_genes))
  colnames(mat) <- c(paste0("ROI_", seq_len(n_roi)), paste0("CTRL_", seq_len(n_ctrl)))
  list(
    mat = Matrix::Matrix(mat, sparse = TRUE),
    roi_ids = paste0("ROI_", seq_len(n_roi)),
    ctrl_ids = paste0("CTRL_", seq_len(n_ctrl))
  )
}

test_that("spatial_roi_de_bpcells rejects invalid or under-populated ROI inputs", {
  fx <- .mk_roi_fixture()

  err1 <- tryCatch(spatial_roi_de_bpcells(fx$mat, character(0)), error = function(e) e)
  expect_s3_class(err1, "spatial_roi_de_error")

  err2 <- tryCatch(spatial_roi_de_bpcells(fx$mat, c("ROI_1", "ROI_2")), error = function(e) e)
  expect_s3_class(err2, "spatial_roi_de_error")

  err3 <- tryCatch(spatial_roi_de_bpcells(fx$mat, c("MISSING_1", "MISSING_2", "MISSING_3")), error = function(e) e)
  expect_s3_class(err3, "spatial_roi_de_error")
})

test_that("spatial_roi_de_bpcells detects upregulated ROI markers with BH adjustment", {
  fx <- .mk_roi_fixture()
  res <- spatial_roi_de_bpcells(fx$mat, roi_ids = fx$roi_ids, min_pct = 0.1, logfc_threshold = 0.25)

  expect_identical(res$backend, "Matrix")
  expect_identical(res$n_roi_total, 40L)
  expect_identical(res$n_ctrl_total, 80L)
  expect_true(is.data.frame(res$markers))
  expect_true(all(c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj") %in% names(res$markers)))

  sig <- res$markers[res$markers$p_val_adj < 0.05 & res$markers$avg_log2FC > 1, ]
  expect_true(all(paste0("Gene", 1:5) %in% sig$gene))
})

test_that("spatial_roi_de_bpcells enforces deterministic max_cells_per_group subsampling", {
  fx <- .mk_roi_fixture(n_roi = 60L, n_ctrl = 100L)
  r1 <- spatial_roi_de_bpcells(fx$mat, roi_ids = fx$roi_ids, max_cells_per_group = 25L, seed = 123L)
  r2 <- spatial_roi_de_bpcells(fx$mat, roi_ids = fx$roi_ids, max_cells_per_group = 25L, seed = 123L)

  expect_identical(r1$n_roi_used, 25L)
  expect_identical(r1$n_ctrl_used, 25L)
  expect_identical(r1$markers, r2$markers)
})

test_that("spatial_roi_de_bpcells supports BPCells on-disk IterableMatrix", {
  skip_if_not_installed("BPCells")
  fx <- .mk_roi_fixture()
  bp_dir <- tempfile("bp_roi_de_")
  on.exit(unlink(bp_dir, recursive = TRUE), add = TRUE)
  BPCells::write_matrix_dir(fx$mat, dir = bp_dir)
  bp_mat <- BPCells::open_matrix_dir(bp_dir)

  res <- spatial_roi_de_bpcells(bp_mat, roi_ids = fx$roi_ids, min_pct = 0.1, logfc_threshold = 0.25)
  expect_identical(res$backend, "BPCells")
  sig <- res$markers[res$markers$p_val_adj < 0.05 & res$markers$avg_log2FC > 1, ]
  expect_true(all(paste0("Gene", 1:5) %in% sig$gene))
})
