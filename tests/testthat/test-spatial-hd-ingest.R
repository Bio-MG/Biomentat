# ==============================================================================
# test-spatial-hd-ingest.R — BK-10 MERSCOPE/MERFISH & Visium HD Ingestion
# ==============================================================================

skip_if_not_installed("Matrix")

.hd_test_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("detect_spatial_hd_format identifies MERSCOPE/MERFISH, Visium HD, and unknown", {
  source(.hd_test_root("tools/pkg_skeleton/R/spatial_hd_ingest.R"), local = globalenv())

  d_mer <- tempfile("mer_")
  dir.create(d_mer, recursive = TRUE)
  on.exit(unlink(d_mer, recursive = TRUE), add = TRUE)
  file.create(file.path(d_mer, "cell_by_gene.csv"))
  file.create(file.path(d_mer, "cell_metadata.csv"))
  det_mer <- detect_spatial_hd_format(d_mer)
  expect_equal(det_mer$platform, "merscope_merfish")

  d_vhd <- tempfile("vhd_")
  dir.create(file.path(d_vhd, "binned_outputs", "square_002um"), recursive = TRUE)
  dir.create(file.path(d_vhd, "binned_outputs", "square_008um"), recursive = TRUE)
  on.exit(unlink(d_vhd, recursive = TRUE), add = TRUE)
  det_vhd <- detect_spatial_hd_format(d_vhd)
  expect_equal(det_vhd$platform, "visium_hd")
  expect_equal(det_vhd$bin_sizes, c("square_002um", "square_008um"))

  d_unk <- tempfile("unk_")
  dir.create(d_unk, recursive = TRUE)
  on.exit(unlink(d_unk, recursive = TRUE), add = TRUE)
  expect_equal(detect_spatial_hd_format(d_unk)$platform, "unknown")
})

test_that("read_merscope_merfish ingests counts, drops blanks, filters cells, and validates errors", {
  source(.hd_test_root("tools/pkg_skeleton/R/spatial_hd_ingest.R"), local = globalenv())

  cnt_f <- tempfile(fileext = ".csv")
  meta_f <- tempfile(fileext = ".csv")
  bad_meta_f <- tempfile(fileext = ".csv")
  on.exit(unlink(c(cnt_f, meta_f, bad_meta_f)), add = TRUE)

  write.csv(data.frame(
    EntityID = c("C1", "C2", "C3"),
    EPCAM = c(5, 0, 12),
    CD3D = c(2, 0, 8),
    `Blank-1` = c(1, 0, 0),
    blank_2 = c(0, 0, 1),
    check.names = FALSE
  ), cnt_f, row.names = FALSE)

  write.csv(data.frame(
    EntityID = c("C1", "C2", "C3"),
    center_x = c(10.5, 20.0, 30.5),
    center_y = c(100.0, 110.0, 120.0),
    fov = c(1L, 1L, 2L)
  ), meta_f, row.names = FALSE)

  res <- read_merscope_merfish(cnt_f, meta_f, drop_blanks = TRUE, min_counts = 1L)
  expect_s4_class(res$counts, "dgCMatrix")
  expect_equal(res$n_cells, 2L)
  expect_equal(res$n_genes, 2L)
  expect_equal(res$n_blanks_removed, 2L)
  expect_equal(colnames(res$counts), c("C1", "C3"))
  expect_equal(res$coords$x, c(10.5, 30.5))

  expect_error(read_merscope_merfish("nonexistent_file.csv", meta_f), class = "spatial_hd_ingest_error")
  write.csv(data.frame(EntityID = c("C1"), fov = 1L), bad_meta_f, row.names = FALSE)
  expect_error(read_merscope_merfish(cnt_f, bad_meta_f), class = "spatial_hd_ingest_error")
})

test_that("aggregate_visium_hd_bins pools 4x4 grid into 4 super-bins preserving total counts", {
  source(.hd_test_root("tools/pkg_skeleton/R/spatial_hd_ingest.R"), local = globalenv())

  grid <- expand.grid(x = c(0, 2, 4, 6), y = c(0, 2, 4, 6))
  bin_ids <- paste0("B", seq_len(nrow(grid)))
  coords <- data.frame(id = bin_ids, x = grid$x, y = grid$y, stringsAsFactors = FALSE)

  set.seed(11)
  cnt <- Matrix::Matrix(
    matrix(rpois(3L * 16L, lambda = 4), nrow = 3L, ncol = 16L,
           dimnames = list(paste0("G", 1:3), bin_ids)),
    sparse = TRUE
  )

  out <- aggregate_visium_hd_bins(cnt, coords, bin_factor = 2L)
  expect_s4_class(out$counts, "dgCMatrix")
  expect_equal(out$n_input_bins, 16L)
  expect_equal(out$n_output_bins, 4L)
  expect_equal(sum(out$counts), sum(cnt))
  expect_equal(sum(out$coords$n_subbins), 16L)
})

test_that("load_spatial_merscope and apply_visium_hd_bin_pooling return valid Seurat objects", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("Matrix")
  source(.hd_test_root("tools/pkg_skeleton/R/spatial_io.R"), local = globalenv())
  source(.hd_test_root("tools/pkg_skeleton/R/spatial_hd_ingest.R"), local = globalenv())

  tmp_dir <- tempfile("merscope_seurat_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  cnt_df <- data.frame(
    EntityID = sprintf("c%d", 1:6),
    Gene1 = c(5, 2, 4, 1, 3, 6),
    Gene2 = c(1, 3, 2, 4, 2, 1),
    Blank_ctrl = c(1, 0, 1, 0, 0, 1),
    stringsAsFactors = FALSE
  )
  meta_df <- data.frame(
    EntityID = sprintf("c%d", 1:6),
    center_x = c(0, 2, 0, 2, 4, 4),
    center_y = c(0, 0, 2, 2, 0, 2),
    stringsAsFactors = FALSE
  )
  utils::write.csv(cnt_df, file.path(tmp_dir, "cell_by_gene.csv"), row.names = FALSE)
  utils::write.csv(meta_df, file.path(tmp_dir, "cell_metadata.csv"), row.names = FALSE)

  obj <- load_spatial_merscope(tmp_dir, sample_name = "MS_Test")
  expect_s4_class(obj, "Seurat")
  expect_equal(ncol(obj), 6L)
  expect_equal(nrow(obj), 2L)
  expect_equal(attr(obj, "ts_spatial_mode"), "merscope")
  expect_equal(attr(obj, "ts_merscope_blanks_removed"), 1L)
  expect_true(all(c("coord_x", "coord_y") %in% colnames(obj@meta.data)))

  same_obj <- apply_visium_hd_bin_pooling(obj, bin_factor = 1L)
  expect_identical(ncol(same_obj), ncol(obj))

  pooled_obj <- apply_visium_hd_bin_pooling(obj, bin_factor = 2L)
  expect_s4_class(pooled_obj, "Seurat")
  expect_true(ncol(pooled_obj) < ncol(obj))
  expect_equal(attr(pooled_obj, "ts_hd_pool_factor"), 2L)
  expect_true("n_subbins" %in% colnames(pooled_obj@meta.data))
})
