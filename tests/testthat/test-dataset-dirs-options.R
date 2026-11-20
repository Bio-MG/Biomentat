# =============================================================================
# test-dataset-dirs-options.R — Saved dataset discovery folders + SC/Spatial/Bulk spread
# =============================================================================

source_project_file("tools/pkg_skeleton/R/session_checkpoint.R")
source_project_file("tools/pkg_skeleton/R/dataset_directory.R")
source_project_file("modules/import/mod_dataset_dirs_options.R")

test_that("ts_save_dataset_dir, ts_saved_dataset_dirs, ts_remove_dataset_dir persist per modality", {
  tmp_cfg <- tempfile(fileext = ".json")
  on.exit(unlink(tmp_cfg), add = TRUE)

  expect_identical(ts_saved_dataset_dirs("sc", config_path = tmp_cfg), character(0))

  d1 <- file.path(tempdir(), "ds_sc_extra1")
  d2 <- file.path(tempdir(), "ds_bulk_extra1")
  ts_save_dataset_dir(d1, modality = "sc", config_path = tmp_cfg)
  ts_save_dataset_dir(d2, modality = "bulk", config_path = tmp_cfg)

  sc_dirs <- ts_saved_dataset_dirs("sc", config_path = tmp_cfg)
  bulk_dirs <- ts_saved_dataset_dirs("bulk", config_path = tmp_cfg)
  expect_length(sc_dirs, 1L)
  expect_length(bulk_dirs, 1L)
  expect_identical(sc_dirs, normalizePath(d1, winslash = "/", mustWork = FALSE))
  expect_identical(bulk_dirs, normalizePath(d2, winslash = "/", mustWork = FALSE))

  ts_remove_dataset_dir(d1, modality = "sc", config_path = tmp_cfg)
  expect_identical(ts_saved_dataset_dirs("sc", config_path = tmp_cfg), character(0))
  expect_length(ts_saved_dataset_dirs("bulk", config_path = tmp_cfg), 1L)
})

test_that("ts_sc_dataset_dir_scan and ts_spatial_dataset_dir_scan scan multiple saved folders", {
  root1 <- tempfile("sc_root1_")
  root2 <- tempfile("sc_root2_")
  on.exit(unlink(c(root1, root2), recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(file.path(root1, "sample_a"), recursive = TRUE)
  writeLines("x", file.path(root1, "sample_a", "matrix.mtx"))
  writeLines("x", file.path(root1, "sample_a", "barcodes.tsv"))
  writeLines("x", file.path(root1, "sample_a", "features.tsv"))
  dir.create(root2, recursive = TRUE)
  writeLines("rds_stub", file.path(root2, "pbmc_saved.rds"))

  sc_res <- ts_sc_dataset_dir_scan(c(root1, root2))
  expect_true("sample_a" %in% names(sc_res))
  expect_true("pbmc_saved.rds" %in% names(sc_res))

  sp1 <- tempfile("sp_root1_")
  sp2 <- tempfile("sp_root2_")
  on.exit(unlink(c(sp1, sp2), recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(file.path(sp1, "visium_1", "spatial"), recursive = TRUE)
  dir.create(sp2, recursive = TRUE)
  writeLines("rds_stub", file.path(sp2, "brain_saved.rds"))

  sp_res <- ts_spatial_dataset_dir_scan(c(sp1, sp2))
  expect_true("visium_1" %in% names(sp_res))
  expect_true("brain_saved.rds" %in% names(sp_res))
})

test_that("dataset_dirs_options_ui exposes folder save controls", {
  skip_if_not_installed("shiny")
  suppressPackageStartupMessages(library(shiny))
  if (!exists(".tr_plain", envir = globalenv(), inherits = FALSE)) {
    assign(".tr_plain", function(x) x, envir = globalenv())
    on.exit(rm(".tr_plain", envir = globalenv()), add = TRUE)
  }
  ui <- dataset_dirs_options_ui(NS("t"), list(t = function(x) x))
  html <- as.character(ui)
  expect_match(html, "t-custom_dataset_dir")
  expect_match(html, "t-btn_save_dataset_dir")
  expect_match(html, "t-saved_dataset_dirs_ui")
})
