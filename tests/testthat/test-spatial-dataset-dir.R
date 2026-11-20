# =============================================================================
# test-spatial-dataset-dir.R — Découverte de datasets Spatial + MCP modality
# =============================================================================
# Complète la découverte « Datasets découverts » (déjà présente en Bulk et SC)
# pour le domaine Spatial (`data_spatial/datasets`, `ts_spatial_dataset_dir_scan`,
# `mod_import_spatial.R`) et expose le filtre `modality` dans `transcripto_drive_datasets`.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/session_checkpoint.R")
source_project_file("tools/pkg_skeleton/R/dataset_directory.R")

test_that("ts_spatial_dataset_dir_scan : dossier absent ou vide retourne list()", {
  expect_identical(ts_spatial_dataset_dir_scan(file.path(tempdir(), "no_such_spatial_dir_xyz")), list())
  empty_dir <- tempfile("sp_empty_")
  dir.create(empty_dir, recursive = TRUE)
  on.exit(unlink(empty_dir, recursive = TRUE), add = TRUE)
  expect_identical(ts_spatial_dataset_dir_scan(empty_dir), list())
})

test_that("ts_spatial_dataset_dir_scan : détecte Visium, fichiers .rds/.h5ad et rejette vides", {
  root <- tempfile("sp_scan_")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  vis <- file.path(root, "visium_brain")
  dir.create(file.path(vis, "spatial"), recursive = TRUE)
  writeLines("barcode,in_tissue", file.path(vis, "spatial", "tissue_positions_list.csv"))
  writeBin(as.raw(c(1, 2, 3, 4)), file.path(vis, "filtered_feature_bc_matrix.h5"))

  bad_dir <- file.path(root, "empty_folder")
  dir.create(bad_dir)

  h5ad_file <- file.path(root, "sample.h5ad")
  writeBin(as.raw(c(1, 2, 3)), h5ad_file)

  zero_file <- file.path(root, "empty.rds")
  file.create(zero_file)

  res <- ts_spatial_dataset_dir_scan(root)
  expect_true("visium_brain" %in% names(res))
  expect_identical(res[["visium_brain"]]$status, "valid")
  expect_identical(res[["visium_brain"]]$format, "visium")

  expect_true("empty_folder" %in% names(res))
  expect_identical(res[["empty_folder"]]$status, "invalid")

  expect_true("sample.h5ad" %in% names(res))
  expect_identical(res[["sample.h5ad"]]$status, "valid")
  expect_identical(res[["sample.h5ad"]]$format, "h5ad")

  expect_true("empty.rds" %in% names(res))
  expect_identical(res[["empty.rds"]]$status, "invalid")
})

test_that("mod_import_spatial.R expose Datasets découverts et ts_spatial_dataset_dir_scan", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "modules/import/mod_import_spatial.R"),
                         warn = FALSE), collapse = "\n")
  expect_match(src, "Datasets découverts")
  expect_match(src, 'actionButton\\(ns\\("btn_scan_datasets"\\)')
  expect_match(src, "ts_spatial_dataset_dir_scan\\(")
})

test_that("mcp_server.R : transcripto_drive_datasets accepte modality = sc/bulk/spatial/all", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "scripts/mcp_server.R"),
                         warn = FALSE), collapse = "\n")
  expect_match(src, "\\.ts_tool_datasets <- function\\(modality = NULL\\)")
})
