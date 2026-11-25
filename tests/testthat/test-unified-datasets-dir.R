# Tests for W4.2 unified datasets/data_{bulk,sc,spatial} discovery

source_project_file("tools/pkg_skeleton/R/session_checkpoint.R")
source_project_file("tools/pkg_skeleton/R/dataset_directory.R")

test_that(".ts_unified_datasets_candidates discovers existing unified modality folders", {
  root <- tempfile("unified_ds_")
  dir.create(file.path(root, "datasets", "data_bulk"), recursive = TRUE)
  dir.create(file.path(root, "datasets", "data_sc", "datasets"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  bulk_cands <- .ts_unified_datasets_candidates("bulk", project_root = root)
  expect_true(any(grepl("/datasets/data_bulk$", bulk_cands)))

  sc_cands <- .ts_unified_datasets_candidates("sc", project_root = root)
  expect_true(any(grepl("/datasets/data_sc$", sc_cands)))
  expect_true(any(grepl("/datasets/data_sc/datasets$", sc_cands)))

  sp_cands <- .ts_unified_datasets_candidates("spatial", project_root = root)
  expect_identical(sp_cands, character(0))
})

test_that(".ts_unified_datasets_candidates respects external TRANSCRIPTO_DATASETS_DIR", {
  ext <- tempfile("ext_ds_")
  dir.create(file.path(ext, "data_spatial"), recursive = TRUE)
  on.exit(unlink(ext, recursive = TRUE), add = TRUE)

  old_env <- Sys.getenv("TRANSCRIPTO_DATASETS_DIR", unset = NA_character_)
  Sys.setenv(TRANSCRIPTO_DATASETS_DIR = ext)
  on.exit({
    if (is.na(old_env)) Sys.unsetenv("TRANSCRIPTO_DATASETS_DIR")
    else Sys.setenv(TRANSCRIPTO_DATASETS_DIR = old_env)
  }, add = TRUE)

  sp_cands <- .ts_unified_datasets_candidates("spatial", project_root = "")
  expect_true(any(grepl("/data_spatial$", sp_cands)))
  dirs <- ts_datasets_dirs(modality = "spatial")
  expect_true(any(grepl("/data_spatial$", dirs)))
})
