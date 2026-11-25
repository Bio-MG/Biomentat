# Tests for W4.2 harmonized visual status emoji badges in Datasets decouverts

source_project_file("modules/import/mod_import_bulk.R")
source_project_file("modules/import/mod_import_sc.R")

test_that("Bulk and SC status badge helpers render 4-state emoji badges", {
  for (fn in list(.bulk_format_ds_status_badge, .sc_format_ds_status_badge)) {
    expect_identical(fn(NULL), "")
    expect_identical(fn(NA_character_), "")
    expect_identical(fn(""), "")

    b_valid <- as.character(fn("valid"))
    expect_true(grepl("bg-success", b_valid, fixed = TRUE))
    expect_true(grepl("\U0001F7E2 valid", b_valid, fixed = TRUE))

    b_warn <- as.character(fn("warning"))
    expect_true(grepl("bg-warning", b_warn, fixed = TRUE))
    expect_true(grepl("\U0001F7E1 warning", b_warn, fixed = TRUE))

    b_partial <- as.character(fn("partial"))
    expect_true(grepl("bg-warning", b_partial, fixed = TRUE))
    expect_true(grepl("\U0001F7E1 partial", b_partial, fixed = TRUE))

    b_empty <- as.character(fn("empty"))
    expect_true(grepl("bg-secondary", b_empty, fixed = TRUE))
    expect_true(grepl("\u26AA empty", b_empty, fixed = TRUE))

    b_inv <- as.character(fn("invalid"))
    expect_true(grepl("bg-danger", b_inv, fixed = TRUE))
    expect_true(grepl("\U0001F534 invalid", b_inv, fixed = TRUE))
  }
})

test_that(".bulk_datasets_table_rows and .sc_datasets_table_rows embed emoji badges", {
  scan_bulk <- list(
    ds1 = list(status = "valid", raison = "", manifest = list(species = "human", matrix_type = "counts"),
               checks = list(experimental_design = list(n_samples = 6L))),
    ds2 = list(status = "invalid", raison = "missing", manifest = list(species = "human", matrix_type = "counts"),
               checks = list())
  )
  rows_b <- .bulk_datasets_table_rows(scan_bulk, ns = function(x) x, tr = function(x) x)
  expect_true(grepl("\U0001F7E2 valid", as.character(rows_b$ds1$row), fixed = TRUE))
  expect_true(grepl("\U0001F534 invalid", as.character(rows_b$ds2$row), fixed = TRUE))

  scan_sc <- list(
    sc1 = list(status = "valid", format = "10x", n_samples = 2L, raison = ""),
    sc2 = list(status = "empty", format = "unknown", n_samples = 0L, raison = "empty")
  )
  rows_s <- .sc_datasets_table_rows(scan_sc, ns = function(x) x, tr = function(x) x)
  expect_true(grepl("\U0001F7E2 valid", as.character(rows_s$sc1$row), fixed = TRUE))
  expect_true(grepl("\u26AA empty", as.character(rows_s$sc2$row), fixed = TRUE))
})

test_that("mod_import_spatial.R renders harmonized 4-state emoji badges", {
  src <- paste(readLines(testthat::test_path("../..", "modules/import/mod_import_spatial.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("bg-success", src, fixed = TRUE))
  expect_true(grepl("bg-warning", src, fixed = TRUE))
  expect_true(grepl("bg-danger", src, fixed = TRUE))
  expect_true(grepl("\\U0001F7E2 valid", src, fixed = TRUE))
  expect_true(grepl("\\U0001F534", src, fixed = TRUE))
})
