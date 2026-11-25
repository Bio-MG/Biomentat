# Tests for BK-32 1-click offline demo button on Bulk import

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
})

source_project_file("tools/pkg_skeleton/R/rdata_io.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("modules/import/mod_rdata_picker.R")
source_project_file("tools/pkg_skeleton/R/bulk_import_engine.R")
source_project_file("tools/pkg_skeleton/R/bulk_multi.R")
source_project_file("tools/pkg_skeleton/R/demo_loader.R")
source_project_file("modules/import/mod_import_bulk.R")

test_that("mod_import_bulk_ui renders btn_load_demo", {
  i18n_stub <- list(t = function(k) k)
  assign("i18n", i18n_stub, envir = .GlobalEnv)
  if (!exists(".tr_plain", envir = .GlobalEnv, inherits = FALSE)) {
    assign(".tr_plain", function(k) k, envir = .GlobalEnv)
    on.exit(rm(".tr_plain", envir = .GlobalEnv), add = TRUE)
  }
  on.exit(rm("i18n", envir = .GlobalEnv), add = TRUE)

  ui_html <- as.character(mod_import_bulk_ui("import_bulk"))
  expect_true(grepl("import_bulk-btn_load_demo", ui_html, fixed = TRUE))
})

test_that("clicking btn_load_demo loads 200x8 Bulk demo and updates business_state", {
  gd <- reactiveValues(
    language = "fr",
    bulk_obj = NULL,
    bulk_obj_epoch = 0L,
    business_state = NULL,
    session = NULL,
    i18n = list(t = function(k) k)
  )

  testServer(mod_import_bulk_server, args = list(global_data = gd), {
    session$setInputs(min_counts = 10, btn_load_demo = 1L)
    expect_false(is.null(gd$bulk_obj))
    expect_identical(nrow(gd$bulk_obj$counts), 200L)
    expect_identical(ncol(gd$bulk_obj$counts), 8L)
    expect_identical(gd$bulk_obj$project, "Bulk_Demo")
    expect_true(isTRUE(gd$business_state$has_data))
    expect_equal(gd$business_state$n_elements, 8)
  })
})
