# Tests for BK-32 1-click offline demo buttons on SC and Spatial import

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(shinyFiles)
  library(Seurat)
})

source_project_file("tools/pkg_skeleton/R/rdata_io.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("modules/import/mod_rdata_picker.R")
source_project_file("tools/pkg_skeleton/R/sc_multi.R")
source_project_file("tools/pkg_skeleton/R/demo_loader.R")
source_project_file("modules/import/mod_import_sc.R")
source_project_file("modules/import/mod_import_spatial.R")

test_that("mod_import_sc_ui and mod_import_spatial_ui render btn_load_demo", {
  i18n_stub <- list(t = function(k) k)
  assign("i18n", i18n_stub, envir = .GlobalEnv)
  if (!exists(".tr_plain", envir = .GlobalEnv, inherits = FALSE)) {
    assign(".tr_plain", function(k) k, envir = .GlobalEnv)
    on.exit(rm(".tr_plain", envir = .GlobalEnv), add = TRUE)
  }
  on.exit(rm("i18n", envir = .GlobalEnv), add = TRUE)

  sc_html <- as.character(mod_import_sc_ui("import_sc"))
  expect_true(grepl("import_sc-btn_load_demo", sc_html, fixed = TRUE))

  sp_html <- as.character(mod_import_spatial_ui("import_spatial"))
  expect_true(grepl("import_spatial-btn_load_demo", sp_html, fixed = TRUE))
})

test_that("clicking btn_load_demo in mod_import_sc_server loads 1500x400 real SC demo", {
  gd <- reactiveValues(
    language = "fr",
    sc_obj = NULL,
    sc_obj_epoch = 0L,
    sc_datasets = list(),
    business_state = NULL,
    session = NULL,
    i18n = list(t = function(k) k)
  )

  testServer(mod_import_sc_server, args = list(global_data = gd), {
    session$setInputs(btn_load_demo = 1L)
    expect_false(is.null(gd$sc_obj))
    expect_equal(nrow(gd$sc_obj), 1500)
    expect_equal(ncol(gd$sc_obj), 400)
    expect_true(isTRUE(gd$business_state$has_data))
    expect_equal(gd$business_state$n_elements, 400)

    # Downstream Seurat workflow smoke check
    so <- Seurat::NormalizeData(gd$sc_obj, verbose = FALSE)
    so <- Seurat::FindVariableFeatures(so, nfeatures = 300, verbose = FALSE)
    so <- Seurat::ScaleData(so, verbose = FALSE)
    so <- Seurat::RunPCA(so, npcs = 10, verbose = FALSE)
    expect_true("pca" %in% names(so@reductions))
  })
})

test_that("clicking btn_load_demo in mod_import_spatial_server loads 1500x400 real Spatial demo", {
  gd <- reactiveValues(
    language = "fr",
    spatial_datasets = list(),
    active_spatial_dataset = NULL,
    spatial_reimport_signal = NULL,
    spatial_results_cache = list(),
    business_state = NULL,
    session = NULL,
    i18n = list(t = function(k) k)
  )

  testServer(mod_import_spatial_server, args = list(global_data = gd), {
    session$setInputs(btn_load_demo = 1L)
    expect_identical(gd$active_spatial_dataset, "Spatial_Demo")
    pkg <- gd$spatial_datasets[["Spatial_Demo"]]
    expect_false(is.null(pkg))
    expect_equal(nrow(pkg$sketch), 1500)
    expect_equal(ncol(pkg$sketch), 400)
    expect_equal(nrow(pkg$coords), 400)
    expect_true(isTRUE(gd$business_state$has_data))
    expect_equal(gd$business_state$n_elements, 400)
    data_mat <- SeuratObject::LayerData(pkg$sketch, assay = SeuratObject::DefaultAssay(pkg$sketch), layer = "data")
    expect_identical(dim(data_mat), c(1500L, 400L))
  })
})
