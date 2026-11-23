# ==============================================================================
# test-spatial-p24-wiring.R — W2.3 P2.4 Spatial UI Wiring & Async Cancel Helper
# ==============================================================================

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
suppressPackageStartupMessages(library(DT))

.sp24_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("APP_DAEMON_SOURCE_FILES preloads P2.4 spatial engines and cancel helper is safe", {
  source(.sp24_root("tools/pkg_skeleton/R/spatial_async.R"), local = globalenv())

  expect_true("tools/pkg_skeleton/R/spatial_banksy.R" %in% APP_DAEMON_SOURCE_FILES)
  expect_true("tools/pkg_skeleton/R/spatial_roi_de.R" %in% APP_DAEMON_SOURCE_FILES)
  expect_true("tools/pkg_skeleton/R/spatial_milo.R" %in% APP_DAEMON_SOURCE_FILES)

  expect_false(ts_async_cancel_mirai(NULL))
  expect_false(ts_async_cancel_mirai(list(a = 1)))
})

test_that("mod_spatial_cluster_controls_ui exposes use_agf checkbox", {
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)
  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }

  source(.sp24_root("modules/spatial/mod_spatial_cluster.R"), local = globalenv())
  ui_html <- as.character(mod_spatial_cluster_controls_ui("sp_cl"))
  expect_true(grepl("sp_cl-use_agf", ui_html, fixed = TRUE))
})

test_that("mod_spatial_multi_ui exposes spatial_milo tab and controls", {
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)
  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }

  source(.sp24_root("modules/spatial/mod_spatial_multi.R"), local = globalenv())
  ui_html <- as.character(mod_spatial_multi_ui("sp_m"))
  expect_true(grepl("spatial_milo", ui_html, fixed = TRUE))
  expect_true(grepl("sp_m-milo_ref_datasets", ui_html, fixed = TRUE))
  expect_true(grepl("sp_m-milo_target_datasets", ui_html, fixed = TRUE))
  expect_true(grepl("sp_m-btn_spatial_milo", ui_html, fixed = TRUE))
  expect_true(grepl("sp_m-spatial_milo_table", ui_html, fixed = TRUE))
})

test_that("mod_spatial_multi_server runs spatial_milo_da via testServer", {
  skip_if_not_installed("RANN")
  skip_if_not_installed("edgeR")

  source(.sp24_root("tools/pkg_skeleton/R/defaults.R"), local = globalenv())
  source(.sp24_root("tools/pkg_skeleton/R/spatial_async.R"), local = globalenv())
  source(.sp24_root("tools/pkg_skeleton/R/spatial_milo.R"), local = globalenv())
  source(.sp24_root("modules/spatial/mod_spatial_multi.R"), local = globalenv())

  if (!exists("ts_datatable", envir = globalenv())) {
    assign("ts_datatable", function(df, ...) DT::datatable(df), envir = globalenv())
  }
  if (!exists(".strip_i18n_html", envir = globalenv())) {
    assign(".strip_i18n_html", function(x) x, envir = globalenv())
  }
  if (!exists("spatial_log_path", envir = globalenv())) {
    assign("spatial_log_path", function(session, tag) tempfile(fileext = ".log"), envir = globalenv())
  }
  if (!exists("spatial_discrete_colors", envir = globalenv())) {
    assign("spatial_discrete_colors", function(v, shared_rv = NULL, kind = "cluster") {
      lv <- unique(as.character(v))
      stats::setNames(grDevices::hcl.colors(max(length(lv), 1L)), lv)
    }, envir = globalenv())
  }

  set.seed(2026)
  n_per <- 25L
  ds_names <- c("R1", "R2", "T1", "T2")
  emb <- data.frame(
    dim1 = c(rnorm(2L * n_per, -0.3, 1.0), rnorm(2L * n_per, 0.3, 1.0)),
    dim2 = rnorm(4L * n_per, 0, 1.0),
    dataset = rep(ds_names, each = n_per),
    cluster = rep(c("C0", "C1"), length.out = 4L * n_per),
    original_id = paste0("spot_", seq_len(4L * n_per)),
    stringsAsFactors = FALSE
  )

  global_data <- shiny::reactiveValues(
    language = "fr",
    i18n = list(t = function(k) k),
    spatial_datasets = list(),
    spatial_multi_integration = list(
      embeddings = emb,
      n_per_dataset = table(emb$dataset),
      reduction_used = "pca"
    )
  )
  shared_rv <- shiny::reactiveValues()

  shiny::testServer(
    mod_spatial_multi_server,
    args = list(global_data = global_data, shared_rv = shared_rv),
    {
      session$setInputs(
        milo_ref_datasets = c("R1", "R2"),
        milo_target_datasets = c("T1", "T2")
      )
      session$setInputs(btn_spatial_milo = 1L)
      st <- output$spatial_milo_status
      expect_true(grepl("Milo spatial :", st, fixed = TRUE))
      expect_true(grepl("voisinages testés", st, fixed = TRUE))
      expect_false(is.null(output$spatial_milo_table))
    }
  )
})
