# ==============================================================================
# test-mod-import-spatial-hd-ui.R — BK-10 UI/server wiring in mod_import_spatial
# ==============================================================================

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))

.hd_ui_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("mod_import_spatial_ui exposes MERSCOPE technology and controls", {
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)
  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }

  source(.hd_ui_root("modules/import/mod_import_spatial.R"), local = globalenv())
  ui_html <- as.character(mod_import_spatial_ui("sp_hd_test"))

  expect_true(grepl('value="merscope"', ui_html, fixed = TRUE))
  expect_true(grepl("sp_hd_test-merscope_drop_blanks", ui_html, fixed = TRUE))
  expect_true(grepl("sp_hd_test-min_counts_ms", ui_html, fixed = TRUE))
})

test_that("run_spatial_import imports MERSCOPE directory end-to-end", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("BPCells")

  source_project_file("tools/pkg_skeleton/R/error_state.R")
  source_project_file("tools/pkg_skeleton/R/state.R")
  source_project_file("tools/pkg_skeleton/R/io_helpers.R")
  source_project_file("tools/pkg_skeleton/R/spatial_io.R")
  source_project_file("tools/pkg_skeleton/R/spatial_hd_ingest.R")
  source_project_file("modules/import/mod_import_spatial.R")

  tmp_dir <- tempfile("merscope_ui_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  cnt_df <- data.frame(
    EntityID = sprintf("cell_%02d", 1:12),
    GeneA = rep(5, 12),
    GeneB = rep(3, 12),
    Blank_1 = rep(1, 12),
    stringsAsFactors = FALSE
  )
  meta_df <- data.frame(
    EntityID = sprintf("cell_%02d", 1:12),
    center_x = seq(10, 120, by = 10),
    center_y = seq(20, 130, by = 10),
    stringsAsFactors = FALSE
  )
  utils::write.csv(cnt_df, file.path(tmp_dir, "cell_by_gene.csv"), row.names = FALSE)
  utils::write.csv(meta_df, file.path(tmp_dir, "cell_metadata.csv"), row.names = FALSE)

  e <- new.env(parent = globalenv())
  e$.tr <- function(s, ...) s
  e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
  e$global_data <- new.env(parent = emptyenv())
  e$global_data$spatial_datasets <- list()
  e$global_data$active_spatial_dataset <- NULL
  e$global_data$spatial_results_cache <- list()
  e$add_log <- function(msg) invisible(NULL)
  e$showNotification <- function(msg, ...) invisible(NULL)
  e$withProgress <- function(expr, ...) expr
  e$incProgress <- function(...) invisible(NULL)
  e$input <- list(
    technology = "merscope", sample_name = "MERSCOPE_UI_Test",
    min_counts = 0, min_features = 0, min_counts_ss = 0, min_features_ss = 0,
    max_sketch = 5000, norm_method = "lognorm", simplify_tol = 20,
    load_raw_also = FALSE, orient_swap_xy = FALSE,
    orient_flip_x = FALSE, orient_flip_y = FALSE,
    hd_bin_size = NULL, hd_pool_factor = "1",
    merscope_drop_blanks = TRUE, min_counts_ms = 1
  )

  writer <- ts_ast_assignment("modules/import/mod_import_spatial.R", "run_spatial_import", e)
  res <- writer(tmp_dir, "MERSCOPE_UI_Test", "merscope")
  expect_true(isTRUE(res$ok))
  expect_equal(res$technology, "merscope")
  expect_equal(res$n_total, 12L)
  ds <- e$global_data$spatial_datasets[["MERSCOPE_UI_Test"]]
  expect_false(is.null(ds))
  expect_equal(ds$technology, "merscope")
  expect_equal(nrow(ds$sketch), 2L)
})
