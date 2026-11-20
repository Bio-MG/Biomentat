# =============================================================================
# test-mod-import-bulk-generic-omics.R — BK-03 generic omics wiring
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/rdata_io.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("tools/pkg_skeleton/R/provenance.R")
source_project_file("tools/pkg_skeleton/R/bulk_import_engine.R")
source_project_file("modules/import/mod_rdata_picker.R")
source_project_file("modules/import/mod_import_bulk.R")

test_that("mod_import_bulk_ui exposes generic omics data type and preprocessing controls (BK-03)", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  suppressPackageStartupMessages(library(shiny))
  suppressPackageStartupMessages(library(bslib))
  if (!exists("i18n", envir = globalenv(), inherits = FALSE)) {
    assign("i18n", list(t = function(x) x), envir = globalenv())
    on.exit(rm("i18n", envir = globalenv()), add = TRUE)
  }
  if (!exists(".tr_plain", envir = globalenv(), inherits = FALSE)) {
    assign(".tr_plain", function(x) x, envir = globalenv())
    on.exit(rm(".tr_plain", envir = globalenv()), add = TRUE)
  }
  ui <- mod_import_bulk_ui("t")
  html <- as.character(ui)
  expect_match(html, "t-bulk_data_type")
  expect_match(html, "t-generic_modality")
  expect_match(html, "t-generic_log_transform")
  expect_match(html, "t-generic_impute")
  expect_match(html, "t-generic_normalize")
})

test_that("mod_import_bulk_server dispatches generic_omics mode through prepare_generic_omics_matrix (BK-03)", {
  tmp_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp_csv), add = TRUE)
  mat_df <- data.frame(
    feature = c("P1", "P2", "P3"),
    S1 = c(10.5, 20.1, NA_real_),
    S2 = c(15.2, 25.4, 30.0),
    S3 = c(12.0, 22.8, 28.5),
    stringsAsFactors = FALSE
  )
  utils::write.csv(mat_df, tmp_csv, row.names = FALSE)

  global_data <- shiny::reactiveValues(
    bulk_obj = NULL,
    bulk_obj_epoch = 0L,
    bulk_datasets = list()
  )

  shiny::testServer(
    mod_import_bulk_server,
    args = list(global_data = global_data),
    {
      session$setInputs(
        bulk_import_mode = "merged_matrix",
        counts_file = list(name = "proteomics.csv", datapath = tmp_csv, size = file.info(tmp_csv)$size, type = "text/csv"),
        counts_format = "rows",
        counts_has_header = TRUE,
        counts_has_rownames = TRUE,
        bulk_data_type = "generic_omics",
        generic_modality = "proteomics",
        generic_log_transform = TRUE,
        generic_impute = "min_half",
        generic_normalize = "median_center",
        min_counts = 10,
        project_name = "Prot_Test",
        multi_label = "",
        btn_load = 1L
      )
      session$flushReact()

      obj <- global_data$bulk_obj
      expect_false(is.null(obj))
      expect_identical(obj$type, "bulk")
      expect_identical(obj$data_type, "generic_omics")
      expect_identical(obj$modality, "proteomics")
      expect_false(isTRUE(obj$is_count_matrix))
      expect_true(is.list(obj$generic_omics))
      expect_identical(obj$generic_omics$type, "generic_omics_matrix")
      expect_identical(sum(is.na(obj$counts)), 0L)
      expect_identical(nrow(obj$counts), 3L)
      expect_identical(ncol(obj$counts), 3L)
    }
  )
})
