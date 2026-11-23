# ==============================================================================
# test-mod-report-attrition-ui.R — BK-11 cohort attrition UI/server wiring
# ==============================================================================

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))

.attr_test_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("mod_sc_report_consolidated UI exposes attrition outputs without navset", {
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)

  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }

  source(.attr_test_root("modules/sc/mod_sc_report_consolidated.R"), local = globalenv())
  ui_left <- as.character(mod_sc_report_consolidated_ui("sc_rep_test"))
  ui_main <- as.character(mod_sc_report_consolidated_output_ui("sc_rep_test"))

  expect_true(grepl("sc_rep_test-dl_attrition_csv", ui_left, fixed = TRUE))
  expect_true(grepl("sc_rep_test-rep_attrition_plot", ui_main, fixed = TRUE))
  expect_true(grepl("sc_rep_test-rep_attrition_mermaid", ui_main, fixed = TRUE))
  expect_false(grepl("navset_tab(", ui_main, fixed = TRUE))
})

test_that("mod_sc_report_consolidated_server renders Mermaid and exports attrition CSV", {
  source(.attr_test_root("tools/pkg_skeleton/R/state.R"), local = globalenv())
  source(.attr_test_root("tools/pkg_skeleton/R/cohort_attrition.R"), local = globalenv())
  source(.attr_test_root("tools/pkg_skeleton/R/report_collector.R"), local = globalenv())
  source(.attr_test_root("tools/pkg_skeleton/R/report_validator.R"), local = globalenv())
  source(.attr_test_root("tools/pkg_skeleton/R/report_render.R"), local = globalenv())
  source(.attr_test_root("tools/pkg_skeleton/R/report_bundle.R"), local = globalenv())
  source(.attr_test_root("modules/sc/mod_sc_report_consolidated.R"), local = globalenv())

  mock_sc <- matrix(1, nrow = 5, ncol = 120)
  shared_rv <- shiny::reactiveValues(
    n_cells_raw = 200L,
    n_cells_post_qc = 150L,
    n_cells_post_doublets = 130L
  )
  global_data <- shiny::reactiveValues(
    language = "fr",
    i18n = list(t = function(k) k),
    sc_obj = mock_sc
  )

  shiny::testServer(
    mod_sc_report_consolidated_server,
    args = list(global_data = global_data, shared_rv = shared_rv),
    {
      mermaid_txt <- output$rep_attrition_mermaid
      expect_true(is.character(mermaid_txt) && grepl("flowchart TD", mermaid_txt, fixed = TRUE))
      expect_true(grepl("n = 200", mermaid_txt, fixed = TRUE))
      expect_true(grepl("n = 120", mermaid_txt, fixed = TRUE))

      csv_path <- output$dl_attrition_csv
      expect_true(is.character(csv_path) && file.exists(csv_path))
      df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
      expect_true(all(c("step_index", "step", "unit", "n_in", "n_excluded",
                        "n_remaining", "pct_cumulative_retained") %in% colnames(df)))
      expect_equal(df$n_remaining, c(200L, 150L, 130L, 120L))
    }
  )
})

test_that("mod_bulk_report UI and server expose dl_bulk_attrition_csv", {
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)
  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }
  if (!exists("i18n_plain", envir = globalenv())) {
    assign("i18n_plain", function(k) k, envir = globalenv())
  }

  source(.attr_test_root("tools/pkg_skeleton/R/cohort_attrition.R"), local = globalenv())
  source(.attr_test_root("modules/bulk/mod_bulk_report.R"), local = globalenv())

  ui_html <- as.character(mod_bulk_report_ui("bulk_rep_test"))
  expect_true(grepl("bulk_rep_test-dl_bulk_attrition_csv", ui_html, fixed = TRUE))

  shared_rv <- shiny::reactiveValues(
    vst_mat = matrix(1, nrow = 60, ncol = 4),
    filtered_counts = matrix(1, nrow = 80, ncol = 4),
    contrasts = list(),
    active_contrast = NULL,
    lfc_thresh = 1,
    padj_thresh = 0.05,
    last_report = NULL
  )
  global_data <- shiny::reactiveValues(
    language = "fr",
    i18n = list(t = function(k) k),
    bulk_obj = list(counts = matrix(1, nrow = 100, ncol = 4), metadata = data.frame(x = 1:4))
  )

  shiny::testServer(
    mod_bulk_report_server,
    args = list(global_data = global_data, shared_rv = shared_rv),
    {
      csv_path <- output$dl_bulk_attrition_csv
      expect_true(is.character(csv_path) && file.exists(csv_path))
      df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
      expect_equal(df$n_remaining, c(100L, 80L, 60L))
      expect_equal(df$unit[1L], "genes")
      expect_equal(shared_rv$last_report$kind, "Attrition cohorte")
    }
  )
})
