# ==============================================================================
# test-mod-bulk-report-qmd.R — BK-12 Quarto audit notebook export UI/server
# ==============================================================================

skip_if_not_installed("shiny")
suppressPackageStartupMessages(library(shiny))

.qmd_test_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("mod_bulk_report_ui exposes dl_audit_qmd downloadButton", {
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)

  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }
  if (!exists("i18n_plain", envir = globalenv())) {
    assign("i18n_plain", function(k) k, envir = globalenv())
  }

  source(.qmd_test_root("modules/bulk/mod_bulk_report.R"), local = globalenv())
  ui_html <- as.character(mod_bulk_report_ui("rep_test"))
  expect_true(grepl("rep_test-dl_audit_qmd", ui_html, fixed = TRUE))
})

test_that("mod_bulk_report_server registers dl_audit_qmd and exports valid .qmd", {
  if (!exists("TS_DRIVE_DIR", envir = globalenv())) {
    assign("TS_DRIVE_DIR", "tools/_drive", envir = globalenv())
  }
  source(.qmd_test_root("tools/pkg_skeleton/R/drive_watcher.R"), local = globalenv())
  source(.qmd_test_root("modules/bulk/mod_bulk_report.R"), local = globalenv())

  shared_rv <- shiny::reactiveValues(
    vst_mat = NULL,
    filtered_counts = NULL,
    contrasts = list(),
    active_contrast = NULL,
    lfc_thresh = 1,
    padj_thresh = 0.05,
    last_report = NULL
  )
  global_data <- shiny::reactiveValues(
    language = "fr",
    i18n = list(t = function(k) k),
    bulk_obj = list(counts = NULL, metadata = NULL)
  )

  shiny::testServer(
    mod_bulk_report_server,
    args = list(global_data = global_data, shared_rv = shared_rv),
    {
      dl_path <- output$dl_audit_qmd
      expect_true(is.character(dl_path) && file.exists(dl_path))
      lines <- readLines(dl_path, warn = FALSE)
      expect_true(any(grepl("^format: html", lines)))
      expect_true(any(grepl("## Audit Events", lines, fixed = TRUE)))
      expect_equal(shared_rv$last_report$kind, "Audit Quarto")
    }
  )
})
