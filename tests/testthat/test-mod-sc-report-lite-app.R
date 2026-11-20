# =============================================================================
# test-mod-sc-report-lite-app.R — W2.3 BK-01 UI wiring in mod_sc_report_consolidated.R
# =============================================================================
# `generate_biomentat_lite_app()` (report_bundle.R:259) génère un explorateur
# Shiny autonome (app.R + lite_bundle.rds + README.txt), mais n'avait aucun
# point d'entrée UI dans le rapport consolidé.
# =============================================================================

test_that("W2.3 BK-01 : mod_sc_report_consolidated.R expose dl_lite_explorer et appelle generate_biomentat_lite_app", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "modules/sc/mod_sc_report_consolidated.R"),
                         warn = FALSE), collapse = "\n")
  expect_match(src, 'downloadButton\\(ns\\("dl_lite_explorer"\\)')
  expect_match(src, 'output\\$dl_lite_explorer\\s*<-\\s*downloadHandler')
  expect_match(src, 'generate_biomentat_lite_app\\(')
})
