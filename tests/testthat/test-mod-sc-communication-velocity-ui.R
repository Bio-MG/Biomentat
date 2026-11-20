# =============================================================================
# test-mod-sc-communication-velocity-ui.R — W2.2 §4.3 : surface UI de
# session$userData$vel_directionality (synthèse d'alignement, repli NULL)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_communication_velocity.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_velocity_directional.R")
source_project_file("modules/sc/mod_sc_communication_velocity.R")

.mvu_src <- function() {
  paste(readLines(file.path(ts_project_root(), "modules/sc/mod_sc_communication_velocity.R"),
                  encoding = "UTF-8", warn = FALSE), collapse = "\n")
}

test_that("UI expose la sortie vel_dir_summary", {
  expect_true(grepl('textOutput(ns("vel_dir_summary"))', .mvu_src(), fixed = TRUE))
})

test_that("serveur rend vel_dir_summary depuis vel_directionality() avec repli NULL", {
  src <- .mvu_src()
  expect_true(grepl("output$vel_dir_summary <- renderText", src, fixed = TRUE))
  expect_true(grepl("is.null(d)", src, fixed = TRUE))
  expect_true(grepl("transition_aligned", src, fixed = TRUE))
  expect_true(grepl('return("")', src, fixed = TRUE))
})

test_that("aucune nouvelle clé i18n et surface publique inchangée", {
  src <- .mvu_src()
  expect_false(grepl("(\\.tr|i18n\\$t)\\(\"[^\"]*(T_ij|rho|aligned)", src))
  expect_true(exists("mod_sc_communication_velocity_server", mode = "function"))
  expect_true(exists("mod_sc_communication_velocity_ui", mode = "function"))
})
