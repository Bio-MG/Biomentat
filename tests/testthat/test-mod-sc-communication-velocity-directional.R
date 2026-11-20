# =============================================================================
# test-mod-sc-communication-velocity-directional.R — W2.2 §4.3 : câblage de la
# directionnalité vélocité (T_ij + corrélation vitesse-récepteur) dans le module
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_communication_velocity.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_velocity_directional.R")
source_project_file("modules/sc/mod_sc_communication_velocity.R")

.mvd_src <- function() {
  paste(readLines(file.path(ts_project_root(), "modules/sc/mod_sc_communication_velocity.R"),
                  encoding = "UTF-8", warn = FALSE), collapse = "\n")
}

test_that("module appelle compute_communication_velocity_directionality sous garde exists()", {
  src <- .mvd_src()
  expect_true(grepl("compute_communication_velocity_directionality", src, fixed = TRUE))
  expect_true(grepl('exists\\("compute_communication_velocity_directionality",\\s*mode = "function"', src))
  expect_true(grepl("tryCatch", src, fixed = TRUE))
})

test_that("T_ij et vecteurs de vitesse absents => NULL, jamais fabriqués", {
  src <- .mvd_src()
  expect_true(grepl("transition_matrix", src, fixed = TRUE))
  expect_true(grepl("vel_directionality", src, fixed = TRUE))
  expect_false(grepl("velocity_graph|scv\\$", src))
})

test_that("surface top-level du module inchangée et aucune nouvelle clé i18n directionnelle", {
  src <- .mvd_src()
  expect_false(grepl("(\\.tr|i18n\\$t)\\(\"[^\"]*[Dd]irectionnalit", src))
  expect_true(exists("mod_sc_communication_velocity_server", mode = "function"))
  expect_true(exists("mod_sc_communication_velocity_ui", mode = "function"))
})
