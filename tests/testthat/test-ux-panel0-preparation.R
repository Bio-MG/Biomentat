# =============================================================================
# test-ux-panel0-preparation.R — step-0 preparation marker harmonised across
# Bulk & Spatial (VELOCITY_IMPROVEMENTS_PROPOSAL_2026-11-17 #8).
# Static source checks (self-contained; no app sourcing).
# =============================================================================

.p0_root <- function() {
  for (p in c(".", "..", "../..")) {
    if (file.exists(file.path(p, "modules", "bulk", "mod_bulk.R"))) return(p)
  }
  skip("repo root not found")
}
.p0_src <- function(rel) paste(readLines(file.path(.p0_root(), rel), warn = FALSE,
                                         encoding = "UTF-8"), collapse = "\n")

test_that("bulk panel 0 carries the shared preparation marker, value preserved", {
  s <- .p0_src("modules/bulk/mod_bulk.R")
  expect_match(s, 'panel_0_preparation"', fixed = TRUE)
  expect_match(s, "grp_0_preparation", fixed = TRUE)
  expect_match(s, 'value = "0_autopipeline"', fixed = TRUE)
  expect_match(s, "BULK-9", fixed = TRUE)
})

test_that("spatial panel 0 carries the shared preparation marker, value preserved", {
  s <- .p0_src("modules/spatial/mod_spatial.R")
  expect_match(s, 'panel_0_preparation"', fixed = TRUE)
  expect_match(s, "grp_0_preparation", fixed = TRUE)
  expect_match(s, 'value = "panel_pipeline"', fixed = TRUE)
})
