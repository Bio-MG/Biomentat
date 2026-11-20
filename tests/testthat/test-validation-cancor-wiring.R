# =============================================================================
# test-validation-cancor-wiring.R — W2.3 PA-5 check_design_cancor() dans validate_bulk_design()
# =============================================================================
# `check_design_cancor()` (validation.R:69) détecte la confusion numérique
# (`|cor| > 0.9`), mais `validate_bulk_design()` n'appelait que
# `check_design_confounding()` (catégoriel uniquement).
# =============================================================================

source_project_file("tools/pkg_skeleton/R/validation.R")

test_that("W2.3 PA-5 : validate_bulk_design() signale une covariable numérique confondue", {
  meta <- data.frame(
    condition = factor(c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt")),
    age = c(20, 21, 22, 60, 61, 62),
    rnd = c(10, 50, 20, 12, 48, 22),
    stringsAsFactors = FALSE
  )
  probs_hit <- validate_bulk_design(meta, "condition", covariates = "age")
  expect_true(length(probs_hit) >= 1L)
  expect_match(paste(probs_hit, collapse = "\n"), "age")

  probs_ok <- validate_bulk_design(meta, "condition", covariates = "rnd")
  expect_length(probs_ok, 0L)
})
