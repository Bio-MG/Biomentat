# ==============================================================================
# test-cohort-attrition.R — BK-11 OMOP / Cohort Attrition Engine
# ==============================================================================

.att_test_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("build_attrition_table computes monotonic counts and retention percentages", {
  source(.att_test_root("tools/pkg_skeleton/R/cohort_attrition.R"), local = globalenv())

  steps <- data.frame(
    step = c("Raw", "QC", "Final"),
    n_remaining = c(1000L, 750L, 600L),
    reason = c("Initial", "Failed QC", "Missing metadata"),
    stringsAsFactors = FALSE
  )
  tbl <- build_attrition_table(steps, unit = "patients")
  expect_s3_class(tbl, "biomentat_attrition")
  expect_equal(nrow(tbl), 3L)
  expect_equal(tbl$n_in, c(1000L, 1000L, 750L))
  expect_equal(tbl$n_excluded, c(0L, 250L, 150L))
  expect_equal(tbl$pct_excluded, c(0, 25, 20))
  expect_equal(tbl$pct_cumulative_retained, c(100, 75, 60))
  expect_equal(tbl$unit, rep("patients", 3L))
})

test_that("build_attrition_table rejects non-monotonic or invalid inputs", {
  source(.att_test_root("tools/pkg_skeleton/R/cohort_attrition.R"), local = globalenv())

  bad_steps <- data.frame(
    step = c("A", "B"),
    n_remaining = c(100L, 120L),
    stringsAsFactors = FALSE
  )
  expect_error(build_attrition_table(bad_steps), class = "cohort_attrition_error")
  expect_error(build_attrition_table(bad_steps, unit = ""), class = "cohort_attrition_error")
})

test_that("OMOP and SC convenience builders + Mermaid + ggplot2 diagram work", {
  source(.att_test_root("tools/pkg_skeleton/R/cohort_attrition.R"), local = globalenv())

  omop <- attrition_from_omop_cohort(5000, 1200, 950, 800, 780)
  expect_equal(nrow(omop), 5L)
  expect_equal(omop$n_excluded, c(0L, 3800L, 250L, 150L, 20L))

  sc <- attrition_from_sc_qc(29065, 28894, 27500, 27500)
  expect_equal(nrow(sc), 4L)
  expect_equal(sc$n_excluded[2L], 171L)

  mmd <- format_attrition_mermaid(omop)
  expect_true(grepl("^flowchart TD", mmd))
  expect_true(grepl("Source Population", mmd, fixed = TRUE))
  expect_true(grepl("-3800 patients", mmd, fixed = TRUE))

  skip_if_not_installed("ggplot2")
  p <- plot_attrition_diagram(omop, title = "OMOP Cohort")
  expect_s3_class(p, "ggplot")
})
