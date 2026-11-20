# =============================================================================
# test-mod-sc-pipeline-qc-label.R — W2.4 #2 : Honest QC label in mod_sc_pipeline.R
# =============================================================================

source_project_file("modules/sc/mod_sc_pipeline.R")

test_that("QC summary helper returns none when no qc_snapshot is available", {
  res_null <- .sc_pipeline_qc_summary(NULL, params = list(qc_min_gene = 100, qc_max_gene = 8000, qc_mt = 20))
  expect_identical(res_null$status, "none")
  expect_null(res_null$text)

  empty_df <- data.frame(
    Echantillon = character(0),
    Cellules_avant = integer(0),
    Cellules_apres = integer(0),
    stringsAsFactors = FALSE
  )
  res_empty <- .sc_pipeline_qc_summary(empty_df)
  expect_identical(res_empty$status, "none")
  expect_null(res_empty$text)
})

test_that("QC summary accurately reports retained/removed cells, features, and applied thresholds", {
  snap <- data.frame(
    Echantillon = c("S1", "S2"),
    Cellules_avant = c(120L, 80L),
    Cellules_apres = c(100L, 65L),
    Pct_Mito_Median = c(4.2, 6.1),
    Pct_Conserve = c(83.3, 81.2),
    stringsAsFactors = FALSE
  )
  obj <- matrix(0, nrow = 1450L, ncol = 165L)
  params <- list(qc_min_gene = 150, qc_max_gene = 6500, qc_mt = 15, doublets_enabled = FALSE)

  out <- .sc_pipeline_qc_summary(snap, params = params, obj = obj)
  expect_identical(out$status, "applied")
  expect_identical(out$n_before, 200L)
  expect_identical(out$n_retained, 165L)
  expect_identical(out$n_removed, 35L)
  expect_identical(out$n_features, 1450L)
  expect_equal(out$pct_retained, 82.5)

  expect_match(out$counts_text, "165", fixed = TRUE)
  expect_match(out$counts_text, "200", fixed = TRUE)
  expect_match(out$counts_text, "35", fixed = TRUE)
  expect_match(out$counts_text, "82.5%", fixed = TRUE)
  expect_match(out$counts_text, "1 450|1450", fixed = FALSE)

  expect_match(out$thresholds_text, "nFeature_RNA > 150", fixed = TRUE)
  expect_match(out$thresholds_text, "nFeature_RNA < 6500", fixed = TRUE)
  expect_match(out$thresholds_text, "% Mito < 15%", fixed = TRUE)
  expect_null(out$doublets_text)
})

test_that("QC summary accounts for excluded doublets and never fabricates doublet filtering when skipped", {
  snap_dbl <- data.frame(
    Echantillon = c("S1", "S2"),
    Cellules_avant = c(100L, 100L),
    Cellules_apres = c(90L, 80L),
    Doublets_exclus = c(5L, 3L),
    stringsAsFactors = FALSE
  )
  obj_dbl <- matrix(0, nrow = 1200L, ncol = 162L)
  params_dbl <- list(
    qc_min_gene = 100, qc_max_gene = 8000, qc_mt = 20,
    doublets_enabled = TRUE, doublet_score_cutoff = 0.3
  )

  out_dbl <- .sc_pipeline_qc_summary(snap_dbl, params = params_dbl, obj = obj_dbl)
  expect_identical(out_dbl$n_before, 200L)
  expect_identical(out_dbl$n_retained, 162L)
  expect_identical(out_dbl$n_removed, 38L)
  expect_identical(out_dbl$n_doublets, 8L)
  expect_match(out_dbl$doublets_text, "8", fixed = TRUE)
  expect_match(out_dbl$doublets_text, "0.3", fixed = TRUE)

  # Even if doublets_enabled = TRUE in UI, if Doublets_exclus is absent (skipped),
  # the label must NOT claim doublets were filtered.
  snap_no_dbl <- snap_dbl[, c("Echantillon", "Cellules_avant", "Cellules_apres")]
  obj_no_dbl <- matrix(0, nrow = 1200L, ncol = 170L)
  out_no_dbl <- .sc_pipeline_qc_summary(snap_no_dbl, params = params_dbl, obj = obj_no_dbl)
  expect_identical(out_no_dbl$n_retained, 170L)
  expect_null(out_no_dbl$doublets_text)
})

test_that("mod_sc_pipeline.R wires uiOutput(ns('qc_summary_label')) and captures launch params", {
  src <- paste(
    readLines(file.path(ts_project_root(), "modules/sc/mod_sc_pipeline.R"),
              warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  expect_match(src, 'uiOutput(ns("qc_summary_label"))', fixed = TRUE)
  expect_match(src, "output$qc_summary_label <- renderUI(", fixed = TRUE)
  expect_match(src, ".sc_pipeline_qc_summary(", fixed = TRUE)
})
