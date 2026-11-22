# =============================================================================
# tests/testthat/test-bulk-mnar-de.R — Unit tests for MNAR/MAR DE engine (BK-05)
# =============================================================================
if (!exists("run_bulk_mnar_de", mode = "function")) {
  if (exists("source_project_file", mode = "function")) {
    source_project_file("tools/pkg_skeleton/R/bulk_mnar_de.R")
  } else {
    source(file.path("..", "..", "tools", "pkg_skeleton", "R", "bulk_mnar_de.R"))
  }
}

test_that("classify_missingness_regime classifies complete, MNAR, and MAR features accurately", {
  groups <- factor(c(rep("ctrl", 5L), rep("case", 5L)), levels = c("ctrl", "case"))

  mat <- matrix(
    10.0,
    nrow = 4L,
    ncol = 10L,
    dimnames = list(
      c("feat_complete", "feat_mnar_fisher", "feat_mnar_half", "feat_mar"),
      paste0("S", 1:10)
    )
  )

  # feat_complete: 0 NAs -> "complete"
  # feat_mnar_fisher: all NA in ctrl (5/5 NA), all observed in case (0/5 NA) -> Fisher p < 0.05 -> "MNAR"
  mat["feat_mnar_fisher", 1:5] <- NA_real_

  # feat_mnar_half: all NA in ctrl (5/5 NA), 3/5 (60% >= 50%) observed in case -> "MNAR"
  mat["feat_mnar_half", 1:5] <- NA_real_
  mat["feat_mnar_half", 6:7] <- NA_real_

  # feat_mar: 1 NA in ctrl, 1 NA in case -> balanced missingness -> "MAR"
  mat["feat_mar", c(1L, 6L)] <- NA_real_

  cls <- classify_missingness_regime(mat, groups, mnar_p_cutoff = 0.05)

  expect_s3_class(cls, "data.frame")
  expect_true(all(c("feature", "pct_na_g1", "pct_na_g2", "pvalue_missingness", "regime") %in% colnames(cls)))
  expect_equal(cls["feat_complete", "regime"], "complete")
  expect_equal(cls["feat_complete", "pct_na_g1"], 0)
  expect_equal(cls["feat_complete", "pct_na_g2"], 0)

  expect_equal(cls["feat_mnar_fisher", "regime"], "MNAR")
  expect_equal(cls["feat_mnar_fisher", "pct_na_g1"], 1)
  expect_equal(cls["feat_mnar_fisher", "pct_na_g2"], 0)
  expect_lt(cls["feat_mnar_fisher", "pvalue_missingness"], 0.05)

  expect_equal(cls["feat_mnar_half", "regime"], "MNAR")
  expect_equal(cls["feat_mnar_half", "pct_na_g1"], 1)
  expect_equal(cls["feat_mnar_half", "pct_na_g2"], 0.4)

  expect_equal(cls["feat_mar", "regime"], "MAR")
  expect_equal(cls["feat_mar", "pct_na_g1"], 0.2)
  expect_equal(cls["feat_mar", "pct_na_g2"], 0.2)
  expect_gt(cls["feat_mar", "pvalue_missingness"], 0.5)
})

test_that("run_bulk_mnar_de detects pure-MNAR, complete DE, and MAR features and preserves .Random.seed", {
  set.seed(12345L)
  groups <- c(rep("ctrl", 5L), rep("case", 5L))
  n_bg <- 20L
  feat_ids <- c("PROT_COMPLETE_DE", "PROT_PURE_MNAR", "PROT_MAR_NULL", "PROT_TOO_FEW", paste0("BG_", seq_len(n_bg)))

  mat <- matrix(
    stats::rnorm(length(feat_ids) * 10L, mean = 20, sd = 0.4),
    nrow = length(feat_ids),
    ncol = 10L,
    dimnames = list(feat_ids, paste0("S", 1:10))
  )

  # 1. Complete DE feature: strong upregulation in case
  mat["PROT_COMPLETE_DE", 1:5] <- c(18.0, 18.2, 17.9, 18.1, 18.0)
  mat["PROT_COMPLETE_DE", 6:10] <- c(23.1, 22.9, 23.0, 23.2, 22.8)

  # 2. Pure-MNAR feature: present in case (~22, sd ~0.5), all-NA in ctrl
  mat["PROT_PURE_MNAR", 1:5] <- NA_real_
  mat["PROT_PURE_MNAR", 6:10] <- c(22.0, 22.7, 21.4, 22.3, 21.6)

  # 3. MAR feature: 1 NA in each group, no true difference
  mat["PROT_MAR_NULL", 1:5] <- c(NA_real_, 20.0, 20.1, 19.9, 20.0)
  mat["PROT_MAR_NULL", 6:10] <- c(20.0, NA_real_, 19.9, 20.1, 20.0)

  # 4. Feature with < min_valid_total (only 1 non-NA value) -> filtered out
  mat["PROT_TOO_FEW", ] <- NA_real_
  mat["PROT_TOO_FEW", 6L] <- 21.0

  set.seed(999L)
  seed_before <- get(".Random.seed", envir = .GlobalEnv)

  res <- run_bulk_mnar_de(
    mat = mat,
    groups = groups,
    case_level = "case",
    ctrl_level = "ctrl",
    mnar_shift_sd = 1.8,
    mnar_width_sd = 0.3,
    min_valid_total = 2L,
    seed = 42L,
    fdr_method = "BH"
  )

  seed_after <- get(".Random.seed", envir = .GlobalEnv)
  expect_identical(seed_before, seed_after)

  expect_type(res, "list")
  expect_true(all(c("table", "n_features", "n_mnar", "n_mar", "n_complete", "provenance") %in% names(res)))
  expect_equal(res$n_features, nrow(mat) - 1L)
  expect_false("PROT_TOO_FEW" %in% res$table$feature)
  expect_equal(res$n_mnar, 1L)
  expect_equal(res$n_mar, 1L)
  expect_equal(res$n_complete, n_bg + 1L)
  expect_equal(res$provenance$engine, "biomentat_mnar_de")
  expect_equal(res$provenance$case_level, "case")
  expect_equal(res$provenance$ctrl_level, "ctrl")

  # Sorted by FDR, pvalue_combined
  expect_false(is.unsorted(res$table$FDR))

  # Pure-MNAR row check
  mnar_row <- res$table[res$table$feature == "PROT_PURE_MNAR", , drop = FALSE]
  expect_equal(mnar_row$regime, "MNAR")
  expect_gt(mnar_row$log2FC, 0.5)
  expect_lt(mnar_row$pvalue_missingness, 0.05)
  expect_lt(mnar_row$pvalue_combined, mnar_row$pvalue_intensity)
  expect_lt(mnar_row$FDR, 0.05)

  # Complete DE row check
  comp_row <- res$table[res$table$feature == "PROT_COMPLETE_DE", , drop = FALSE]
  expect_equal(comp_row$regime, "complete")
  expect_gt(comp_row$log2FC, 4.0)
  expect_equal(comp_row$pvalue_combined, comp_row$pvalue_intensity)
  expect_lt(comp_row$FDR, 0.05)

  # MAR row check
  mar_row <- res$table[res$table$feature == "PROT_MAR_NULL", , drop = FALSE]
  expect_equal(mar_row$regime, "MAR")
  expect_lt(abs(mar_row$log2FC), 0.2)
  expect_gt(mar_row$FDR, 0.2)
})

test_that("classify_missingness_regime and run_bulk_mnar_de reject invalid inputs cleanly", {
  mat_ok <- matrix(
    c(10, 11, 12, 13, 14, 15),
    nrow = 2L,
    ncol = 3L,
    dimnames = list(c("F1", "F2"), c("S1", "S2", "S3"))
  )
  grp_ok <- c("ctrl", "ctrl", "case")

  expect_error(
    classify_missingness_regime(NULL, grp_ok),
    class = "bulk_mnar_error"
  )
  expect_error(
    classify_missingness_regime(unname(mat_ok), grp_ok),
    class = "bulk_mnar_error"
  )
  expect_error(
    classify_missingness_regime(mat_ok, c("ctrl", "case")),
    class = "bulk_mnar_error"
  )
  expect_error(
    classify_missingness_regime(mat_ok, c("ctrl", "ctrl", "ctrl")),
    class = "bulk_mnar_error"
  )
  expect_error(
    run_bulk_mnar_de(mat_ok, grp_ok, case_level = "nonexistent"),
    class = "bulk_mnar_error"
  )

  mat_all_na <- mat_ok
  mat_all_na[] <- NA_real_
  expect_error(
    run_bulk_mnar_de(mat_all_na, grp_ok, min_valid_total = 2L),
    class = "bulk_mnar_error"
  )
})
