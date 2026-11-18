# Tests W3.4 Quick-Wins : BK-03 (generic omics), BK-02 (CMAP offline), BK-01 (biomentat-lite)

test_that("BK-03 : prepare_generic_omics_matrix & run_generic_omics_de gèrent flottants/négatifs/NA", {
  set.seed(42)
  mat <- matrix(rnorm(60, mean = 5, sd = 1), nrow = 10, ncol = 6)
  rownames(mat) <- paste0("PROT_", seq_len(10))
  colnames(mat) <- paste0("S", seq_len(6))
  mat[1, 4:6] <- mat[1, 4:6] + 4.0
  mat[2, 1] <- NA_real_
  mat[3, 1:5] <- NA_real_ # < 50% présent -> filtré

  prep <- prepare_generic_omics_matrix(
    mat, modality = "proteomics", log_transform = FALSE,
    impute = "median", normalize = "median_center", min_present_frac = 0.5
  )
  expect_identical(prep$type, "generic_omics_matrix")
  expect_identical(prep$modality, "proteomics")
  expect_false(isTRUE(prep$qc$is_count_matrix))
  expect_identical(prep$qc$n_features_in, 10L)
  expect_identical(prep$qc$n_features_kept, 9L)
  expect_identical(prep$qc$n_na_after, 0L)

  de <- run_generic_omics_de(prep, group = c("ctrl", "ctrl", "ctrl", "case", "case", "case"),
                             case_level = "case", control_level = "ctrl")
  expect_identical(de$type, "generic_omics_de")
  expect_true(all(c("feature", "logFC", "P.Value", "adj.P.Val", "modality_mode") %in% colnames(de$table)))
  expect_true(all(de$table$modality_mode == "generic_quantitative"))
  expect_identical(de$table$feature[1L], "PROT_1")
  expect_true(de$table$logFC[1L] > 2)
})

test_that("BK-02 : score_cmap_connectivity distingue reversal et mimic hors-ligne", {
  genes <- paste0("G", seq_len(20))
  de_stats <- stats::setNames(seq(10, -10, length.out = 20), genes)

  pack <- list(
    list(id = "DRUG_REVERSAL", target = "KIN1", moa = "inhibitor",
         up = c("G18", "G19", "G20"), down = c("G1", "G2", "G3")),
    list(id = "DRUG_MIMIC", target = "KIN2", moa = "activator",
         up = c("G1", "G2", "G3"), down = c("G18", "G19", "G20"))
  )
  res <- score_cmap_connectivity(de_stats, pack)
  expect_identical(res$type, "cmap_offline_result")
  expect_true(isTRUE(res$qc$offline))
  tbl <- res$table
  rev_row <- tbl[tbl$id == "DRUG_REVERSAL", ]
  mim_row <- tbl[tbl$id == "DRUG_MIMIC", ]
  expect_true(rev_row$connectivity_score < -0.5)
  expect_identical(rev_row$direction, "reverses")
  expect_true(mim_row$connectivity_score > 0.5)
  expect_identical(mim_row$direction, "mimics")
})

test_that("BK-01 : generate_biomentat_lite_app génère un bundle autonome valide", {
  out_dir <- tempfile("biomentat_lite_")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  df1 <- data.frame(gene = c("TP53", "EGFR"), logFC = c(2.1, -1.5), stringsAsFactors = FALSE)
  res <- generate_biomentat_lite_app(
    out_dir = out_dir,
    title = "Test Lite",
    tables = list(de_results = df1),
    metadata = list(source = "unit_test")
  )
  expect_true(file.exists(res$app_path))
  expect_true(file.exists(res$bundle_path))
  expect_true(file.exists(res$readme_path))
  expect_identical(res$n_tables, 1L)
  expect_no_error(parse(file = res$app_path))
  b <- readRDS(res$bundle_path)
  expect_identical(b$title, "Test Lite")
  expect_identical(nrow(b$tables$de_results), 2L)
})
