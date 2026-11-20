# =============================================================================
# test-sc-merge-gene-report.R — W2.4 SC merge gene-loss report across all paths
# =============================================================================

test_that("mod_import_sc uses .sc_log_merge_gene_report across all 3 merge paths", {
  repo_root <- if (file.exists("modules/import/mod_import_sc.R")) "." else "../.."
  src <- paste(readLines(file.path(repo_root, "modules/import/mod_import_sc.R"), warn = FALSE),
               collapse = "\n")

  expect_match(src, "\\.sc_log_merge_gene_report\\s*<-\\s*function")
  hits <- gregexpr("\\.sc_log_merge_gene_report\\(", src)[[1]]
  # 3 call sites (design multi-sample, multi-file, 10x_multi)
  expect_gte(length(hits[hits > 0L]), 3L)
})
