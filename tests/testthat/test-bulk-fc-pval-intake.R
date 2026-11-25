# Tests for BK-30 pre-computed log2FC / P-value intake lane

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
})

source_project_file("tools/pkg_skeleton/R/rdata_io.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("modules/import/mod_rdata_picker.R")
source_project_file("tools/pkg_skeleton/R/bulk_import_engine.R")
source_project_file("tools/pkg_skeleton/R/bulk_multi.R")
source_project_file("modules/import/mod_import_bulk.R")

test_that(".bulk_parse_fc_pval_table parses DESeq2/edgeR/limma and LFC-only tables", {
  df <- data.frame(
    gene = c("TP53", "MYC", "EGFR", "BRCA1"),
    logFC = c(2.1, -1.8, 0.5, NA_real_),
    P.Value = c(1e-5, 2e-4, 0.3, 0.9),
    adj.P.Val = c(1e-4, 1e-3, 0.6, 0.99),
    stringsAsFactors = FALSE
  )
  parsed <- .bulk_parse_fc_pval_table(df)
  expect_identical(nrow(parsed$table), 3L)
  expect_identical(parsed$lfc_col, "logFC")
  expect_identical(parsed$pval_col, "P.Value")
  expect_identical(parsed$padj_col, "adj.P.Val")
  expect_identical(rownames(parsed$lfc_matrix), c("TP53", "MYC", "EGFR"))
  expect_equal(as.numeric(parsed$lfc_matrix[, 1L]), c(2.1, -1.8, 0.5))

  # LFC-only table with rownames
  mat <- matrix(c(1.2, -0.9, 0.3), ncol = 1L, dimnames = list(c("G1", "G2", "G3"), "lfc"))
  p2 <- .bulk_parse_fc_pval_table(mat)
  expect_identical(nrow(p2$table), 3L)
  expect_true(all(is.na(p2$table$pvalue)))
  expect_true(all(is.na(p2$table$padj)))
})

test_that(".bulk_parse_fc_pval_table rejects invalid inputs with bulk_import_error", {
  expect_error(.bulk_parse_fc_pval_table(NULL), class = "bulk_import_error")
  expect_error(.bulk_parse_fc_pval_table(data.frame(a = c("x", "y"), stringsAsFactors = FALSE)),
               class = "bulk_import_error")
})

test_that("mod_import_bulk_ui exposes fc_pval data type choice", {
  i18n_stub <- list(t = function(k) k)
  assign("i18n", i18n_stub, envir = .GlobalEnv)
  if (!exists(".tr_plain", envir = .GlobalEnv, inherits = FALSE)) {
    assign(".tr_plain", function(k) k, envir = .GlobalEnv)
    on.exit(rm(".tr_plain", envir = .GlobalEnv), add = TRUE)
  }
  on.exit(rm("i18n", envir = .GlobalEnv), add = TRUE)

  html <- as.character(mod_import_bulk_ui("import_bulk"))
  expect_true(grepl("fc_pval", html, fixed = TRUE))
})
