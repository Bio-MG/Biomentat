# Tests for BK-31 PCA loadings per-PC pathway enrichment (global BH PK-3)

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
})

source_project_file("modules/bulk/mod_bulk_filter.R")

test_that(".bulk_pca_loadings_enrichment detects pathway driving PC1 with global BH", {
  set.seed(42L)
  ng <- 60L; ns <- 6L
  mat <- matrix(stats::rnorm(ng * ns, mean = 0, sd = 0.2), nrow = ng, ncol = ns)
  rownames(mat) <- paste0("G", seq_len(ng))
  colnames(mat) <- paste0("S", seq_len(ns))
  # Make G1..G10 strongly separate S1..3 vs S4..6 (driving PC1)
  mat[1:10, 1:3] <- mat[1:10, 1:3] + 4
  mat[1:10, 4:6] <- mat[1:10, 4:6] - 4

  gene_sets <- list(
    PC1_Driver = paste0("G", 1:10),
    Random_Set = paste0("G", 31:40)
  )
  res <- .bulk_pca_loadings_enrichment(mat, gene_sets, max_pcs = 3L, min_set_size = 5L)
  expect_true(is.data.frame(res))
  expect_true(all(c("pc", "pathway", "set_size", "z_score", "direction", "pvalue", "padj") %in% names(res)))
  pc1_drv <- res[res$pc == "PC1" & res$pathway == "PC1_Driver", , drop = FALSE]
  expect_identical(nrow(pc1_drv), 1L)
  expect_true(abs(pc1_drv$z_score) > 3)
  expect_true(pc1_drv$padj < 0.01)
  expect_equal(res$padj, stats::p.adjust(res$pvalue, method = "BH"))
})

test_that(".bulk_pca_loadings_enrichment handles invalid/small inputs cleanly", {
  expect_identical(nrow(.bulk_pca_loadings_enrichment(NULL, list(A = letters[1:6]))), 0L)
  small <- matrix(1:4, nrow = 2, dimnames = list(c("G1", "G2"), c("S1", "S2")))
  expect_identical(nrow(.bulk_pca_loadings_enrichment(small, list(A = c("G1", "G2")))), 0L)
})

test_that("mod_bulk_filter_pca_ui exposes btn_pca_enrich and pca_enrich_table", {
  i18n_stub <- list(t = function(k) k)
  assign("i18n", i18n_stub, envir = .GlobalEnv)
  assign(".tr_placeholder", function() "Select...", envir = .GlobalEnv)
  on.exit({
    rm("i18n", envir = .GlobalEnv)
    rm(".tr_placeholder", envir = .GlobalEnv)
  }, add = TRUE)

  html <- as.character(mod_bulk_filter_pca_ui("bulk_filter"))
  expect_true(grepl("bulk_filter-btn_pca_enrich", html, fixed = TRUE))
  expect_true(grepl("bulk_filter-pca_enrich_table", html, fixed = TRUE))
})
