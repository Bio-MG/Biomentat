# Tests for BK-33 per-cluster k-means ORA enrichment hook (PK-3 global BH)

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
})

source_project_file("modules/bulk/mod_bulk_pattern.R")

test_that(".bulk_pattern_cluster_ora computes hypergeometric ORA and global BH across clusters", {
  univ <- paste0("G", seq_len(100))
  clusters_df <- data.frame(
    gene = c(paste0("G", 1:10), paste0("G", 11:20)),
    cluster = c(rep("Cluster_1", 10), rep("Cluster_2", 10)),
    stringsAsFactors = FALSE
  )
  gene_sets <- list(
    Set_A = paste0("G", 1:10),
    Set_B = paste0("G", 11:20),
    Set_C = c("G1", "G11", "G50", "G51")
  )

  res <- .bulk_pattern_cluster_ora(clusters_df, gene_sets, universe = univ, min_overlap = 2L)
  expect_true(is.data.frame(res))
  expect_true(all(c("cluster", "pathway", "overlap", "cluster_size", "set_size",
                    "pvalue", "padj", "genes") %in% names(res)))
  expect_identical(nrow(res), 2L)
  expect_identical(res$pathway, c("Set_A", "Set_B"))
  expect_identical(res$overlap, c(10L, 10L))
  # Global BH across both rows
  expect_equal(res$padj, stats::p.adjust(res$pvalue, method = "BH"))
  expect_true(all(res$padj < 1e-6))
})

test_that(".bulk_pattern_cluster_ora handles empty/invalid inputs cleanly", {
  res_empty <- .bulk_pattern_cluster_ora(NULL, list(A = c("G1", "G2")), universe = c("G1", "G2"))
  expect_identical(nrow(res_empty), 0L)

  clusters_df <- data.frame(gene = c("G1", "G2"), cluster = c("C1", "C1"), stringsAsFactors = FALSE)
  res_no_hit <- .bulk_pattern_cluster_ora(clusters_df, list(A = c("G9", "G10")),
                                          universe = paste0("G", 1:20), min_overlap = 2L)
  expect_identical(nrow(res_no_hit), 0L)
})

test_that("mod_bulk_pattern_output_ui includes Enrichissement par cluster tab", {
  i18n_stub <- list(t = function(k) k)
  assign("i18n", i18n_stub, envir = .GlobalEnv)
  on.exit(rm("i18n", envir = .GlobalEnv), add = TRUE)

  html <- as.character(mod_bulk_pattern_output_ui("bulk_pattern"))
  expect_true(grepl("bulk_pattern-btn_enrich_clusters", html, fixed = TRUE))
  expect_true(grepl("bulk_pattern-pattern_enrich_table", html, fixed = TRUE))
  expect_true(grepl("bulk_pattern-enrich_pack", html, fixed = TRUE))
})
