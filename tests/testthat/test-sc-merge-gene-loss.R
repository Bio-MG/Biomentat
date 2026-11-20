# =============================================================================
# test-sc-merge-gene-loss.R — W2.4 #2: SC merge gene-loss summary helper
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_helpers.R")
source_project_file("tools/pkg_skeleton/R/sc_pipeline.R")

test_that("summarize_merge_gene_loss computes intersection defaults when merged_genes is NULL", {
  gene_lists <- list(
    sampleA = c("G1", "G2", "G3", "G4"),
    sampleB = c("G2", "G3", "G5")
  )

  res <- summarize_merge_gene_loss(gene_lists)

  expect_setequal(res$union_genes, c("G1", "G2", "G3", "G4", "G5"))
  expect_setequal(res$common_genes, c("G2", "G3"))
  expect_identical(res$total_union, 5L)
  expect_identical(res$total_retained, 2L)
  expect_identical(res$total_lost_vs_union, 3L)

  expect_s3_class(res$per_sample, "data.frame")
  expect_identical(res$per_sample$sample, c("sampleA", "sampleB"))
  expect_identical(res$per_sample$n_input, c(4L, 3L))
  expect_identical(res$per_sample$n_retained, c(2L, 2L))
  expect_identical(res$per_sample$n_lost, c(2L, 1L))
  expect_equal(res$per_sample$pct_lost, c(50, 100 / 3), tolerance = 1e-6)

  expect_setequal(res$sample_lost_genes$sampleA, c("G1", "G4"))
  expect_setequal(res$sample_lost_genes$sampleB, "G5")
})

test_that("summarize_merge_gene_loss compares against explicit merged_genes (union & custom)", {
  gene_lists <- list(
    s1 = c("A", "B", "C"),
    s2 = c("B", "C", "D")
  )

  # Union merge (zero-filled): 0 genes lost
  res_union <- summarize_merge_gene_loss(gene_lists, merged_genes = c("A", "B", "C", "D"))
  expect_identical(res_union$total_union, 4L)
  expect_identical(res_union$total_retained, 4L)
  expect_identical(res_union$total_lost_vs_union, 0L)
  expect_identical(res_union$per_sample$n_lost, c(0L, 0L))
  expect_equal(res_union$per_sample$pct_lost, c(0, 0))
  expect_identical(res_union$sample_lost_genes$s1, character(0))
  expect_identical(res_union$sample_lost_genes$s2, character(0))

  # Filtered/custom merge: only A and B retained
  res_sub <- summarize_merge_gene_loss(gene_lists, merged_genes = c("A", "B"))
  expect_identical(res_sub$total_retained, 2L)
  expect_identical(res_sub$total_lost_vs_union, 2L)
  expect_identical(res_sub$per_sample$n_retained, c(2L, 1L))
  expect_identical(res_sub$per_sample$n_lost, c(1L, 2L))
  expect_setequal(res_sub$sample_lost_genes$s1, "C")
  expect_setequal(res_sub$sample_lost_genes$s2, c("C", "D"))
})

test_that("summarize_merge_gene_loss handles single sample, empty vectors, duplicates, and unnamed lists", {
  # Single sample
  res_one <- summarize_merge_gene_loss(list(only = c("X", "Y", "X", NA, "")))
  expect_identical(res_one$total_union, 2L)
  expect_identical(res_one$total_retained, 2L)
  expect_identical(res_one$total_lost_vs_union, 0L)
  expect_identical(res_one$per_sample$n_input, 2L)
  expect_identical(res_one$per_sample$n_lost, 0L)

  # Empty gene vector in one dataset -> pct_lost must be 0, not NaN
  res_empty <- summarize_merge_gene_loss(list(c("A", "B"), character(0)))
  expect_identical(res_empty$per_sample$sample, c("sample_1", "sample_2"))
  expect_identical(res_empty$common_genes, character(0))
  expect_identical(res_empty$total_retained, 0L)
  expect_identical(res_empty$total_lost_vs_union, 2L)
  expect_equal(res_empty$per_sample$pct_lost, c(100, 0))
})

test_that("summarize_merge_gene_loss rejects invalid inputs with sc_pipeline_error", {
  expect_error(summarize_merge_gene_loss(list()), class = "sc_pipeline_error")
  expect_error(summarize_merge_gene_loss(NULL), class = "sc_pipeline_error")
  expect_error(summarize_merge_gene_loss(c("A", "B")), class = "sc_pipeline_error")
  expect_error(
    summarize_merge_gene_loss(list(s1 = "A"), merged_genes = list("A")),
    class = "sc_pipeline_error"
  )
})
