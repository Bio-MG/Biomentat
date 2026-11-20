# =============================================================================
# test-presto-markers.R — W2.3 presto::wilcoxauc() acceleration for FindAllMarkers
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_helpers.R")

test_that("sc_fast_find_all_markers uses presto::wilcoxauc for wilcox and returns Seurat schema", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("presto")

  set.seed(42)
  counts <- matrix(rpois(20 * 24, lambda = 2), nrow = 20, ncol = 24)
  counts[1:5, 1:12] <- counts[1:5, 1:12] + 30L
  counts[6:10, 13:24] <- counts[6:10, 13:24] + 30L
  rownames(counts) <- paste0("G", seq_len(20))
  colnames(counts) <- paste0("C", seq_len(24))

  obj <- SeuratObject::CreateSeuratObject(counts = counts)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  SeuratObject::Idents(obj) <- factor(c(rep("A", 12), rep("B", 12)))

  res_presto <- sc_fast_find_all_markers(
    obj,
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    only.pos = TRUE,
    verbose = FALSE
  )
  expect_true(is.data.frame(res_presto))
  expect_identical(attr(res_presto, "engine"), "presto")
  expect_true(all(c("gene", "cluster", "p_val", "p_val_adj", "avg_log2FC", "pct.1", "pct.2") %in% colnames(res_presto)))
  expect_true(nrow(res_presto) > 0L)
  expect_true(all(res_presto$avg_log2FC >= 0.25))

  res_seurat <- sc_fast_find_all_markers(
    obj,
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0.25,
    only.pos = TRUE,
    verbose = FALSE,
    use_presto = FALSE
  )
  expect_true(is.data.frame(res_seurat))
  expect_identical(attr(res_seurat, "engine"), "seurat")
})

test_that("SC marker sites route through sc_fast_find_all_markers", {
  for (f in c("modules/sc/mod_sc_markers.R", "modules/sc/mod_sc.R", "tools/pkg_skeleton/R/sc_pipeline.R")) {
    p <- if (file.exists(f)) f else file.path("../..", f)
    txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
    expect_match(txt, "sc_fast_find_all_markers")
  }
})
