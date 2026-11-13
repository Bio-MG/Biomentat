# =============================================================================
# test-packs-mouse-orthology.R — P2.4 (D4) : orthologie souris -> humain
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/pathway_helpers.R")

test_that("P2.4 D4 : packs_map_mouse_to_human_symbols() mappe deterministiquement hors-ligne", {
  skip_if_not_installed("babelgene")
  withr::local_options(list(biomentat.offline_strict = TRUE))

  res <- packs_map_mouse_to_human_symbols(c("Cd4", "Cd8a", "Trp53", "Actb", "FakeMouseGene999"))
  expect_true(is.list(res))
  expect_true(all(c("CD4", "CD8A", "TP53", "ACTB") %in% res$human_genes))
  expect_true(is.data.frame(res$mapping))
  expect_true(all(c("mouse_symbol", "human_symbol", "support_n") %in% colnames(res$mapping)))
  expect_identical(res$stats$method, "babelgene")
  expect_identical(as.integer(res$stats$n_input), 5L)
  expect_gte(as.integer(res$stats$n_mapped), 4L)
  expect_gte(as.integer(res$stats$n_unmapped), 1L)
})

test_that("P2.4 D4 : refus honnete pathway_error sur entree invalide ou vide", {
  expect_error(packs_map_mouse_to_human_symbols(character(0)), class = "pathway_error")
  expect_error(packs_map_mouse_to_human_symbols(c("", NA_character_)), class = "pathway_error")
})

test_that("P2.4 D4 : run_packs_enrichment(..., organism = 'mouse') enrichit un GMT humain", {
  skip_if_not_installed("babelgene")
  skip_if_not_installed("clusterProfiler")

  t2g <- data.frame(
    term = c(rep("T_CELL_ACTIVATION", 10L), rep("OTHER_PATHWAY", 10L)),
    gene = c("CD4", "CD8A", "LCK", "ZAP70", paste0("TGENE", 1:6), paste0("OGENE", 1:10)),
    stringsAsFactors = FALSE
  )
  res <- run_packs_enrichment(
    genes = c("Cd4", "Cd8a", "Lck", "Zap70"),
    term2gene = t2g,
    universe = c("CD4", "CD8A", "LCK", "ZAP70", paste0("TGENE", 1:6), paste0("OGENE", 1:10), "TP53", "ACTB"),
    pval_cutoff = 1,
    qval_cutoff = 1,
    organism = "mouse"
  )
  expect_true(is.data.frame(res))
  expect_true("T_CELL_ACTIVATION" %in% res$ID)
  st <- attr(res, "orthology_stats", exact = TRUE)
  expect_true(is.list(st))
  expect_identical(st$method, "babelgene")
  expect_gte(as.integer(st$n_mapped), 4L)
})
