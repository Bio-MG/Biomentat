# =============================================================================
# tests/testthat/test-bulk-dtu.R — Unit tests for pure-R DTU engine
# =============================================================================
if (!exists("run_bulk_dtu", mode = "function")) {
  if (exists("source_project_file", mode = "function")) {
    source_project_file("tools/pkg_skeleton/R/bulk_dtu.R")
  } else {
    source(file.path("..", "..", "tools", "pkg_skeleton", "R", "bulk_dtu.R"))
  }
}

test_that("run_bulk_dtu detects strong isoform switching and opposite logFC signs", {
  skip_if_not_installed("edgeR")

  set.seed(101)
  n_samples <- 6L
  groups <- factor(c("A", "A", "A", "B", "B", "B"))

  # 5 multi-isoform genes (2 tx each = 10 tx) + 1 single-isoform gene (1 tx)
  tx_ids <- c(
    "tx_sw_1", "tx_sw_2",
    "tx_bg1_1", "tx_bg1_2",
    "tx_bg2_1", "tx_bg2_2",
    "tx_bg3_1", "tx_bg3_2",
    "tx_bg4_1", "tx_bg4_2",
    "tx_single_1"
  )
  gene_ids <- c(
    "GENE_SWITCH", "GENE_SWITCH",
    "GENE_BG1", "GENE_BG1",
    "GENE_BG2", "GENE_BG2",
    "GENE_BG3", "GENE_BG3",
    "GENE_BG4", "GENE_BG4",
    "GENE_SINGLE"
  )
  tx2gene <- data.frame(
    transcript_id = tx_ids,
    gene_id = gene_ids,
    stringsAsFactors = FALSE
  )

  counts <- matrix(
    stats::rpois(length(tx_ids) * n_samples, lambda = 120),
    nrow = length(tx_ids),
    ncol = n_samples,
    dimnames = list(tx_ids, paste0("S", seq_len(n_samples)))
  )
  # Strong isoform switch in GENE_SWITCH: tx_sw_1 high in A, tx_sw_2 high in B
  counts["tx_sw_1", 1:3] <- c(800L, 820L, 790L)
  counts["tx_sw_1", 4:6] <- c(80L, 75L, 85L)
  counts["tx_sw_2", 1:3] <- c(80L, 85L, 75L)
  counts["tx_sw_2", 4:6] <- c(810L, 795L, 825L)

  res <- run_bulk_dtu(counts, tx2gene, groups, test_method = "Simes")

  expect_type(res, "list")
  expect_true(all(c("gene_table", "transcript_table", "n_genes_tested",
                    "n_transcripts_tested", "test_method", "provenance") %in% names(res)))
  expect_equal(res$n_genes_tested, 5L)
  expect_equal(res$n_transcripts_tested, 10L)
  expect_equal(res$test_method, "Simes")
  expect_equal(res$provenance$engine, "edgeR::diffSpliceDGE")
  expect_true(nzchar(res$provenance$edgeR_version))

  sw_row <- res$gene_table[res$gene_table$GeneID == "GENE_SWITCH", , drop = FALSE]
  expect_equal(nrow(sw_row), 1L)
  expect_lt(sw_row$FDR, 0.05)

  sw_tx <- res$transcript_table[res$transcript_table$GeneID == "GENE_SWITCH", , drop = FALSE]
  expect_equal(nrow(sw_tx), 2L)
  expect_lt(prod(sw_tx$logFC), 0)
})

test_that("summarize_dtu_gene returns sorted transcript rows for a gene", {
  skip_if_not_installed("edgeR")

  set.seed(202)
  groups <- factor(c("ctrl", "ctrl", "ctrl", "treat", "treat", "treat"))
  tx_ids <- c("t1", "t2", "t3", "t4", "t5", "t6")
  gene_ids <- c("G1", "G1", "G2", "G2", "G3", "G3")
  tx2gene <- data.frame(tx_id = tx_ids, gene_name = gene_ids, stringsAsFactors = FALSE)

  counts <- matrix(
    stats::rpois(36L, lambda = 150),
    nrow = 6L,
    ncol = 6L,
    dimnames = list(tx_ids, paste0("sample_", 1:6))
  )
  counts["t1", 1:3] <- 900L
  counts["t1", 4:6] <- 90L
  counts["t2", 1:3] <- 90L
  counts["t2", 4:6] <- 900L

  design <- stats::model.matrix(~ groups)
  res <- run_bulk_dtu(counts, tx2gene, design, test_method = "F")

  expect_equal(res$test_method, "F")
  expect_true("gene.F" %in% colnames(res$gene_table) || "gene.LR" %in% colnames(res$gene_table))

  sub <- summarize_dtu_gene(res, "G1")
  expect_s3_class(sub, "data.frame")
  expect_equal(nrow(sub), 2L)
  expect_true(all(sub$GeneID == "G1"))
  expect_true(sub$FDR[1] <= sub$FDR[2])
  expect_lt(prod(sub$logFC), 0)

  empty_sub <- summarize_dtu_gene(res, "NONEXISTENT_GENE")
  expect_equal(nrow(empty_sub), 0L)
})

test_that("run_bulk_dtu validates inputs cleanly", {
  skip_if_not_installed("edgeR")

  counts <- matrix(
    100L,
    nrow = 4L,
    ncol = 4L,
    dimnames = list(paste0("tx", 1:4), paste0("s", 1:4))
  )
  groups_ok <- factor(c("A", "A", "B", "B"))
  tx2gene_ok <- data.frame(
    transcript_id = paste0("tx", 1:4),
    gene_id = c("g1", "g1", "g2", "g2"),
    stringsAsFactors = FALSE
  )

  # Missing columns in tx2gene
  bad_tx2gene <- data.frame(foo = paste0("tx", 1:4), bar = c("g1", "g1", "g2", "g2"))
  expect_error(
    run_bulk_dtu(counts, bad_tx2gene, groups_ok),
    class = "bulk_dtu_error"
  )

  # Single-level group vector
  single_grp <- factor(c("A", "A", "A", "A"))
  expect_error(
    run_bulk_dtu(counts, tx2gene_ok, single_grp),
    class = "bulk_dtu_error"
  )

  # No multi-isoform genes (1 transcript per gene)
  mono_tx2gene <- data.frame(
    transcript_id = paste0("tx", 1:4),
    gene_id = c("g1", "g2", "g3", "g4"),
    stringsAsFactors = FALSE
  )
  expect_error(
    run_bulk_dtu(counts, mono_tx2gene, groups_ok),
    class = "bulk_dtu_error"
  )

  # Counts without rownames
  counts_no_rn <- counts
  rownames(counts_no_rn) <- NULL
  expect_error(
    run_bulk_dtu(counts_no_rn, tx2gene_ok, groups_ok),
    class = "bulk_dtu_error"
  )
})
