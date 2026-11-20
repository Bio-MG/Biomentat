# =============================================================================
# test-ccc-consensus.R — W2.2 §11.8 Gem #2 : consensus multi-methodes CCC
# par rangs intra-methode (CellChat x LIANA)
# =============================================================================

.consensus_load_views <- function() {
  repo_root <- if (file.exists("tools/pkg_skeleton/R/sc_communication_views.R")) "." else "../.."
  if (!exists("finalize_communication_result", mode = "function")) {
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication.R"), local = FALSE)
  }
  if (!exists("plot_communication_dotplot", mode = "function")) {
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication_views.R"), local = FALSE)
  }
}

.consensus_make_score_result <- function(scores = c(0.9, 0.8, 0.5, 0.4, 0.2, 0.1),
                                         senders = c("CD4 T", "B", "CD4 T", "NK", "B", "NK"),
                                         receivers = c("B", "CD8 T", "CD8 T", "B", "NK", "CD4 T"),
                                         ligands = paste0("L", seq_along(scores)),
                                         receptors = paste0("R", seq_along(scores)),
                                         sender_mapped = NULL,
                                         receiver_mapped = NULL,
                                         source_method = "cellchat") {
  n <- length(scores)
  tbl <- data.frame(
    sender = senders,
    receiver = receivers,
    ligand = ligands,
    receptor = receptors,
    interaction = paste0(ligands, "_", receptors),
    pathway = paste0("P", seq_len(n)),
    score = as.numeric(scores),
    p_value = rep(0.01, n),
    p_adjusted = rep(0.02, n),
    source_method = rep(source_method, n),
    source_file = rep("mock_score.csv", n),
    source_cell_identity_level = rep("celltype", n),
    stringsAsFactors = FALSE
  )
  if (!is.null(sender_mapped)) tbl$sender_mapped <- sender_mapped
  if (!is.null(receiver_mapped)) tbl$receiver_mapped <- receiver_mapped
  finalize_communication_result(canonical_table = tbl, source_method = source_method)
}

.consensus_make_rank_result <- function(ranks = c(1, 2, 3, 4, 5, 6),
                                        senders = c("CD4 T", "B", "CD4 T", "NK", "B", "NK"),
                                        receivers = c("B", "CD8 T", "CD8 T", "B", "NK", "CD4 T"),
                                        ligands = paste0("L", seq_along(ranks)),
                                        receptors = paste0("R", seq_along(ranks)),
                                        sender_mapped = NULL,
                                        receiver_mapped = NULL,
                                        source_method = "liana") {
  n <- length(ranks)
  tbl <- data.frame(
    sender = senders,
    receiver = receivers,
    ligand = ligands,
    receptor = receptors,
    interaction = paste0(ligands, "_", receptors),
    pathway = rep(NA_character_, n),
    score = rep(NA_real_, n),
    p_value = rep(0.005, n),
    p_adjusted = rep(NA_real_, n),
    source_method = rep(source_method, n),
    source_file = rep("mock_rank.csv", n),
    source_cell_identity_level = rep("celltype", n),
    rank = as.numeric(ranks),
    rank_direction = rep("lower_is_better", n),
    rank_aggregation_mode = rep("specificity", n),
    stringsAsFactors = FALSE
  )
  if (!is.null(sender_mapped)) tbl$sender_mapped <- sender_mapped
  if (!is.null(receiver_mapped)) tbl$receiver_mapped <- receiver_mapped
  finalize_communication_result(
    canonical_table = tbl,
    source_method = source_method,
    external_consensus = TRUE
  )
}

test_that("build_communication_consensus validates inputs and raises communication_context_error", {
  .consensus_load_views()
  r_cc <- .consensus_make_score_result()
  r_li <- .consensus_make_rank_result()

  expect_error(
    build_communication_consensus(NULL),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_consensus(list(r_cc)),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_consensus(list(r_cc, list(type = "invalid"))),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_consensus(list(r_cc, r_li), top_pct = 0),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_consensus(list(r_cc, r_li), top_pct = 1.1),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_consensus(list(r_cc, r_li), top_pct = NA_real_),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_consensus(list(r_cc, r_li), min_methods = 1L),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_consensus(list(r_cc, r_li), min_methods = NA_integer_),
    class = "communication_context_error"
  )
})

test_that("build_communication_consensus computes within-method percentiles for CellChat x LIANA without averaging raw scores", {
  .consensus_load_views()
  # 5 pairs in CellChat (scores: 100, 80, 50, 20, 10 -> best = 100 -> pct 0.2, 0.4, 0.6, 0.8, 1.0)
  r_cc <- .consensus_make_score_result(
    scores = c(100, 80, 50, 20, 10),
    senders = c("CD4 T", "B", "CD4 T", "NK", "B"),
    receivers = c("B", "CD8 T", "CD8 T", "B", "NK"),
    ligands = c("L1", "L2", "L3", "L4", "L5"),
    receptors = c("R1", "R2", "R3", "R4", "R5")
  )
  # 5 pairs in LIANA (ranks: 1, 2, 10, 20, 50 -> best = 1 -> pct 0.2, 0.4, 0.6, 0.8, 1.0)
  # Notice L1_R1 and L2_R2 are top 2 in both (pct <= 0.4), while L6_R6 is only in LIANA
  r_li <- .consensus_make_rank_result(
    ranks = c(1, 2, 10, 20, 50),
    senders = c("CD4 T", "B", "CD4 T", "NK", "Mono"),
    receivers = c("B", "CD8 T", "CD8 T", "B", "B"),
    ligands = c("L1", "L2", "L3", "L4", "L6"),
    receptors = c("R1", "R2", "R3", "R4", "R6")
  )

  res <- build_communication_consensus(
    list(cellchat = r_cc, liana = r_li),
    top_pct = 0.40,
    min_methods = 2L
  )

  expect_true(is.list(res))
  expect_identical(
    res$epistemic_note,
    "Consensus par rangs intra-méthode — scores bruts jamais moyennés entre outils"
  )
  expect_identical(res$methods, c("cellchat", "liana"))
  expect_equal(res$params$top_pct, 0.40)
  expect_identical(res$params$min_methods, 2L)
  expect_identical(res$qc$n_methods, 2L)
  expect_identical(res$qc$n_total_pairs, 6L)
  expect_identical(res$qc$n_consensus, 2L)

  tbl <- res$table
  expect_true(is.data.frame(tbl))
  req_cols <- c(
    "sender_node", "receiver_node", "ligand", "receptor",
    "n_methods_present", "n_methods_top", "is_consensus",
    "mean_within_rank_pct", "methods_present",
    "pct_cellchat", "pct_liana"
  )
  expect_true(all(req_cols %in% colnames(tbl)))
  # Raw scores are NEVER averaged into a consensus score column
  expect_false(any(c("mean_score", "consensus_score", "raw_score_mean") %in% colnames(tbl)))

  # Check ordering: consensus first (L1_R1 pct 0.2, then L2_R2 pct 0.4)
  expect_identical(tbl$ligand[1:2], c("L1", "L2"))
  expect_identical(tbl$is_consensus[1:2], c(TRUE, TRUE))
  expect_identical(tbl$n_methods_top[1:2], c(2L, 2L))
  expect_equal(tbl$pct_cellchat[1:2], c(0.2, 0.4))
  expect_equal(tbl$pct_liana[1:2], c(0.2, 0.4))
  expect_equal(tbl$mean_within_rank_pct[1:2], c(0.2, 0.4))
  expect_identical(tbl$methods_present[1:2], c("cellchat, liana", "cellchat, liana"))

  # Non-consensus pairs
  expect_true(all(!tbl$is_consensus[3:6]))
  l6_row <- tbl[tbl$ligand == "L6", ]
  expect_identical(l6_row$n_methods_present, 1L)
  expect_true(is.na(l6_row$pct_cellchat))
  expect_equal(l6_row$pct_liana, 1.0)
  expect_identical(l6_row$methods_present, "liana")
})

test_that("build_communication_consensus handles unnamed lists, duplicate source_methods, and 3+ methods", {
  .consensus_load_views()
  r1 <- .consensus_make_score_result(
    scores = c(0.9, 0.5, 0.1),
    senders = c("A", "B", "C"),
    receivers = c("B", "C", "A"),
    ligands = c("L1", "L2", "L3"),
    receptors = c("R1", "R2", "R3"),
    source_method = "cellchat"
  )
  r2 <- .consensus_make_rank_result(
    ranks = c(1, 2, 3),
    senders = c("A", "B", "C"),
    receivers = c("B", "C", "A"),
    ligands = c("L1", "L2", "L3"),
    receptors = c("R1", "R2", "R3"),
    source_method = "liana"
  )
  r3 <- .consensus_make_score_result(
    scores = c(10, 2, 8),
    senders = c("A", "B", "C"),
    receivers = c("B", "C", "A"),
    ligands = c("L1", "L2", "L3"),
    receptors = c("R1", "R2", "R3"),
    source_method = "cellphonedb"
  )

  # Unnamed list derives names from source_method
  res3 <- build_communication_consensus(
    list(r1, r2, r3),
    top_pct = 0.40,
    min_methods = 3L
  )
  expect_identical(res3$methods, c("cellchat", "liana", "cellphonedb"))
  expect_true(all(c("pct_cellchat", "pct_liana", "pct_cellphonedb") %in% colnames(res3$table)))
  # L1_R1 is rank 1/3 (pct = 1/3 <= 0.40) in all 3 methods
  expect_identical(res3$qc$n_consensus, 1L)
  expect_identical(res3$table$ligand[1L], "L1")
  expect_identical(res3$table$n_methods_top[1L], 3L)

  # Duplicate source_method in unnamed list gets disambiguated deterministically
  res_dup <- build_communication_consensus(list(r1, r1), top_pct = 0.5)
  expect_equal(length(unique(res_dup$methods)), 2L)
})

test_that("build_communication_consensus harmonizes node keys via .communication_node_keys and handles ties deterministically", {
  .consensus_load_views()
  # Method 1 uses raw labels "T_CD4" mapped to "CD4 T", with tied scores
  r1 <- .consensus_make_score_result(
    scores = c(0.8, 0.8, 0.2, 0.1),
    senders = c("T_CD4", "T_CD4", "B_cell", "NK"),
    receivers = c("B_cell", "NK", "T_CD4", "B_cell"),
    ligands = c("L1", "L2", "L3", "L4"),
    receptors = c("R1", "R2", "R3", "R4"),
    sender_mapped = c("CD4 T", "CD4 T", "B", "NK"),
    receiver_mapped = c("B", "NK", "CD4 T", "B")
  )
  # Method 2 uses already-harmonized labels "CD4 T" and "B", with tied ranks
  r2 <- .consensus_make_rank_result(
    ranks = c(1, 1, 3, 4),
    senders = c("CD4 T", "CD4 T", "B", "NK"),
    receivers = c("B", "NK", "CD4 T", "B"),
    ligands = c("L1", "L2", "L3", "L4"),
    receptors = c("R1", "R2", "R3", "R4")
  )

  res <- build_communication_consensus(list(cellchat = r1, liana = r2), top_pct = 0.30)
  expect_identical(res$qc$n_total_pairs, 4L)
  expect_identical(res$qc$n_consensus, 2L)
  # Ties on best rank/score both get min percentile 1/4 = 0.25 <= 0.30
  top_rows <- res$table[res$table$is_consensus, ]
  expect_equal(top_rows$pct_cellchat, c(0.25, 0.25))
  expect_equal(top_rows$pct_liana, c(0.25, 0.25))
  expect_identical(top_rows$sender_node, c("CD4 T", "CD4 T"))
  expect_identical(top_rows$receiver_node, c("B", "NK"))
})

test_that("build_communication_consensus and plot_communication_consensus handle empty overlap or 0 consensus gracefully", {
  .consensus_load_views()
  r1 <- .consensus_make_score_result(
    scores = c(0.9, 0.4),
    senders = c("A", "B"),
    receivers = c("B", "A"),
    ligands = c("L1", "L2"),
    receptors = c("R1", "R2")
  )
  r2 <- .consensus_make_rank_result(
    ranks = c(1, 2),
    senders = c("C", "D"),
    receivers = c("D", "C"),
    ligands = c("L9", "L10"),
    receptors = c("R9", "R10")
  )

  res_disjoint <- build_communication_consensus(list(cellchat = r1, liana = r2), top_pct = 0.5)
  expect_identical(res_disjoint$qc$n_total_pairs, 4L)
  expect_identical(res_disjoint$qc$n_consensus, 0L)
  expect_true(all(!res_disjoint$table$is_consensus))
  expect_true(all(res_disjoint$table$n_methods_present == 1L))

  # Empty state plot when 0 consensus pairs
  p_empty <- plot_communication_consensus(res_disjoint)
  expect_s3_class(p_empty, "ggplot")
  expect_match(p_empty$labels$subtitle, res_disjoint$epistemic_note, fixed = TRUE)
})

test_that("plot_communication_consensus renders top consensus interactions and validates input", {
  .consensus_load_views()
  expect_error(
    plot_communication_consensus(list(bad = TRUE)),
    class = "communication_context_error"
  )

  r_cc <- .consensus_make_score_result(
    scores = c(0.95, 0.85, 0.75, 0.40, 0.10),
    senders = c("CD4 T", "B", "NK", "CD4 T", "B"),
    receivers = c("B", "CD8 T", "B", "NK", "NK"),
    ligands = c("L1", "L2", "L3", "L4", "L5"),
    receptors = c("R1", "R2", "R3", "R4", "R5")
  )
  r_li <- .consensus_make_rank_result(
    ranks = c(1, 2, 3, 4, 5),
    senders = c("CD4 T", "B", "NK", "CD4 T", "B"),
    receivers = c("B", "CD8 T", "B", "NK", "NK"),
    ligands = c("L1", "L2", "L3", "L4", "L5"),
    receptors = c("R1", "R2", "R3", "R4", "R5")
  )

  cons <- build_communication_consensus(
    list(cellchat = r_cc, liana = r_li),
    top_pct = 0.60,
    min_methods = 2L
  )
  expect_identical(cons$qc$n_consensus, 3L)

  p <- plot_communication_consensus(cons, top_n = 2L)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 2L)
  expect_match(p$labels$subtitle, cons$epistemic_note, fixed = TRUE)
})
