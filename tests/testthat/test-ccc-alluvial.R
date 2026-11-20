# =============================================================================
# test-ccc-alluvial.R — W2.2 §11.8 Gem #4 : vue alluviale 3 tiers
# (Sender -> Pathway/Ligand -> Receiver) en ggplot2 pur
# =============================================================================

.alluvial_make_score_result <- function(n_paths = 4L) {
  senders <- c("CD4 T", "B", "NK", "Monocyte")
  receivers <- c("B", "CD8 T", "CD4 T", "NK")
  rows <- lapply(seq_len(n_paths), function(i) {
    data.frame(
      sender = senders[((i - 1L) %% length(senders)) + 1L],
      receiver = receivers[((i - 1L) %% length(receivers)) + 1L],
      ligand = paste0("LIG", i),
      receptor = paste0("REC", i),
      interaction = paste0("LIG", i, " -> REC", i),
      pathway = paste0("PW", i),
      score = (n_paths - i + 1L) * 0.1,
      p_value = 0.01,
      p_adjusted = 0.02,
      source_method = "cellchat",
      source_file = "mock_alluvial.csv",
      source_cell_identity_level = "celltype",
      stringsAsFactors = FALSE
    )
  })
  tbl <- do.call(rbind, rows)
  finalize_communication_result(canonical_table = tbl, source_method = "cellchat")
}

.alluvial_make_rank_result <- function() {
  tbl <- data.frame(
    sender = c("CD4 T", "B", "NK"),
    receiver = c("B", "CD8 T", "CD4 T"),
    ligand = c("IL7", "CD40LG", "CCL5"),
    receptor = c("IL7R", "CD40", "CCR5"),
    interaction = c("IL7 -> IL7R", "CD40LG -> CD40", "CCL5 -> CCR5"),
    pathway = c("IL7 signaling", "CD40 signaling", "CCL signaling"),
    score = c(NA_real_, NA_real_, NA_real_),
    p_value = c(0.01, 0.02, 0.03),
    p_adjusted = c(NA_real_, NA_real_, NA_real_),
    source_method = rep("liana", 3L),
    source_file = rep("mock_liana.csv", 3L),
    source_cell_identity_level = rep("celltype", 3L),
    rank = c(1, 2, 3),
    rank_direction = rep("lower_is_better", 3L),
    rank_aggregation_mode = rep("specificity", 3L),
    stringsAsFactors = FALSE
  )
  finalize_communication_result(
    canonical_table = tbl,
    source_method = "liana",
    external_consensus = TRUE
  )
}

test_that("plot_communication_alluvial is registered in communication_views_public_api()", {
  expect_true("plot_communication_alluvial" %in% communication_views_public_api())
})

test_that("plot_communication_alluvial validates invalid inputs with communication_import_error", {
  expect_error(
    plot_communication_alluvial(NULL),
    class = "communication_import_error"
  )
  expect_error(
    plot_communication_alluvial(list(type = "wrong")),
    class = "communication_import_error"
  )
  expect_error(
    plot_communication_alluvial(.alluvial_make_score_result(), middle_by = "invalid"),
    regexp = "pathway|ligand"
  )
})

test_that("plot_communication_alluvial renders pathway and ligand middle tiers in pure ggplot2", {
  res <- .alluvial_make_score_result(n_paths = 4L)

  p_pw <- plot_communication_alluvial(res, middle_by = "pathway")
  expect_s3_class(p_pw, "ggplot")
  expect_true(length(p_pw$layers) >= 3L)
  expect_match(p_pw$labels$subtitle, "cellchat", ignore.case = TRUE)
  expect_match(p_pw$labels$subtitle, "aucune causalite", ignore.case = TRUE)

  p_lig <- plot_communication_alluvial(res, middle_by = "ligand")
  expect_s3_class(p_lig, "ggplot")
  expect_true(length(p_lig$layers) >= 3L)

  # Extract node labels to confirm pathway vs ligand nodes are plotted
  pw_labels <- p_pw$layers[[3]]$data$node
  lig_labels <- p_lig$layers[[3]]$data$node
  expect_true(any(grepl("^PW1", pw_labels)))
  expect_true(any(grepl("^LIG1", lig_labels)))
  expect_false(any(grepl("^LIG1", pw_labels)))
})

test_that("plot_communication_alluvial enforces top_n_pathways and top_n_flows guardrails", {
  res <- .alluvial_make_score_result(n_paths = 12L)

  p_top_pw <- plot_communication_alluvial(
    res,
    top_n_pathways = 3L,
    middle_by = "pathway",
    top_n_flows = 30L
  )
  expect_s3_class(p_top_pw, "ggplot")
  pw_labels <- p_top_pw$layers[[3]]$data$node
  # Top 3 by score are PW1, PW2, PW3 (scores 1.2, 1.1, 1.0)
  expect_true(all(c("PW1", "PW2", "PW3") %in% pw_labels))
  expect_false("PW12" %in% pw_labels)

  p_top_flows <- plot_communication_alluvial(
    res,
    top_n_pathways = 12L,
    middle_by = "pathway",
    top_n_flows = 2L
  )
  expect_s3_class(p_top_flows, "ggplot")
  flow_labels <- p_top_flows$layers[[3]]$data$node
  expect_true(all(c("PW1", "PW2") %in% flow_labels))
  expect_false("PW5" %in% flow_labels)
})

test_that("plot_communication_alluvial returns empty view plot on 0-row or zero-score inputs", {
  res <- .alluvial_make_score_result(n_paths = 3L)
  empty_tbl <- res$canonical_table[0, , drop = FALSE]

  p_empty <- plot_communication_alluvial(res, filtered_table = empty_tbl)
  expect_s3_class(p_empty, "ggplot")
  expect_match(p_empty$labels$subtitle, "Aucun", ignore.case = TRUE)
  expect_match(p_empty$labels$subtitle, "cellchat", ignore.case = TRUE)

  # Zero-score rows should also be filtered out and fall back to empty view plot
  zero_tbl <- res$canonical_table
  zero_tbl$score <- 0
  p_zero <- plot_communication_alluvial(res, filtered_table = zero_tbl)
  expect_s3_class(p_zero, "ggplot")
  expect_match(p_zero$labels$subtitle, "Aucun", ignore.case = TRUE)

  # Missing pathway values when middle_by = "pathway" should fall back cleanly
  no_pw_tbl <- res$canonical_table
  no_pw_tbl$pathway <- NA_character_
  p_no_pw <- plot_communication_alluvial(res, filtered_table = no_pw_tbl, middle_by = "pathway")
  expect_s3_class(p_no_pw, "ggplot")
  expect_match(p_no_pw$labels$subtitle, "Aucun", ignore.case = TRUE)
})

test_that("plot_communication_alluvial supports score-less rank-only input with epistemic subtitle", {
  res_rank <- .alluvial_make_rank_result()
  p_rank <- plot_communication_alluvial(res_rank, middle_by = "pathway")
  expect_s3_class(p_rank, "ggplot")
  expect_true(length(p_rank$layers) >= 3L)
  expect_match(p_rank$labels$subtitle, "liana", ignore.case = TRUE)
  expect_match(p_rank$labels$subtitle, "rang importe", ignore.case = TRUE)
  expect_match(p_rank$labels$subtitle, "aucune causalite", ignore.case = TRUE)
})

test_that("plot_communication_alluvial harmonizes node names via sender_mapped/receiver_mapped", {
  res <- .alluvial_make_score_result(n_paths = 2L)
  tbl <- res$canonical_table
  tbl$sender_mapped <- c("Helper_T", "B_Cell")
  tbl$receiver_mapped <- c("B_Cell", "Cytotoxic_T")
  res_mapped <- finalize_communication_result(canonical_table = tbl, source_method = "cellchat")

  p <- plot_communication_alluvial(res_mapped, middle_by = "pathway")
  expect_s3_class(p, "ggplot")
  lbls <- p$layers[[3]]$data$node
  expect_true(all(c("Helper_T", "B_Cell", "Cytotoxic_T") %in% lbls))
})
