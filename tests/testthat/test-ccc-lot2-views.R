# =============================================================================
# test-ccc-lot2-views.R — W1.3 Lot 2 §11 view-layer polish tests
# =============================================================================

.lot2_load_comm_views <- function() {
  repo_root <- if (file.exists("tools/pkg_skeleton/R/sc_communication_views.R")) "." else "../.."
  if (!exists("finalize_communication_result", mode = "function")) {
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication.R"), local = FALSE)
  }
  if (!exists("plot_communication_dotplot", mode = "function")) {
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication_views.R"), local = FALSE)
  }
}

.lot2_make_comm_result <- function(with_rank = FALSE) {
  src <- if (with_rank) "liana" else "cellchat"
  tbl <- data.frame(
    sender = c("A", "A", "B", "C", "D"),
    receiver = c("B", "B", "C", "D", "A"),
    ligand = c("L1", "L2", "L3", "L4", "L5"),
    receptor = c("R1", "R2", "R3", "R4", "R5"),
    interaction = paste0("L", 1:5, "_R", 1:5),
    pathway = c("P1", "P1", "P2", "P3", "P4"),
    score = if (with_rank) rep(NA_real_, 5L) else c(0.6, 0.5, 0.8, 0.2, 0.1),
    p_value = c(0.01, 0.02, 0.01, 0.04, 0.05),
    p_adjusted = c(0.02, 0.03, 0.02, 0.05, 0.06),
    source_method = rep(src, 5L),
    source_file = rep("mock.csv", 5L),
    source_cell_identity_level = rep("celltype", 5L),
    stringsAsFactors = FALSE
  )
  if (with_rank) {
    tbl$rank <- c(1, 2, 1, 5, 10)
  }
  finalize_communication_result(
    canonical_table = tbl,
    source_method = src
  )
}

test_that("plot_communication_dotplot supports top_n_pairs guardrail", {
  .lot2_load_comm_views()
  r <- .lot2_make_comm_result(with_rank = FALSE)

  p_all <- plot_communication_dotplot(r)
  expect_s3_class(p_all, "ggplot")
  expect_equal(nrow(p_all$data), 4L)

  p_top2 <- plot_communication_dotplot(r, top_n_pairs = 2L)
  expect_s3_class(p_top2, "ggplot")
  expect_equal(nrow(p_top2$data), 2L)
  # Pairs by sum_score: A->B (1.1), B->C (0.8), C->D (0.2), D->A (0.1)
  pairs_kept <- paste(p_top2$data$sender_node, p_top2$data$receiver_node, sep = "->")
  expect_setequal(pairs_kept, c("A->B", "B->C"))

  # Non-positive / non-finite values are no-ops
  expect_equal(nrow(plot_communication_dotplot(r, top_n_pairs = 0L)$data), 4L)
  expect_equal(nrow(plot_communication_dotplot(r, top_n_pairs = NA_real_)$data), 4L)
  expect_equal(nrow(plot_communication_dotplot(r, top_n_pairs = Inf)$data), 4L)

  # Rank-based result falls back to interaction count
  r_rk <- .lot2_make_comm_result(with_rank = TRUE)
  p_rk_top1 <- plot_communication_dotplot(r_rk, top_n_pairs = 1L)
  expect_equal(nrow(p_rk_top1$data), 1L)
  expect_equal(p_rk_top1$data$sender_node, "A")
  expect_equal(p_rk_top1$data$receiver_node, "B")
})

test_that("plot_communication_circle supports top_n_edges guardrail", {
  .lot2_load_comm_views()
  r <- .lot2_make_comm_result(with_rank = FALSE)

  pc_all <- plot_communication_circle(r)
  expect_s3_class(pc_all, "ggplot")
  expect_equal(length(unique(pc_all$layers[[1L]]$data$edge_id)), 4L)

  pc_top2 <- plot_communication_circle(r, top_n_edges = 2L)
  expect_s3_class(pc_top2, "ggplot")
  edge_df <- pc_top2$layers[[1L]]$data
  expect_equal(length(unique(edge_df$edge_id)), 2L)
  expect_true(all(unique(edge_df$weight) %in% c(1.1, 0.8)))
  expect_match(pc_top2$labels$subtitle, "2 arete\\(s\\)")

  # Non-positive / non-finite values are no-ops
  expect_equal(
    length(unique(plot_communication_circle(r, top_n_edges = -2L)$layers[[1L]]$data$edge_id)),
    4L
  )
  expect_equal(
    length(unique(plot_communication_circle(r, top_n_edges = Inf)$layers[[1L]]$data$edge_id)),
    4L
  )
})

test_that("plot_communication_centrality_quadrant validates, aggregates, and plots", {
  .lot2_load_comm_views()

  expect_error(
    plot_communication_centrality_quadrant(data.frame(bad = 1)),
    class = "communication_import_error"
  )
  expect_true(
    "plot_communication_centrality_quadrant" %in% communication_views_public_api()
  )

  roles <- data.frame(
    pathway = c("P1", "P1", "P1", "P2", "P2", "P2"),
    group = c("A", "B", "C", "A", "B", "C"),
    outdeg = c(10, 2, 5, 6, 1, 3),
    indeg = c(1, 12, 6, 2, 8, 4),
    role = c("sender", "receiver", "mediator", "sender", "receiver", "mediator"),
    stringsAsFactors = FALSE
  )

  # Aggregated across all pathways
  p_agg <- plot_communication_centrality_quadrant(roles)
  expect_s3_class(p_agg, "ggplot")
  expect_equal(nrow(p_agg$data), 3L)
  a_row <- p_agg$data[p_agg$data$group == "A", ]
  expect_equal(a_row$outdeg, 16)
  expect_equal(a_row$indeg, 3)

  layer_classes <- vapply(p_agg$layers, function(l) class(l$geom)[1L], character(1L))
  expect_true("GeomVline" %in% layer_classes)
  expect_true("GeomHline" %in% layer_classes)
  expect_true("GeomPoint" %in% layer_classes)
  expect_true("GeomText" %in% layer_classes)

  # Single pathway filter
  p_p1 <- plot_communication_centrality_quadrant(roles, pathway = "P1")
  expect_s3_class(p_p1, "ggplot")
  expect_equal(nrow(p_p1$data), 3L)
  expect_equal(p_p1$data$outdeg[p_p1$data$group == "A"], 10)
  expect_match(p_p1$labels$subtitle, "P1", fixed = TRUE)

  # Empty / unknown pathway returns graceful empty plot
  p_empty <- plot_communication_centrality_quadrant(roles, pathway = "UNKNOWN")
  expect_s3_class(p_empty, "ggplot")
})
