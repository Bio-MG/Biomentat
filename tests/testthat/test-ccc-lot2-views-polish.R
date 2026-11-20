# =============================================================================
# test-ccc-lot2-views-polish.R — W1.3 Lot 2 §11 heatmap role/top_n + circle self-loops
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

.polish_result <- function() {
  .lot2_load_comm_views()
  res <- .lot2_make_comm_result()
  tb <- res$canonical_table
  tb$receiver[5L] <- tb$sender[5L]   # D -> D autocrine
  res$canonical_table <- tb
  res
}

test_that("heatmap role outgoing/incoming aggregates with Score cumul\u00e9 legend", {
  r <- .polish_result()
  p_both <- plot_communication_pathway_heatmap(r)
  expect_false(identical(p_both$scales$get_scales("fill")$name, "Score cumul\u00e9"))
  p_out <- plot_communication_pathway_heatmap(r, role = "outgoing")
  expect_identical(p_out$scales$get_scales("fill")$name, "Score cumul\u00e9")
  expect_setequal(unique(p_out$data$pair_label), c("A", "B", "C", "D"))
  a_p1 <- p_out$data[p_out$data$pair_label == "A" & p_out$data$pathway == "P1", ]
  expect_equal(a_p1$mean_score, 1.1)
  p_in <- plot_communication_pathway_heatmap(r, role = "incoming")
  expect_setequal(unique(p_in$data$pair_label), c("B", "C", "D"))
  expect_error(plot_communication_pathway_heatmap(r, role = "bogus"))
})

test_that("heatmap top_n caps pathways and aliases top_n_pathways", {
  r <- .polish_result()
  p <- plot_communication_pathway_heatmap(r, top_n = 2L)
  expect_lte(length(unique(p$data$pathway)), 2L)
  expect_true("P1" %in% p$data$pathway)
  p2 <- plot_communication_pathway_heatmap(r, top_n_pathways = 1L, top_n = 3L)
  expect_identical(length(unique(p2$data$pathway)), 1L)
  expect_gt(length(unique(plot_communication_pathway_heatmap(r)$data$pathway)), 2L)
})

test_that("circle draws autocrine self-loops with sender-colored edges", {
  r <- .polish_result()
  p <- plot_communication_circle(r)
  expect_s3_class(p, "ggplot")
  # layer 1 stays the non-self edge layer (frozen contract).
  expect_equal(length(unique(p$layers[[1L]]$data$edge_id)), 3L)
  expect_true("sender" %in% names(p$layers[[1L]]$data))
  loops <- p$layers[[2L]]$data
  expect_identical(unique(loops$sender), "D")
  expect_identical(length(unique(loops$loop_id)), 1L)
  expect_match(p$labels$subtitle, "1 auto-interaction(s) dessinee(s) en boucle", fixed = TRUE)
  expect_identical(p$scales$get_scales("colour")$name, "Sender")
})

test_that("circle with only self-interactions draws loops, not an empty plot", {
  r <- .polish_result()
  tb <- r$canonical_table
  tb$receiver <- tb$sender
  r$canonical_table <- tb
  p <- plot_communication_circle(r)
  expect_s3_class(p, "ggplot")
  expect_identical(nrow(p$layers[[1L]]$data), 0L)
  expect_identical(length(unique(p$layers[[2L]]$data$loop_id)), 4L)
})

test_that("public api names unchanged by polish", {
  .lot2_load_comm_views()
  api <- communication_views_public_api()
  expect_true(all(c("plot_communication_pathway_heatmap", "plot_communication_circle") %in% api))
})
