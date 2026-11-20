# =============================================================================
# test-ccc-lot2-ui.R — W1.3 Lot 2 §11 CCC UI & readability enhancements
# =============================================================================

test_that("plot_communication_roles_dot supports hide_minor and top_n_pathways", {
  repo_root <- if (file.exists("tools/pkg_skeleton/R/sc_communication_views.R")) "." else "../.."
  if (!exists("plot_communication_roles_dot", mode = "function")) {
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication_views.R"), local = FALSE)
  }

  roles <- data.frame(
    pathway = c("P1", "P1", "P2", "P2", "P3", "P3"),
    group = c("A", "B", "A", "B", "A", "B"),
    outdeg = c(0.8, 0.0, 0.4, 0.0, 0.1, 0.0),
    indeg = c(0.0, 0.6, 0.0, 0.2, 0.0, 0.0),
    role = c("sender", "receiver", "sender", "minor", "minor", "minor"),
    stringsAsFactors = FALSE
  )

  p_filtered <- plot_communication_roles_dot(roles, hide_minor = TRUE, top_n_pathways = 1L)
  expect_s3_class(p_filtered, "ggplot")
  expect_false("minor" %in% as.character(p_filtered$data$role))
  expect_equal(length(unique(as.character(p_filtered$data$pathway))), 1L)
  expect_equal(unique(as.character(p_filtered$data$pathway)), "P1")
  expect_equal(p_filtered$theme$axis.text.x$angle, 35)
})

test_that("plot_communication_pathway_heatmap supports top_n_pathways", {
  repo_root <- if (file.exists("tools/pkg_skeleton/R/sc_communication_views.R")) "." else "../.."
  if (!exists("finalize_communication_result", mode = "function")) {
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication.R"), local = FALSE)
  }
  if (!exists("plot_communication_pathway_heatmap", mode = "function")) {
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication_views.R"), local = FALSE)
  }

  tbl <- data.frame(
    sender = c("A", "A", "B", "B", "C", "C"),
    receiver = c("B", "C", "A", "C", "A", "B"),
    ligand = c("L1", "L2", "L3", "L4", "L5", "L6"),
    receptor = c("R1", "R2", "R3", "R4", "R5", "R6"),
    interaction = paste0("L", 1:6, "_", "R", 1:6),
    pathway = c("PW1", "PW1", "PW2", "PW2", "PW3", "PW3"),
    score = c(0.9, 0.8, 0.5, 0.4, 0.1, 0.05),
    p_value = rep(0.01, 6),
    p_adjusted = rep(0.02, 6),
    source_method = rep("cellchat", 6),
    source_file = rep("mock.csv", 6),
    source_cell_identity_level = rep("celltype", 6),
    stringsAsFactors = FALSE
  )
  comm_res <- finalize_communication_result(
    canonical_table = tbl,
    source_method = "cellchat"
  )

  p_top2 <- plot_communication_pathway_heatmap(comm_res, top_n_pathways = 2L)
  expect_s3_class(p_top2, "ggplot")
  kept <- unique(as.character(p_top2$data$pathway))
  expect_equal(length(kept), 2L)
  expect_setequal(kept, c("PW1", "PW2"))
})

test_that("mod_sc_communication wires Lot 2 controls and inactive filter disable", {
  repo_root <- if (file.exists("modules/sc/mod_sc_communication.R")) "." else "../.."
  src <- paste(readLines(file.path(repo_root, "modules/sc/mod_sc_communication.R"), warn = FALSE),
               collapse = "\n")

  expect_match(src, 'comm_centrality_hide_minor', fixed = TRUE)
  expect_match(src, 'comm_centrality_top_n', fixed = TRUE)
  expect_match(src, 'comm_heatmap_top_n', fixed = TRUE)
  expect_match(src, 'comm_dotplot_top_n', fixed = TRUE)
  expect_match(src, 'top_n_pairs', fixed = TRUE)
  expect_match(src, 'comm_circle_top_n', fixed = TRUE)
  expect_match(src, 'top_n_edges', fixed = TRUE)
  expect_match(src, 'comm_centrality_quadrant_plot', fixed = TRUE)
  expect_match(src, 'plot_communication_centrality_quadrant', fixed = TRUE)
  expect_lt(
    regexpr('DT::dataTableOutput(ns("comm_centrality")', src, fixed = TRUE)[1],
    regexpr('plotOutput(ns("comm_centrality_quadrant_plot")', src, fixed = TRUE)[1]
  )
  expect_match(src, '.toggle_filter', fixed = TRUE)
  expect_match(src, 'shinyjs::disable', fixed = TRUE)
})
