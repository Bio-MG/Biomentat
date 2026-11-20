# =============================================================================
# test-sc-velocity-trajectory-wiring.R — W1.1 items #10 & #12
# =============================================================================
# VELOCITY_IMPROVEMENTS_PROPOSAL_2026-11-17.md:
#   - Item #10: move the velocity log box (`velocity_log`) into the workflow
#     sidebar (`mod_sc_velocity_ui`) directly under the scVelo opt-in button
#     with bounded scrollable height, instead of inside the results card
#     (`mod_sc_velocity_output_ui`).
#   - Item #12: wire auto-pipeline / metadata-persisted trajectory results
#     (`obj@meta.data$pseudotime` + `traj_*`) into `mod_sc_trajectory_server`
#     so "Plot Trajectoire" (`output$trajectory_plot`) and its export render
#     without requiring a manual click on "Calculer Trajectoire".
# =============================================================================

source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/palettes.R")
source_project_file("tools/pkg_skeleton/R/error_log.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("tools/pkg_skeleton/R/sc_trajectory.R")
source_project_file("modules/sc/mod_sc_velocity.R")
source_project_file("modules/sc/mod_sc_trajectory.R")

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
})

# Minimal i18n / .tr_plain shim for standalone UI rendering tests
.w11_with_ui_shims <- function(expr) {
  had_i18n <- exists("i18n", envir = globalenv(), inherits = FALSE)
  old_i18n <- if (had_i18n) get("i18n", envir = globalenv(), inherits = FALSE) else NULL
  had_trp  <- exists(".tr_plain", envir = globalenv(), inherits = FALSE)
  old_trp  <- if (had_trp) get(".tr_plain", envir = globalenv(), inherits = FALSE) else NULL

  assign("i18n", list(t = function(x) x), envir = globalenv())
  assign(".tr_plain", function(x) x, envir = globalenv())
  on.exit({
    if (had_i18n) assign("i18n", old_i18n, envir = globalenv())
    else if (exists("i18n", envir = globalenv(), inherits = FALSE)) rm("i18n", envir = globalenv())
    if (had_trp) assign(".tr_plain", old_trp, envir = globalenv())
    else if (exists(".tr_plain", envir = globalenv(), inherits = FALSE)) rm(".tr_plain", envir = globalenv())
  }, add = TRUE)
  force(expr)
}

# Build a minimal Seurat object mimicking SC auto-pipeline output
.w11_make_seurat <- function(with_pseudotime = TRUE, reductions = c("pca", "umap"),
                             traj_red = "umap", root_cell = 3L,
                             traj_method = "exploratory_knn") {
  set.seed(42)
  n_genes <- 20L
  n_cells <- 16L
  mat <- matrix(rpois(n_genes * n_cells, lambda = 4), nrow = n_genes, ncol = n_cells,
                dimnames = list(paste0("G", seq_len(n_genes)),
                                paste0("C", seq_len(n_cells))))
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = mat))
  seu@meta.data$seurat_clusters <- factor(rep(c("0", "1"), length.out = n_cells))

  if ("pca" %in% reductions) {
    pca_emb <- matrix(rnorm(n_cells * 2L), nrow = n_cells, ncol = 2L,
                      dimnames = list(colnames(seu), c("PC_1", "PC_2")))
    seu[["pca"]] <- Seurat::CreateDimReducObject(embeddings = pca_emb, key = "PC_", assay = "RNA")
  }
  if ("umap" %in% reductions) {
    umap_emb <- matrix(runif(n_cells * 2L), nrow = n_cells, ncol = 2L,
                       dimnames = list(colnames(seu), c("UMAP_1", "UMAP_2")))
    seu[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap_emb, key = "UMAP_", assay = "RNA")
  }

  if (isTRUE(with_pseudotime)) {
    pt <- seq(0, 1, length.out = n_cells)
    seu@meta.data$pseudotime                 <- pt
    seu@meta.data$traj_in_component          <- rep(TRUE, n_cells)
    seu@meta.data$traj_method                <- rep(traj_method, n_cells)
    seu@meta.data$traj_computation_reduction <- rep(traj_red, n_cells)
    seu@meta.data$traj_root_method           <- rep("diameter", n_cells)
    seu@meta.data$traj_root_cell             <- rep(root_cell, n_cells)
    seu@meta.data$traj_root_cluster          <- rep(NA_character_, n_cells)
    seu@meta.data$traj_root_component_size   <- rep(n_cells, n_cells)
  }
  seu
}

# =============================================================================
# Item #10 — Velocity log box ergonomics (sidebar under opt-in button)
# =============================================================================

test_that("Item #10: velocity_log lives in sidebar under velocity_scvelo_run, not in output card", {
  skip_if_not_installed("shinyWidgets")

  ui_rendered <- .w11_with_ui_shims({
    list(
      sidebar = as.character(mod_sc_velocity_ui("vel")),
      output  = as.character(mod_sc_velocity_output_ui("vel"))
    )
  })

  # Sidebar contains the log box with bounded scrollable height
  expect_match(ui_rendered$sidebar, 'id="vel-velocity_log"', fixed = TRUE)
  expect_match(ui_rendered$sidebar, "max-height:140px;overflow-y:auto;", fixed = TRUE)

  # Positioned AFTER the scVelo opt-in button and BEFORE the phase portrait gene selector
  pos_optin <- regexpr('id="vel-velocity_scvelo_run"', ui_rendered$sidebar, fixed = TRUE)[1]
  pos_log   <- regexpr('id="vel-velocity_log"', ui_rendered$sidebar, fixed = TRUE)[1]
  pos_gene  <- regexpr('id="vel-velocity_gene"', ui_rendered$sidebar, fixed = TRUE)[1]
  expect_true(pos_optin > 0L && pos_log > pos_optin && pos_gene > pos_log)

  # Output results card no longer holds velocity_log, but keeps drawer + 4 tabs
  expect_false(grepl('id="vel-velocity_log"', ui_rendered$output, fixed = TRUE))
  expect_match(ui_rendered$output, 'id="vel-velocity_drawer"', fixed = TRUE)
  expect_match(ui_rendered$output, 'id="vel-velocity_out_tabs"', fixed = TRUE)
})

# =============================================================================
# Item #12 — Auto-pipeline trajectory -> Plot Trajectoire wiring
# =============================================================================

test_that("Item #12: .traj_result_from_meta hydrates metadata provenance safely", {
  skip_if_not_installed("Seurat")

  expect_null(.traj_result_from_meta(NULL))
  seu_no_pt <- .w11_make_seurat(with_pseudotime = FALSE)
  expect_null(.traj_result_from_meta(seu_no_pt))

  seu_ap <- .w11_make_seurat(with_pseudotime = TRUE, traj_red = "umap", root_cell = 5L)
  res <- .traj_result_from_meta(seu_ap)
  expect_true(is.list(res))
  expect_identical(res$source, "metadata")
  expect_identical(res$method_type, "exploratory_knn")
  expect_identical(res$computation_reduction, "umap")
  expect_identical(res$root_cell, 5L)
  expect_null(res$graph)
  expect_identical(res$n_cells, as.integer(ncol(seu_ap)))
  expect_identical(res$root_component_size, as.integer(ncol(seu_ap)))
  expect_equal(res$pseudotime, seu_ap@meta.data$pseudotime)

  # Out-of-range or NA root_cell must sanitize to NULL (never crash plot_trajectory)
  seu_bad_root <- .w11_make_seurat(with_pseudotime = TRUE, root_cell = 9999L)
  expect_null(.traj_result_from_meta(seu_bad_root)$root_cell)
  seu_na_root <- .w11_make_seurat(with_pseudotime = TRUE, root_cell = NA_integer_)
  expect_null(.traj_result_from_meta(seu_na_root)$root_cell)
})

test_that("Item #12: .resolve_traj_result preserves matching manual result and re-hydrates stale/NULL", {
  skip_if_not_installed("Seurat")

  seu_ap <- .w11_make_seurat(with_pseudotime = TRUE, traj_red = "umap", root_cell = 4L)
  # When current_res is NULL (auto-pipeline case), hydrates from meta.data
  hydrated <- .resolve_traj_result(seu_ap, NULL)
  expect_identical(hydrated$source, "metadata")
  expect_identical(hydrated$root_cell, 4L)

  # When current_res matches (manual run with graph), keeps manual result intact
  manual_res <- list(
    pseudotime = seu_ap@meta.data$pseudotime,
    graph = structure(list(mock = TRUE), class = "mock_graph"),
    root_cell = 4L,
    n_cells = ncol(seu_ap),
    method_type = "exploratory_knn"
  )
  kept <- .resolve_traj_result(seu_ap, manual_res)
  expect_identical(kept$graph, manual_res$graph)
  expect_null(kept$source)

  # When pseudotime changes (new pipeline run), stale manual_res is replaced
  seu_ap2 <- seu_ap
  seu_ap2@meta.data$pseudotime <- rev(seu_ap@meta.data$pseudotime)
  replaced <- .resolve_traj_result(seu_ap2, manual_res)
  expect_identical(replaced$source, "metadata")
  expect_equal(replaced$pseudotime, seu_ap2@meta.data$pseudotime)

  # When new object lacks pseudotime, returns NULL
  seu_no_pt <- .w11_make_seurat(with_pseudotime = FALSE)
  expect_null(.resolve_traj_result(seu_no_pt, manual_res))
})

test_that("Item #12: .resolve_traj_display_reduction picks valid reduction with clean fallback", {
  skip_if_not_installed("Seurat")

  seu_both <- .w11_make_seurat(with_pseudotime = TRUE, reductions = c("pca", "umap"))
  expect_identical(.resolve_traj_display_reduction(seu_both, "umap", NULL, "pca"), "umap")
  expect_identical(.resolve_traj_display_reduction(seu_both, "pca", NULL, "umap"), "pca")
  expect_identical(.resolve_traj_display_reduction(seu_both, NULL, list(computation_reduction = "umap"), NULL), "umap")

  seu_umap_only <- .w11_make_seurat(with_pseudotime = TRUE, reductions = "umap")
  expect_identical(.resolve_traj_display_reduction(seu_umap_only, "pca", list(computation_reduction = "umap"), "umap"), "umap")
})

test_that("Item #12: testServer renders Plot Trajectoire directly after auto-pipeline without manual calc", {
  skip_if_not_installed("Seurat")

  seu_ap <- .w11_make_seurat(with_pseudotime = TRUE, traj_red = "umap", root_cell = 3L)
  seu_plain <- .w11_make_seurat(with_pseudotime = FALSE)

  gd <- shiny::reactiveValues(sc_obj = NULL, i18n = NULL, language = "fr")
  sr <- shiny::reactiveValues(
    traj_reduction     = NULL,
    traj_method        = NULL,
    traj_genes         = character(0),
    active_tab         = NULL,
    sc_palette         = "default",
    sc_manual_gradient = NULL
  )

  testServer(mod_sc_trajectory_server,
             args = list(id = "trajectory", global_data = gd, shared_rv = sr), {
    # Simulate SC auto-pipeline finishing (sets sc_obj + shared_rv$traj_reduction/method,
    # NEVER clicks input$calc_trajectory)
    gd$sc_obj <- seu_ap
    sr$traj_reduction <- "umap"
    sr$traj_method <- "exploratory_knn"
    session$setInputs(traj_reduction = "umap", traj_color = "pseudotime", traj_show_edges = FALSE)
    session$flushReact()

    # Local traj_result is hydrated from auto-pipeline metadata
    res <- traj_result()
    expect_false(is.null(res))
    expect_identical(res$source, "metadata")
    expect_identical(res$root_cell, 3L)
    expect_identical(res$computation_reduction, "umap")
    expect_match(traj_status_rv(), "OK Pseudotemps", fixed = TRUE)

    # Plot Trajectoire renders cleanly (non-NULL PNG payload from renderPlot)
    plt_out <- output$trajectory_plot
    expect_false(is.null(plt_out))
    expect_true(is.list(plt_out) && !is.null(plt_out$src))

    # Replacing sc_obj with a dataset lacking pseudotime clears traj_result
    gd$sc_obj <- seu_plain
    session$flushReact()
    expect_null(traj_result())
  })
})
