# =============================================================================
# test-mod-sc-communication-spatial-2d-ui.R — W1.3 Lot 1 : Spatial 2D UI polish
# =============================================================================
# Self-sourcing pour isolation par shard :
# - grisement de spc_radius / spc_nperm en mode 2D (avec garde testServer)
# - masquage des 4 colonnes 100%-NA dans spc_table
# - affichage de l'étendue (extent) et des unités des coordonnées dans l'UI
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_communication_spatial.R")

.spc_2d_ui_env <- function() {
  skip_if_not_installed("shiny")
  skip_if_not_installed("DT")
  env <- new.env(parent = asNamespace("shiny"))
  env$i18n <- list(t = function(x) x)
  env$.tr_plain <- function(x) x
  env$.strip_i18n_html <- function(x) x
  env$`%||%` <- function(a, b) if (is.null(a)) b else a
  env$build_communication_spatial_context <- build_communication_spatial_context
  env$plot_communication_spatial_edges <- plot_communication_spatial_edges
  env$plot_communication_spatial_distance_summary <- plot_communication_spatial_distance_summary
  env$build_communication_spatial_export <- build_communication_spatial_export
  env$velocity_object_fingerprint <- function(obj) "fp_test"
  env$state_get <- function(rv, key) {
    if (is.null(rv)) return(NULL)
    rv[[key]]
  }
  env$provenance_append <- function(rv, prov) invisible(NULL)
  env$showNotification <- function(...) invisible(NULL)
  env$ts_datatable <- function(data, ...) {
    attr(data, "ts_dt_colnames") <- colnames(data)
    DT::datatable(data)
  }
  env$ts_export_plot <- function(...) invisible(NULL)
  sys.source(
    file.path(ts_project_root(), "modules", "sc", "mod_sc_communication_spatial.R"),
    envir = env
  )
  env
}

.spc_2d_mock_sc_obj <- function() {
  skip_if_not_installed("SeuratObject")
  coords <- .ccx_coords()
  ids <- .ccx_identities()
  counts <- matrix(
    rpois(length(ids) * 5L, lambda = 3),
    nrow = 5L,
    ncol = length(ids),
    dimnames = list(paste0("G", 1:5), names(ids))
  )
  obj <- SeuratObject::CreateSeuratObject(counts = counts)
  obj$cell_type <- unname(ids[colnames(obj)])
  emb <- as.matrix(coords[match(colnames(obj), rownames(coords)), c("x", "y")])
  rownames(emb) <- colnames(obj)
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_")
  obj
}

test_that("UI spatial expose les contrôles 2D, le statut et la table DT", {
  env <- .spc_2d_ui_env()
  ui <- env$mod_sc_communication_spatial_ui("spc")
  html <- as.character(shiny::tagList(ui))

  expect_match(html, 'id="spc-spc_source"', fixed = TRUE)
  expect_match(html, 'id="spc-spc_radius"', fixed = TRUE)
  expect_match(html, 'id="spc-spc_nperm"', fixed = TRUE)
  expect_match(html, 'id="spc-spc_status"', fixed = TRUE)
  expect_match(html, 'id="spc-spc_spatial_info"', fixed = TRUE)
  expect_match(html, 'id="spc-spc_table"', fixed = TRUE)
})

test_that("W1.3 : shinyjs::disable/enable est gardé pour compatibilité testServer", {
  src <- paste(
    readLines(file.path(ts_project_root(), "modules", "sc", "mod_sc_communication_spatial.R"),
              warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  expect_match(src, 'shinyjs::disable\\("spc_radius"\\)')
  expect_match(src, 'shinyjs::disable\\("spc_nperm"\\)')
  expect_match(src, 'shinyjs::enable\\("spc_radius"\\)')
  expect_match(src, 'shinyjs::enable\\("spc_nperm"\\)')
  expect_match(src, "getDefaultReactiveDomain\\(\\)|tryCatch")
})

test_that("server testServer mode 2D : étendue/unités dans spc_status et 4 colonnes 100%-NA masquées", {
  env <- .spc_2d_ui_env()
  sc_obj <- .spc_2d_mock_sc_obj()
  res <- .ccx_result()

  comm_state <- shiny::reactiveValues(
    result = res,
    object_fingerprint = "fp_test",
    context_cells = NULL,
    ctx_results = list(),
    ctx_view = list()
  )
  coords_mat <- .ccx_coords()
  coords_df <- data.frame(
    id = rownames(coords_mat),
    x = coords_mat[, "x"],
    y = coords_mat[, "y"],
    stringsAsFactors = FALSE
  )
  global_data <- shiny::reactiveValues(
    sc_obj = sc_obj,
    spatial_obj = list(coords = coords_df),
    i18n = list(t = function(x) x)
  )
  shared_rv <- shiny::reactiveValues(
    comm_tabs_current = "ccc_tab_spatial"
  )

  shiny::testServer(
    env$mod_sc_communication_spatial_server,
    args = list(comm_state = comm_state, global_data = global_data, shared_rv = shared_rv),
    {
      session$setInputs(
        spc_source = "reduction",
        spc_reduction = "umap",
        spc_radius = NA_real_,
        spc_nperm = 0L
      )
      session$flushReact()

      status_txt <- output$spc_status
      expect_true(is.character(status_txt) && nzchar(status_txt))
      # Étendue X [-2, 10] / Y [-1, 1] et mention projection 2D (non physique)
      expect_match(status_txt, "X \\[", perl = TRUE)
      expect_match(status_txt, "Y \\[", perl = TRUE)
      expect_match(status_txt, "projection 2D", fixed = TRUE)

      # Les 4 colonnes 100%-NA sont absentes de la table filtrée pour DT
      ctx <- comm_state$ctx_results$spatial
      expect_false(is.null(ctx))
      dropped <- .spc_drop_all_na_cols(ctx$pair_table)
      expect_false(any(c(
        "frac_receiver_within_radius",
        "frac_sender_within_radius",
        "z_score",
        "p_perm"
      ) %in% colnames(dropped)))
      expect_true(all(c("sender_node", "receiver_node", "centroid_distance") %in% colnames(dropped)))

      # Passage en mode spatial_object : spc_spatial_info et spc_status exposent étendue + unités natives
      session$setInputs(spc_source = "spatial_object", spc_radius = 3, spc_nperm = 5L)
      session$flushReact()

      info_txt <- output$spc_spatial_info
      expect_match(info_txt, "X \\[", perl = TRUE)
      expect_match(info_txt, "Y \\[", perl = TRUE)
      expect_match(info_txt, "unit", fixed = TRUE)

      session$setInputs(spc_compute = 1L)
      session$flushReact()
      status_spatial <- output$spc_status
      expect_match(status_spatial, "unit", fixed = TRUE)
      expect_match(status_spatial, "X \\[", perl = TRUE)

      ctx_sp <- comm_state$ctx_results$spatial
      kept <- .spc_drop_all_na_cols(ctx_sp$pair_table)
      expect_true(all(c(
        "frac_receiver_within_radius",
        "frac_sender_within_radius",
        "z_score",
        "p_perm"
      ) %in% colnames(kept)))
    }
  )
})
