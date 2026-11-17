# =============================================================================
# test-mod-spatial-padlock.R -- W1.3 : Spatial padlock & disabled gating
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/spatial_async.R")
source_project_file("tools/pkg_skeleton/R/spatial_deconv_prep.R")
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(bsicons)
})
source_project_file("modules/spatial/mod_spatial_pipeline.R")
source_project_file("modules/spatial/mod_spatial_cluster.R")
source_project_file("modules/spatial/deconv/mod_spatial_deconv_ui.R")
source_project_file("modules/spatial/deconv/mod_spatial_deconv_reference.R")
source_project_file("modules/spatial/deconv/mod_spatial_deconv_refviz.R")
source_project_file("modules/spatial/deconv/mod_spatial_deconv_outputs.R")
source_project_file("modules/spatial/deconv/mod_spatial_deconv.R")
source_project_file("modules/spatial/mod_spatial_viz.R")
source_project_file("modules/spatial/mod_spatial_multi.R")
source_project_file("modules/spatial/mod_spatial_niche.R")
source_project_file("modules/spatial/mod_spatial_export.R")
source_project_file("modules/spatial/mod_spatial_report.R")
source_project_file("modules/spatial/mod_spatial_qc.R")
source_project_file("modules/spatial/mod_spatial.R")

.msp_stub_symbols <- function(stubs) {
  saved <- lapply(names(stubs), function(nm) get0(nm, envir = globalenv(), inherits = FALSE))
  names(saved) <- names(stubs)
  for (nm in names(stubs)) assign(nm, stubs[[nm]], envir = globalenv())
  saved
}

.msp_restore_symbols <- function(saved) {
  for (nm in names(saved)) {
    if (is.null(saved[[nm]])) {
      rm(list = intersect(nm, ls(globalenv())), envir = globalenv())
    } else {
      assign(nm, saved[[nm]], envir = globalenv())
    }
  }
}

test_that("mod_spatial_server : padlock titles & breadcrumb lock badge react to spatial_obj, datasets, and language", {
  noop_srv <- function(id, global_data, shared_rv) NULL
  saved <- .msp_stub_symbols(list(
    .strip_i18n_html = function(x) x,
    spatial_daemons_ready = function() TRUE,
    init_spatial_daemons = function(...) TRUE,
    spatial_daemon_status = function() "ready",
    mod_spatial_qc_server = noop_srv,
    mod_spatial_cluster_server = noop_srv,
    mod_spatial_deconv_server = noop_srv,
    mod_spatial_viz_server = noop_srv,
    mod_spatial_multi_server = noop_srv,
    mod_spatial_niche_server = noop_srv,
    mod_spatial_pipeline_server = noop_srv,
    mod_spatial_export_server = noop_srv,
    mod_spatial_report_server = noop_srv
  ))
  on.exit(.msp_restore_symbols(saved), add = TRUE)

  gd <- shiny::reactiveValues(
    language = "fr",
    i18n = list(t = function(k) {
      if (identical(gd$language, "en") && identical(k, "Importez d'abord un echantillon spatial.")) {
        "Import a spatial sample first."
      } else k
    }),
    spatial_obj = NULL,
    spatial_datasets = list(),
    active_spatial_dataset = NULL,
    spatial_results_cache = list(),
    workflow_display_mode = "specialise",
    wf_mode_version = 0L
  )

  shiny::testServer(mod_spatial_server, args = list(global_data = gd), expr = {
    session$flushReact()

    for (nm in paste0("step_title_", c("qc", "cluster", "deconv", "viz", "multi", "niche", "export"))) {
      html <- paste(as.character(output[[nm]]$html), collapse = "\n")
      expect_match(html, "\U0001F512", fixed = TRUE, info = paste("locked padlock expected on", nm))
      expect_match(html, "#8a8f98", fixed = TRUE, info = paste("muted color expected on", nm))
    }
    expect_match(paste(as.character(output$step_title_qc$html), collapse = "\n"),
                 "fa-circle-info|circle-info", perl = TRUE)

    # Step to panel_qc -> breadcrumb shows lock badge while btn_next_step stays clickable
    session$setInputs(steps = "panel_qc")
    session$flushReact()
    bc_locked <- paste(as.character(output$active_step_breadcrumb_ui$html), collapse = "\n")
    expect_match(bc_locked, "\U0001F512", fixed = TRUE)
    expect_false(grepl('btn_next_step"[^>]*disabled', bc_locked))

    # Switch language -> breadcrumb re-renders in English
    gd$language <- "en"
    session$flushReact()
    bc_en <- paste(as.character(output$active_step_breadcrumb_ui$html), collapse = "\n")
    expect_match(bc_en, "Import a spatial sample first.", fixed = TRUE)

    # Load 1 spatial sample -> qc..niche, export unlock; multi stays locked
    gd$spatial_obj <- list(bpcells_dir = "bp", coords = data.frame(id = "s1", x = 1, y = 1))
    gd$spatial_datasets <- list(sample1 = gd$spatial_obj)
    session$flushReact()

    for (nm in paste0("step_title_", c("qc", "cluster", "deconv", "viz", "niche", "export"))) {
      html <- paste(as.character(output[[nm]]$html), collapse = "\n")
      expect_false(grepl("\U0001F512", html, fixed = TRUE), info = paste("expected unlocked:", nm))
    }
    expect_match(paste(as.character(output$step_title_multi$html), collapse = "\n"),
                 "\U0001F512", fixed = TRUE)

    # Load 2nd spatial sample -> multi unlocks
    gd$spatial_datasets <- list(sample1 = gd$spatial_obj, sample2 = gd$spatial_obj)
    session$flushReact()
    expect_false(grepl("\U0001F512",
                       paste(as.character(output$step_title_multi$html), collapse = "\n"),
                       fixed = TRUE))
  })
})

test_that("mod_spatial_cluster : btn_find_markers starts disabled and toggleState observer runs cleanly", {
  saved <- .msp_stub_symbols(list(
    i18n = list(t = function(x) x),
    .tr_plain = function(x) x,
    i18n_plain = function(x) x
  ))
  on.exit(.msp_restore_symbols(saved), add = TRUE)

  ui_html <- paste(as.character(mod_spatial_cluster_controls_ui("spatial-cluster")), collapse = "\n")
  expect_match(ui_html, 'disabled[^>]*id="spatial-cluster-btn_find_markers"', perl = TRUE)

  gd <- shiny::reactiveValues(language = "fr", i18n = NULL, spatial_obj = NULL)
  srv <- shiny::reactiveValues(cluster_labels = NULL, qc_pass_idx = NULL)
  shiny::testServer(mod_spatial_cluster_server, args = list(global_data = gd, shared_rv = srv), expr = {
    session$flushReact()
    gd$spatial_obj <- list(bpcells_dir = "bp", coords = data.frame(id = "c1", x = 1, y = 1))
    session$flushReact()
    srv$cluster_labels <- c(c1 = "1")
    session$flushReact()
    expect_true(TRUE)
  })
})

test_that("mod_spatial_deconv : btn_ref_viz/btn_deconv start disabled and toggleState observer runs cleanly", {
  saved <- .msp_stub_symbols(list(
    i18n = list(t = function(x) x),
    .tr_plain = function(x) x,
    i18n_plain = function(x) x
  ))
  on.exit(.msp_restore_symbols(saved), add = TRUE)

  ui_html <- paste(as.character(mod_spatial_deconv_controls_ui("spatial-deconv")), collapse = "\n")
  expect_match(ui_html, 'disabled[^>]*id="spatial-deconv-btn_ref_viz"', perl = TRUE)
  expect_match(ui_html, 'disabled[^>]*id="spatial-deconv-btn_deconv"', perl = TRUE)

  gd <- shiny::reactiveValues(language = "fr", i18n = NULL, spatial_obj = NULL, spatial_reference = NULL)
  srv <- shiny::reactiveValues(deconv_props = NULL, deconv_params = NULL)
  shiny::testServer(mod_spatial_deconv_server, args = list(global_data = gd, shared_rv = srv), expr = {
    session$setInputs(mode = "rctd", ref_source = "local", ref_celltype_col = "")
    session$flushReact()
    session$setInputs(mode = "stdeconvolve")
    session$flushReact()
    gd$spatial_obj <- list(bpcells_dir = "bp", coords = data.frame(id = "c1", x = 1, y = 1))
    gd$spatial_reference <- list(path = "ref.rds", n_cells = 10L)
    session$setInputs(mode = "rctd", ref_source = "shared")
    session$flushReact()
    expect_true(TRUE)
  })
})
