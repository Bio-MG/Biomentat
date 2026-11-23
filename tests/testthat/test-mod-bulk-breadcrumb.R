# =============================================================================
# test-mod-bulk-breadcrumb.R — Bulk prev/next step breadcrumb bar (W1.2 Part D)
# =============================================================================

suppressPackageStartupMessages(library(shiny))

.load_bulk_breadcrumb_env <- function() {
  e <- new.env(parent = globalenv())
  e$i18n <- list(t = function(x) x)
  e$.tr_plain <- function(x) x
  e$.tr <- function(x) x
  e$.strip_i18n_html <- function(x) x
  e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
  e$layout_sidebar <- bslib::layout_sidebar
  e$sidebar <- bslib::sidebar
  e$accordion <- bslib::accordion
  e$accordion_panel <- bslib::accordion_panel
  e$navset_card_underline <- bslib::navset_card_underline
  e$nav_menu <- bslib::nav_menu
  e$nav_panel <- bslib::nav_panel
  e$nav_select <- function(...) invisible(NULL)
  e$card <- bslib::card
  e$card_header <- bslib::card_header
  e$ts_drive_feed_business_state <- function(...) invisible(NULL)
  e$bulk_purge_shared_results <- function(...) invisible(NULL)
  e$.ts_nul_free_text <- function(x) x

  child_uis <- c(
    "mod_bulk_mapping_ui", "mod_bulk_filter_ui", "mod_bulk_de_ui",
    "mod_bulk_pathways_ui", "mod_bulk_signatures_ui", "mod_bulk_wgcna_ui",
    "mod_bulk_survival_ui", "mod_bulk_pattern_ui", "mod_bulk_dose_response_ui",
    "mod_bulk_network_ui", "mod_bulk_report_ui", "mod_bulk_datasets_ui",
    "mod_bulk_merge_ui", "mod_bulk_de_contrast_bar_ui",
    "mod_bulk_filter_pca_ui", "mod_bulk_filter_qc_ui", "mod_bulk_filter_batch_ui",
    "mod_bulk_de_volcano_ui", "mod_bulk_de_ma_ui", "mod_bulk_de_heatmap_ui",
    "mod_bulk_de_table_ui", "mod_bulk_de_summary_ui", "mod_bulk_de_dtu_ui",
    "mod_bulk_de_multimethod_ui", "mod_bulk_de_venn_ui",
    "mod_bulk_pathways_output_ui", "mod_bulk_signatures_output_ui",
    "mod_bulk_wgcna_output_ui", "mod_bulk_pattern_output_ui",
    "mod_bulk_dose_response_output_ui", "mod_bulk_survival_output_ui",
    "mod_bulk_network_output_ui", "mod_bulk_multi_ui",
    "mod_bulk_merge_output_ui"
  )
  for (fn in child_uis) {
    e[[fn]] <- function(id) tags$div(id = id)
  }

  child_servers <- c(
    "mod_bulk_mapping_server", "mod_bulk_filter_server", "mod_bulk_de_server",
    "mod_bulk_pathways_server", "mod_bulk_signatures_server",
    "mod_bulk_wgcna_server", "mod_bulk_pattern_server",
    "mod_bulk_dose_response_server", "mod_bulk_survival_server",
    "mod_bulk_datasets_server", "mod_bulk_multi_server",
    "mod_bulk_merge_server", "mod_bulk_network_server", "mod_bulk_report_server"
  )
  for (fn in child_servers) {
    e[[fn]] <- function(...) invisible(NULL)
  }

  sys.source(file.path(ts_project_root(), "modules/bulk/mod_bulk.R"), envir = e)
  e
}

test_that("mod_bulk_ui mounts active_step_breadcrumb_ui and ts-bulk-step-bar CSS", {
  e <- .load_bulk_breadcrumb_env()
  html <- as.character(e$mod_bulk_ui("bulk"))
  expect_true(grepl("bulk-active_step_breadcrumb_ui", html, fixed = TRUE))
  expect_true(grepl(".ts-bulk-step-bar", html, fixed = TRUE))
})

test_that("mod_bulk.R implements breadcrumb prev/next controls matching spatial pattern", {
  src <- paste(
    readLines(file.path(ts_project_root(), "modules/bulk/mod_bulk.R"),
              warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  expect_match(src, "active_step_breadcrumb_ui", fixed = TRUE)
  expect_match(src, 'session$ns("btn_prev_step")', fixed = TRUE)
  expect_match(src, 'session$ns("btn_next_step")', fixed = TRUE)
  expect_match(src, "d-flex align-items-center", fixed = TRUE)
  expect_match(src, 'sprintf("%d / %d", idx, length(.step_order))', fixed = TRUE)
  expect_match(src, "observeEvent(input$btn_prev_step", fixed = TRUE)
  expect_match(src, "observeEvent(input$btn_next_step", fixed = TRUE)
})

test_that("mod_bulk_server breadcrumb renders step counter and navigates prev/next", {
  e <- .load_bulk_breadcrumb_env()
  gd <- shiny::reactiveValues(
    language = "fr",
    i18n = list(t = function(x) x),
    workflow_display_mode = "specialise",
    wf_mode_version = 0L,
    bulk_obj = NULL,
    bulk_datasets = list()
  )

  shiny::testServer(e$mod_bulk_server, args = list(global_data = gd), {
    session$setInputs(main_tabs = "tab_pca")

    ui1 <- as.character(output$active_step_breadcrumb_ui$html)
    expect_true(grepl("Bulk RNA", ui1, fixed = TRUE))
    expect_true(grepl("1 / 12", ui1, fixed = TRUE))
    expect_true(grepl("btn_prev_step", ui1, fixed = TRUE))
    expect_true(grepl("btn_next_step", ui1, fixed = TRUE))
    expect_true(grepl("btn-sm btn-outline-secondary py-0 px-2 disabled", ui1, fixed = TRUE))

    # Clicking prev at step 1 stays at step 1
    session$setInputs(btn_prev_step = 1L)
    expect_identical(current_step_rv(), "panel_filter")

    # Clicking next advances to step 2 (panel_de -> tab_volcano)
    session$setInputs(btn_next_step = 1L)
    expect_identical(current_step_rv(), "panel_de")
    expect_identical(shared_rv$active_tab, "tab_volcano")
    ui2 <- as.character(output$active_step_breadcrumb_ui$html)
    expect_true(grepl("2 / 12", ui2, fixed = TRUE))

    # Simulate client nav_select updating input$main_tabs, then click prev
    session$setInputs(main_tabs = "tab_volcano")
    session$setInputs(btn_prev_step = 2L)
    expect_identical(current_step_rv(), "panel_filter")
    expect_identical(shared_rv$active_tab, "tab_pca")
    ui3 <- as.character(output$active_step_breadcrumb_ui$html)
    expect_true(grepl("1 / 12", ui3, fixed = TRUE))
  })
})
