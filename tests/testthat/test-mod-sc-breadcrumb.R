# =============================================================================
# test-mod-sc-breadcrumb.R — W1.2 Part D: SC prev/next step breadcrumb bar
# =============================================================================

suppressPackageStartupMessages(library(shiny))

.sc_bc_root <- function() {
  if (exists("ts_project_root", mode = "function")) {
    return(ts_project_root())
  }
  normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
}

.load_sc_bc_env <- function() {
  e <- new.env(parent = globalenv())
  e$i18n <- list(t = function(x) x)
  e$.tr_plain <- function(x) x
  e$.tr <- function(x) x
  e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
  e$layout_sidebar <- bslib::layout_sidebar
  e$sidebar <- bslib::sidebar
  e$accordion <- bslib::accordion
  e$accordion_panel <- bslib::accordion_panel
  e$navset_card_underline <- bslib::navset_card_underline
  e$navset_pill <- bslib::navset_pill
  e$nav_menu <- bslib::nav_menu
  e$nav_panel <- bslib::nav_panel
  e$card <- bslib::card
  e$card_header <- bslib::card_header

  child_uis <- c(
    "mod_sc_annotation_output_ui", "mod_sc_annotation_ui",
    "mod_sc_communication_output_ui", "mod_sc_communication_ui",
    "mod_sc_corr_output_ui", "mod_sc_corr_ui",
    "mod_sc_da_cross_output_ui", "mod_sc_da_cross_ui",
    "mod_sc_da_design_output_ui", "mod_sc_da_design_ui",
    "mod_sc_da_milo_output_ui", "mod_sc_da_milo_ui",
    "mod_sc_da_sccoda_output_ui", "mod_sc_da_sccoda_ui",
    "mod_sc_datasets_ui", "mod_sc_mapping_ui",
    "mod_sc_markers_output_ui", "mod_sc_markers_ui",
    "mod_sc_metadata_ui", "mod_sc_pathways_output_ui",
    "mod_sc_pathways_ui", "mod_sc_pipeline_ui",
    "mod_sc_pseudobulk_output_ui", "mod_sc_pseudobulk_ui",
    "mod_sc_rarity_output_ui", "mod_sc_rarity_ui",
    "mod_sc_report_consolidated_output_ui", "mod_sc_report_consolidated_ui",
    "mod_sc_trajectory_output_ui", "mod_sc_trajectory_ui",
    "mod_sc_velocity_output_ui", "mod_sc_velocity_ui",
    "mod_sc_viz_output_ui", "mod_sc_viz_ui"
  )
  for (fn in child_uis) {
    e[[fn]] <- function(id) tags$div(id = id)
  }

  sys.source(file.path(.sc_bc_root(), "modules/sc/mod_sc.R"), envir = e)
  e
}

test_that("mod_sc_ui mounts active_step_breadcrumb_ui in the results header", {
  e <- .load_sc_bc_env()
  ui_html <- as.character(e$mod_sc_ui("sc"))

  expect_match(ui_html, 'id="sc-active_step_breadcrumb_ui"', fixed = TRUE)
  expect_match(ui_html, 'id="sc-sc_step_sync_bar"', fixed = TRUE)
})

test_that(".sc_breadcrumb_ui_tag renders prev/next buttons, d / N counter, and boundary disabled states", {
  e <- .load_sc_bc_env()
  ns <- NS("sc")
  steps <- e$.SC_STEP_ORDER
  n_steps <- length(steps)
  expect_true(n_steps >= 10L)

  has_disabled <- function(html, btn_id) {
    m <- regmatches(html, regexpr(sprintf('<button[^>]*id="%s"[^>]*>', btn_id), html, perl = TRUE))
    length(m) == 1L && nzchar(m) && grepl("\\bdisabled\\b", m, perl = TRUE)
  }

  # First step: 1 / N, btn_prev_step disabled, btn_next_step enabled
  html_first <- as.character(e$.sc_breadcrumb_ui_tag(ns, steps[1L]))
  expect_match(html_first, "ts-sc-step-sync-bar", fixed = TRUE)
  expect_match(html_first, "d-flex align-items-center", fixed = TRUE)
  expect_match(html_first, 'id="sc-btn_prev_step"', fixed = TRUE)
  expect_match(html_first, 'id="sc-btn_next_step"', fixed = TRUE)
  expect_match(html_first, sprintf("1 / %d", n_steps), fixed = TRUE)
  expect_true(has_disabled(html_first, "sc-btn_prev_step"))
  expect_false(has_disabled(html_first, "sc-btn_next_step"))

  # Middle step (1_pipeline): 2 / N, both buttons enabled
  html_mid <- as.character(e$.sc_breadcrumb_ui_tag(ns, "1_pipeline"))
  expect_match(html_mid, sprintf("2 / %d", n_steps), fixed = TRUE)
  expect_false(has_disabled(html_mid, "sc-btn_prev_step"))
  expect_false(has_disabled(html_mid, "sc-btn_next_step"))
  expect_match(html_mid, "Pipeline (QC", fixed = TRUE)

  # Last step: N / N, btn_prev_step enabled, btn_next_step disabled
  html_last <- as.character(e$.sc_breadcrumb_ui_tag(ns, steps[n_steps]))
  expect_match(html_last, sprintf("%d / %d", n_steps, n_steps), fixed = TRUE)
  expect_false(has_disabled(html_last, "sc-btn_prev_step"))
  expect_true(has_disabled(html_last, "sc-btn_next_step"))
})

test_that(".sc_breadcrumb_adjacent_step navigates prev/next across .SC_STEP_ORDER safely", {
  e <- .load_sc_bc_env()
  steps <- e$.SC_STEP_ORDER
  n_steps <- length(steps)

  expect_identical(e$.sc_breadcrumb_step_index(steps[1L]), 1L)
  expect_identical(e$.sc_breadcrumb_step_index("1_pipeline"), 2L)
  expect_identical(e$.sc_breadcrumb_step_index(steps[n_steps]), n_steps)
  expect_identical(e$.sc_breadcrumb_step_index("unknown_step"), 1L)
  expect_identical(e$.sc_breadcrumb_step_index(NULL), 1L)

  # Boundary guards
  expect_null(e$.sc_breadcrumb_adjacent_step(steps[1L], "prev"))
  expect_identical(e$.sc_breadcrumb_adjacent_step(steps[1L], "next"), steps[2L])
  expect_identical(e$.sc_breadcrumb_adjacent_step(steps[n_steps], "prev"), steps[n_steps - 1L])
  expect_null(e$.sc_breadcrumb_adjacent_step(steps[n_steps], "next"))

  # Full forward walk reaches every step in order
  cur <- steps[1L]
  visited <- cur
  repeat {
    nxt <- e$.sc_breadcrumb_adjacent_step(cur, "next")
    if (is.null(nxt)) break
    visited <- c(visited, nxt)
    cur <- nxt
  }
  expect_identical(visited, steps)

  # Full backward walk returns to step 1
  back <- cur
  repeat {
    prv <- e$.sc_breadcrumb_adjacent_step(cur, "prev")
    if (is.null(prv)) break
    back <- c(back, prv)
    cur <- prv
  }
  expect_identical(back, rev(steps))

  # Unknown / NULL step returns NULL
  expect_null(e$.sc_breadcrumb_adjacent_step("unknown_step", "next"))
  expect_null(e$.sc_breadcrumb_adjacent_step(NULL, "prev"))
})

test_that("mod_sc.R wires active_step_breadcrumb_ui and btn_prev_step / btn_next_step observers", {
  src <- paste(readLines(file.path(.sc_bc_root(), "modules/sc/mod_sc.R"), warn = FALSE),
               collapse = "\n")

  expect_match(src, 'uiOutput(ns("active_step_breadcrumb_ui"))', fixed = TRUE)
  expect_match(src, "output$active_step_breadcrumb_ui <- renderUI({", fixed = TRUE)
  expect_match(src, "observeEvent(input$btn_prev_step, {", fixed = TRUE)
  expect_match(src, "observeEvent(input$btn_next_step, {", fixed = TRUE)
  expect_match(src, 'actionButton(ns("btn_prev_step")', fixed = TRUE)
  expect_match(src, 'actionButton(ns("btn_next_step")', fixed = TRUE)
  expect_match(src, 'sprintf("%d / %d", idx, length(step_order))', fixed = TRUE)
})
