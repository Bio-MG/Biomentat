# =============================================================================
# test-mod-sc-ui-pills-prep.R — Tests for SC UI fixes:
#   1. Annotation card-header-pills wrapping inside #main_tabs
#   2. Specialized workflow keeping grp_prep (0_autopipeline & 0_mapping) visible
# =============================================================================

suppressPackageStartupMessages(library(shiny))

.sc_pills_prep_root <- function() {
  if (exists("ts_project_root", mode = "function")) {
    return(ts_project_root())
  }
  normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
}

.load_sc_pills_prep_env <- function() {
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

  sys.source(file.path(.sc_pills_prep_root(), "modules/sc/mod_sc.R"), envir = e)
  e
}

test_that("mod_sc_ui emits wrapping CSS for .card-header:has(.card-header-pills) and .card-header-pills", {
  e <- .load_sc_pills_prep_env()
  ui_html <- as.character(e$mod_sc_ui("sc"))

  expect_match(
    ui_html,
    "\\.bslib-card:has\\(#sc-main_tabs\\)\\s+\\.card-header:has\\(\\.card-header-pills\\)\\s*\\{[^}]*flex-wrap:\\s*wrap;[^}]*gap:\\s*0\\.4rem;",
    perl = TRUE
  )
  expect_match(
    ui_html,
    "\\.bslib-card:has\\(#sc-main_tabs\\)\\s+\\.card-header-pills\\s*\\{[^}]*flex-wrap:\\s*wrap;[^}]*max-width:\\s*100%;[^}]*margin-left:\\s*0;",
    perl = TRUE
  )
})

test_that(".sc_wf_mode_css keeps grp_prep visible and does not hide acc_prep siblings", {
  e <- .load_sc_pills_prep_env()
  expect_true(exists(".sc_wf_mode_css", envir = e, mode = "function"))

  ns <- NS("sc")

  # When active step is in acc_prep (e.g. 1_pipeline or 0_autopipeline),
  # grp_prep stays visible and acc_prep siblings (0_autopipeline, 0_mapping, 1_pipeline) are NOT hidden.
  for (prep_step in c("0_autopipeline", "1_pipeline")) {
    css_prep <- e$.sc_wf_mode_css(ns, prep_step, e$.SC_STEP_TAB_MAP[[prep_step]])
    expect_match(
      css_prep,
      "#sc-acc_workflow > .accordion-item:not([data-value='grp_prep'])",
      fixed = TRUE
    )
    expect_false(grepl("#sc-acc_prep > .accordion-item:not(", css_prep, fixed = TRUE))
  }

  # When active step is in another group (e.g. 2_annotation in grp_analyse),
  # BOTH grp_analyse and grp_prep stay visible in #sc-acc_workflow,
  # and only non-active siblings in #sc-acc_analyse are hidden.
  css_annot <- e$.sc_wf_mode_css(ns, "2_annotation", e$.SC_STEP_TAB_MAP[["2_annotation"]])
  expect_match(
    css_annot,
    "#sc-acc_workflow > .accordion-item:not([data-value='grp_analyse']):not([data-value='grp_prep']) { display: none; }",
    fixed = TRUE
  )
  expect_match(
    css_annot,
    "#sc-acc_analyse > .accordion-item:not([data-value='2_annotation']) { display: none; }",
    fixed = TRUE
  )

  # NULL info returns NULL
  expect_null(e$.sc_wf_mode_css(ns, "unknown_step", NULL))
})

test_that("output$wf_mode_style delegates to .sc_wf_mode_css inside tags$style(HTML(...))", {
  src <- paste(
    readLines(file.path(.sc_pills_prep_root(), "modules/sc/mod_sc.R"), warn = FALSE),
    collapse = "\n"
  )
  expect_match(src, ".sc_wf_mode_css(ns, step, info)", fixed = TRUE)
  expect_match(src, "tags$style(HTML(css))", fixed = TRUE)
})
