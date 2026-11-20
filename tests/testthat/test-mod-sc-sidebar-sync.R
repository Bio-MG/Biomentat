# =============================================================================
# test-mod-sc-sidebar-sync.R — W1.2 #3b: SC sidebar/accordion <-> nav_menu
# dropdown synchronization without reintroducing the H1 0.6s ping-pong storm.
# =============================================================================

suppressPackageStartupMessages(library(shiny))

.sc_sync_root <- function() {
  if (exists("ts_project_root", mode = "function")) {
    return(ts_project_root())
  }
  normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
}

.load_sc_sync_env <- function() {
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

  sys.source(file.path(.sc_sync_root(), "modules/sc/mod_sc.R"), envir = e)
  e
}

test_that(".sc_tab_sync_plan resyncs accordion even when active_sc_step_rv already matches", {
  e <- .load_sc_sync_env()
  expect_true(exists(".sc_tab_sync_plan", envir = e, mode = "function"))

  # Case 1: active_sc_step_rv() was already updated (e.g., via shared_rv$active_tab
  # or prior navigation), but acc_workflow / acc_analyse are still on grp_prep / NULL.
  plan <- e$.sc_tab_sync_plan(
    tab = "tab_viz",
    current_step = "3_viz",
    open_workflow = "grp_prep",
    open_acc_panel = NULL
  )
  expect_false(is.null(plan))
  expect_identical(plan$info$group, "grp_analyse")
  expect_identical(plan$info$acc, "acc_analyse")
  expect_identical(plan$info$panel, "3_viz")
  expect_false(plan$need_step)
  expect_true(plan$need_workflow)
  expect_true(plan$need_acc)
  expect_true(plan$need_guard)

  # Case 2: acc_workflow was switched to grp_prep while active_sc_step_rv() stayed on 6_pathway
  plan2 <- e$.sc_tab_sync_plan(
    tab = "tab_pathway",
    current_step = "6_pathway",
    open_workflow = "grp_prep",
    open_acc_panel = "6_pathway"
  )
  expect_false(plan2$need_step)
  expect_true(plan2$need_workflow)
  expect_false(plan2$need_acc)
})

test_that(".sc_tab_sync_plan handles nav_menu dropdown navigation across DA sub-tabs", {
  e <- .load_sc_sync_env()

  # Navigating via nav_menu from tab_da_design to tab_da_sccoda keeps panel = 8_da_merged
  # but must still sync da_tabs -> da_methods and da_methods_tabs -> da_sccoda.
  plan <- e$.sc_tab_sync_plan(
    tab = "tab_da_sccoda",
    current_step = "8_da_merged",
    open_workflow = "grp_abondance",
    open_acc_panel = "8_da_merged",
    da_tab = "da_design",
    da_method_tab = "da_milo"
  )
  expect_false(plan$need_step)
  expect_false(plan$need_workflow)
  expect_false(plan$need_acc)
  expect_true(plan$need_da_tab)
  expect_identical(plan$target_da_tab, "da_methods")
  expect_true(plan$need_da_method_tab)
  expect_identical(plan$target_da_method_tab, "da_sccoda")
  expect_true(plan$need_guard)

  # Navigating to tab_da_cross syncs da_tabs -> da_cross
  plan_cross <- e$.sc_tab_sync_plan(
    tab = "tab_da_cross",
    current_step = "8_da_merged",
    open_workflow = "grp_abondance",
    open_acc_panel = "8_da_merged",
    da_tab = "da_methods",
    da_method_tab = "da_milo"
  )
  expect_true(plan_cross$need_da_tab)
  expect_identical(plan_cross$target_da_tab, "da_cross")
  expect_false(plan_cross$need_da_method_tab)
  expect_true(plan_cross$need_guard)
})

test_that(".sc_tab_sync_plan avoids redundant updates when already synchronized", {
  e <- .load_sc_sync_env()

  plan <- e$.sc_tab_sync_plan(
    tab = "tab_communication",
    current_step = "8b_communication",
    open_workflow = c("grp_prep", "grp_analyse"),
    open_acc_panel = c("2_annotation", "8b_communication")
  )
  expect_false(plan$need_step)
  expect_false(plan$need_workflow)
  expect_false(plan$need_acc)
  expect_false(plan$need_da_tab)
  expect_false(plan$need_da_method_tab)
  expect_false(plan$need_guard)

  expect_null(e$.sc_tab_sync_plan(NULL, "1_pipeline"))
  expect_null(e$.sc_tab_sync_plan("", "1_pipeline"))
  expect_null(e$.sc_tab_sync_plan("tab_unknown", "1_pipeline"))
})

test_that(".sc_sync_is_echo preserves the 0.6s H1 anti-ping-pong guard", {
  e <- .load_sc_sync_env()
  expect_true(exists(".sc_sync_is_echo", envir = e, mode = "function"))

  # Same step within 0.6s -> echo suppressed
  expect_true(e$.sc_sync_is_echo("8_da_merged", "8_da_merged", guard_ts = 10.0, now_ts = 10.2))
  # Same step at >= 0.6s -> not an echo
  expect_false(e$.sc_sync_is_echo("8_da_merged", "8_da_merged", guard_ts = 10.0, now_ts = 10.65))
  # Different step within 0.6s -> genuine user click, not an echo
  expect_false(e$.sc_sync_is_echo("4b_pseudobulk", "8_da_merged", guard_ts = 10.0, now_ts = 10.1))
  # NULL guard step -> not an echo
  expect_false(e$.sc_sync_is_echo("8_da_merged", NULL, guard_ts = 10.0, now_ts = 10.1))
})

test_that("mod_sc.R wires .sc_tab_sync_plan in input$main_tabs without breaking H1 guard or mode gating", {
  src <- paste(readLines(file.path(.sc_sync_root(), "modules/sc/mod_sc.R"), warn = FALSE),
               collapse = "\n")

  expect_match(src, ".sc_tab_sync_plan(", fixed = TRUE)
  expect_match(src, ".sc_sync_is_echo(", fixed = TRUE)
  expect_match(src, 'identical(global_data$workflow_display_mode %||% "specialise", "specialise")', fixed = TRUE)
  expect_match(src, ".sc_sync_guard$from_tab_step <- info$panel", fixed = TRUE)
  expect_match(src, "wf_mode_version", fixed = TRUE)
})
