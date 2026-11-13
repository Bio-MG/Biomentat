# =============================================================================
# test-mod-bulk-padlock-stepsync.R — Bulk padlocks (3e, 3g) & step-sync (§3.1 & §3.2)
# =============================================================================

suppressPackageStartupMessages(library(shiny))

.load_mod_bulk_env <- function() {
  e <- new.env(parent = globalenv())
  e$i18n <- list(t = function(x) x)
  e$.tr_plain <- function(x) x
  e$.tr <- function(x) x
  e$layout_sidebar <- bslib::layout_sidebar
  e$sidebar <- bslib::sidebar
  e$accordion <- bslib::accordion
  e$accordion_panel <- bslib::accordion_panel
  e$navset_card_underline <- bslib::navset_card_underline
  e$nav_menu <- bslib::nav_menu
  e$nav_panel <- bslib::nav_panel
  e$card <- bslib::card
  e$card_header <- bslib::card_header

  # Stub child UI functions so mod_bulk_ui() can be evaluated in isolation
  child_uis <- c(
    "mod_bulk_mapping_ui", "mod_bulk_filter_ui", "mod_bulk_de_ui",
    "mod_bulk_pathways_ui", "mod_bulk_signatures_ui", "mod_bulk_wgcna_ui",
    "mod_bulk_survival_ui", "mod_bulk_pattern_ui", "mod_bulk_dose_response_ui",
    "mod_bulk_network_ui", "mod_bulk_report_ui", "mod_bulk_datasets_ui",
    "mod_bulk_merge_ui", "mod_bulk_de_contrast_bar_ui",
    "mod_bulk_filter_pca_ui", "mod_bulk_filter_qc_ui", "mod_bulk_filter_batch_ui",
    "mod_bulk_de_volcano_ui", "mod_bulk_de_ma_ui", "mod_bulk_de_heatmap_ui",
    "mod_bulk_de_table_ui", "mod_bulk_de_summary_ui",
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

  sys.source(file.path(ts_project_root(), "modules/bulk/mod_bulk.R"), envir = e)
  e
}

test_that("sidebar panels 3e and 3g have greyed span + padlock gating matching 3c/3d/3f", {
  e <- .load_mod_bulk_env()
  html <- as.character(e$mod_bulk_ui("bulk"))

  expected_locked <- c(
    "3c. Réseau de co-expression (WGCNA)",
    "3d. Survie & association clinique",
    "3e. Clustering de profils",
    "3f. Dose–réponse / time-course",
    "3g. Réseau PCSF (interactome)"
  )
  for (lbl in expected_locked) {
    pat <- paste0('<span style="color:#8a8f98">', lbl, "</span> 🔒")
    expect_true(
      grepl(pat, html, fixed = TRUE),
      info = paste("Expected greyed + padlock label for:", lbl)
    )
  }
})

test_that("mod_bulk_ui mounts step_sync_bar and signatures_guard_ui in header", {
  e <- .load_mod_bulk_env()
  html <- as.character(e$mod_bulk_ui("bulk"))
  expect_true(grepl("bulk-bulk_step_sync_bar", html, fixed = TRUE))
  expect_true(grepl("bulk-bulk_signatures_guard_ui", html, fixed = TRUE))
})

test_that(".bulk_tab_step_meta maps 4 · Signatures & profils and 5 · Réseau PCSF to sidebar steps 3b-3g", {
  e <- .load_mod_bulk_env()
  expect_true(exists(".bulk_tab_step_meta", envir = e, inherits = FALSE))

  expect_identical(e$.bulk_tab_step_meta("tab_signatures")$step, "3b")
  expect_identical(e$.bulk_tab_step_meta("tab_wgcna")$step, "3c")
  expect_identical(e$.bulk_tab_step_meta("tab_survival")$step, "3d")
  expect_identical(e$.bulk_tab_step_meta("tab_pattern")$step, "3e")
  expect_identical(e$.bulk_tab_step_meta("tab_dose")$step, "3f")
  expect_identical(e$.bulk_tab_step_meta("tab_bulk_network")$step, "3g")
  expect_identical(e$.bulk_tab_step_meta("tab_bulk_network")$panel, "panel_network")
})

test_that(".bulk_signatures_ready requires both vst_mat and active_contrast", {
  e <- .load_mod_bulk_env()
  expect_true(exists(".bulk_signatures_ready", envir = e, inherits = FALSE))

  expect_false(e$.bulk_signatures_ready(list(vst_mat = NULL, active_contrast = NULL))$ready)
  expect_false(e$.bulk_signatures_ready(list(vst_mat = matrix(1, 2, 2), active_contrast = NULL))$ready)
  expect_false(e$.bulk_signatures_ready(list(vst_mat = NULL, active_contrast = "A_vs_B"))$ready)
  expect_true(e$.bulk_signatures_ready(list(vst_mat = matrix(1, 2, 2), active_contrast = "A_vs_B"))$ready)
})
