# =============================================================================
# test-ux-spread-bulk-sc.R — Sticky side-head & compact progress ribbon for
# Bulk & SC (PROPOSAL_UX_SPREAD_SPATIAL_SLIDE_PANEL_2026-11-19.md Parts A & B)
# =============================================================================

suppressPackageStartupMessages(library(shiny))

.load_bulk_ui_env <- function() {
  e <- new.env(parent = globalenv())
  e$i18n <- list(t = function(x) x)
  e$.tr_plain <- function(x) x
  e$.tr <- function(x) x
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

.load_sc_ui_env <- function() {
  e <- new.env(parent = globalenv())
  e$i18n <- list(t = function(x) x)
  e$.tr_plain <- function(x) x
  e$.tr <- function(x) x
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

  sys.source(file.path(ts_project_root(), "modules/sc/mod_sc.R"), envir = e)
  e
}

# ── Bulk pure helper tests ──────────────────────────────────────────────────

test_that(".bulk_pipeline_progress_items handles NULL, empty, and empty-container states", {
  e <- .load_bulk_ui_env()
  expected_codes <- c("map", "filter", "de", "pathway", "report")

  for (st in list(
    NULL,
    list(),
    list(
      mapping_applied = FALSE,
      filtered_counts = matrix(numeric(0), nrow = 0, ncol = 0),
      contrasts       = list(),
      pathway_results = data.frame(),
      vst_mat         = matrix(numeric(0), nrow = 0, ncol = 0)
    )
  )) {
    items <- e$.bulk_pipeline_progress_items(st)
    expect_type(items, "list")
    expect_length(items, 5L)

    codes <- vapply(items, `[[`, character(1), "code")
    expect_identical(codes, expected_codes)

    dones <- vapply(items, `[[`, logical(1), "done")
    expect_true(all(!dones))

    by_code <- stats::setNames(items, codes)
    expect_identical(by_code$map$symbol, "\u26aa")
    expect_identical(by_code$filter$symbol, "\u26aa")
    expect_identical(by_code$de$symbol, "\U0001f512")
    expect_identical(by_code$pathway$symbol, "\U0001f512")
    expect_identical(by_code$report$symbol, "\U0001f512")
  }
})

test_that(".bulk_pipeline_progress_items reflects partial and full progress accurately", {
  e <- .load_bulk_ui_env()

  partial <- e$.bulk_pipeline_progress_items(list(
    mapping_applied = TRUE,
    filtered_counts = matrix(c(10, 20), nrow = 2, ncol = 1)
  ))
  by_code <- stats::setNames(partial, vapply(partial, `[[`, character(1), "code"))
  expect_true(by_code$map$done)
  expect_identical(by_code$map$symbol, "\u2705")
  expect_true(by_code$filter$done)
  expect_identical(by_code$filter$symbol, "\u2705")
  expect_false(by_code$de$done)
  expect_false(by_code$de$locked)
  expect_identical(by_code$de$symbol, "\u26aa")
  expect_true(by_code$pathway$locked)
  expect_identical(by_code$pathway$symbol, "\U0001f512")
  expect_true(by_code$report$locked)
  expect_identical(by_code$report$symbol, "\U0001f512")

  full <- e$.bulk_pipeline_progress_items(list(
    mapping_applied = TRUE,
    filtered_counts = matrix(1, nrow = 2, ncol = 2),
    contrasts       = list(A_vs_B = data.frame(log2FoldChange = 1)),
    pathway_results = list(GO = data.frame(ID = "GO:0001")),
    vst_mat         = matrix(5, nrow = 2, ncol = 2)
  ))
  by_full <- stats::setNames(full, vapply(full, `[[`, character(1), "code"))
  expect_true(all(vapply(full[1:4], `[[`, logical(1), "done")))
  expect_true(all(vapply(full[1:4], `[[`, character(1), "symbol") == "\u2705"))
  expect_false(by_full$report$locked)
  expect_identical(by_full$report$symbol, "\u26aa")
})

# ── SC pure helper tests ────────────────────────────────────────────────────

test_that(".sc_pipeline_progress_items handles NULL, empty, and empty-container states", {
  e <- .load_sc_ui_env()
  expected_codes <- c("qc", "norm", "cluster", "umap", "annot", "markers")

  for (st in list(
    NULL,
    list(),
    list(
      meta_cols        = character(0),
      has_var_features = FALSE,
      reductions       = character(0),
      has_annot        = FALSE,
      markers_data     = data.frame()
    )
  )) {
    items <- e$.sc_pipeline_progress_items(st)
    expect_type(items, "list")
    expect_length(items, 6L)

    codes <- vapply(items, `[[`, character(1), "code")
    expect_identical(codes, expected_codes)

    dones <- vapply(items, `[[`, logical(1), "done")
    expect_true(all(!dones))

    syms <- vapply(items, `[[`, character(1), "symbol")
    expect_true(all(syms == "\u26aa"))
  }
})

test_that(".sc_pipeline_progress_items reflects partial and full progress accurately", {
  e <- .load_sc_ui_env()

  partial <- e$.sc_pipeline_progress_items(list(
    meta_cols        = c("nCount_RNA", "percent.mt", "seurat_clusters"),
    has_var_features = TRUE,
    reductions       = c("pca")
  ))
  by_code <- stats::setNames(partial, vapply(partial, `[[`, character(1), "code"))
  expect_true(by_code$qc$done)
  expect_identical(by_code$qc$symbol, "\u2705")
  expect_true(by_code$norm$done)
  expect_true(by_code$cluster$done)
  expect_false(by_code$umap$done)
  expect_identical(by_code$umap$symbol, "\u26aa")
  expect_false(by_code$annot$done)
  expect_false(by_code$markers$done)

  full <- e$.sc_pipeline_progress_items(list(
    meta_cols        = c("percent.mt", "seurat_clusters"),
    has_var_features = TRUE,
    reductions       = c("pca", "umap"),
    has_annot        = TRUE,
    markers_data     = data.frame(gene = "CD3D", cluster = "0")
  ))
  expect_true(all(vapply(full, `[[`, logical(1), "done")))
  expect_true(all(vapply(full, `[[`, character(1), "symbol") == "\u2705"))
})

# ── Structural UI tests + BULK-9 guard ──────────────────────────────────────

test_that("mod_bulk_ui renders .ts-bulk-side-head + bulk_progress_ribbon_ui without fillable (BULK-9)", {
  e <- .load_bulk_ui_env()
  ui_obj <- e$mod_bulk_ui("bulk")
  html <- as.character(ui_obj)

  expect_match(html, "ts-bulk-side-head", fixed = TRUE)
  expect_match(html, "bulk-bulk_progress_ribbon_ui", fixed = TRUE)
  expect_match(html, "#bulk-acc_bulk .accordion-body", fixed = TRUE)

  # BULK-9 guard: mod_bulk_ui must NOT introduce fillable = TRUE
  bulk_src <- paste(readLines(file.path(ts_project_root(), "modules/bulk/mod_bulk.R"),
                              warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  ui_fn_start <- regexpr("mod_bulk_ui\\s*<-\\s*function", bulk_src)
  ui_fn_end <- regexpr("mod_bulk_server\\s*<-\\s*function", bulk_src)
  ui_src <- substr(bulk_src, ui_fn_start, ui_fn_end)
  expect_false(grepl("fillable\\s*=\\s*TRUE", ui_src),
               info = "BULK-9: mod_bulk_ui must NOT set fillable = TRUE")
  app_src <- paste(readLines(file.path(ts_project_root(), "app.R"),
                             warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_false(grepl('fillable\\s*=[^\\n]*tab_bulk', app_src),
               info = "BULK-9: app.R must NOT opt tab_bulk into fillable")
})

test_that("mod_sc_ui renders .ts-sc-side-head + sc_progress_ribbon_ui", {
  e <- .load_sc_ui_env()
  ui_obj <- e$mod_sc_ui("sc")
  html <- as.character(ui_obj)

  expect_match(html, "ts-sc-side-head", fixed = TRUE)
  expect_match(html, "sc-sc_progress_ribbon_ui", fixed = TRUE)
  expect_match(html, "#sc-acc_workflow .accordion-body", fixed = TRUE)
})
