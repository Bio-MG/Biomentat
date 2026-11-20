# =============================================================================
# test-spatial-sidebar-option-a.R — option A du proposal sidebar Spatial.
# =============================================================================
# WHY THIS FILE EXISTS. The operator reported the "Workflow Spatial" sidebar
# (340 px, 9 accordion panels) feeling like an "infinite scroll" — the sidebar
# scrolled INSIDE the page which itself scrolled, with no visual anchor. The
# mini-proposal (docs/proposals/
# PROPOSAL_SPATIAL_WORKFLOW_SIDEBAR_LENGTH_2026-11-14.md, LOCAL) picked
# option A (quick win CSS) before option B (finish the V1.x split, effort M).
#
# Option A is implemented as three measured layers, each pinned here:
#
#   A1 (structural): app.R passes `fillable = "tab_spatial"` to page_navbar().
#   bslib then marks ONLY that pane html-fill-item/html-fill-container (verified
#   empirically), the parent layout_sidebar inherits a DEFINITE height and the
#   NATIVE sidebar overflow:auto engages (bslib components/scss/sidebar.scss:
#   `.bslib-sidebar-layout > .sidebar { overflow: auto }`) — the sidebar becomes
#   an independent scroll region. Pinned at source level (rendering the whole
#   page_navbar here would drag the entire app UI into this suite).
#
#   A2 (anchor): mod_spatial_ui() wraps the always-visible header block (mirai
#   badge + dataset selector + daemon status + reset button) in a
#   `ts-spatial-side-head` div pinned sticky INSIDE that scroll region.
#
#   A3 (density): accordion body/button paddings reduced via the scoped
#   #spatial-steps stylesheet, and the 3 long pedagogical alerts of panel 6
#   (niches) are folded into <details> — content WORD-FOR-WORD unchanged,
#   container only.
#
# Garde-fou (same as the proposal): NO control id is renamed by option A — the
# drive allowlist pins (test-drive-allowlist-spatial.R) must stay green WITHOUT
# touching them; a rename belongs to option B and must update allowlist + pins
# in the SAME commit.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/defaults.R")
suppressPackageStartupMessages(library(shiny))

source_project_file("tools/pkg_skeleton/R/io_helpers.R")            # %||%
source_project_file("tools/pkg_skeleton/R/error_state.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("tools/pkg_skeleton/R/spatial_stats.R")
load_project_pkg()
source_project_file("tools/pkg_skeleton/R/spatial_async.R")
source_project_file("modules/spatial/mod_spatial_qc.R")

if (!exists(".tr", envir = globalenv()))
  assign(".tr", function(key) key, envir = globalenv())
if (!exists(".t_fmt", envir = globalenv())) {
  assign(".t_fmt", function(template, ...) {
    vals <- list(...)
    for (nm in names(vals)) {
      template <- gsub(paste0("{", nm, "}"), format(vals[[nm]]), template, fixed = TRUE)
    }
    template
  }, envir = globalenv())
}
if (!exists(".tr_plain", envir = globalenv()))
  assign(".tr_plain", function(key) key, envir = globalenv())

# Same stub/restore seam as test-mod-spatial-qc-tabs.R (i18n / .tr_plain at
# call time from globalenv).
.tssoa_stub_ui_symbols <- function() {
  stubbed <- c("i18n", ".tr_plain")
  before <- stats::setNames(
    vapply(stubbed, exists, logical(1), envir = globalenv(), inherits = FALSE),
    stubbed)
  old <- lapply(stubbed, function(nm) {
    if (isTRUE(before[[nm]])) get(nm, envir = globalenv(), inherits = FALSE) else NULL
  })
  names(old) <- stubbed
  assign("i18n", list(t = function(x) x), envir = globalenv())
  assign(".tr_plain", function(x) x, envir = globalenv())
  list(stubbed = stubbed, before = before, old = old)
}

.tssoa_restore_ui_symbols <- function(saved) {
  for (nm in saved$stubbed) {
    if (isTRUE(saved$before[[nm]])) assign(nm, saved$old[[nm]], envir = globalenv())
    else if (exists(nm, envir = globalenv(), inherits = FALSE)) rm(list = nm, envir = globalenv())
  }
  invisible(NULL)
}

.tssoa_ascii <- function(html) {
  out <- iconv(html, from = "UTF-8", to = "ASCII", sub = "?")
  if (is.na(out)) html else out
}

# -----------------------------------------------------------------------------
# A1 — the page-level fillable opt-in (source contract on app.R).
# -----------------------------------------------------------------------------
test_that("A1: app.R opts ONLY the spatial pane into the page fill chain", {
  app_src <- paste(readLines(file.path("..", "..", "app.R"), warn = FALSE),
                   collapse = "\n")
  expect_true(grepl('fillable = "tab_spatial"', app_src, fixed = TRUE),
              info = "page_navbar must declare fillable = 'tab_spatial' (option A structural layer)")
  # The other analysis tabs must NOT be dragged into the fill chain (BULK-9
  # history: the fill chain is a sensitive seam — only tab_spatial is listed).
  expect_false(grepl('fillable = c(', app_src, fixed = TRUE),
               info = "option A must stay a single-pane opt-in (no multi-pane fillable vector)")
  expect_false(grepl("fillable = TRUE", app_src, fixed = TRUE),
               info = "option A must not flip the whole page fillable (BULK-9 regression guard)")
})

# -----------------------------------------------------------------------------
# A2 + A3 — the rendered sidebar: sticky header wrapper + scoped accordion CSS.
# -----------------------------------------------------------------------------
test_that("A2/A3: rendered sidebar carries sticky head + density CSS, ids unchanged", {
  library(bslib)
  saved <- .tssoa_stub_ui_symbols()
  on.exit(.tssoa_restore_ui_symbols(saved), add = TRUE)

  for (f in c("tools/pkg_skeleton/R/spatial_deconv_prep.R",
              "modules/spatial/mod_spatial_pipeline.R",
              "modules/spatial/mod_spatial_cluster.R",
              "modules/spatial/deconv/mod_spatial_deconv_ui.R",
              "modules/spatial/deconv/mod_spatial_deconv.R",
              "modules/spatial/mod_spatial_viz.R",
              "modules/spatial/mod_spatial_multi.R",
              "modules/spatial/mod_spatial_niche.R",
              "modules/spatial/mod_spatial_export.R",
              "modules/spatial/mod_spatial_report.R",
              "modules/spatial/mod_spatial_qc.R",
              "modules/spatial/mod_spatial.R")) {
    source_project_file(f)
  }

  html <- .tssoa_ascii(as.character(mod_spatial_ui("spatial")))

  # A2: the always-visible header block is wrapped for stickiness.
  expect_true(grepl('class="ts-spatial-side-head"', html, fixed = TRUE),
              info = "the sticky wrapper div must be rendered")
  expect_true(grepl(".ts-spatial-side-head { position: sticky", html, fixed = TRUE),
              info = "the sticky rule must ship in the scoped stylesheet")

  # A3: density rules are scoped to the workflow accordion id.
  expect_true(grepl("#spatial-steps .accordion-body { padding: 0.5rem 0.65rem; }",
                    html, fixed = TRUE),
              info = "accordion body padding reduction must be scoped to #spatial-steps")
  expect_true(grepl("#spatial-steps .accordion-button { padding: 0.45rem 0.75rem; }",
                    html, fixed = TRUE),
              info = "accordion button padding reduction must be scoped to #spatial-steps")

  # Garde-fou: NO control id was renamed by option A — the drive-plan seams
  # (dataset selector, daemon status, reset button) are still the same ns ids.
  for (id in c("spatial-active_dataset_ui", "spatial-daemon_status_ui",
               "spatial-btn_reset_daemons", "spatial-steps")) {
    expect_true(grepl(paste0('id="', id, '"'), html, fixed = TRUE),
                info = paste("option A must not rename drive-addressed ids —", id))
  }
})

# -----------------------------------------------------------------------------
# A3 (panel 6) — the 3 theory alerts are folded into <details>, word-for-word.
# -----------------------------------------------------------------------------
test_that("A3: niche panel folds the 3 pedagogical alerts, content unchanged", {
  saved <- .tssoa_stub_ui_symbols()
  on.exit(.tssoa_restore_ui_symbols(saved), add = TRUE)
  source_project_file("modules/spatial/mod_spatial_niche.R")

  html <- .tssoa_ascii(as.character(
    mod_spatial_niche_controls_ui("spatial-niche")))

  # Three folded <details> (one per theory block).
  expect_equal(lengths(regmatches(html, gregexpr("<details", html))), 4L,
               info = "3 new theory summaries + the pre-existing 'Options avancees' details = 4 <details>")

  # Each theory summary is present with its one-click reveal label.
  for (lbl in c("Qu'est-ce qu'une niche spatiale ?",
                "Comment lire le z-score de co-occurrence ?",
                "Qu'est-ce que l'etiquetage aleatoire (Ripley's K) ?")) {
    lbl_ascii <- .tssoa_ascii(lbl)
    expect_true(grepl(lbl_ascii, html, fixed = TRUE),
                info = paste("folded summary label must render:", lbl_ascii))
  }

  # The ORIGINAL explanation strings survive word-for-word (container-only
  # change — the i18n keys are untouched).
  for (key in c(
    "Regroupe les spots/cellules selon la COMPOSITION de leur voisinage spatial",
    "Pour chaque paire de labels (cluster, type cellulaire dominant ou niche)",
    "Pour UN label cible, compare son agregation spatiale observee")) {
    expect_true(grepl(.tssoa_ascii(key), html, fixed = TRUE),
                info = paste("original explanation text must be unchanged:", key))
  }

  # And the operational control ids of the niche module are intact.
  for (id in c("spatial-niche-k_neighbors", "spatial-niche-n_niches",
               "spatial-niche-btn_niches", "spatial-niche-btn_enrichment",
               "spatial-niche-k_neighbors_enrich", "spatial-niche-n_perm_enrich")) {
    expect_true(grepl(paste0('id="', id, '"'), html, fixed = TRUE),
                info = paste("option A must not rename niche control ids —", id))
  }
})
