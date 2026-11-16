# =============================================================================
# test-p1-workflow-display-mode.R — Étape 3 (P1) du proposal
# PROPOSAL_UX_UI_HARMONISATION_WORKFLOW_2026-11-15 : mode d'affichage du
# workflow { Spécialisé | Complet }, défaut SPÉCIALISÉ (arbitrage opérateur).
# =============================================================================
# WHY SOURCE-LEVEL. The mechanism is (a) ONE page-level radio feeding
# global_data$workflow_display_mode and (b) per-module shinyjs visibility
# observers toggling `.accordion-item[data-value]` CSS selectors — the
# behaviour lives in observers, which only a live browser can exercise (done
# on the boot of the day, handoff §U). These pins freeze the CONTRACT: the
# parameter name, the DEFAULT (Spécialisé), the drive-safe mechanism (CSS
# selectors on data-value, NO value/id renames) and the spécialisé-only
# results→sidebar auto-sync for Bulk/SC (Spatial already ships its
# bidirectional sync in both modes).
# =============================================================================

.src_of <- function(...) {
  paste(readLines(file.path("..", "..", ...), warn = FALSE), collapse = "\n")
}

test_that("P1: the page parameter exists, defaults to Spécialisé, feeds global_data", {
  src <- .src_of("app.R")
  expect_match(src, 'radioButtons("workflow_display_mode"', fixed = TRUE,
               info = "the ⚙️ Paramètres section must carry the workflow display radio")
  expect_match(src, 'selected = "specialise"', fixed = TRUE,
               info = "the DEFAULT must be Spécialisé (operator arbitrage 2026-11-15)")
  # P2 init déplacé DANS la création de global_data (avant les mod_*_server :
  # une lecture de champ NULL ne crée PAS de dépendance — MEASURÉ boot33d).
  expect_match(src, 'workflow_display_mode = "specialise"', fixed = TRUE,
               info = "global_data must be initialised to Spécialisé")
  expect_match(src, 'wf_mode_version = 0L', fixed = TRUE,
               info = "the version counter must exist BEFORE the module servers")
  # Both directions wired: radio -> global_data, and restored session -> radio.
  expect_match(src, 'global_data$workflow_display_mode <- input$workflow_display_mode', fixed = TRUE)
  expect_match(src, 'updateRadioButtons(session, "workflow_display_mode"', fixed = TRUE)
})

test_that("P1: the mechanism is drive-safe (CSS on data-value, no value renames)", {
  for (f in c(file.path("modules", "spatial", "mod_spatial.R"),
              file.path("modules", "bulk", "mod_bulk.R"),
              file.path("modules", "sc", "mod_sc.R"))) {
    src <- .src_of(f)
    expect_match(src, ".accordion-item:not([data-value=", fixed = TRUE,
                 info = paste("visibility must hide via .accordion-item[data-value] selectors —", f))
    # per-line check: TRE's `.` crosses newlines, and the guard must reject
    # only a SAME-LINE mixing of an accordion definition with the mode.
    lines <- readLines(file.path("..", "..", f), warn = FALSE)
    expect_false(any(grepl("accordion_panel\\(.*workflow_display_mode", lines)),
                 info = paste("the mode must NOT alter accordion definitions —", f))
  }
  # The pinned workflow values survive untouched (drive allowlist seams).
  sc <- .src_of("modules", "sc", "mod_sc.R")
  for (v in c('"2_annotation"', '"8b_communication"', '"8_da_merged"', '"9b_report_consolide"')) {
    expect_match(sc, paste0("value = ", v), fixed = TRUE,
                 info = paste("SC accordion value must be untouched:", v))
  }
  bulk <- .src_of("modules", "bulk", "mod_bulk.R")
  for (v in c('"panel_filter"', '"panel_network"', '"panel_report"')) {
    expect_match(bulk, paste0("value = ", v), fixed = TRUE,
                 info = paste("Bulk accordion value must be untouched:", v))
  }
})

test_that("P1: visibility is a reactive STYLE (survives reconnect); auto-sync Spécialisé-only", {
  # MEASURÉ boot33 : les messages shinyjs envoyés avant la connexion client
  # sont PERDUS — la visibilité voyage donc par un style CSS renderUI
  # (re-envoyé au client à chaque connect), pas par shinyjs.
  sp <- .src_of("modules", "spatial", "mod_spatial.R")
  expect_match(sp, "output$wf_mode_style <- renderUI({", fixed = TRUE,
               info = "Spatial visibility must be a reactive style output")
  expect_match(sp, ".ts_wf_last_panel", fixed = TRUE,
               info = "Spécialisé must always keep one step visible (reopen last)")
  expect_false(grepl("shinyjs::hide", sp),
               info = "no message-loss channel (shinyjs) may drive the visibility")

  bulk <- .src_of("modules", "bulk", "mod_bulk.R")
  expect_match(bulk, "output$wf_mode_style <- renderUI({", fixed = TRUE)
  expect_match(bulk, 'identical(global_data$workflow_display_mode %||% "specialise", "specialise")',
               fixed = TRUE,
               info = "the Bulk results→sidebar auto-sync must run in Spécialisé ONLY")

  sc <- .src_of("modules", "sc", "mod_sc.R")
  expect_match(sc, "output$wf_mode_style <- renderUI({", fixed = TRUE)
  expect_match(sc, 'identical(global_data$workflow_display_mode %||% "specialise", "specialise")',
               fixed = TRUE,
               info = "the SC results→group/panel auto-sync must run in Spécialisé ONLY")
})
