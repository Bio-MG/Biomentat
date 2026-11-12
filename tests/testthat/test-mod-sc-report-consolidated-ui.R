# =============================================================================
# test-mod-sc-report-consolidated-ui.R — modules/sc/mod_sc_report_consolidated.R
# (tiroir « État du rapport »)
# =============================================================================
# CCC_NAV_DRAWER clôture C5 (2026-11-12, GO opérateur « C5 + C6 drawers ») —
# candidat FAIBLE assumé du packet gate 2 : le panneau 9b est PLAT (0 navset,
# 1 carte verdicts) — il n'existe AUCUNE cible de saut. Le tiroir porte donc
# la SEULE moitié applicable du pattern (5f0a906 + bc9c09a) : le RÉSUMÉ d'état
# vivant (verdicts compacts + sections refusées + provenance) — PAS de nav,
# PAS de navset, PAS de marqueur actif (déviation documentée ici et dans la
# fiche promue ; pas de clonage cargo-cult du pattern multi-panneaux).
# =============================================================================

.REPC_UI_FILE <- "modules/sc/mod_sc_report_consolidated.R"

.n_hits_repc <- function(src, pattern) {
  lengths(regmatches(src, gregexpr(pattern, src, fixed = TRUE)))[1]
}

test_that("C5 : le panneau reste PLAT (0 navset) — pas de fausse navigation clonée", {
  src <- paste(readLines(file.path(ts_project_root(), .REPC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Gel du contrat plat : aucune structure de navigation ne doit apparaître
  # (le drawer C5 n'a PAS de cible de saut — déviation assumée) :
  expect_identical(.n_hits_repc(src, "navset_tab("), 0L,
                   info = "un navset est apparu sur le panneau plat 9b")
  expect_identical(.n_hits_repc(src, "nav_select("), 0L,
                   info = "un nav_select sans cible est apparu")
})

test_that("C5 : tiroir état — pastille, .offcanvas, aria-live, Fermer .tr_plain", {
  src <- paste(readLines(file.path(ts_project_root(), .REPC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, '`data-bs-toggle` = "offcanvas"', fixed = TRUE,
               info = "pastille d'ouverture du tiroir absente")
  expect_match(src, "d-flex justify-content-between align-items-center gap-2", fixed = TRUE)
  expect_match(src, '`aria-controls` = ns("rep_drawer")', fixed = TRUE)
  expect_match(src, 'id = ns("rep_drawer")', fixed = TRUE)
  expect_match(src, 'class = "offcanvas offcanvas-end"', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("rep_drawer_summary"))', fixed = TRUE,
               info = "le résumé d'état du tiroir est absent")
  expect_match(src, '`aria-live` = "polite"', fixed = TRUE)
  expect_match(src, "output$rep_drawer_summary <- renderUI(", fixed = TRUE,
               info = "le résumé n'est pas rendu côté serveur")
  expect_match(src, '`aria-label` = .tr_plain("Fermer")', fixed = TRUE,
               info = "aria-label du btn-close sans .tr_plain")
  expect_identical(.n_hits_repc(src, '`aria-label` = i18n$t('), 0L,
                   info = "un aria-label consomme le $t() brut (span i18n dans l'attribut)")
  expect_false(grepl('"Close"', src, fixed = TRUE),
               info = "aria-label anglais hardcoded « Close »")
})

test_that("C5 : résumé honnête — verdicts du contrat, sections refusées, provenance", {
  src <- paste(readLines(file.path(ts_project_root(), .REPC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Le résumé consomme la validation EXISTANTE (zéro recalcul, zéro métrique
  # inventée) :
  expect_match(src, ".collect_and_validate()", fixed = TRUE,
               info = "le résumé ne passe pas par la collecte existante")
  # Verdicts compacts par section (style P3 — liste 2 colonnes, badges d'état) :
  expect_match(src, "column-count:2", fixed = TRUE,
               info = "la liste compacte 2 colonnes est absente")
  expect_match(src, "state_label_color(v$state)", fixed = TRUE,
               info = "les badges d'état ne suivent pas la palette existante")
  # Sections refusées + compteurs — MÊMES clés que le statut existant :
  expect_match(src, '"\u26a0 sections refus\u00e9es (provenance absente) :"', fixed = TRUE)
  expect_match(src, '"section(s) avec r\u00e9sultat."', fixed = TRUE)
  expect_match(src, "cv$ri$provenance_df", fixed = TRUE,
               info = "le compte de provenance n'est pas consommé du contrat")
  # État PENDING : pas de crash sur objet SC absent :
  expect_match(src, "if (is.null(cv))", fixed = TRUE,
               info = "pas de branche NULL (crash tiroir avant tout traitement)")
})
