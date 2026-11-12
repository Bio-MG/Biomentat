# =============================================================================
# test-mod-sc-viz-drawer-ui.R — modules/sc/mod_sc_viz.R (tiroir « État de la
# vue »)
# =============================================================================
# CCC_NAV_DRAWER clôture C6 (2026-11-12, GO opérateur « C5 + C6 drawers ») —
# candidat MOYEN du packet gate 2, re-mesuré : le panneau Viz est un CONTENEUR
# UNIQUE (1 plot_container scrollable, 0 navset) — pas de cible de saut non
# plus. Le tiroir porte la moitié applicable du pattern : le RÉSUMÉ d'état
# vivant (catégorie courante, palette, thème/dimensions/format d'export — des
# réglages sinon cachés dans le popover) — PAS de nav (déviation documentée).
#
# La carte des catégories (.viz_category_map, inline dans
# output$viz_category_badge) est FACTORISÉE en helper local UNIQUE — le tiroir
# consomme le même helper : zéro divergence possible sur le libellé.
# =============================================================================

.VIZ_UI_FILE <- "modules/sc/mod_sc_viz.R"

.n_hits_viz <- function(src, pattern) {
  lengths(regmatches(src, gregexpr(pattern, src, fixed = TRUE)))[1]
}

test_that("C6 : le panneau reste un conteneur unique (0 navset) — pas de fausse navigation clonée", {
  src <- paste(readLines(file.path(ts_project_root(), .VIZ_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_identical(.n_hits_viz(src, "navset_tab("), 0L,
                   info = "un navset est apparu sur le conteneur unique Viz")
  expect_identical(.n_hits_viz(src, "nav_select("), 0L,
                   info = "un nav_select sans cible est apparu")
})

test_that("C6 : tiroir état — pastille, .offcanvas, aria-live, Fermer .tr_plain", {
  src <- paste(readLines(file.path(ts_project_root(), .VIZ_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, '`data-bs-toggle` = "offcanvas"', fixed = TRUE,
               info = "pastille d'ouverture du tiroir absente")
  expect_match(src, '`aria-controls` = ns("viz_drawer")', fixed = TRUE)
  expect_match(src, 'id = ns("viz_drawer")', fixed = TRUE)
  expect_match(src, 'class = "offcanvas offcanvas-end"', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("viz_drawer_summary"))', fixed = TRUE,
               info = "le résumé d'état du tiroir est absent")
  expect_match(src, '`aria-live` = "polite"', fixed = TRUE)
  expect_match(src, "output$viz_drawer_summary <- renderUI(", fixed = TRUE,
               info = "le résumé n'est pas rendu côté serveur")
  expect_match(src, '`aria-label` = .tr_plain("Fermer")', fixed = TRUE,
               info = "aria-label du btn-close sans .tr_plain")
  expect_identical(.n_hits_viz(src, '`aria-label` = i18n$t('), 0L,
                   info = "un aria-label consomme le $t() brut (span i18n dans l'attribut)")
  expect_false(grepl('"Close"', src, fixed = TRUE),
               info = "aria-label anglais hardcoded « Close »")
})

test_that("C6 : résumé honnête — catégorie factorisée, palette, réglages d'export", {
  src <- paste(readLines(file.path(ts_project_root(), .VIZ_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # La carte des catégories est factorisée en helper UNIQUE (consommé par le
  # badge ET le tiroir — zéro divergence possible) :
  expect_identical(.n_hits_viz(src, ".viz_category_label <- function"), 1L,
                   info = "le helper de catégorie doit exister UNE fois")
  expect_match(src, "output$viz_category_badge <- renderUI({\n      req(input$viz_type)\n      global_data$language\n      cat_label <- .viz_category_label(input$viz_type)",
               fixed = TRUE,
               info = "le badge n'utilise pas le helper factorisé")
  expect_match(src, ".viz_category_label(input$viz_type)", fixed = TRUE)
  # Réglages d'export consommés de l'état (sinon cachés dans le popover) :
  for (pat in c("input$plot_theme", "input$plot_width", "input$plot_height",
                "input$export_format", "input$sc_palette")) {
    expect_match(src, pat, fixed = TRUE,
                 info = paste0("réglage non consommé : ", pat))
  }
})
