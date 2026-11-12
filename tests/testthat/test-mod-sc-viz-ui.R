# =============================================================================
# test-mod-sc-viz-ui.R — modules/sc/mod_sc_viz.R (UI sidebar)
# =============================================================================
# Item 2E du packet round-7 (O28 batch B3).
#
# Constat 2E : à HEAD, les bannières « Bloc 2 — Régulation & Fonction (à
# venir) » et « Bloc 3 » surchargeaient la sidebar en permanence (preuve
# mod_sc_viz.R:227-232 à HEAD 88f9a43 ; spéc §2.10). Correctif : encapsulation
# dans un accordion FERMÉ par défaut « Fonctionnalités avancées / Roadmap » —
# textes des deux blocs INCHANGÉS, espace vertical libéré immédiatement.
#
# Preuve comportementale : mod_sc_viz_ui() est rendue headlessly (i18n et
# .tr_plain shimmés) ; bslib marque l'état des panneaux via aria-expanded
# (mesuré : closed -> aria-expanded="false" 1x / "true" 0x ; open -> inverse).
# =============================================================================

source_project_file("modules/sc/mod_sc_viz.R")

suppressPackageStartupMessages(library(shiny))  # NS(), tags, radioButtons...
suppressPackageStartupMessages(library(bslib))  # accordion / accordion_panel

# Shims de construction UI (résolus depuis globalenv à l'appel) : identité —
# sémantique exacte de global.R:345/:369 quand I18N_AVAILABLE n'est pas TRUE.
.vizui_shim <- function(code) {
  old_i18n <- get0("i18n", envir = globalenv())
  old_trp  <- get0(".tr_plain", envir = globalenv())
  old_tr   <- get0(".tr", envir = globalenv())
  assign("i18n", list(t = function(x) x), envir = globalenv())
  assign(".tr_plain", function(x) x, envir = globalenv())
  assign(".tr", function(x) x, envir = globalenv())
  on.exit({
    if (is.null(old_i18n)) rm(list = "i18n", envir = globalenv())
    else assign("i18n", old_i18n, envir = globalenv())
    if (is.null(old_trp)) {
      if (exists(".tr_plain", envir = globalenv()))
        rm(list = ".tr_plain", envir = globalenv())
    } else {
      assign(".tr_plain", old_trp, envir = globalenv())
    }
    if (is.null(old_tr)) {
      if (exists(".tr", envir = globalenv()))
        rm(list = ".tr", envir = globalenv())
    } else {
      assign(".tr", old_tr, envir = globalenv())
    }
  }, add = TRUE)
  force(code)
}

test_that("2E : la sidebar Viz ne contient plus l'accordéon roadmap obsolète", {
  suppressPackageStartupMessages(library(bslib))
  ui <- .vizui_shim(mod_sc_viz_ui("viz"))
  # tagList() de la sidebar => shiny.tag.list (et non shiny.tag simple) :
  expect_true(inherits(ui, "shiny.tag.list") || inherits(ui, "shiny.tag"))
  html <- paste(as.character(ui), collapse = " ")
  # Le panneau de roadmap obsolète a été retiré (demande opérateur 2026-11-12) :
  expect_false(grepl("Fonctionnalit\u00e9s avanc\u00e9es / Roadmap", html, fixed = TRUE))
  expect_false(grepl("Bloc 2 \u2014 R\u00e9gulation", html, fixed = TRUE))
  expect_false(grepl("Bloc 3 \u2014 Dynamique", html, fixed = TRUE))
})
