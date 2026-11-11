# =============================================================================
# test-mod-sc-da-cross-ui.R — modules/sc/mod_sc_da_cross.R (tiroir
# « Navigation »)
# =============================================================================
# CCC_NAV_DRAWER vague B jalon J4 (2026-11-12, GO opérateur « GO J3 + J4 »
# après vérité-terrain live du tiroir CCC chargé) — généralisation du pattern
# tiroir optimisé au panneau Vues croisées (8f) — candidat C4 du packet gate
# 2 (le pont 8d/8e tire sa valeur du fait que Milo + scCODA sont déjà tirés).
# Ordre : J4 AVANT J3 (mod_sc_velocity.R EN VOL lane scvelo — staging par
# hunks au J3, hunks étrangers jamais commités).
#
# ⚠️ Limite documentée (motif test-mod-sc-communication-ui.R) : verrou source
# structurel + identité AST des jumelles (rendu headless complet non réalisable
# à coût raisonnable).
#
# ⚠️ test-da-cross-views-contract-freeze.R gèle les top-level de
# mod_sc_da_cross.R à exactement 3 fonctions — la carte vit EN LOCAL de
# fonction, en DEUX exemplaires jumelles.
# =============================================================================

.CROSS_UI_FILE <- "modules/sc/mod_sc_da_cross.R"

.n_hits_cross <- function(src, pattern) {
  lengths(regmatches(src, gregexpr(pattern, src, fixed = TRUE)))[1]
}

.CROSS_DRAWER_VALUES <- c("cross_tab_concordance", "cross_tab_map",
                          "cross_tab_compo", "cross_tab_effects",
                          "cross_tab_limites")

test_that("cross drawer : navset id + 5 valeurs ASCII explicites (cibles du saut)", {
  src <- paste(readLines(file.path(ts_project_root(), .CROSS_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_identical(.n_hits_cross(src, "navset_tab("), 1L)
  expect_match(src, 'id = ns("cross_out_tabs")', fixed = TRUE,
               info = "le navset Vues croisées n'a pas d'id — nav_select ne peut pas cibler")
  for (v in .CROSS_DRAWER_VALUES) {
    expect_match(src, sprintf('value = "%s"', v), fixed = TRUE, info = v)
  }
})

# Extrait TOUTES les affectations `.cross_drawer_groups` du module, À TOUTE
# profondeur (gel test-da-cross-views-contract-freeze.R — la carte vit en
# jumelles locales, jamais en constante top-level).
.cross_extract_drawer_groups <- function() {
  p <- ts_ast_parse(.CROSS_UI_FILE)
  hits <- ts_ast_find_all(p, function(x) {
    if (!is.call(x) || !identical(x[[1]], as.name("<-"))) return(FALSE)
    lhs <- tryCatch(ts_ast_deparse(x[[2]]), error = function(e) "")
    identical(lhs, ".cross_drawer_groups")
  })
  lapply(hits, function(h) eval(h[[3]], envir = new.env(parent = baseenv())))
}

test_that("cross drawer : .cross_drawer_groups (2 jumelles) = 2 groupes, 5 panneaux, valeurs uniques", {
  twins <- .cross_extract_drawer_groups()
  expect_identical(length(twins), 2L,
                   info = "attendu : jumelle output_ui + jumelle serveur, rien d'autre")
  expect_false(is.null(twins[[1]]), info = "jumelle introuvable dans le module")
  expect_identical(twins[[1]], twins[[2]],
                   info = "les deux cartes .cross_drawer_groups divergent")
  grp <- twins[[1]]
  expect_identical(length(grp), 2L)
  vals <- unlist(lapply(grp, function(g) vapply(g, function(p) p$value, character(1))))
  labs <- unlist(lapply(grp, function(g) vapply(g, function(p) p$label, character(1))))
  expect_identical(length(vals), 5L)
  expect_identical(anyDuplicated(vals), 0L)
  expect_setequal(vals, .CROSS_DRAWER_VALUES)
  src <- paste(readLines(file.path(ts_project_root(), .CROSS_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  for (lb in labs) {
    expect_match(src, sprintf('nav_panel(i18n$t("%s")', lb), fixed = TRUE, info = lb)
  }
})

test_that("cross drawer : pastille, .offcanvas, sauts + marqueur actif, compteurs, aria-live, Fermer .tr_plain", {
  src <- paste(readLines(file.path(ts_project_root(), .CROSS_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, '`data-bs-toggle` = "offcanvas"', fixed = TRUE,
               info = "pastille d'ouverture du tiroir absente")
  expect_match(src, "d-flex justify-content-between align-items-center gap-2", fixed = TRUE)
  expect_match(src, '`aria-controls` = ns("cross_drawer")', fixed = TRUE)
  expect_match(src, 'id = ns("cross_drawer")', fixed = TRUE)
  expect_match(src, 'class = "offcanvas offcanvas-end"', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("cross_drawer_summary"))', fixed = TRUE,
               info = "le résumé d'état du tiroir est absent")
  expect_match(src, '`aria-live` = "polite"', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("cross_drawer_nav"))', fixed = TRUE)
  expect_match(src, "output$cross_drawer_nav <- renderUI(", fixed = TRUE,
               info = "le serveur ne rend pas le nav tiroir")
  expect_match(src, "identical(sel, pn$value)", fixed = TRUE,
               info = "le marqueur du panneau courant est absent")
  expect_match(src, "fw-bold", fixed = TRUE)
  expect_match(src, "badge text-bg-secondary", fixed = TRUE,
               info = "le compteur de panneaux par groupe est absent")
  expect_match(src, 'observeEvent(input[[paste0("cross_jump_", pn$value)]], {', fixed = TRUE)
  expect_match(src, 'nav_select("cross_out_tabs", selected = pn$value)', fixed = TRUE)
  expect_match(src, '`aria-label` = .tr_plain("Fermer")', fixed = TRUE,
               info = "aria-label du btn-close sans .tr_plain")
  expect_identical(.n_hits_cross(src, '`aria-label` = i18n$t('), 0L,
                   info = "un aria-label consomme le $t() brut (span i18n dans l'attribut)")
  expect_false(grepl('"Close"', src, fixed = TRUE),
               info = "aria-label anglais hardcoded « Close »")
})

test_that("cross drawer : résumé d'état honnête — disponibilité 8d/8e + comparabilité du contrat", {
  src <- paste(readLines(file.path(ts_project_root(), .CROSS_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, "output$cross_drawer_summary <- renderUI(", fixed = TRUE,
               info = "le résumé d'état n'est pas rendu côté serveur")
  # Disponibilité des DEUX résultats canoniques (consommés, jamais recalculés) :
  expect_match(src, ".has_both()", fixed = TRUE,
               info = "la disponibilité 8d/8e ne passe pas par .has_both()")
  expect_match(src, "shared_rv$da_milo_result", fixed = TRUE,
               info = "disponibilité Milo non consommée")
  expect_match(src, "shared_rv$da_sccoda_result", fixed = TRUE,
               info = "disponibilité scCODA non consommée")
  # Comparabilité DU CONTRAT (build_da_cross_method_summary / .cross_summary) :
  expect_match(src, ".cross_summary()", fixed = TRUE,
               info = "la comparabilité ne consomme pas le résumé du contrat")
  expect_match(src, "fully_comparable", fixed = TRUE,
               info = "le verdict de comparabilité du contrat n'est pas consommé")
  # État PENDING : message d'attente honnête quand un résultat manque (jamais
  # de crash) — même clé que la boîte de comparabilité du panneau :
  expect_match(src, "Lancez Milo (8d) ET scCODA (8e) pour activer la lecture croisée.",
               fixed = TRUE,
               info = "pas de branche attente (crash tiroir avant les 2 résultats)")
})
