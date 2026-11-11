# =============================================================================
# test-mod-sc-da-sccoda-ui.R — modules/sc/mod_sc_da_sccoda.R (tiroir
# « Navigation »)
# =============================================================================
# CCC_NAV_DRAWER vague A jalon 2 (2026-11-11, GO opérateur gate 2 : « aria fix
# + P3 + Vague A ») — généralisation du pattern tiroir optimisé au panneau
# scCODA (8e) — candidat C3 du packet d'arbitrage gate 2 (même famille DA que
# Milo, jalon 1 commité 6c8eb84).
#
# ⚠️ Limite documentée (motif test-mod-sc-communication-ui.R) : le rendu
# headless complet n'est pas réalisable à coût raisonnable ; la preuve est un
# verrou source structurel + identité AST des jumelles.
#
# ⚠️ test-sccoda-contract-freeze.R gèle les top-level de mod_sc_da_sccoda.R à
# exactement 3 fonctions — la carte .sccoda_drawer_groups vit EN LOCAL de
# fonction, en DEUX exemplaires jumelles, comme .milo_drawer_groups.
# =============================================================================

.SCCODA_UI_FILE <- "modules/sc/mod_sc_da_sccoda.R"

.n_hits_sccoda <- function(src, pattern) {
  lengths(regmatches(src, gregexpr(pattern, src, fixed = TRUE)))[1]
}

.SCCODA_DRAWER_VALUES <- c("sccoda_tab_effects", "sccoda_tab_compo",
                           "sccoda_tab_hdi", "sccoda_tab_limites")

test_that("scCODA drawer : navset id + 4 valeurs ASCII explicites (cibles du saut)", {
  src <- paste(readLines(file.path(ts_project_root(), .SCCODA_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Un seul navset_tab, avec un id (cible nav_select du tiroir — motif 2L CCC) :
  expect_identical(.n_hits_sccoda(src, "navset_tab("), 1L)
  expect_match(src, 'id = ns("sccoda_out_tabs")', fixed = TRUE,
               info = "le navset scCODA n'a pas d'id — nav_select ne peut pas cibler")
  for (v in .SCCODA_DRAWER_VALUES) {
    expect_match(src, sprintf('value = "%s"', v), fixed = TRUE, info = v)
  }
})

# Extrait TOUTES les affectations `.sccoda_drawer_groups` du module, À TOUTE
# profondeur (gel test-sccoda-contract-freeze.R — la carte vit en jumelles
# locales, jamais en constante top-level).
.sccoda_extract_drawer_groups <- function() {
  p <- ts_ast_parse(.SCCODA_UI_FILE)
  hits <- ts_ast_find_all(p, function(x) {
    if (!is.call(x) || !identical(x[[1]], as.name("<-"))) return(FALSE)
    lhs <- tryCatch(ts_ast_deparse(x[[2]]), error = function(e) "")
    identical(lhs, ".sccoda_drawer_groups")
  })
  lapply(hits, function(h) eval(h[[3]], envir = new.env(parent = baseenv())))
}

test_that("scCODA drawer : .sccoda_drawer_groups (2 jumelles) = 2 groupes, 4 panneaux, valeurs uniques", {
  twins <- .sccoda_extract_drawer_groups()
  expect_identical(length(twins), 2L,
                   info = "attendu : jumelle output_ui + jumelle serveur, rien d'autre")
  expect_false(is.null(twins[[1]]), info = "jumelle introuvable dans le module")
  # Identité stricte des deux jumelles (aucun drift possible sans échec) :
  expect_identical(twins[[1]], twins[[2]],
                   info = "les deux cartes .sccoda_drawer_groups divergent")
  grp <- twins[[1]]
  expect_identical(length(grp), 2L)
  vals <- unlist(lapply(grp, function(g) vapply(g, function(p) p$value, character(1))))
  labs <- unlist(lapply(grp, function(g) vapply(g, function(p) p$label, character(1))))
  expect_identical(length(vals), 4L)
  expect_identical(anyDuplicated(vals), 0L)
  expect_setequal(vals, .SCCODA_DRAWER_VALUES)
  src <- paste(readLines(file.path(ts_project_root(), .SCCODA_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  for (lb in labs) {
    expect_match(src, sprintf('nav_panel(i18n$t("%s")', lb), fixed = TRUE, info = lb)
  }
})

test_that("scCODA drawer : pastille, .offcanvas, sauts + marqueur actif, compteurs, aria-live, Fermer .tr_plain", {
  src <- paste(readLines(file.path(ts_project_root(), .SCCODA_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, '`data-bs-toggle` = "offcanvas"', fixed = TRUE,
               info = "pastille d'ouverture du tiroir absente")
  expect_match(src, "d-flex justify-content-between align-items-center gap-2", fixed = TRUE)
  expect_match(src, '`aria-controls` = ns("sccoda_drawer")', fixed = TRUE)
  expect_match(src, 'id = ns("sccoda_drawer")', fixed = TRUE)
  expect_match(src, 'class = "offcanvas offcanvas-end"', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("sccoda_drawer_summary"))', fixed = TRUE,
               info = "le résumé d'état du tiroir est absent")
  expect_match(src, '`aria-live` = "polite"', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("sccoda_drawer_nav"))', fixed = TRUE)
  expect_match(src, "output$sccoda_drawer_nav <- renderUI(", fixed = TRUE,
               info = "le serveur ne rend pas le nav tiroir")
  expect_match(src, "identical(sel, pn$value)", fixed = TRUE,
               info = "le marqueur du panneau courant est absent")
  expect_match(src, "fw-bold", fixed = TRUE)
  expect_match(src, "badge text-bg-secondary", fixed = TRUE,
               info = "le compteur de panneaux par groupe est absent")
  expect_match(src, 'observeEvent(input[[paste0("sccoda_jump_", pn$value)]], {', fixed = TRUE)
  expect_match(src, 'nav_select("sccoda_out_tabs", selected = pn$value)', fixed = TRUE)
  # A11y : .tr_plain, JAMAIS de $t() brut en contexte attribut (motif bc9c09a) :
  expect_match(src, '`aria-label` = .tr_plain("Fermer")', fixed = TRUE,
               info = "aria-label du btn-close sans .tr_plain")
  expect_identical(.n_hits_sccoda(src, '`aria-label` = i18n$t('), 0L,
                   info = "un aria-label consomme le $t() brut (span i18n dans l'attribut)")
  expect_false(grepl('"Close"', src, fixed = TRUE),
               info = "aria-label anglais hardcoded « Close »")
})

test_that("scCODA drawer : résumé d'état honnête — statut, contraste, effets crédibles (convention module), référence", {
  src <- paste(readLines(file.path(ts_project_root(), .SCCODA_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, "output$sccoda_drawer_summary <- renderUI(", fixed = TRUE,
               info = "le résumé d'état n'est pas rendu côté serveur")
  # Les métriques DU CONTRAT — le compteur d'effets crédibles suit LA
  # CONVENTION DU MODULE (length(credible_effects), motif du message de
  # verdict sccoda_state) — jamais un seuil réinventé :
  expect_match(src, "length(r$credible_effects", fixed = TRUE,
               info = "le compteur d'effets crédibles n'utilise pas la convention du module")
  expect_match(src, "reference_identity", fixed = TRUE,
               info = "référence d'identité non consommée du contrat")
  expect_match(src, "parameters$target_condition", fixed = TRUE,
               info = "contraste (cible) non consommé du contrat")
  expect_match(src, ".status_label(", fixed = TRUE,
               info = "statut non libellé via .status_label()")
  # État PENDING : le tiroir affiche le statut d'attente (pas de crash) :
  expect_match(src, "if (is.null(r))", fixed = TRUE,
               info = "pas de branche résultat NULL (crash tiroir avant 1er run)")
})
