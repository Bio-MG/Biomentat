# =============================================================================
# test-mod-sc-da-milo-ui.R — modules/sc/mod_sc_da_milo.R (tiroir « Navigation »)
# =============================================================================
# CCC_NAV_DRAWER vague A jalon 1 (2026-11-11, GO opérateur gate 2 : « aria fix
# + P3 + Vague A ») — généralisation du pattern optimisé du tiroir CCC
# (5f0a906 P1 marqueur actif + P2 compteurs + P4 aria-live ; suite bc9c09a fix
# aria-label + P3 liste compacte) au panneau Milo (8d) — candidat C2 du packet
# d'arbitrage gate 2 (docs/ideation/CCC_NAV_DRAWER_GATE2_ARBITRATION_PACKET_2026-11-11.md :
# profondeur maximale mesurée, même famille DA que scCODA).
#
# ⚠️ Limite documentée (motif test-mod-sc-communication-ui.R) : le rendu
# headless complet n'est pas réalisable à coût raisonnable ; la preuve est un
# verrou source structurel + identité AST des jumelles.
#
# ⚠️ test-milo-contract-freeze.R:111 gèle les top-level de mod_sc_da_milo.R à
# exactement 3 fonctions — la carte .milo_drawer_groups vit EN LOCAL de
# fonction, en DEUX exemplaires jumelles (rendu tiroir output_ui + nav serveur),
# comme .comm_drawer_groups.
# =============================================================================

.MILO_UI_FILE <- "modules/sc/mod_sc_da_milo.R"

.n_hits_milo <- function(src, pattern) {
  lengths(regmatches(src, gregexpr(pattern, src, fixed = TRUE)))[1]
}

.MILO_DRAWER_VALUES <- c("milo_tab_embed", "milo_tab_logfc", "milo_tab_compo",
                         "milo_tab_table", "milo_tab_limites")

test_that("Milo drawer : navset id + 5 valeurs ASCII explicites (cibles du saut)", {
  src <- paste(readLines(file.path(ts_project_root(), .MILO_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Un seul navset_tab, avec un id (cible nav_select du tiroir — motif 2L CCC) :
  expect_identical(.n_hits_milo(src, "navset_tab("), 1L)
  expect_match(src, 'id = ns("milo_out_tabs")', fixed = TRUE,
               info = "le navset Milo n'a pas d'id — nav_select ne peut pas cibler")
  # Les 5 panneaux portent des `value` ASCII explicites :
  for (v in .MILO_DRAWER_VALUES) {
    expect_match(src, sprintf('value = "%s"', v), fixed = TRUE, info = v)
  }
})

# Extrait TOUTES les affectations `.milo_drawer_groups` du module, À TOUTE
# profondeur (gel test-milo-contract-freeze.R:111 — la carte vit en jumelles
# locales, jamais en constante top-level).
.milo_extract_drawer_groups <- function() {
  p <- ts_ast_parse(.MILO_UI_FILE)
  hits <- ts_ast_find_all(p, function(x) {
    if (!is.call(x) || !identical(x[[1]], as.name("<-"))) return(FALSE)
    lhs <- tryCatch(ts_ast_deparse(x[[2]]), error = function(e) "")
    identical(lhs, ".milo_drawer_groups")
  })
  lapply(hits, function(h) eval(h[[3]], envir = new.env(parent = baseenv())))
}

test_that("Milo drawer : .milo_drawer_groups (2 jumelles) = 2 groupes, 5 panneaux, valeurs uniques", {
  twins <- .milo_extract_drawer_groups()
  expect_identical(length(twins), 2L,
                   info = "attendu : jumelle output_ui + jumelle serveur, rien d'autre")
  expect_false(is.null(twins[[1]]), info = "jumelle introuvable dans le module")
  # Identité stricte des deux jumelles (aucun drift possible sans échec) :
  expect_identical(twins[[1]], twins[[2]],
                   info = "les deux cartes .milo_drawer_groups divergent")
  grp <- twins[[1]]
  expect_identical(length(grp), 2L)
  vals <- unlist(lapply(grp, function(g) vapply(g, function(p) p$value, character(1))))
  labs <- unlist(lapply(grp, function(g) vapply(g, function(q) q$label, character(1))))
  expect_identical(length(vals), 5L)
  expect_identical(anyDuplicated(vals), 0L)
  expect_setequal(vals, .MILO_DRAWER_VALUES)
  # Chaque label a une occurrence nav_panel i18n$t("<label>") dans le navset :
  src <- paste(readLines(file.path(ts_project_root(), .MILO_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  for (lb in labs) {
    expect_match(src, sprintf('nav_panel(i18n$t("%s")', lb), fixed = TRUE, info = lb)
  }
})

test_that("Milo drawer : pastille, .offcanvas, sauts + marqueur actif, compteurs, aria-live, Fermer .tr_plain", {
  src <- paste(readLines(file.path(ts_project_root(), .MILO_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Déclencheur (pastille DANS l'en-tête de carte — motif 2L CCC, aucune
  # rangée verticale ajoutée) :
  expect_match(src, '`data-bs-toggle` = "offcanvas"', fixed = TRUE,
               info = "pastille d'ouverture du tiroir absente")
  expect_match(src, "d-flex justify-content-between align-items-center gap-2", fixed = TRUE)
  expect_match(src, '`aria-controls` = ns("milo_drawer")', fixed = TRUE)
  # Le tiroir Bootstrap 5 natif, à droite :
  expect_match(src, 'id = ns("milo_drawer")', fixed = TRUE)
  expect_match(src, 'class = "offcanvas offcanvas-end"', fixed = TRUE)
  # Résumé d'état (uiOutput, aria-live sur le CONTAINER) :
  expect_match(src, 'uiOutput(ns("milo_drawer_summary"))', fixed = TRUE,
               info = "le résumé d'état du tiroir est absent")
  expect_match(src, '`aria-live` = "polite"', fixed = TRUE)
  # Navigation dynamique (marqueur « où suis-je ? » + compteurs) :
  expect_match(src, 'uiOutput(ns("milo_drawer_nav"))', fixed = TRUE)
  expect_match(src, "output$milo_drawer_nav <- renderUI(", fixed = TRUE,
               info = "le serveur ne rend pas le nav tiroir")
  expect_match(src, "identical(sel, pn$value)", fixed = TRUE,
               info = "le marqueur du panneau courant est absent")
  expect_match(src, "fw-bold", fixed = TRUE)
  expect_match(src, "badge text-bg-secondary", fixed = TRUE,
               info = "le compteur de panneaux par groupe est absent")
  # Sauts directs + tiroir reste ouvert (motif 2L CCC verbatim) :
  expect_match(src, 'observeEvent(input[[paste0("milo_jump_", pn$value)]], {', fixed = TRUE)
  expect_match(src, 'nav_select("milo_out_tabs", selected = pn$value)', fixed = TRUE)
  # A11y : aria-label du btn-close passe par .tr_plain (JAMAIS de $t() brut en
  # contexte attribut — BUG A11Y « Fermer » fixé par bc9c09a sur Communication) :
  expect_match(src, '`aria-label` = .tr_plain("Fermer")', fixed = TRUE,
               info = "aria-label du btn-close sans .tr_plain")
  expect_identical(.n_hits_milo(src, '`aria-label` = i18n$t('), 0L,
                   info = "un aria-label consomme le $t() brut (span i18n dans l'attribut)")
  expect_false(grepl('"Close"', src, fixed = TRUE),
               info = "aria-label anglais hardcoded « Close »")
})

test_that("Milo drawer : résumé d'état honnête — statut, contraste, voisinages, couverture", {
  src <- paste(readLines(file.path(ts_project_root(), .MILO_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Le résumé est rendu côté serveur :
  expect_match(src, "output$milo_drawer_summary <- renderUI(", fixed = TRUE,
               info = "le résumé d'état n'est pas rendu côté serveur")
  # Les 4 métriques du contrat (assert_milo_result) — JAMAIS de métrique
  # inventée (pas de seuil de significativité fabriqué) :
  expect_match(src, "neighbourhood_summary", fixed = TRUE,
               info = "n voisinages non consommé du contrat")
  expect_match(src, "fraction_cells_in_nhoods", fixed = TRUE,
               info = "couverture cellules non consommée du contrat")
  expect_match(src, "tested_contrast", fixed = TRUE,
               info = "contraste testé non consommé du contrat")
  expect_match(src, ".status_label(", fixed = TRUE,
               info = "statut non libellé via .status_label()")
  # État PENDING : le tiroir affiche le statut d'attente (pas de crash sur
  # résultat NULL) :
  expect_match(src, "if (is.null(r))", fixed = TRUE,
               info = "pas de branche résultat NULL (crash tiroir avant 1er run)")
})
