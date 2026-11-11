# =============================================================================
# test-mod-sc-velocity-ui.R — modules/sc/mod_sc_velocity.R (tiroir
# « Navigation »)
# =============================================================================
# CCC_NAV_DRAWER vague B jalon J3 (2026-11-12, GO opérateur « GO J3 + J4 »)
# — généralisation du pattern tiroir optimisé au panneau Velocity (C1 du
# packet gate 2, 4 vues nettes). J4 commité AVANT (997bc14) car ce fichier
# est EN VOL lane scvelo (fix opt-in loom, session concurrente) — le commit
# J3 stage SES SEULS HUNKS (le hunk étranger loom reste NON commité, jamais
# touché).
#
# ⚠️ Limite documentée (motif test-mod-sc-communication-ui.R) : verrou source
# structurel + identité AST des jumelles.
#
# ⚠️ test-velocity-contract-freeze.R:152 gèle les top-level de
# mod_sc_velocity.R à exactement 3 fonctions — la carte vit EN LOCAL de
# fonction, en DEUX exemplaires jumelles.
# =============================================================================

.VELO_UI_FILE <- "modules/sc/mod_sc_velocity.R"

.n_hits_velo <- function(src, pattern) {
  lengths(regmatches(src, gregexpr(pattern, src, fixed = TRUE)))[1]
}

.VELO_DRAWER_VALUES <- c("velocity_tab_phase", "velocity_tab_vectors",
                         "velocity_tab_coverage", "velocity_tab_qc")

test_that("velocity drawer : navset id + 4 valeurs ASCII explicites (cibles du saut)", {
  src <- paste(readLines(file.path(ts_project_root(), .VELO_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_identical(.n_hits_velo(src, "navset_tab("), 1L)
  expect_match(src, 'id = ns("velocity_out_tabs")', fixed = TRUE,
               info = "le navset Velocity n'a pas d'id — nav_select ne peut pas cibler")
  for (v in .VELO_DRAWER_VALUES) {
    expect_match(src, sprintf('value = "%s"', v), fixed = TRUE, info = v)
  }
})

# Extrait TOUTES les affectations `.velocity_drawer_groups` du module, À TOUTE
# profondeur (gel test-velocity-contract-freeze.R — la carte vit en jumelles
# locales, jamais en constante top-level).
.velo_extract_drawer_groups <- function() {
  p <- ts_ast_parse(.VELO_UI_FILE)
  hits <- ts_ast_find_all(p, function(x) {
    if (!is.call(x) || !identical(x[[1]], as.name("<-"))) return(FALSE)
    lhs <- tryCatch(ts_ast_deparse(x[[2]]), error = function(e) "")
    identical(lhs, ".velocity_drawer_groups")
  })
  lapply(hits, function(h) eval(h[[3]], envir = new.env(parent = baseenv())))
}

test_that("velocity drawer : .velocity_drawer_groups (2 jumelles) = 1 groupe, 4 panneaux, valeurs uniques", {
  twins <- .velo_extract_drawer_groups()
  expect_identical(length(twins), 2L,
                   info = "attendu : jumelle output_ui + jumelle serveur, rien d'autre")
  expect_false(is.null(twins[[1]]), info = "jumelle introuvable dans le module")
  expect_identical(twins[[1]], twins[[2]],
                   info = "les deux cartes .velocity_drawer_groups divergent")
  grp <- twins[[1]]
  expect_identical(length(grp), 1L)
  vals <- unlist(lapply(grp, function(g) vapply(g, function(p) p$value, character(1))))
  labs <- unlist(lapply(grp, function(g) vapply(g, function(p) p$label, character(1))))
  expect_identical(length(vals), 4L)
  expect_identical(anyDuplicated(vals), 0L)
  expect_setequal(vals, .VELO_DRAWER_VALUES)
  src <- paste(readLines(file.path(ts_project_root(), .VELO_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  for (lb in labs) {
    expect_match(src, sprintf('nav_panel(i18n$t("%s")', lb), fixed = TRUE, info = lb)
  }
})

test_that("velocity drawer : pastille, .offcanvas, sauts + marqueur actif, compteurs, aria-live, Fermer .tr_plain", {
  src <- paste(readLines(file.path(ts_project_root(), .VELO_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, '`data-bs-toggle` = "offcanvas"', fixed = TRUE,
               info = "pastille d'ouverture du tiroir absente")
  expect_match(src, "d-flex justify-content-between align-items-center gap-2", fixed = TRUE)
  expect_match(src, '`aria-controls` = ns("velocity_drawer")', fixed = TRUE)
  expect_match(src, 'id = ns("velocity_drawer")', fixed = TRUE)
  expect_match(src, 'class = "offcanvas offcanvas-end"', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("velocity_drawer_summary"))', fixed = TRUE,
               info = "le résumé d'état du tiroir est absent")
  expect_match(src, '`aria-live` = "polite"', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("velocity_drawer_nav"))', fixed = TRUE)
  expect_match(src, "output$velocity_drawer_nav <- renderUI(", fixed = TRUE,
               info = "le serveur ne rend pas le nav tiroir")
  expect_match(src, "identical(sel, pn$value)", fixed = TRUE,
               info = "le marqueur du panneau courant est absent")
  expect_match(src, "fw-bold", fixed = TRUE)
  expect_match(src, "badge text-bg-secondary", fixed = TRUE,
               info = "le compteur de panneaux par groupe est absent")
  expect_match(src, 'observeEvent(input[[paste0("velocity_jump_", pn$value)]], {', fixed = TRUE)
  expect_match(src, 'nav_select("velocity_out_tabs", selected = pn$value)', fixed = TRUE)
  expect_match(src, '`aria-label` = .tr_plain("Fermer")', fixed = TRUE,
               info = "aria-label du btn-close sans .tr_plain")
  expect_identical(.n_hits_velo(src, '`aria-label` = i18n$t('), 0L,
                   info = "un aria-label consomme le $t() brut (span i18n dans l'attribut)")
  expect_false(grepl('"Close"', src, fixed = TRUE),
               info = "aria-label anglais hardcoded « Close »")
})

test_that("velocity drawer : résumé d'état honnête — statut vivant + vecteurs opt-in", {
  src <- paste(readLines(file.path(ts_project_root(), .VELO_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, "output$velocity_drawer_summary <- renderUI(", fixed = TRUE,
               info = "le résumé d'état n'est pas rendu côté serveur")
  # Statut VIVANT (même source que la ligne de statut du panneau — zéro
  # divergence) :
  expect_match(src, "velocity_status_rv()", fixed = TRUE,
               info = "le statut du tiroir n'est pas la source vivante du panneau")
  # Disponibilité des vecteurs opt-in consommée de l'état (jamais inventée) :
  expect_match(src, "velocity_state$umap_vectors", fixed = TRUE,
               info = "la disponibilité vecteurs UMAP opt-in n'est pas consommée de l'état")
})
