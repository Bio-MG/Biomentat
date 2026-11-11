# =============================================================================
# test-mod-sc-communication-ui.R — modules/sc/mod_sc_communication.R (onglet
# « Limites »)
# =============================================================================
# Item 2G du packet round-7 (O28 batch B3).
#
# Constat 2G : à HEAD, l'onglet « Limites » empilait 5 bannières `alert-*`
# pleine largeur (1 warning, 2 info, 2 light ; preuve
# mod_sc_communication.R:247-256 à HEAD 88f9a43 ; spéc §2.15). Correctif : UNE
# carte synthétique à puces thématiques, hiérarchie de sévérité conservée via
# la couleur des puces (warning > info > light).
#
# ⚠️ Coordination 3D (packet row 2G) : le TEXTE de la 5ᵉ bannière (objet RDS
# CellChat) est l'objet du constat 3D (batch B2). Ce lot le reproduit VERBATIM
# et le verrou ci-dessous garantit qu'il n'a PAS été reformulé ici.
#
# ⚠️ Limite documentée : mod_sc_communication_ui() embarque les UI de 4
# sous-modules (spatial/trajectory/velocity_ctx/perturbation) — le rendu
# headless complet n'est pas réalisable à coût raisonnable ; la preuve est un
# verrou source structurel (pattern test-mod-sc-velocity.R).
# =============================================================================

.MSCC_UI_FILE <- "modules/sc/mod_sc_communication.R"

.n_hits <- function(src, pattern) {
  lengths(regmatches(src, gregexpr(pattern, src, fixed = TRUE)))[1]
}

test_that("2G : l'onglet Limites est une carte synthétique (header + 5 puces)", {
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Carte synthétique à en-tête (spéc §2.15) :
  expect_match(src, "Points d'attention \u2014 lecture des r\u00e9sultats", fixed = TRUE)
  # Les 5 puces thématiques :
  for (chip in c("Validit\u00e9 technique vs biologique",
                 "Comparabilit\u00e9 des scores",
                 "Causalit\u00e9",
                 "Labels cellulaires",
                 "Objets RDS CellChat")) {
    expect_match(src, chip, fixed = TRUE, info = chip)
  }
})

test_that("2G : plus aucune bannière alert-* empilée — la seule restante est l'en-tête du module", {
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # À HEAD : 6 `class = "alert alert-` (1 bannière d'en-tête L26 + 5 dans
  # l'onglet Limites). Après 2G : il ne reste que la bannière d'en-tête.
  expect_identical(.n_hits(src, "class = \"alert alert-"), 1L,
                   info = "des bannières alert-* subsistent dans l'onglet Limites (constat 2G)")
  # Hiérarchie de sévérité conservée sur les puces (1 warning > 2 info > 2 light) :
  expect_identical(.n_hits(src, "badge bg-warning"), 1L)
  expect_identical(.n_hits(src, "badge bg-info"), 2L)
  expect_identical(.n_hits(src, "badge bg-secondary"), 2L)
})

test_that("2G : les textes d'origine sont préservés VERBATIM (coordination 3D — B2 a corrigé la 5ᵉ)", {
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, "Validite technique uniquement : un import valide n'implique aucune validite biologique.",
               fixed = TRUE)
  expect_match(src, "Les scores ne sont pas compares entre sources differentes (CellChat vs CellPhoneDB).",
               fixed = TRUE)
  expect_match(src, "Une association ligand-receptor n'etablit aucune causalite.",
               fixed = TRUE)
  expect_match(src, "Les labels sans correspondance exacte dans la colonne choisie sont conserves tels quels et listes — jamais renommes.",
               fixed = TRUE)
  # 5ᵉ texte = objet du constat 3D — corrigé par B2 (O28, 2026-10-31) : le
  # verrou B3 « verbatim jamais reformulé » s'appliquait TANT QUE B2 n'était
  # pas passé ; le texte est désormais celui aligné sur le moteur (test 3D
  # ci-dessous), l'ancien texte faux ne doit plus exister nulle part.
  expect_match(src, "appari\u00e9es EXACTEMENT sur @LR$LRsig", fixed = TRUE)
  expect_identical(.n_hits(src, "le mapping ligand-receptor vers pathway n'est pas reconstitue"), 0L,
                   info = "l'ancien texte (drift UI/code) subsiste")
})

# =============================================================================
# B4 (O28, 2026-10-30) — items 2F / 2I / 2L sur le même fichier.
# Preuves mesurées : re-sweep 1920×1080 (RE_SWEEP_1920_2026-10-30.md) —
# 14 sous-onglets en 3 lignes (124 px), zone filtres statique 276 px (L4),
# bas de pile 754 px ; cible B4 ≈ 430 px (mesuré après : 392 px).
# =============================================================================

.B4_CCC_VALUES <- c("ccc_tab_interactions", "ccc_tab_table_filtres", "ccc_tab_dotplot",
                    "ccc_tab_heatmap", "ccc_tab_reseau", "ccc_tab_centralite",
  "ccc_tab_comparaison",
                    "ccc_tab_spatial", "ccc_tab_trajectoire", "ccc_tab_velocite",
                    "ccc_tab_perturbation", "ccc_tab_schema", "ccc_tab_harmonisation",
                    "ccc_tab_qc_import", "ccc_tab_limites")

test_that("2F : les 14 sous-onglets CCC sont regroupés en 4 nav_menu thématiques (spéc §2.14)", {
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # 4 groupes thématiques, libellés de la spéc §2.14 :
  for (grp in c("Exploration des interactions", "Contextes multi-modaux",
                "Traçabilité & qualité", "Documentation & méthodes")) {
    expect_match(src, sprintf('nav_menu(i18n$t("%s")', grp), fixed = TRUE, info = grp)
  }
  # Les 14 panneaux portent désormais des `value` ASCII explicites (cibles du
  # saut 2L, indépendantes de la traduction) :
  for (v in .B4_CCC_VALUES) {
    expect_match(src, sprintf('value = "%s"', v), fixed = TRUE, info = v)
  }
  # Un seul navset_tab, avec un id (cible nav_select du tiroir 2L) :
  expect_identical(.n_hits(src, "navset_tab("), 1L)
  expect_match(src, 'id = ns("comm_tabs")', fixed = TRUE)
  # Témoin de non-régression 2G : la carte « Limites » reste dans son panel
  # (déplacé, jamais reformulé — les textes verbatim sont verrouillés plus haut).
})

# Extrait TOUTES les affectations `.comm_drawer_groups` du module, À TOUTE
# profondeur (le gel communication interdit toute affectation top-level hors
# les 3 fonctions d'orchestration — la carte vit donc EN LOCAL de fonction, en
# deux exemplaires jumelles : rendu tiroir dans output_ui, observateurs de
# saut dans le serveur).
.B4_extract_drawer_groups <- function() {
  p <- ts_ast_parse(.MSCC_UI_FILE)
  hits <- ts_ast_find_all(p, function(x) {
    if (!is.call(x) || !identical(x[[1]], as.name("<-"))) return(FALSE)
    lhs <- tryCatch(ts_ast_deparse(x[[2]]), error = function(e) "")
    identical(lhs, ".comm_drawer_groups")
  })
  lapply(hits, function(h) eval(h[[3]], envir = new.env(parent = baseenv())))
}

test_that("2F/2L : .comm_drawer_groups (2 jumelles) = 4 groupes, 14 panneaux, valeurs uniques", {
  twins <- .B4_extract_drawer_groups()
  # Deux exemplaires (UI + serveur), aucune constante top-level :
  expect_identical(length(twins), 2L,
                   info = "attendu : jumelle output_ui + jumelle serveur, rien d'autre")
  expect_false(is.null(twins[[1]]), info = "jumelle introuvable dans le module")
  # Identité stricte des deux jumelles (aucun drift possible sans échec) :
  expect_identical(twins[[1]], twins[[2]],
                   info = "les deux cartes .comm_drawer_groups divergent")
  grp <- twins[[1]]
  expect_identical(length(grp), 4L)
  vals <- unlist(lapply(grp, function(g) vapply(g, function(p) p$value, character(1))))
  labs <- unlist(lapply(grp, function(g) vapply(g, function(p) p$label, character(1))))
  expect_identical(length(vals), 15L)  # CC-2 (2026-11-09) : onglet Comparaison
  expect_identical(anyDuplicated(vals), 0L)
  expect_setequal(vals, .B4_CCC_VALUES)
  # Chaque label a une occurrence nav_panel i18n$t("<label>") dans le navset :
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  for (lb in labs) {
    expect_match(src, sprintf('nav_panel(i18n$t("%s")', lb), fixed = TRUE, info = lb)
  }
})

test_that("B4/freeze : le module ne définit toujours que ses 3 fonctions top-level", {
  # Garde de COORDINATION B4 (verrou source : test-communication-contract-freeze.R:140-145,
  # helper :23-35 — aucune affectation top-level hors les 3 fonctions
  # d'orchestration ; la carte 2F/2L vit en local de fonction, en jumelles).
  exprs <- parse(file.path(ts_project_root(), .MSCC_UI_FILE), keep.source = FALSE)
  nms <- character(0)
  for (e in exprs) {
    if (is.call(e) && identical(deparse(e[[1L]]), "<-") && is.symbol(e[[2L]])) {
      nms <- c(nms, as.character(e[[2L]]))
    }
  }
  expect_setequal(unique(nms),
                  c("mod_sc_communication_ui", "mod_sc_communication_output_ui",
                    "mod_sc_communication_server"))
})

test_that("2I : la zone de filtres statique est un accordéon FERMÉ par défaut (spéc §2.20)", {
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # L'accordéon existe, fermé au boot (L4 : 276 px → ~54 px mesurés après) :
  expect_match(src, 'id = ns("comm_filters_accordion")', fixed = TRUE)
  expect_match(src, "open = FALSE", fixed = TRUE)
  expect_match(src, "i18n$t(\"Filtres d'interaction\")", fixed = TRUE)
  # L'ancien bloc statique a disparu :
  expect_identical(.n_hits(src, 'class = "border rounded p-2 m-2"'), 0L,
                   info = "la zone de filtres statique (L4) subsiste")
  # Les 7 inputs de filtres sont TOUJOURS câblés (mêmes ids, une instance) :
  for (id in c("comm_f_score_min", "comm_f_pmax", "comm_f_pathways", "comm_f_rank_max",
               "comm_f_senders", "comm_f_receivers", "comm_f_include_self")) {
    expect_identical(.n_hits(src, paste0('ns("', id, '")')), 1L,
                     info = paste0("input filtre ", id, " absent ou dupliqué"))
  }
  # Le résumé de filtrage est factorisé : une source réactive, deux sorties
  # (ligne sous les filtres 2I = phrase + tiroir 2L = liste compacte P3,
  # CCC-DRAW-OPT suite 2026-11-11 — la phrase n'a PLUS qu'UNE sortie).
  expect_match(src, ".filter_summary_txt <- reactive(", fixed = TRUE)
  expect_identical(.n_hits(src, "renderText({ .filter_summary_txt() })"), 1L)
})

test_that("2L : tiroir offcanvas secondaire — pastille, .offcanvas, sauts, zéro duplication (spéc §2.23)", {
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Déclencheur (pastille DANS l'en-tête de carte — aucune rangée verticale
  # ajoutée, le bas de pile 2I ne doit pas regonfler) :
  expect_match(src, '`data-bs-toggle` = "offcanvas"', fixed = TRUE)
  expect_match(src, "d-flex justify-content-between align-items-center gap-2", fixed = TRUE)
  # Le tiroir Bootstrap 5 natif, à droite :
  expect_match(src, 'id = ns("comm_drawer")', fixed = TRUE)
  expect_match(src, 'class = "offcanvas offcanvas-end"', fixed = TRUE)
  # Sauts directs : un observer par panneau, ciblant les `value` posés par 2F ;
  # le tiroir reste ouvert après le saut (accès aux panneaux sans replier le
  # workflow — verdict re-sweep 1920).
  expect_match(src, 'observeEvent(input[[paste0("comm_jump_", pn$value)]], {', fixed = TRUE)
  expect_match(src, 'nav_select("comm_tabs", selected = pn$value)', fixed = TRUE)
  # ZÉRO duplication des filtres (arbitrage 2I : instance unique dans
  # l'accordéon ; le tiroir ne porte que le résumé + la navigation). Compte
  # mesuré : 1 UI (accordéon) + 2 updateNumericInput (reset état / échantillon)
  # + 1 input$ (comm_filtered) + 4 updateNumericInput 3F (O28, 2026-10-31 :
  # normalisation contextuelle score/rang par source — spéc §2.18 point 4 ;
  # branches NULL, LIANA, CellPhoneDB, défaut CellChat) — coordonné B4↔B2 :
  # chaque occurrence sert UN des deux rôles (rendu / reset / consommation /
  # adaptation par source), aucune n'en duplique un autre.
  expect_identical(.n_hits(src, 'comm_f_score_min'), 8L,
                   info = "occurrences comm_f_score_min (1 UI + 2 resets + 1 input$ + 4 updates 3F)")
})

# =============================================================================
# CCC-DRAW-OPT (2026-11-11, boost mode) — optimisation du tiroir AVANT toute
# généralisation (gate 1 de CCC_NAV_DRAWER_GENERALIZATION, auto-revue datée
# docs/ideation/CCC_NAV_DRAWER_AUTOREVIEW_2026-11-11.md) : P1 marqueur du
# panneau courant, P2 compteur par groupe, P4 aria-live + libellé Fermer i18n.
# P3 (résumé compact) DIFFÉRÉ — noté dans l'auto-revue.
# =============================================================================

test_that("CCC-DRAW-OPT : tiroir optimisé — nav dynamique (marqueur actif), compteurs, aria-live, Fermer i18n", {
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # P1 : la navigation du tiroir est rendue DYNAMIQUEMENT (uiOutput + renderUI)
  # et marque le panneau courant (bold) — l'utilisateur sait où il est.
  expect_match(src, 'uiOutput(ns("comm_drawer_nav"))', fixed = TRUE,
               info = "P1 : le rendu dynamique du nav tiroir est absent")
  expect_match(src, "output$comm_drawer_nav <- renderUI(", fixed = TRUE,
               info = "P1 : le serveur ne rend pas le nav tiroir")
  expect_match(src, "identical(sel, pn$value)", fixed = TRUE,
               info = "P1 : le marqueur du panneau courant est absent")
  expect_match(src, "fw-bold", fixed = TRUE,
               info = "P1 : le style du marqueur actif est absent")
  # P2 : compteur de panneaux par groupe (badge)
  # P2 : compteur de panneaux par groupe — style BS 5.3 text-bg-* DELIBÉRÉMENT
  # distinct des badges 2G (« badge bg-secondary », compteur épingle à 2 au
  # bloc 2G) : deux rôles, deux classes, zéro collision de pin.
  expect_match(src, "badge text-bg-secondary", fixed = TRUE,
               info = "P2 : le compteur de panneaux par groupe est absent")
  # P4 : aria-live sur le résumé du tiroir (SR annoncent les changements de
  # filtrage) + libellé Fermer i18n. RE-ÉPINGLÉ CCC-DRAW-OPT suite (2026-11-11)
  # : le raw $t() en contexte ATTRIBUT posait le span i18n brut dans
  # l'attribut (BUG A11Y reporté STATUS §0 CQF-V ; même famille que le
  # boot-fail CQF-3 173c43e) — le libellé passe par .tr_plain (strip span).
  expect_match(src, 'aria-live', fixed = TRUE,
               info = "P4 : aria-live absent du tiroir")
  expect_match(src, '`aria-label` = .tr_plain("Fermer")', fixed = TRUE,
               info = "le aria-label du btn-close ne passe pas par .tr_plain")
  expect_false(grepl('"Close"', src, fixed = TRUE),
               info = "P4 : le aria-label anglais hardcoded « Close » persiste")
})

# =============================================================================
# CCC-DRAW-OPT SUITE (2026-11-11, GO opérateur gate 2 : « aria fix + P3 +
# Vague A ») — FIX BUG A11Y Fermer + P3 résumé compact. Fix : dans un
# ATTRIBUT, le $t() BUILD-TIME renvoie un shiny.tag span (le shim JS i18n ne
# passe pas sur les attributs) → l'attribut portait le HTML brut
# <span class="i18n" ...> au lieu du texte « Fermer ». P3 (auto-revue O3) :
# le résumé du tiroir passe de la phrase sprintf dense (8 métriques en une
# ligne, wrap dense à 400 px) à une liste compacte 2 colonnes — la ligne sous
# les filtres 2I garde la phrase (sortie inchangée, P4 aria-live survit sur
# le CONTAINER du tiroir).
# =============================================================================
test_that("CCC-DRAW-OPT suite : aria-label Fermer texte plat (.tr_plain) + P3 liste compacte 2 colonnes", {
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # FIX A11Y : le btn-close consomme .tr_plain (strip span, motif global.R /
  # 173c43e) — plus JAMAIS de $t() brut en contexte attribut dans ce fichier :
  expect_identical(.n_hits(src, '`aria-label` = i18n$t('), 0L,
                   info = "un aria-label consomme encore le $t() brut (span i18n dans l'attribut)")
  # P3 : le tiroir rend une LISTE compacte (uiOutput + renderUI serveur,
  # tags$ul 2 colonnes) ; la phrase sprintf reste UNIQUEMENT sur la ligne 2I.
  expect_match(src, 'uiOutput(ns("comm_filter_summary_drawer"))', fixed = TRUE,
               info = "P3 : le tiroir ne rend pas le résumé en uiOutput")
  expect_match(src, "output$comm_filter_summary_drawer <- renderUI(", fixed = TRUE,
               info = "P3 : le résumé compact du tiroir n'est pas rendu côté serveur")
  expect_match(src, 'class = "list-unstyled', fixed = TRUE,
               info = "P3 : la liste compacte est absente")
  expect_match(src, "column-count:2", fixed = TRUE,
               info = "P3 : les 2 colonnes compactes sont absentes")
  # Les 8 métriques du résumé sont TOUJOURS portées (labels i18n listés) :
  for (lab in c("Interactions affich\u00e9es", "Retir\u00e9es par score",
                "Retir\u00e9es par rang", "Retir\u00e9es par p-value",
                "Retir\u00e9es par pathway", "Retir\u00e9es par sender",
                "Retir\u00e9es par receiver", "Retir\u00e9es \u2014 auto-interactions")) {
    expect_match(src, paste0('.tr("', lab, '")'), fixed = TRUE,
                 info = paste0("P3 : label résumé absent — ", lab))
  }
  # La ligne 2I garde la phrase (sortie inchangée) :
  expect_match(src, "output$comm_filter_summary <- renderText({ .filter_summary_txt() })",
               fixed = TRUE)
  # aria-live SURVIT à P3 (container du tiroir) :
  expect_match(src, '`aria-live` = "polite"', fixed = TRUE,
               info = "P4 : aria-live perdu au passage P3")
})


# =============================================================================
# B2 (O28, 2026-10-31) — item 3D : la bannière « Limites » est alignée sur le
# moteur. Constat 3D (packet round-7 ; correctif 3D RÉFUTÉ round 2 jamais
# posé) : le texte affirmait « le mapping ligand-receptor vers pathway n'est
# pas reconstitue » alors que le moteur le résout par appariement EXACT sur
# @LR$LRsig (sc_communication_engine.R:260-268) et remonte un engine_failure
# explicite quand @LR$LRsig est absent (sc_communication_engine.R:238-243).
# Formulation de la spéc §2.16 (tranchée à l'implémentation, comme autorisé).
# ⚠️ B4 (pastille offcanvas de l'en-tête) et B5 (surface « Harmonisation
# identites ») restent intacts — la correction ne touche que le TEXTE de la
# 5ᵉ puce de l'onglet « Limites » (ccc_tab_limites).
# =============================================================================
test_that("3D : la bannière RDS CellChat est alignée sur le moteur (exact match + engine_failure)", {
  src <- paste(readLines(file.path(ts_project_root(), .MSCC_UI_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Le moteur résout ligand/receptor/pathway par appariement EXACT :
  expect_match(src, "appari\u00e9es EXACTEMENT sur @LR$LRsig (ligand, receptor, pathway_name)",
               fixed = TRUE)
  expect_match(src, "aucune reconstruction heuristique", fixed = TRUE)
  # Le comportement d'échec documenté du moteur est ÉNONCÉ (pas nié) :
  expect_match(src, "engine_failure", fixed = TRUE)
  # Le texte faux (prémisse réfutée round 2) a disparu de TOUT le fichier :
  expect_identical(.n_hits(src, "n'est pas reconstitue"), 0L,
                   info = "l'ancien texte contredisant le moteur subsiste")
})

# ---------------------------------------------------------------------------
# P1-6b (2026-11-07, plan post-E2E #6) : LIANA par échantillon → worker mirai
# (mesuré E2E : >2h15min synchrone, session bloquée). Garde structurelle du
# câblage — même famille que les gardes P1-6a du pipeline manuel.
# ---------------------------------------------------------------------------
test_that("P1-6b : LIANA délégué à un worker mirai, repli synchrone tracé", {
  src <- paste(readLines(file.path(ts_project_root(), "modules", "sc",
                                   "mod_sc_communication.R"),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # Tâche async + pool unique + bouton lié :
  expect_match(src, "ext_task_cls_liana$new(function(obj, args, deps = liana_core_deps)",
               fixed = TRUE, info = "tâche ExtendedTask LIANA absente")
  expect_match(src, "ensure_app_daemons()", fixed = TRUE,
               info = "le pool applicatif unique n'est pas armé")
  expect_match(src, "run_liana_by_sample(\n              obj,", fixed = TRUE,
               info = "le moteur pur ne tourne pas en daemon")
  expect_match(src, 'bslib::bind_task_button(t, "comm_compute_liana")', fixed = TRUE,
               info = "bouton compute LIANA non lié à la tâche")
  expect_match(src, 'identical(liana_task$status(), "running")', fixed = TRUE,
               info = "pas de garde re-entrance sur le calcul LIANA")
  # Déps daemon embarquées avec la tâche (moteur NON dans le preload) :
  expect_match(src, '"sc_communication_engine.R"', fixed = TRUE)
  # REPLI SYNCHRONE TRACÉ (jamais muet) :
  expect_match(src, "Worker mirai indisponible", fixed = TRUE)
  # L'application du résultat est FACTORISÉE (async et sync partagent) :
  expect_match(src, ".apply_liana_result <- function(collection, obj)", fixed = TRUE)
  expect_match(src, ".apply_liana_result(liana_task$result(), global_data$sc_obj)",
               fixed = TRUE, info = "l'achèvement async n'applique pas le résultat")
  # Le moteur n'est PAS re-câblé en dur dans l'observer sync principal :
  expect_match(src, "liana_task$invoke(obj, args)", fixed = TRUE)
})
