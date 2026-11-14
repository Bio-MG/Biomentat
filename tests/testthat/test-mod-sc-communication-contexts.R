# =============================================================================
# test-mod-sc-communication-contexts.R — les 3 modules « contextes » CCC
# (mod_sc_communication_velocity.R / _trajectory.R / _perturbation.R) + la
# normalisation des filtres score/rang du parent.
# =============================================================================
# Item 3F du packet round-7 (O28 batch B2, 2026-10-31 — spéc §2.18).
#
# Constat 3F : gardes RÉACTIVES existantes (bandeau vel_info :58-61, stop :84)
# mais AUCUNE garde PRÉVENTIVE (disable/bandeau avant clic) et 0 withProgress
# dans les trois modules (git grep mesuré à HEAD 7550e39, packet row 3F).
#
# ⚠️ Limite documentée (même motif que test-mod-sc-communication-ui.R) : les
# gardes vivent DANS les moduleServer — le rendu headless complet des
# sous-modules n'est pas réalisable à coût raisonnable (comm_state, objet
# Seurat, collection LIANA). La preuve est donc un verrou source structurel :
# chaque garde est verrouillée par sa signature de code, et la sémantique
# (états utilisables du contrat velocity) est citée depuis la source de
# vérité (sc_velocity.R:141-151).
# =============================================================================

.VEL_FILE  <- "modules/sc/mod_sc_communication_velocity.R"
.TRJ_FILE  <- "modules/sc/mod_sc_communication_trajectory.R"
.SPATIAL_FILE <- "modules/sc/mod_sc_communication_spatial.R"
.PERT_FILE <- "modules/sc/mod_sc_communication_perturbation.R"
.PAR_FILE  <- "modules/sc/mod_sc_communication.R"

.read <- function(rel) {
  paste(readLines(file.path(ts_project_root(), rel),
                  warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}
.n_hits <- function(src, pattern) {
  lengths(regmatches(src, gregexpr(pattern, src, fixed = TRUE)))[1]
}

# ── Vélocité : garde préventive (spéc §2.18 point 1) ─────────────────────────
test_that("3F vélocité : bouton désactivé préventivement + bandeau — états du contrat consommés", {
  src <- .read(.VEL_FILE)
  # Le bouton est piloté par updateActionButton(disabled=) sur le résultat
  # velocity (le stop post-clic :84 reste en seconde ligne) :
  expect_match(src, 'updateActionButton(session, "vel_compute", disabled',
               fixed = TRUE)
  # États utilisables = ceux du CONTRAT (source de vérité sc_velocity.R) —
  # valid_no_vectors et les invalid_* restent bloqués :
  expect_match(src, 'st %in% c("valid", "valid_partial_embedding")', fixed = TRUE)
  # Bandeau préventif dédié (l'infobulle vel_info seule ne désactive rien) :
  expect_match(src, 'output$vel_guard <- renderUI', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("vel_guard"))', fixed = TRUE)
  # Texte de la spéc (valid_no_vectors / absent) :
  expect_match(src,
    "Le calcul de contexte n\u00e9cessite des vecteurs de v\u00e9locit\u00e9 pr\u00e9calcul\u00e9s complets (absents dans l'import actuel).",
    fixed = TRUE)
})

test_that("3F vélocité : jauge withProgress autour du calcul de contexte", {
  src <- .read(.VEL_FILE)
  expect_match(src, 'withProgress(message = .tr("Calcul du contexte v\u00e9locit\u00e9...")',
               fixed = TRUE)
  # Le calcul est BIEN encapsulé (build dans le bloc withProgress) :
  expect_match(src, "build_communication_velocity_context(", fixed = TRUE)
})

# ── Trajectoire : garde préventive (spéc §2.18 point 2) ──────────────────────
test_that("3F trajectoire : bascule auto vers global sans slingshot + bouton désactivé sans pseudo-temps", {
  src <- .read(.TRJ_FILE)
  # Détecteur slingshot + bascule automatique (option documentée de la spéc) :
  expect_match(src, 'grepl("^pseudotime_slingshot_", meta_cols)', fixed = TRUE)
  expect_match(src, 'updateRadioButtons(session, "trj_lineage_mode", selected = "global")',
               fixed = TRUE)
  # Bouton désactivé tant qu'aucune colonne de pseudo-temps n'existe :
  expect_match(src, 'updateActionButton(session, "trj_compute", disabled',
               fixed = TRUE)
  # Bandeau préventif dédié :
  expect_match(src, 'output$trj_guard <- renderUI', fixed = TRUE)
  expect_match(src, 'uiOutput(ns("trj_guard"))', fixed = TRUE)
  expect_match(src,
    "Aucune colonne de pseudo-temps dans l'objet courant : lancez d'abord Trajectoire (8.4)",
    fixed = TRUE)
})

test_that("3F trajectoire : jauge withProgress autour de l'extraction + du calcul", {
  src <- .read(.TRJ_FILE)
  expect_match(src, 'withProgress(message = .tr("Calcul du contexte trajectoire...")',
               fixed = TRUE)
  expect_match(src, "communication_fetch_expression_matrix(", fixed = TRUE)
  expect_match(src, "build_communication_trajectory_context(", fixed = TRUE)
})

# ── Perturbation : jauge (spéc §2.18 point 3) + garde préventive ─────────────
test_that("3F perturbation : simulation encapsulée dans withProgress + bouton gardé sans résultat", {
  src <- .read(.PERT_FILE)
  # Texte exact de la spéc :
  expect_match(src, 'withProgress(message = .tr("Simulation de perturbation in silico...")',
               fixed = TRUE)
  expect_match(src, "build_communication_perturbation(", fixed = TRUE)
  # Garde préventive : sans résultat importé, la simulation n'a AUCUNE table :
  expect_match(src, 'updateActionButton(session, "pert_compute"', fixed = TRUE)
  expect_match(src, "disabled = is.null(comm_state$result)", fixed = TRUE)
})

# ── Parent : normalisation contextuelle des filtres score / rang (point 4) ───
test_that("3F filtres : intitulés et bornes suivent la source (un rang n'est pas un score)", {
  src <- .read(.PAR_FILE)
  # L'adaptation réagit au résultat déposé :
  expect_match(src, "observeEvent(comm_state$result, {\n      r <- comm_state$result",
               fixed = TRUE)
  # CellChat — probabilité [0.0 - 1.0] (spéc) :
  expect_match(src, "Probabilit\u00e9 de communication [0.0 - 1.0]", fixed = TRUE)
  # CellPhoneDB — moyenne d'expression, borne haute = max de la TABLE :
  expect_match(src, "Moyenne d'expression [0.0 - Max]", fixed = TRUE)
  # LIANA — rang agrégé (1 = plus significatif), score sans objet :
  expect_match(src, "Rang agr\u00e9g\u00e9 [1 - N] (1 = plus significatif)", fixed = TRUE)
  expect_match(src,
    "Score minimum \u2014 sans objet pour LIANA (un rang n'est pas un score : utiliser Rang maximum)",
    fixed = TRUE)
  # Aucune borne fabriquée : les max viennent de la table canonique :
  expect_match(src, 'max(t$score, na.rm = TRUE)', fixed = TRUE)
  expect_match(src, 'max(t$rank, na.rm = TRUE)', fixed = TRUE)
})

# ── Lane 3b (2026-11-14, retour opérateur) : AUTO-CALCUL des contextes ──────
# Demande opérateur : « do we need button suppression in Contextes
# multi-modaux, and active data instead? » — les 3 contextes sans paramètre
# utilisateur obligatoire (spatial, vélocité, trajectoire) se calculent
# SEULS dès que leurs entrées sont réunies et qu'aucun contexte n'est établi
# (plus de clic requis après un reload). Le bouton RESTE (recalcul explicite
# après changement de paramètre). La perturbation garde son bouton SANS auto :
# c'est une simulation pilotée par l'opérateur (cible/mode à choisir).
# Gardes 3F IDENTIQUES (ready/disabled) — l'auto ne court-circuite aucune
# garde, et un échec d'auto reste classé (notification, jamais muet).
test_that("lane3b vélocité : calcul factorisé + auto-calcul sans clic", {
  src <- .read(.VEL_FILE)
  expect_match(src, ".vel_compute_ctx <- function", fixed = TRUE,
               info = "calcul vélocité non factorisé")
  expect_match(src, "observeEvent(input$vel_compute, { .vel_compute_ctx() })",
               fixed = TRUE, info = "le bouton vel_compute ne délègue pas")
  expect_match(src,
               "observeEvent(c(comm_state$result, .vel_vr(), global_data$sc_obj), {",
               fixed = TRUE, info = "auto-calcul vélocité absent")
  expect_match(src, "if (!is.null(vel_ctx())) return()", fixed = TRUE,
               info = "garde anti-recalcul absente (vélocité)")
  # garde 3F consommée PAR l'auto (jamais court-circuitée) :
  expect_match(src, "if (!.vel_ready(.vel_vr())) return()", fixed = TRUE,
               info = "l'auto vélocité ne respecte pas la garde 3F")
})

test_that("lane3b spatial : calcul factorisé + auto-calcul quand la réduction est prête", {
  src <- .read(.SPATIAL_FILE)
  expect_match(src, ".spc_compute_ctx <- function", fixed = TRUE,
               info = "calcul spatial non factorisé")
  expect_match(src, "observeEvent(input$spc_compute, { .spc_compute_ctx() })",
               fixed = TRUE, info = "le bouton spc_compute ne délègue pas")
  expect_match(src,
               "observeEvent(c(comm_state$result, global_data$sc_obj, input$spc_reduction), {",
               fixed = TRUE, info = "auto-calcul spatial absent")
  expect_match(src, "if (!is.null(spc_ctx())) return()", fixed = TRUE,
               info = "garde anti-recalcul absente (spatial)")
})

test_that("lane3b trajectoire : calcul factorisé + auto-calcul garde pseudo-temps", {
  src <- .read(.TRJ_FILE)
  expect_match(src, ".trj_compute_ctx <- function", fixed = TRUE,
               info = "calcul trajectoire non factorisé")
  expect_match(src, "observeEvent(input$trj_compute, { .trj_compute_ctx() })",
               fixed = TRUE, info = "le bouton trj_compute ne délègue pas")
  # Fix live 2026-11-14 : input$trj_pt_col DOIT être un déclencheur — après un
  # reload, sa valeur n'arrive qu'au round-trip client du updateSelectInput
  # (sélection « pseudotime » posée par l'observer d'objet), APRÈS le dernier
  # changement des autres déclencheurs ; sans cet abonnement l'auto rendait la
  # main sur l'input NULL et ne se ré-exécutait jamais (mesuré : repro
  # testServer snapshot réel POST-roitest, scratch_repro_contexts.R).
  expect_match(src,
               paste0("observeEvent(c(comm_state$result, .trj_meta_cols(), ",
                      "global_data$sc_obj,\n                   input$trj_pt_col), {"),
               fixed = TRUE,
               info = paste("auto-calcul trajectoire sans abonnement input$trj_pt_col",
                            "(l'auto ne se ré-exécute jamais au round-trip client)"))
  expect_match(src, "if (!.trj_has_pseudotime(.trj_meta_cols())) return()", fixed = TRUE,
               info = "l'auto trajectoire ne respecte pas la garde pseudo-temps")
})

test_that("lane3b perturbation : PAS d'auto-calcul (simulation pilotée opérateur)", {
  src <- .read(.PERT_FILE)
  expect_match(src, "observeEvent(input$pert_compute, {", fixed = TRUE,
               info = "le bouton perturbation a disparu")
  expect_identical(.n_hits(src, "pert_compute_ctx"), 0L,
                   info = "la perturbation ne doit PAS s'auto-calculer")
})

# ── Lane 3c (2026-11-14, retour live « white card ») : reset conditionnel par identité ──
# Mesuré en live (snapshot POST-roitest, parcours opérateur propre, log DIAG
# /tmp/app_boot8.log) : au reload, le parent fait clear (result <- NULL) PUIS
# le hook restaure le résultat ; l'auto-calcul s'exécute ENTRE les deux (il
# voit la valeur COURANTE) et calcule le contexte ; le reset — dont le 1er
# run est sauté par ignoreInit — arrive ensuite sur l'assignation restaurée
# et efface un contexte calculé SUR le résultat courant ; le garde
# anti-recalcul de l'auto (vel_ctx n'est pas son propre déclencheur) ne
# re-déclenche pas → panneau blanc jusqu'au recalcul manuel. Ordre mesuré :
# AUTO FIRE, AUTO FIRE, RESET(result non-null), RENDER ctx NULL.
# Fix : le reset n'invalide que si le résultat courant diffère (identical)
# du résultat source du contexte. Contrat : les 3 modules portent le garde.
test_that("lane3c white card : reset des contextes conditionné par l'identité du résultat", {
  for (f in c(.VEL_FILE, .TRJ_FILE, .SPATIAL_FILE)) {
    src <- .read(f)
    expect_match(src, "if (identical(comm_state$result, ", fixed = TRUE,
                 info = paste("garde d'identité absente du reset dans", f))
    expect_match(src, "_src(comm_state$result)", fixed = TRUE,
                 info = paste("résultat source non enregistré au calcul dans", f))
    expect_match(src, "_src(NULL)", fixed = TRUE,
                 info = paste("résultat source non purgé à l'erreur dans", f))
  }
})

test_that("lane3c perturbation : switch de cible gardé contre l'input non round-trippé", {
  # Mesuré au repro testServer (cascade clear->restore, scratch_repro_contexts.R
  # phase 3) : observeEvent(comm_state$result) -> .update_value_choices() ->
  # .target_col() -> switch(input$pert_target) plantait (« EXPR must be a
  # length 1 vector ») quand la cible n'a pas encore de valeur client.
  src <- .read(.PERT_FILE)
  expect_match(src, 'switch(input$pert_target %||% "ligand",', fixed = TRUE,
               info = "switch pert_target sans garde NULL")
})

# ── Lane 3c : sentinelle du crash client shiny 1.14 (white card, cause racine) ──
# Mesuré en live (collecteur window.onerror + trace XHR + log websocket) :
# un updateSelect* à choices vides sans `selected` fait planter le CLIENT
# (shiny.min.js s.setValue(u[0].value) sur liste d'options vide ; pour
# server=TRUE le payload perd sa clé value — dropNulls retire le NULL).
# L'exception ABORTE le traitement du lot : le websocket ne délivre plus
# aucune valeur du module (vélocité/trajectoire/spatial calculées mais
# jamais rendues). Le serveur était innocent (repro testServer : sorties OK).
# Défense dans global.R : wrappers interceptant les appels non qualifiés —
# selected = première valeur (comportement conservé) ou character(0) si
# choices vides (payload value: [] -> branche sûre du client).
test_that("lane3c sentinelle client : wrappers updateSelect* dans global.R", {
  src <- .read("global.R")
  expect_match(src, ".ts_update_safe_select_value <- function", fixed = TRUE,
               info = "helper de valeur sûre absent de global.R")
  expect_match(src, "updateSelectInput <- function(session, inputId", fixed = TRUE,
               info = "wrapper updateSelectInput absent (résolution lexicale)")
  expect_match(src, "updateSelectizeInput <- function(session, inputId", fixed = TRUE,
               info = "wrapper updateSelectizeInput absent (résolution lexicale)")
  # comportement : première valeur conservée, vide -> value vide (branche sûre)
  expect_match(src, "if (length(vals)) vals[1] else character(0)", fixed = TRUE,
               info = "la sémantique première-valeur / value-vide n'est pas verrouillée")
  # label-only (choices NULL) et groupes (non atomique) traversés sans toucher
  expect_match(src, "if (!is.atomic(choices)) return(NULL)", fixed = TRUE,
               info = "les choices non atomiques ne sont pas protégés du changement de comportement")
})
