# =============================================================================
# test-sc-async-optin.R — Lane async opt-in (2026-11-14, §M priorité 3)
# =============================================================================
# Retour opérateur : les calculs lourds (scVelo opt-in, CellChat) gèlent la
# session. MESURÉ au source : CellChat est STRUCTURELLEMENT synchrone
# (mod_sc_communication.R, observer comm_compute — jamais migré) et scVelo
# opt-in aussi (mod_sc_velocity.R, observer velocity_scvelo_run — notification
# persistante + tryCatch bloquant 10-90 s).
#
# Le pattern existe déjà et est validé en production : LIANA P1-6b
# (mod_sc_communication.R:1280-1474) — ExtendedTask + pool unique
# (ensure_app_daemons, spatial_async.R) + mirai (deps sourcés idempotents,
# chemins absolus) + .timeout + garde de ré-entrée + repli synchrone TRACÉ
# (jamais muet) + observateur d'achèvement (mêmes surfaces que le chemin
# synchrone : statut, notification classée, purge).
#
# Contrainte moteur (sc_communication_engine.R, INVARIANT) : « future forcé
# sequential à l'INTÉRIEUR du job, mirai reste le seul orchestrateur » — les
# deux moteurs sont déjà worker-safe ; seul le câblage module manque.
#
# Limite documentée (même motif que test-mod-sc-communication-ui.R) : les
# ExtendedTask vivent DANS les moduleServer avec des gardes sur l'objet
# Seurat/LIANA — la preuve est un verrou source structurel (contrats ci-
# dessous), la mécanique ExtendedTask/mirai étant déjà validée E2E par
# LIANA en production.
# =============================================================================

Sys.setenv(OMNIPATHR_LOGFILE = "none")

source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/validation.R")
source_project_file("tools/pkg_skeleton/R/provenance.R")
source_project_file("tools/pkg_skeleton/R/error_state.R")
source_project_file("tools/pkg_skeleton/R/dataset_directory.R")
source_project_file("tools/pkg_skeleton/R/datatable.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
source_project_file("tools/pkg_skeleton/R/sc_communication.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_input.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_engine.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_views.R")
source_project_file("modules/sc/mod_sc_velocity.R")
source_project_file("modules/sc/mod_sc_communication.R")

suppressPackageStartupMessages({
  library(shiny)
})

# Sous-modules contextuels hors périmètre : no-op (cf. B5 / CQF-2).
assign("mod_sc_communication_spatial_server",
       function(...) invisible(NULL), envir = globalenv())
assign("mod_sc_communication_trajectory_server",
       function(...) invisible(NULL), envir = globalenv())
assign("mod_sc_communication_velocity_server",
       function(...) invisible(NULL), envir = globalenv())
assign("mod_sc_communication_perturbation_server",
       function(...) invisible(NULL), envir = globalenv())
assign("showNotification", function(ui, type = "default", duration = 5, ...)
  invisible(NULL), envir = globalenv())
assign("showModal", function(ui, ...) invisible(NULL), envir = globalenv())
assign("removeModal", function(...) invisible(NULL), envir = globalenv())
assign(".sc_project_root", function() ts_project_root(), envir = globalenv())

.vel_file <- "modules/sc/mod_sc_velocity.R"
.comm_file <- "modules/sc/mod_sc_communication.R"
.read_src <- function(rel) {
  paste(readLines(file.path(ts_project_root(), rel),
                  warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}
.n_hits_src <- function(src, pattern) {
  lengths(regmatches(src, gregexpr(pattern, src, fixed = TRUE)))[1]
}

# ── scVelo opt-in : ExtendedTask + mirai (pattern LIANA P1-6b) ───────────────
test_that("lane async : scVelo opt-in part en daemon (ExtendedTask + pool unique)", {
  src <- .read_src(.vel_file)
  expect_match(src, "scvelo_task <- ", fixed = TRUE,
               info = "tâche scVelo absente")
  # résolution de la classe ExtendedTask (shiny/bslib/globalenv, idiom LIANA)
  # + construction $new (l'idiome exact de LIANA : ext_task_cls_*$new)
  expect_match(src, 'exists("ExtendedTask", where = asNamespace("shiny"), inherits = FALSE)',
               fixed = TRUE,
               info = "la classe ExtendedTask n'est pas résolue pour scVelo")
  expect_match(src, "ext_task_cls_vel$new(", fixed = TRUE,
               info = "la tâche scVelo n'utilise pas ExtendedTask$new")
  # le daemon charge le moteur (idempotent, chemins absolus)
  expect_match(src, "sc_velocity.R", fixed = TRUE,
               info = "le moteur scVelo n'est pas sourcé dans le daemon")
  expect_match(src, "error_state.R", fixed = TRUE,
               info = "error_state.R absent des deps daemon scVelo")
  expect_match(src, "ensure_app_daemons()", fixed = TRUE,
               info = "le pool unique n'est pas assuré avant le mirai")
  expect_match(src, ".timeout = TS_SC_AUTO_PIPELINE_TIMEOUT_S * 1000L", fixed = TRUE,
               info = "pas de timeout mirai sur la tâche scVelo")
})

test_that("lane async : scVelo garde de ré-entrée + repli synchrone TRACÉ", {
  src <- .read_src(.vel_file)
  # ré-entrée : un calcul en vol ne relance pas
  expect_match(src, 'identical(scvelo_task$status(), "running")', fixed = TRUE,
               info = "garde de ré-entrée absente sur l'opt-in scVelo")
  # repli synchrone tracé (jamais muet — convention run_job)
  expect_match(src, "exécution synchrone", fixed = TRUE,
               info = "le repli synchrone scVelo n'est pas tracé")
  # bind_task_button tenté (comme LIANA — actionButton toléré via tryCatch)
  expect_match(src, 'bind_task_button(t, "velocity_scvelo_run"', fixed = TRUE,
               info = "bind_task_button absent pour la tâche scVelo")
})

test_that("lane async : achèvement scVelo = MÊMES surfaces que le synchrone", {
  src <- .read_src(.vel_file)
  # l'application du résultat est FACTORISÉE (zéro voie parallèle) : utilisée
  # par le repli synchrone ET l'observateur d'achèvement
  expect_match(src, ".apply_scvelo_result <- function", fixed = TRUE,
               info = "application du résultat scVelo non factorisée")
  n <- .n_hits_src(src, ".apply_scvelo_result(")
  expect_gte(n, 2L)
  # observateur d'achèvement sur le statut de la tâche
  expect_match(src, "observeEvent(scvelo_task$status(), {", fixed = TRUE,
               info = "observateur d'achèvement scVelo absent")
  # échec async = état persistant optin_error (même surface que le synchrone)
  expect_match(src, "velocity_state$optin_error <- list(", fixed = TRUE,
               info = "l'échec async scVelo ne laisse pas d'état persistant")
})

# ── CellChat : ExtendedTask + mirai (même pattern) ───────────────────────────
test_that("lane async : CellChat part en daemon (ExtendedTask + pool unique)", {
  src <- .read_src(.comm_file)
  expect_match(src, "cellchat_task <- ", fixed = TRUE,
               info = "tâche CellChat absente")
  expect_match(src, 'exists("ExtendedTask", where = asNamespace("shiny"), inherits = FALSE)',
               fixed = TRUE,
               info = "la classe ExtendedTask n'est pas résolue pour CellChat")
  expect_match(src, "ext_task_cls_cellchat$new(", fixed = TRUE,
               info = "la tâche CellChat n'utilise pas ExtendedTask$new")
  expect_match(src, "sc_communication_input.R", fixed = TRUE,
               info = "sc_communication_input.R absent des deps daemon CellChat")
  expect_match(src, "sc_communication_engine.R", fixed = TRUE,
               info = "le moteur CellChat n'est pas sourcé dans le daemon")
  expect_match(src, "ensure_app_daemons()", fixed = TRUE,
               info = "le pool unique n'est pas assuré avant le mirai CellChat")
  expect_match(src, ".timeout = TS_SC_AUTO_PIPELINE_TIMEOUT_S * 1000L", fixed = TRUE,
               info = "pas de timeout mirai sur la tâche CellChat")
})

test_that("lane async : CellChat garde de ré-entrée + repli synchrone TRACÉ", {
  src <- .read_src(.comm_file)
  expect_match(src, 'identical(cellchat_task$status(), "running")', fixed = TRUE,
               info = "garde de ré-entrée absente sur le calcul CellChat")
  expect_match(src, "exécution synchrone", fixed = TRUE,
               info = "le repli synchrone CellChat n'est pas tracé")
})

test_that("lane async : achèvement CellChat = .store_result (mêmes surfaces)", {
  src <- .read_src(.comm_file)
  expect_match(src, ".apply_cellchat_result <- function", fixed = TRUE,
               info = "application du résultat CellChat non factorisée")
  n <- .n_hits_src(src, ".apply_cellchat_result(")
  expect_gte(n, 2L)
  expect_match(src, "observeEvent(cellchat_task$status(), {", fixed = TRUE,
               info = "observateur d'achèvement CellChat absent")
  # la voie async applique le résultat PAR LE MÊME helper que le synchrone
  expect_match(src, ".apply_cellchat_result(cellchat_task$result()", fixed = TRUE,
               info = "l'achèvement async CellChat ne passe pas par le helper partagé")
})

# ── R-3bis (idiome repo) : miroir de visibilité du drawer comm_tabs ─────────
# Mesuré en live (E2E async CellChat, boot21) : comm_table livrée mais
# comm_status — invalide au MÊME flush — jamais reçue (websocket $values) ;
# la famille « sortie invalidée pendant cachée ne livre jamais » (mod_bulk.R
# fix 2026-11-06, drive_allowlist.R) frappe le drawer entier. Cure établie :
# le parent publie la feuille active (shared_rv$comm_tabs_current) et les
# renders consomment le miroir — re-exécution à l'ouverture du sous-onglet.
test_that("lane async : miroir R-3bis du drawer consommé par les renders", {
  src <- .read_src(.comm_file)
  expect_match(src, 'observeEvent(input$comm_tabs, {', fixed = TRUE,
               info = "publication du miroir de feuille absente (parent)")
  expect_match(src, 'shared_rv$comm_tabs_current <- input$comm_tabs', fixed = TRUE,
               info = "la feuille active n'est pas publiée dans shared_rv")
  expect_match(src, 'state_get(shared_rv, "comm_tabs_current")', fixed = TRUE,
               info = "comm_status ne consomme pas le miroir R-3bis")
  # les 4 modules de contexte consomment le miroir avec LEUR feuille
  # (correctif 2026-11-14 : le panneau PRINCIPAL mod_sc_velocity.R n'a PAS de
  # garde R-3bis — c'est le module de CONTEXTE vélocité, ci-dessous, qui la
  # porte ; l'ancienne variable indéfinie .VEL_FILE masquait cette erreur)
  for (f in c("modules/sc/mod_sc_communication_velocity.R",
              "modules/sc/mod_sc_communication_trajectory.R",
              "modules/sc/mod_sc_communication_spatial.R",
              "modules/sc/mod_sc_communication_perturbation.R")) {
    s <- .read_src(f)
    expect_match(s, ".ctx_tab_guard <- function", fixed = TRUE,
                 info = paste("helper de garde R-3bis absent de", f))
    expect_match(s, 'state_get(shared_rv, "comm_tabs_current")', fixed = TRUE,
                 info = paste("le miroir R-3bis n'est pas consommé dans", f))
  }
  # chaque garde verrouille SA feuille
  expect_match(.read_src("modules/sc/mod_sc_communication_velocity.R"),
               '"ccc_tab_velocite"', fixed = TRUE)
  expect_match(.read_src("modules/sc/mod_sc_communication_trajectory.R"),
               '"ccc_tab_trajectoire"', fixed = TRUE)
  expect_match(.read_src("modules/sc/mod_sc_communication_spatial.R"),
               '"ccc_tab_spatial"', fixed = TRUE)
  expect_match(.read_src("modules/sc/mod_sc_communication_perturbation.R"),
               '"ccc_tab_perturbation"', fixed = TRUE)
})
