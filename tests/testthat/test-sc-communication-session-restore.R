# =============================================================================
# test-sc-communication-session-restore.R — Lane 3 (2026-11-14)
# =============================================================================
# Sauvegarde de session COMPLÈTE (rapport P0.4 §K proposition 2, demande
# opérateur : « save all session biomentat with cellchat ccc, and velocity »).
#
# Mesuré avant fix (app.R:1530 + §H.7) : le snapshot ne capture que les
# datasets — après « Charger Session », velocity/communication (CellChat)
# étaient silencieusement PERDUS.
#
# Câblage posé :
#   app.R : snapshot porte une sous-liste `results` (velocity_result +
#           communication_result/collection + registre CC-2 + comparatif +
#           provenance) et la restauration pose shared_rv$communication_restore
#           (one-shot) + velocity_result directement ;
#   mod_sc.R : référence partagée publiée à l'app (global_data$sc_shared) ;
#   mod_sc_communication.R : publication du registre CC-2 vers shared_rv
#           (dépôt + comparatif + purge) et hook de réhydratation one-shot
#           qui réinjecte comm_state via le MÊME helper de filtres que le
#           dépôt (zéro voie parallèle, rien n'est ré-inféré).
#
# Ce lot (testServer) : le hook réinjecte l'état, consomme le déclencheur,
# publie vers shared_rv et énonce un statut HONNÊTE au panneau.
# =============================================================================

Sys.setenv(OMNIPATHR_LOGFILE = "none")

source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/validation.R")
source_project_file("tools/pkg_skeleton/R/provenance.R")
source_project_file("tools/pkg_skeleton/R/error_state.R")
source_project_file("tools/pkg_skeleton/R/dataset_directory.R")
source_project_file("tools/pkg_skeleton/R/datatable.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
source_project_file("tools/pkg_skeleton/R/sc_communication.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_input.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_views.R")
source_project_file("modules/sc/mod_sc_communication.R")

suppressPackageStartupMessages({
  library(shiny)
})

# ── Sous-modules contextuels hors périmètre : no-op (cf. B5 / CQF-2) ────────
assign("mod_sc_communication_spatial_server",
       function(...) invisible(NULL), envir = globalenv())
assign("mod_sc_communication_trajectory_server",
       function(...) invisible(NULL), envir = globalenv())
assign("mod_sc_communication_velocity_server",
       function(...) invisible(NULL), envir = globalenv())
assign("mod_sc_communication_perturbation_server",
       function(...) invisible(NULL), envir = globalenv())

# ── Shims : MockShinySession n'implémente pas showNotification etc. ─────────
assign("showNotification", function(ui, type = "default", duration = 5, ...)
  invisible(NULL), envir = globalenv())
assign("showModal", function(ui, ...) invisible(NULL), envir = globalenv())
assign("removeModal", function(...) invisible(NULL), envir = globalenv())
# .sc_project_root est défini DANS app.R (remontée des chemins de projet) —
# le harnais testServer le shimme vers la racine du dépôt (cf. LIANA deps).
assign(".sc_project_root", function() ts_project_root(), envir = globalenv())

# ── Fixture : table canonique minimale (colonnes consommées par les filtres) ─
.l3_canonical <- function() {
  list(
    source_method = "cellchat",
    canonical_table = data.frame(
      pathway = c("MIF", "MIF"), sender = c("A", "B"), receiver = c("C", "D"),
      sender_mapped = NA_character_, receiver_mapped = NA_character_,
      stringsAsFactors = FALSE),
    provenance = list(list(analysis_id = "l3-test", method = "cellchat"))
  )
}

.l3_registry <- function() {
  list("Run 1 (00:00) — 47 voies" = list(label = "Run 1 (00:00) — 47 voies",
                                         object = NULL, ts = Sys.time()))
}

test_that("lane3 : le hook de réhydratation réinjecte l'état Communication (one-shot)", {
  gd <- reactiveValues(language = "fr", i18n = NULL, sc_obj = NULL)
  shared <- reactiveValues()
  testServer(mod_sc_communication_server,
             args = list(id = "communication", global_data = gd,
                         shared_rv = shared),
    {
      # Avant : panneau vide.
      shared$communication_result <- NULL
      # Déclencheur posé par app.R au chargement du snapshot.
      shared$communication_restore <- list(
        result = .l3_canonical(), collection = NULL, active_sample = NULL,
        ccc_registry = .l3_registry(), ccc_comparison = NULL)
      session$flushReact()

      # One-shot : le déclencheur est consommé.
      expect_null(shared$communication_restore)
      # L'état partagé est (re)publié pour les consommateurs (rapport
      # consolidé, contextes).
      expect_identical(shared$communication_result$source_method, "cellchat")
      expect_identical(length(shared$communication_ccc_registry), 1L)
      # Le registre restauré N'EST PAS re-registré (déjà porté par le résultat).
      # Statut HONNÊTE du panneau.
      expect_match(output$comm_status, "restaur")

      # Ré-import impossible : rien ne ré-infère ici — un second dépôt via
      # .store_result re-registrerait, mais le hook, lui, reste one-shot.
      shared$communication_restore <- NULL
      session$flushReact()
      expect_null(shared$communication_restore)
    })
})

test_that("lane3 : le hook ignore proprement un payload invalide (dégradation honnête)", {
  gd <- reactiveValues(language = "fr", i18n = NULL, sc_obj = NULL)
  shared <- reactiveValues()
  testServer(mod_sc_communication_server,
             args = list(id = "communication", global_data = gd,
                         shared_rv = shared),
    {
      # Payload non-liste (snapshot corrompu / format inattendu) : consommé,
      # rien n'est injecté, aucun plantage.
      shared$communication_restore <- "junk"
      session$flushReact()
      expect_null(shared$communication_restore)
      expect_null(shared$communication_result)
      expect_null(shared$communication_ccc_registry)
      # Payload liste SANS résultat : même contrat.
      shared$communication_restore <- list(collection = NULL)
      session$flushReact()
      expect_null(shared$communication_restore)
      expect_null(shared$communication_result)
    })
})
