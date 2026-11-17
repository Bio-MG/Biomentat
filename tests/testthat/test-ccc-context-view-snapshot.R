# =============================================================================
# test-ccc-context-view-snapshot.R — PAS 1 idéation saves (2026-11-17)
# =============================================================================
# Retour opérateur : « even save do not save all plots eg in single cell,
# communication, "context multi-modaux" » (IDEATION_SAVE_PLOTS_2026-11-17.md,
# gap mesuré n°2). Pas 1 (effort S) : l'état LÉGER des contextes CCC —
# paramètres de vue + sélections de paire des tiroirs
# velocity/spatial/trajectory — voyage dans le snapshot (jamais d'objet
# plot : au restore, les tiroirs se re-paramètrent et les contextes se
# recalculent depuis les résultats déjà embarqués, auto-calc Lane 3b).
#
# Câblage posé :
#   mod_sc_communication.R : comm_state$ctx_view (source unique) + miroir
#           snapshot shared_rv$communication_ctx_view + hook de réhydratation
#           (communication_restore$ctx_view, dégradation honnête si absent) ;
#   tiroirs (_velocity/_spatial/_trajectory) : publication AU CALCUL (les
#           params qui ONT PRODUIT le contexte) + ré-application au restore
#           avec gardes de validité ; paires sauvegardées honorées si
#           toujours présentes ;
#   app.R : save (results$communication_ctx_view) + restore (pending payload).
#
# Ce lot (testServer) : contrat parent (miroir + one-shot + purge + payloads
# invalides) et contrat tiroirs (ré-application + garde sans contexte).
# La publication au calcul n'est PAS couverte ici (fixtures de contexte
# complets hors périmètre effort S) — vérifiée en live (§W du rapport).
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
source_project_file("modules/sc/mod_sc_communication_velocity.R")
source_project_file("modules/sc/mod_sc_communication_spatial.R")
source_project_file("modules/sc/mod_sc_communication_trajectory.R")
# 4e tiroir (perturbation) : les VRAIS modules sont instanciés par le parent —
# contrairement au harnais lane 3 (shims no-op), ce lot teste aussi les
# ré-applications de restauration, donc les serveurs réels sont requis.
source_project_file("modules/sc/mod_sc_communication_perturbation.R")

suppressPackageStartupMessages({
  library(shiny)
  library(SeuratObject)
})

# ── Shims : MockShinySession n'implémente pas showNotification etc. ─────────
assign("showNotification", function(ui, type = "default", duration = 5, ...)
  invisible(NULL), envir = globalenv())
assign("showModal", function(ui, ...) invisible(NULL), envir = globalenv())
assign("removeModal", function(...) invisible(NULL), envir = globalenv())
assign(".sc_project_root", function() ts_project_root(), envir = globalenv())

# ── Fixture : table canonique minimale (cf. lane 3) ──────────────────────────
.p1_canonical <- function() {
  list(
    source_method = "cellchat",
    canonical_table = data.frame(
      pathway = c("MIF", "MIF"), sender = c("A", "B"), receiver = c("C", "D"),
      sender_mapped = NA_character_, receiver_mapped = NA_character_,
      stringsAsFactors = FALSE),
    provenance = list(list(analysis_id = "p1-test", method = "cellchat"))
  )
}

# ── Fixture : objet Seurat minimal (déclencheur du purge — assigner NULL sur
# un champ déjà NULL n'invalide PAS les observateurs, mesuré) ─────────────────
.p1_seurat_min <- function() {
  CreateSeuratObject(counts = matrix(1L, nrow = 2, ncol = 3,
                                     dimnames = list(c("g1", "g2"),
                                                     c("c1", "c2", "c3"))))
}

# ══════════════════════════════════════════════════════════════════════════════
# CONTRAT PARENT — miroir snapshot + réhydratation + purge
# ══════════════════════════════════════════════════════════════════════════════

test_that("pas1 : le hook de réhydratation porte l'état léger vers le miroir snapshot", {
  gd <- reactiveValues(language = "fr", i18n = NULL, sc_obj = NULL)
  shared <- reactiveValues()
  testServer(mod_sc_communication_server,
             args = list(id = "communication", global_data = gd,
                         shared_rv = shared),
    {
      ctx <- list(
        velocity = list(min_overlap = 0.35),
        spatial = list(source = "spatial_object", radius = 0.4, nperm = 25L),
        trajectory = list(pt_col = "pseudotime", lineage_mode = "global", nbins = 10L)
      )
      shared$communication_restore <- list(
        result = .p1_canonical(), collection = NULL, active_sample = NULL,
        ccc_registry = NULL, ccc_comparison = NULL, ctx_view = ctx)
      session$flushReact()

      # One-shot consommé, résultat réinjecté (contrat lane 3 inchangé)…
      expect_null(shared$communication_restore)
      expect_identical(shared$communication_result$source_method, "cellchat")
      # …et l'état léger est publié au miroir snapshot TEL QUEL — app.R
      # save_session_btn l'embarquera au prochain save.
      expect_identical(shared$communication_ctx_view$velocity$min_overlap, 0.35)
      expect_identical(shared$communication_ctx_view$spatial$radius, 0.4)
      expect_identical(shared$communication_ctx_view$trajectory$nbins, 10L)
    })
})

test_that("pas1 : snapshot antérieur sans ctx_view = dégradation honnête (miroir NULL)", {
  gd <- reactiveValues(language = "fr", i18n = NULL, sc_obj = NULL)
  shared <- reactiveValues()
  testServer(mod_sc_communication_server,
             args = list(id = "communication", global_data = gd,
                         shared_rv = shared),
    {
      # Résultat présent, état léger absent (format d'avant le PAS 1).
      shared$communication_restore <- list(
        result = .p1_canonical(), collection = NULL, active_sample = NULL,
        ccc_registry = NULL, ccc_comparison = NULL)
      session$flushReact()
      expect_identical(shared$communication_result$source_method, "cellchat")
      expect_null(shared$communication_ctx_view)
    })
})

test_that("pas1 : payload ctx_view non-liste (snapshot corrompu) = NULL, jamais de junk", {
  gd <- reactiveValues(language = "fr", i18n = NULL, sc_obj = NULL)
  shared <- reactiveValues()
  testServer(mod_sc_communication_server,
             args = list(id = "communication", global_data = gd,
                         shared_rv = shared),
    {
      shared$communication_restore <- list(
        result = .p1_canonical(), collection = NULL, active_sample = NULL,
        ccc_registry = NULL, ccc_comparison = NULL, ctx_view = "junk")
      session$flushReact()
      expect_null(shared$communication_restore)
      expect_null(shared$communication_ctx_view)
      # Le résultat, lui, n'est pas affecté par la junk de l'état léger.
      expect_identical(shared$communication_result$source_method, "cellchat")
    })
})

test_that("pas1 : le changement d'objet purge l'état léger ET son miroir", {
  gd <- reactiveValues(language = "fr", i18n = NULL, sc_obj = NULL)
  shared <- reactiveValues()
  testServer(mod_sc_communication_server,
             args = list(id = "communication", global_data = gd,
                         shared_rv = shared),
    {
      shared$communication_restore <- list(
        result = .p1_canonical(), collection = NULL, active_sample = NULL,
        ccc_registry = NULL, ccc_comparison = NULL,
        ctx_view = list(velocity = list(min_overlap = 0.35)))
      session$flushReact()
      expect_identical(shared$communication_ctx_view$velocity$min_overlap, 0.35)

      # Remplacement de jeu (import / activation dataset) : l'observateur
      # d'objet fait .clear_result_state() — l'état léger vit avec le
      # résultat qu'il paramètre, son miroir avec lui.
      gd$sc_obj <- .p1_seurat_min()
      session$flushReact()
      expect_null(shared$communication_result)
      expect_null(shared$communication_ctx_view)
    })
})

# ══════════════════════════════════════════════════════════════════════════════
# CONTRAT TIROIRS — ré-application des params sauvegardés + garde sans contexte
# =============================================================================
# NOTE outillage (mesuré) : MockShinySession n'applique PAS les messages
# update*Input aux inputs mockés (sonde shiny 4.4.2 : input reste NULL après
# update), et aucun miroir inputMessages n'est exposé. La ré-application est
# donc verrouillée par CONTRAT DE SOURCE (idiome repo, cf.
# test-mod-sc-communication-contexts.R) ; les gardes de publication restent
# testées COMPORTEMENT (elles n'exigent aucune réflexion d'input).
# =============================================================================

.read_mod <- function(rel) {
  paste(readLines(file.path(ts_project_root(), rel),
                  warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

test_that("pas1 velocity : contrat de ré-application + publication au calcul", {
  src <- .read_mod("modules/sc/mod_sc_communication_velocity.R")
  # L'observateur de restauration existe, abonné au bon tiroir de l'état.
  expect_match(src, "observeEvent(comm_state$ctx_view$velocity, {", fixed = TRUE,
               info = "observateur de ré-application velocity absent")
  # Mise à jour directe du seuil (garde identique + jamais de NULL).
  expect_match(src, "shiny::updateNumericInput(session, \"vel_min_overlap\", value = mo)",
               fixed = TRUE,
               info = "le seuil sauvegardé n'est pas ré-appliqué")
  # Publication AU CALCUL : les params qui ONT PRODUIT le contexte voyagent.
  expect_match(src, ".ctx_view_set(\"velocity\", list(min_overlap = input$vel_min_overlap %||% NA_real_))",
               fixed = TRUE,
               info = "la publication au calcul du contexte velocity est absente")
  # Aide jumelle : read-modify-write explicite (affectation imbriquée interdite).
  expect_match(src, "if (!is.list(cur)) cur <- list()", fixed = TRUE,
               info = "la garde read-modify-write de .ctx_view_set est absente")
})

test_that("pas1 spatial : contrat de ré-application + gardes de validité + paires", {
  src <- .read_mod("modules/sc/mod_sc_communication_spatial.R")
  expect_match(src, "observeEvent(comm_state$ctx_view$spatial, {", fixed = TRUE,
               info = "observateur de ré-application spatial absent")
  # Radio + numériques ré-appliqués.
  expect_match(src, "shiny::updateRadioButtons(session, \"spc_source\", selected = v$source)",
               fixed = TRUE, info = "la source sauvegardée n'est pas ré-appliquée")
  # La réduction n'est ré-appliquée que si VALIDE sur l'objet courant.
  expect_match(src, "red %in% .spc_valid_reductions()", fixed = TRUE,
               info = "la garde de validité de la réduction restaurée est absente")
  # Garde de publication : SANS contexte établi, rien n'est publié (au restore,
  # les observateurs d'entrée déclenchés par la ré-application ne peuvent pas
  # écraser l'état restauré).
  expect_match(src, "if (is.null(spc_ctx())) return()", fixed = TRUE,
               info = "la garde de publication spatial est absente")
  # La paire sauvegardée est honorée si toujours présente (restauration =
  # même vue qu'avant la sauvegarde).
  expect_match(src, "saved_pair %in% pair_vals", fixed = TRUE,
               info = "l'honoration de la paire sauvegardée est absente")
})

test_that("pas1 trajectory : contrat de ré-application + gardes + paire/lignée", {
  src <- .read_mod("modules/sc/mod_sc_communication_trajectory.R")
  expect_match(src, "observeEvent(comm_state$ctx_view$trajectory, {", fixed = TRUE,
               info = "observateur de ré-application trajectory absent")
  # Colonnes d'objet : ré-appliquées seulement si VALIDES (mêmes gardes que
  # le peuplement).
  expect_match(src, "pt %in% pt_cols", fixed = TRUE,
               info = "la garde de validité de la colonne de pseudo-temps est absente")
  expect_match(src, "asy %in% names(obj@assays)", fixed = TRUE,
               info = "la garde de validité de l'assay est absente")
  # Garde de publication identique au tiroir spatial.
  expect_match(src, "if (is.null(trj_ctx())) return()", fixed = TRUE,
               info = "la garde de publication trajectory est absente")
  # Paire/lignée sauvegardées honorées si toujours présentes.
  expect_match(src, "saved_pair %in% pair_pool", fixed = TRUE,
               info = "l'honoration de la paire sauvegardée est absente")
  expect_match(src, "saved_lin %in% lins", fixed = TRUE,
               info = "l'honoration de la lignée sauvegardée est absente")
})

test_that("pas1 spatial : garde de publication — sans contexte, RIEN n'est publié", {
  gd <- reactiveValues(language = "fr", i18n = NULL, sc_obj = NULL,
                       spatial_obj = NULL)
  shared <- reactiveValues()
  comm_state <- reactiveValues()
  testServer(mod_sc_communication_spatial_server,
             args = list(id = "spatial", comm_state = comm_state,
                         global_data = gd, shared_rv = shared),
    {
      # L'opérateur change des params SANS contexte calculé : aucune
      # publication — au restore, les observateurs d'entrée déclenchés par
      # la ré-application ne peuvent pas écraser l'état restauré.
      session$setInputs(spc_source = "spatial_object")
      session$setInputs(spc_radius = 0.6)
      session$flushReact()
      expect_null(comm_state$ctx_view)
      expect_null(shared$communication_ctx_view)
    })
})

test_that("pas1 trajectory : garde de publication — sans contexte, RIEN n'est publié", {
  gd <- reactiveValues(language = "fr", i18n = NULL, sc_obj = NULL)
  shared <- reactiveValues()
  comm_state <- reactiveValues()
  testServer(mod_sc_communication_trajectory_server,
             args = list(id = "trajectory", comm_state = comm_state,
                         global_data = gd, shared_rv = shared),
    {
      session$setInputs(trj_nbins = 30L)
      session$setInputs(trj_lineage_mode = "global")
      session$flushReact()
      expect_null(comm_state$ctx_view)
      expect_null(shared$communication_ctx_view)
    })
})
