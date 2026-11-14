# =============================================================================
# test-sc-velocity-session-restore.R — Lane 3b (2026-11-14)
# =============================================================================
# Retour d'expérience opérateur sur le premier test save/reload lane 3 :
# « a part of communication is saved (other need re calculus) and velocity
# NOTHING saved ». Bug confirmé au source : la restauration directe
# shared$velocity_result dans le handler de chargement était PURGÉE juste
# après par l'observer d'epoch de mod_sc (F08 : sc_purge_shared_results purge
# « velocity_result ») — la communication survivait (comm_state dans la
# closure + champs communication_* hors liste de purge).
#
# Fix :
#   app.R — PENDING PAYLOAD : le chargement pose global_data
#           $session_results_pending ; l'observer qui réinjecte est créé
#           APRÈS mod_sc_server() (les observers s'exécutent dans l'ordre de
#           création) → il passe APRÈS les purges d'epoch/sc_obj du même
#           cycle de flush ; les hooks one-shot des modules s'exécutent au
#           cycle suivant, plus aucune purge en attente.
#   mod_sc_velocity.R — publication de l'état du PANNEAU
#           (.velocity_publish_state() aux succès validate/opt-in, purge au
#           changement d'objet) + hook one-shot velocity_restore qui
#           réhydrate velocity_state (matrices validées, vecteurs, alignement
#           — JAMAIS les datapath temporaires des fileInput ni les caches).
#
# Ce lot : contrat du hook (one-shot, résultat restauré, panneau réinjecté,
# statut honnête) + payload invalide ignoré.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/validation.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
source_project_file("modules/sc/mod_sc_velocity.R")

# Fixture minimale : objet Seurat + canonique opaque (le hook assigne tel
# quel — la structure du canonique est contractée ailleurs, lane 1).
.l3v_cells <- paste0("cell", seq_len(20))
.l3v_canonical <- function() {
  list(
    status = "valid",
    cell_names = .l3v_cells,
    velocity_vectors = matrix(runif(40), nrow = 20, dimnames = list(.l3v_cells, c("dx", "dy"))),
    params = list(reduction = "umap"),
    provenance = list(list(analysis_id = "l3v-test", method = "velocity"))
  )
}
.l3v_panel <- function() {
  g <- paste0("gene", seq_len(5))
  list(
    result             = .l3v_canonical(),
    validated          = list(
      cell_names = .l3v_cells,
      spliced   = matrix(1, nrow = 5, ncol = 20,
                         dimnames = list(g, .l3v_cells)),
      unspliced = matrix(1, nrow = 5, ncol = 20,
                         dimnames = list(g, .l3v_cells))),
    alignment          = list(cells = .l3v_cells, genes = g),
    umap_vectors       = NULL,
    partial_embedding  = NULL,
    embedding_alignment = NULL,
    vector_validation  = list(field = "scvelo_optin", ok = TRUE),
    object_fingerprint = "snapshot-fp"
  )
}

# Shim : capture du updateSelectizeInput (sous testServer, les update* ne
# déplacent pas input$ — on asserte ce que le SERVEUR envoie).
.l3v_updates <- new.env(parent = emptyenv())
assign("updateSelectizeInput", function(session, inputId, choices = NULL,
                                        selected = NULL, server = FALSE, ...) {
  .l3v_updates[[inputId]] <- list(choices = choices, selected = selected)
  invisible(NULL)
}, envir = globalenv())

test_that("lane3b : le hook velocity réhydrate panneau + résultat (one-shot)", {
  gd <- shiny::reactiveValues(sc_obj = NULL, language = "fr", i18n = NULL)
  shared <- shiny::reactiveValues()
  shiny::testServer(mod_sc_velocity_server,
                    args = list(global_data = gd, shared_rv = shared),
    {
      # Avant restauration : purgé (ce que fait l'epoch/sc_obj au reload).
      shared$velocity_result <- NULL
      shared$velocity_restore <- list(
        velocity_result = .l3v_canonical(),
        panel_state     = .l3v_panel())
      session$flushReact()

      # One-shot : le déclencheur est consommé.
      expect_null(shared$velocity_restore)
      # Le résultat consommé (Communication/viz/rapport) est restauré.
      expect_identical(shared$velocity_result$status, "valid")
      expect_identical(length(shared$velocity_result$velocity_vectors[, 1]), 20L)
      # L'état du panneau est réinjecté — observable via le statut honnête
      # et la vue vecteurs (lisant velocity_state).
      expect_match(output$velocity_status, "restaur")
      expect_match(output$velocity_vector_status, "Statut : ")
    })
})

test_that("lane3b : résultat seul (snapshot lane 3 v1) restaure le consommé sans panneau", {
  gd <- shiny::reactiveValues(sc_obj = NULL, language = "fr", i18n = NULL)
  shared <- shiny::reactiveValues()
  shiny::testServer(mod_sc_velocity_server,
                    args = list(global_data = gd, shared_rv = shared),
    {
      shared$velocity_restore <- list(velocity_result = .l3v_canonical(),
                                      panel_state = NULL)
      session$flushReact()
      expect_null(shared$velocity_restore)
      expect_identical(shared$velocity_result$status, "valid")
      # Pas d'état de panneau : le statut du panneau reste honnête (inchangé).
      expect_false(grepl("restaur", output$velocity_status))
    })
})

test_that("lane3b : la réhydratation peuple le sélecteur de gène (phase portrait)", {
  # Retour opérateur « phase portrait isn't » : sans choices sur
  # velocity_gene, req(input$velocity_gene) laissait le phase portrait vide
  # après reload. Le hook doit émettre le MÊME peuplement que la validation.
  gd <- shiny::reactiveValues(sc_obj = NULL, language = "fr", i18n = NULL)
  shared <- shiny::reactiveValues()
  shiny::testServer(mod_sc_velocity_server,
                    args = list(global_data = gd, shared_rv = shared),
    {
      shared$markers_data <- NULL
      shared$velocity_restore <- list(
        velocity_result = .l3v_canonical(),
        panel_state     = .l3v_panel())
      session$flushReact()
      sent <- .l3v_updates$velocity_gene
      expect_false(is.null(sent),
                   info = "le hook n'a pas émis updateSelectizeInput(velocity_gene)")
      expect_identical(sent$choices, paste0("gene", seq_len(5)))
      expect_false(is.null(sent$selected),
                   info = "aucun gène par défaut sélectionné après restauration")
    })
})

test_that("lane3b : payload invalide ignoré proprement", {
  gd <- shiny::reactiveValues(sc_obj = NULL, language = "fr", i18n = NULL)
  shared <- shiny::reactiveValues()
  shiny::testServer(mod_sc_velocity_server,
                    args = list(global_data = gd, shared_rv = shared),
    {
      shared$velocity_restore <- "junk"
      session$flushReact()
      expect_null(shared$velocity_restore)
      expect_null(shared$velocity_result)
    })
})

# Contrats source : le pending payload d'app.R, la publication panneau et les
# purges cohérentes.
test_that("lane3b contrats source (app.R + mod_sc_velocity.R)", {
  src_app <- paste(readLines(file.path(ts_project_root(), "app.R"),
                             warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src_app, "global_data$session_results_pending <- list(", fixed = TRUE,
               info = "pending payload absent du chargement")
  expect_match(src_app, "observeEvent(global_data$session_results_pending, {",
               fixed = TRUE, info = "observer de ré-injection absent d'app.R")
  expect_match(src_app, "velocity_panel_state        = shared$velocity_panel_state",
               fixed = TRUE, info = "état panneau absent du snapshot")
  src_v <- paste(readLines(file.path(ts_project_root(), "modules", "sc",
                                     "mod_sc_velocity.R"),
                           warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src_v, ".velocity_publish_state <- function", fixed = TRUE,
               info = "helper de publication absent")
  expect_match(src_v, "observeEvent(shared_rv$velocity_restore, {", fixed = TRUE,
               info = "hook de réhydratation velocity absent")
  expect_match(src_v, "État velocity restauré depuis la session sauvegardée",
               fixed = TRUE, info = "statut honnête de restauration absent")
  # cohérence purge : le changement d'objet invalide AUSSI l'état publié
  expect_match(src_v, "shared_rv$velocity_panel_state <- NULL", fixed = TRUE,
               info = "purge de l'état publié absente")
  # zéro voie parallèle : validate et réhydratation partagent le helper gène
  n_uses <- lengths(regmatches(src_v, gregexpr(".velocity_refresh_gene_selector()",
                                               src_v, fixed = TRUE)))[1]
  expect_gte(n_uses, 2L)
})
