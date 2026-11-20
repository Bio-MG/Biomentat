# =============================================================================
# test-mod-sc-communication-perturbation-ui.R — W2.2 §5.4 : UI + serveur
# perturbation multi-cibles, pathway entier et population (sender + receiver)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_communication_perturbation.R")

.pert_ui_env <- function() {
  skip_if_not_installed("shiny")
  skip_if_not_installed("DT")
  env <- new.env(parent = asNamespace("shiny"))
  env$i18n <- list(t = function(x) x)
  env$.tr_plain <- function(x) x
  env$.strip_i18n_html <- function(x) x
  env$communication_perturbation_targets <- communication_perturbation_targets
  env$build_communication_perturbation <- build_communication_perturbation
  env$velocity_object_fingerprint <- function(obj) "fp_test"
  env$state_get <- function(rv, key) {
    if (is.null(rv)) return(NULL)
    rv[[key]]
  }
  env$provenance_append <- function(rv, prov) invisible(NULL)
  env$showNotification <- function(...) invisible(NULL)
  env$withProgress <- function(message, value = NULL, expr) force(expr)
  sys.source(
    file.path(ts_project_root(), "modules", "sc", "mod_sc_communication_perturbation.R"),
    envir = env
  )
  env
}

test_that("UI perturbation exposes pathway, population, multi-select pert_value, and contract guards", {
  env <- .pert_ui_env()
  ui <- env$mod_sc_communication_perturbation_ui("pert")
  html <- as.character(shiny::tagList(ui))

  # Epistemic banner and core controls
  expect_match(html, "IN SILICO PERTURBATION", fixed = TRUE)
  expect_match(html, 'id="pert-pert_only_affected"', fixed = TRUE)
  expect_match(html, 'id="pert-pert_table"', fixed = TRUE)
  expect_match(html, 'id="pert-dl_pert_delta_png"', fixed = TRUE)
  expect_match(html, 'id="pert-dl_pert_delta_pdf"', fixed = TRUE)
  expect_match(html, 'id="pert-dl_pert_nodes_png"', fixed = TRUE)
  expect_match(html, 'id="pert-dl_pert_nodes_pdf"', fixed = TRUE)
  expect_match(html, 'id="pert-dl_pert_csv"', fixed = TRUE)

  # Target choices include all canonical targets + pathway + population
  for (tgt in c(communication_perturbation_targets(), "pathway", "population")) {
    expect_match(html, sprintf('value="%s"', tgt), fixed = TRUE)
  }
  expect_match(html, "Voie de signalisation", fixed = TRUE)
  expect_match(html, "Population cellulaire (sender + receiver)", fixed = TRUE)

  # pert_value is a multiple selectizeInput
  expect_match(html, 'id="pert-pert_value"', fixed = TRUE)
  expect_match(html, 'multiple="multiple"', fixed = TRUE)
})

test_that("server .update_value_choices() and pert_target_info handle pathway (omitting NA) and population", {
  env <- .pert_ui_env()
  updates <- list()
  env$updateSelectizeInput <- function(session, inputId, choices = NULL, selected = NULL, ...) {
    updates[[length(updates) + 1L]] <<- list(
      inputId = inputId, choices = choices, selected = selected
    )
    invisible(NULL)
  }

  r <- .ccx_result()
  # Inject one NA pathway row to verify NA omission in UI choices
  r$canonical_table$pathway[1L] <- NA_character_

  comm_state <- shiny::reactiveValues(
    result = NULL,
    object_fingerprint = "fp_test",
    ctx_results = NULL
  )
  gd <- shiny::reactiveValues(sc_obj = "mock_sc", language = "fr", i18n = NULL)
  srv <- shiny::reactiveValues(comm_tabs_current = "ccc_tab_perturbation")

  shiny::testServer(
    env$mod_sc_communication_perturbation_server,
    args = list(comm_state = comm_state, global_data = gd, shared_rv = srv),
    {
      # Empty result clears choices
      comm_state$result <- r
      session$flushReact()

      # Switch to pathway (P1 is on row 1 which was set to NA, row 2 still has P1; P2, P9 present)
      session$setInputs(pert_target = "pathway")
      session$flushReact()
      last_up <- updates[[length(updates)]]
      expect_identical(last_up$inputId, "pert_value")
      expect_identical(last_up$choices, sort(unique(stats::na.omit(r$canonical_table$pathway))))
      expect_false(anyNA(last_up$choices))
      expect_match(output$pert_target_info, sprintf("^%d valeur\\(s\\)", length(last_up$choices)))

      # Switch to population (union of sender_node and receiver_node: A, B, C, D)
      session$setInputs(pert_target = "population")
      session$flushReact()
      last_up <- updates[[length(updates)]]
      expect_identical(last_up$inputId, "pert_value")
      expect_identical(last_up$choices, c("A", "B", "C", "D"))
      expect_match(output$pert_target_info, "^4 valeur\\(s\\)")

      # Clearing result resets pert_value choices to empty
      comm_state$result <- NULL
      session$flushReact()
      last_up <- updates[[length(updates)]]
      expect_identical(last_up$choices, character(0))
    }
  )
})

test_that("server testServer simulates multi-target, pathway, and population perturbations", {
  env <- .pert_ui_env()
  env$updateSelectizeInput <- function(...) invisible(NULL)
  r <- .ccx_result()

  comm_state <- shiny::reactiveValues(
    result = r,
    object_fingerprint = "fp_test",
    ctx_results = NULL
  )
  gd <- shiny::reactiveValues(sc_obj = "mock_sc", language = "fr", i18n = NULL)
  srv <- shiny::reactiveValues(comm_tabs_current = "ccc_tab_perturbation")

  shiny::testServer(
    env$mod_sc_communication_perturbation_server,
    args = list(comm_state = comm_state, global_data = gd, shared_rv = srv),
    {
      # 1. Multi-target ligand removal
      session$setInputs(
        pert_target = "ligand",
        pert_value = c("L1", "L2"),
        pert_mode = "remove",
        pert_compute = 1L
      )
      ctx_multi <- comm_state$ctx_results$perturbation
      expect_false(is.null(ctx_multi))
      expect_identical(ctx_multi$target, "ligand")
      expect_identical(ctx_multi$value, c("L1", "L2"))
      expect_true(isTRUE(ctx_multi$params$multi_target))
      expect_identical(ctx_multi$qc$n_removed, 5L)
      expect_length(output$pert_recap, 1L)
      expect_match(output$pert_recap, "Cible   : ligand = L1, L2", fixed = TRUE)

      # 2. Whole-pathway removal
      session$setInputs(
        pert_target = "pathway",
        pert_value = "P1",
        pert_mode = "remove",
        pert_compute = 2L
      )
      ctx_pw <- comm_state$ctx_results$perturbation
      expect_identical(ctx_pw$target, "pathway")
      expect_identical(ctx_pw$value, "P1")
      expect_identical(ctx_pw$qc$n_removed, 2L)

      # 3. Whole-population attenuation (sender + receiver roles)
      session$setInputs(
        pert_target = "population",
        pert_value = "A",
        pert_mode = "attenuate",
        pert_factor = 0.5,
        pert_compute = 3L
      )
      ctx_pop <- comm_state$ctx_results$perturbation
      expect_identical(ctx_pop$target, "population")
      expect_identical(ctx_pop$value, "A")
      expect_identical(ctx_pop$qc$n_attenuated, 3L)
    }
  )
})
