# =============================================================================
# test-ccc-lot2-ui-b.R — W1.3 Lot 2 §11 remaining UI polish tests
# =============================================================================

.lot2b_read_comm_src <- function() {
  repo_root <- if (file.exists("modules/sc/mod_sc_communication.R")) "." else "../.."
  paste(
    readLines(file.path(repo_root, "modules/sc/mod_sc_communication.R"),
              warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
}

test_that("mod_sc_communication_ui wraps config controls in comm_config_accordion", {
  src <- .lot2b_read_comm_src()

  expect_match(src, 'id = ns("comm_config_accordion")', fixed = TRUE)
  expect_match(src, 'open = "comm_config_panel"', fixed = TRUE)
  expect_match(src, 'value = "comm_config_panel"', fixed = TRUE)

  pos_acc <- regexpr('id = ns("comm_config_accordion")', src, fixed = TRUE)[1]
  pos_src <- regexpr('radioButtons(ns("comm_source")', src, fixed = TRUE)[1]
  pos_status <- regexpr('textOutput(ns("comm_status"))', src, fixed = TRUE)[1]
  expect_gt(pos_acc, 0L)
  expect_gt(pos_src, pos_acc)
  expect_gt(pos_status, pos_src)
})

test_that("mod_sc_communication_server auto-collapses config panel on Calcul OK / Import OK", {
  src <- .lot2b_read_comm_src()

  expect_match(src, '.collapse_config_panel <- function()', fixed = TRUE)
  expect_match(src, 'bslib::accordion_panel_close(', fixed = TRUE)
  expect_match(src, '"comm_config_accordion"', fixed = TRUE)
  expect_match(src, 'values = "comm_config_panel"', fixed = TRUE)

  n_calls <- lengths(regmatches(src, gregexpr('.collapse_config_panel()', src, fixed = TRUE)))[1]
  # 3 call sites (.apply_cellchat_result, .apply_liana_result, comm_import)
  expect_gte(n_calls, 3L)
})

test_that("mod_sc_communication_server toggles comm_f_pmax and validates comm_centrality_plot", {
  src <- .lot2b_read_comm_src()

  expect_match(src, 'f_pmax_id <- "comm_f_pmax"', fixed = TRUE)
  expect_match(src, '.toggle_filter(f_pmax_id,', fixed = TRUE)

  # comm_centrality_plot uses validate(need(...)) instead of returning a bare string
  pos_cp <- regexpr('output$comm_centrality_plot <- renderPlot({', src, fixed = TRUE)[1]
  pos_cq <- regexpr('output$comm_centrality_quadrant_plot <- renderPlot({', src, fixed = TRUE)[1]
  expect_gt(pos_cp, 0L)
  expect_gt(pos_cq, pos_cp)
  cp_body <- substr(src, pos_cp, pos_cq)
  expect_match(cp_body, 'validate(need(', fixed = TRUE)
})
