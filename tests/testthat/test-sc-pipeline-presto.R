# =============================================================================
# test-sc-pipeline-presto.R — W2.3 wiring of sc_fast_find_all_markers
#   Verifies that tools/pkg_skeleton/R/sc_pipeline.R and modules/sc/mod_sc.R
#   route FindAllMarkers calls through sc_fast_find_all_markers() with a
#   guarded fallback to Seurat::FindAllMarkers / FindAllMarkers.
# =============================================================================

test_that("sc_pipeline.R wires sc_fast_find_all_markers with guarded fallback", {
  p <- ts_ast_parse("tools/pkg_skeleton/R/sc_pipeline.R")

  fast_calls <- ts_ast_find_all(p, function(x) ts_ast_is_call_to(x, "sc_fast_find_all_markers"))
  expect_length(fast_calls, 1L)

  fam_calls <- ts_ast_find_all(p, function(x) {
    ts_ast_is_call_to(x, "FindAllMarkers") || ts_ast_is_call_to(x, "Seurat::FindAllMarkers")
  })
  expect_length(fam_calls, 1L)

  # Verify expected arguments on sc_fast_find_all_markers call
  call_args <- as.list(fast_calls[[1L]])
  expect_identical(ts_ast_deparse(call_args[["only.pos"]]), "TRUE")
  expect_identical(ts_ast_deparse(call_args[["min.pct"]]), "0.1")
  expect_identical(ts_ast_deparse(call_args[["logfc.threshold"]]), "0.25")
  expect_identical(ts_ast_deparse(call_args[["verbose"]]), "FALSE")

  # Verify exists("sc_fast_find_all_markers", mode = "function") guard wraps the call
  if_guards <- ts_ast_find_all(p, function(x) {
    is.call(x) && identical(as.character(x[[1L]])[1L], "if") &&
      grepl('exists("sc_fast_find_all_markers", mode = "function")',
            ts_ast_deparse(x[[2L]]), fixed = TRUE)
  })
  expect_length(if_guards, 1L)
  expect_true(length(if_guards[[1L]]) >= 4L)
})

test_that("mod_sc.R wires sc_fast_find_all_markers with guarded fallback", {
  p <- ts_ast_parse("modules/sc/mod_sc.R")

  fast_calls <- ts_ast_find_all(p, function(x) ts_ast_is_call_to(x, "sc_fast_find_all_markers"))
  expect_length(fast_calls, 1L)

  fam_calls <- ts_ast_find_all(p, function(x) {
    ts_ast_is_call_to(x, "FindAllMarkers") || ts_ast_is_call_to(x, "Seurat::FindAllMarkers")
  })
  expect_length(fam_calls, 1L)

  # Verify input parameters are forwarded to sc_fast_find_all_markers
  call_args <- as.list(fast_calls[[1L]])
  expect_identical(ts_ast_deparse(call_args[["test.use"]]), "input$marker_test")
  expect_identical(ts_ast_deparse(call_args[["min.pct"]]), "input$marker_min_pct")
  expect_identical(ts_ast_deparse(call_args[["logfc.threshold"]]), "input$marker_logfc")
  expect_identical(ts_ast_deparse(call_args[["only.pos"]]), "input$only_pos")
  expect_identical(ts_ast_deparse(call_args[["verbose"]]), "input$verbose")

  # Verify exists("sc_fast_find_all_markers", mode = "function") guard wraps the call
  if_guards <- ts_ast_find_all(p, function(x) {
    is.call(x) && identical(as.character(x[[1L]])[1L], "if") &&
      grepl('exists("sc_fast_find_all_markers", mode = "function")',
            ts_ast_deparse(x[[2L]]), fixed = TRUE)
  })
  expect_length(if_guards, 1L)
  expect_true(length(if_guards[[1L]]) >= 4L)
})
