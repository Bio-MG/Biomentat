# test-app-onstop-daemons.R -- GARDE : app.R enregistre un onStop applicatif qui
# ferme mirai::daemons(0) sans jamais lever et nettoie le pidfile de session.
app_lines <- function() {
  root <- testthat::test_path("..", "..")
  readLines(file.path(root, "app.R"), warn = FALSE, encoding = "UTF-8")
}

test_that("app.R registers a top-level shiny::onStop before shinyApp()", {
  l <- app_lines()
  on_stop <- grep("^shiny::onStop\\(function\\(\\)", l)
  app_call <- grep("^shinyApp\\(ui, server\\)", l)
  expect_length(on_stop, 1L)
  expect_length(app_call, 1L)
  expect_lt(on_stop, app_call)
})

test_that("onStop teardown is guarded and silent", {
  txt <- paste(app_lines(), collapse = "\n")
  expect_match(txt, 'requireNamespace\\("mirai", quietly = TRUE\\)')
  expect_match(txt, "try\\(mirai::daemons\\(0\\), silent = TRUE\\)")
})

test_that("onStop teardown removes the session pidfile if present", {
  txt <- paste(app_lines(), collapse = "\n")
  expect_match(txt, "BIOMENTAT_SESSION_PIDFILE")
  expect_match(txt, "unlink\\(pidfile\\)")
})

test_that("per-session onStop(stop_spatial_daemons) is kept", {
  expect_true(any(grepl("onStop(stop_spatial_daemons)", app_lines(), fixed = TRUE)))
})
