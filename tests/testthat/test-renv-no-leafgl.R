# =============================================================================
# test-renv-no-leafgl.R -- Regression guard: leafgl is absent from renv.lock, app.R, and global.R
# =============================================================================

test_that("renv.lock parses as valid JSON and excludes leafgl", {
  lock_path <- file.path(ts_project_root(), "renv.lock")
  expect_true(file.exists(lock_path))

  parsed <- jsonlite::fromJSON(lock_path, simplifyVector = FALSE)
  expect_type(parsed, "list")
  expect_true("Packages" %in% names(parsed))
  expect_false("leafgl" %in% names(parsed$Packages))
})

test_that("app.R does not reference leafgl", {
  app_path <- file.path(ts_project_root(), "app.R")
  expect_true(file.exists(app_path))

  app_lines <- readLines(app_path, warn = FALSE, encoding = "UTF-8")
  expect_false(any(grepl("leafgl", app_lines, fixed = TRUE)))
})

test_that("global.R does not reference leafgl or has_leafgl", {
  global_path <- file.path(ts_project_root(), "global.R")
  expect_true(file.exists(global_path))

  global_lines <- readLines(global_path, warn = FALSE, encoding = "UTF-8")
  expect_false(any(grepl("leafgl", global_lines, fixed = TRUE)))
  expect_false(any(grepl("has_leafgl", global_lines, fixed = TRUE)))
})
