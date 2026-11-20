# =============================================================================
# test-app-rds-size-wire.R — W1.2 #3: wire ts_estimate_rds_size_mb into app.R
# =============================================================================

app_path <- file.path(ts_project_root(), "app.R")
app_lines <- readLines(app_path, warn = FALSE, encoding = "UTF-8")

extract_save_session_weight_body <- function(lines) {
  idx <- grep("output\\$save_session_weight\\s*<-\\s*renderText\\(", lines)
  if (length(idx) != 1L) {
    stop("Expected exactly one output$save_session_weight <- renderText(...) in app.R")
  }
  expr <- parse(text = lines[idx:length(lines)], n = 1)[[1]]
  # expr is `output$save_session_weight <- renderText({ ... })`
  expr[[3]][[2]]
}

test_that("app.R sources memory.R and wires ts_estimate_rds_size_mb into save_session_weight", {
  expect_true(any(grepl('source\\("tools/pkg_skeleton/R/memory\\.R"\\)', app_lines)))

  body_txt <- paste(deparse(extract_save_session_weight_body(app_lines)), collapse = "\n")
  expect_match(body_txt, "ts_estimate_rds_size_mb\\s*\\(")
  expect_match(body_txt, 'exists\\("ts_estimate_rds_size_mb"')
})

test_that("save_session_weight uses calibrated ts_estimate_rds_size_mb (~0.22x) when available", {
  source_project_file("tools/pkg_skeleton/R/memory.R")
  body_expr <- extract_save_session_weight_body(app_lines)

  # ~8 MB raw numeric vector -> calibrated ~1.76 MB (vs ~8 MB raw)
  big_vec <- numeric(1024 * 1024)
  raw_mb <- as.numeric(object.size(big_vec)) / 1024^2

  env <- new.env(parent = baseenv())
  env$`%||%` <- function(a, b) if (is.null(a)) b else a
  env$object.size <- utils::object.size
  env$ts_estimate_rds_size_mb <- ts_estimate_rds_size_mb
  env$global_data <- list(sc_obj = big_vec, sc_shared = NULL)

  out <- eval(body_expr, envir = env)
  expect_type(out, "character")
  expect_match(out, "Snapshot estimé")

  est_shown <- as.numeric(sub(".*≈\\s*([0-9]+)\\s*Mo.*", "\\1", out))
  expect_true(is.finite(est_shown))
  expect_lt(est_shown, raw_mb * 0.4)
})

test_that("save_session_weight falls back safely when ts_estimate_rds_size_mb is unavailable or returns NA", {
  body_expr <- extract_save_session_weight_body(app_lines)
  big_vec <- numeric(1024 * 1024)
  raw_mb <- as.numeric(object.size(big_vec)) / 1024^2

  # Case 1: function absent from environment
  env_missing <- new.env(parent = baseenv())
  env_missing$`%||%` <- function(a, b) if (is.null(a)) b else a
  env_missing$object.size <- utils::object.size
  env_missing$global_data <- list(sc_obj = big_vec, sc_shared = NULL)

  out_missing <- eval(body_expr, envir = env_missing)
  expect_match(out_missing, sprintf("%.0f Mo", raw_mb))

  # Case 2: function returns NA_real_
  env_na <- new.env(parent = baseenv())
  env_na$`%||%` <- function(a, b) if (is.null(a)) b else a
  env_na$object.size <- utils::object.size
  env_na$ts_estimate_rds_size_mb <- function(...) NA_real_
  env_na$global_data <- list(sc_obj = big_vec, sc_shared = NULL)

  out_na <- eval(body_expr, envir = env_na)
  expect_match(out_na, sprintf("%.0f Mo", raw_mb))
})

test_that("save_session_weight reports empty snapshot when no data or results exist", {
  source_project_file("tools/pkg_skeleton/R/memory.R")
  body_expr <- extract_save_session_weight_body(app_lines)

  env <- new.env(parent = baseenv())
  env$`%||%` <- function(a, b) if (is.null(a)) b else a
  env$object.size <- utils::object.size
  env$ts_estimate_rds_size_mb <- ts_estimate_rds_size_mb
  env$global_data <- list(sc_shared = NULL)

  out <- eval(body_expr, envir = env)
  expect_match(out, "Snapshot : vide")
})
