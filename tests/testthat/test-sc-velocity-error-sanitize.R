# =============================================================================
# test-sc-velocity-error-sanitize.R — W1.1 #1 : assainissement de la branche
# error de scvelo_task (NUL / UTF-8 invalide dans conditionMessage(err))
# =============================================================================
# Cible : modules/sc/mod_sc_velocity.R, observateur `observeEvent(scvelo_task$status(), ...)`,
# branche `else if (identical(st, "error"))`.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")

.vel_err_extract_helpers_env <- function() {
  src_lines <- readLines(
    file.path(ts_project_root(), "modules/sc/mod_sc_velocity.R"),
    warn = FALSE, encoding = "UTF-8"
  )
  idx_start <- grep("^\\s*\\.vel_sanitize_text\\s*<-\\s*function", src_lines)[1]
  idx_end <- grep("^\\s*\\.vel_scvelo_fail\\s*<-\\s*function", src_lines)[1]
  stopifnot(!is.na(idx_start), !is.na(idx_end), idx_end > idx_start)

  env <- new.env(parent = globalenv())
  eval(parse(text = src_lines[idx_start:(idx_end - 1L)]), envir = env)
  env
}

.vel_err_extract_status_handler_expr <- function() {
  exprs <- parse(
    file.path(ts_project_root(), "modules/sc/mod_sc_velocity.R"),
    encoding = "UTF-8"
  )
  found <- NULL
  walk <- function(node) {
    if (!is.null(found) || !is.call(node)) return(invisible(NULL))
    if (identical(node[[1L]], as.name("observeEvent")) && length(node) >= 3L) {
      dep <- paste(deparse(node[[2L]]), collapse = "")
      if (grepl("scvelo_task$status()", dep, fixed = TRUE)) {
        found <<- node[[3L]]
        return(invisible(NULL))
      }
    }
    for (i in seq_along(node)) walk(node[[i]])
  }
  for (e in exprs) walk(e)
  stopifnot(!is.null(found))
  found
}

test_that("W1.1 #1 : .vel_sanitize_text nettoie les octets UTF-8 invalides", {
  env <- .vel_err_extract_helpers_env()

  bad_utf8 <- rawToChar(as.raw(c(0x50, 0x79, 0x20, 0xff, 0xfe, 0x20, 0x65, 0x72, 0x72)))
  Encoding(bad_utf8) <- "UTF-8"
  expect_false(validUTF8(bad_utf8))

  cleaned <- env$.vel_sanitize_text(bad_utf8)
  expect_true(is.character(cleaned))
  expect_length(cleaned, 1L)
  expect_true(validUTF8(cleaned))
  expect_match(cleaned, "Py ", fixed = TRUE)
  expect_match(cleaned, " err", fixed = TRUE)
})

test_that("W1.1 #1 : la branche error de scvelo_task assainit conditionMessage(err) avant .vel_scvelo_fail", {
  env <- .vel_err_extract_helpers_env()
  handler_expr <- .vel_err_extract_status_handler_expr()

  bad_msg <- rawToChar(as.raw(c(0x73, 0x63, 0x56, 0x65, 0x6c, 0x6f, 0x20, 0xff, 0xfe)))
  Encoding(bad_msg) <- "UTF-8"
  expect_false(validUTF8(bad_msg))

  captured <- NULL
  env$.vel_scvelo_fail <- function(e) {
    captured <<- e
    invisible(NULL)
  }
  env$scvelo_task <- list(
    status = function() "error",
    result = function() stop(errorCondition(bad_msg, class = c("scvelo_py_error", "error", "condition")))
  )

  eval(handler_expr, envir = env)

  expect_true(inherits(captured, "error"))
  expect_true(inherits(captured, "scvelo_py_error"))
  expect_true(validUTF8(conditionMessage(captured)))
  expect_match(conditionMessage(captured), "scVelo", fixed = TRUE)
})

test_that("W1.1 #1 : la branche error de scvelo_task assainit aussi une valeur non-condition", {
  env <- .vel_err_extract_helpers_env()
  handler_expr <- .vel_err_extract_status_handler_expr()

  bad_val <- rawToChar(as.raw(c(0x66, 0x61, 0x69, 0x6c, 0x20, 0xff, 0xfe)))
  Encoding(bad_val) <- "UTF-8"

  captured <- NULL
  env$.vel_scvelo_fail <- function(e) {
    captured <<- e
    invisible(NULL)
  }
  env$scvelo_task <- list(
    status = function() "error",
    result = function() list(detail = bad_val)
  )

  eval(handler_expr, envir = env)

  expect_true(inherits(captured, "error"))
  expect_true(validUTF8(conditionMessage(captured)))
  expect_match(conditionMessage(captured), "fail", fixed = TRUE)
})
