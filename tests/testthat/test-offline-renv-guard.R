# =============================================================================
# test-offline-renv-guard.R -- W3.2 : garde renv / CRAN en mode Offline strict
# =============================================================================

source_project_file("tools/pkg_skeleton/R/defaults.R")

.w32_env <- new.env(parent = globalenv())
.w32_env$ts_is_offline_strict <- ts_is_offline_strict
.w32_env$ts_require_online <- ts_require_online

.w32_enforce_guard <- ts_ast_assignment(
  "global.R", "ts_enforce_renv_offline_guard", eval_env = .w32_env
)
.w32_restore_guarded <- ts_ast_assignment(
  "global.R", "ts_renv_restore_guarded", eval_env = .w32_env
)

test_that("W3.2 : ts_require_online('renv_cran') leve biomentat_offline_error en mode strict", {
  withr::local_options(list(biomentat.offline_strict = NULL))
  withr::local_envvar(c(BIOMENTAT_OFFLINE_STRICT = ""))
  expect_true(isTRUE(ts_require_online("renv_cran")))

  withr::with_options(list(biomentat.offline_strict = TRUE), {
    err <- expect_error(
      ts_require_online("renv_cran", "renv::restore() / CRAN egress"),
      class = "biomentat_offline_error"
    )
    expect_identical(err$egress_id, "renv_cran")
  })

  withr::with_envvar(c(BIOMENTAT_OFFLINE_STRICT = "1"), {
    err <- expect_error(
      ts_require_online("renv_cran"),
      class = "biomentat_offline_error"
    )
    expect_identical(err$egress_id, "renv_cran")
  })
})

test_that("W3.2 : ts_enforce_renv_offline_guard() configure renv uniquement si strict offline est actif", {
  withr::local_options(list(
    biomentat.offline_strict = NULL,
    renv.config.auto.snapshot = NULL,
    renv.config.synchronized.check = NULL,
    renv.config.repos.override = NULL,
    repos = getOption("repos")
  ))
  withr::local_envvar(c(BIOMENTAT_OFFLINE_STRICT = ""))

  expect_false(isTRUE(.w32_enforce_guard()))
  expect_null(getOption("renv.config.auto.snapshot"))
  expect_null(getOption("renv.config.synchronized.check"))
  expect_null(getOption("renv.config.repos.override"))

  withr::with_options(list(biomentat.offline_strict = TRUE), {
    expect_true(isTRUE(.w32_enforce_guard()))
    expect_identical(getOption("renv.config.auto.snapshot"), FALSE)
    expect_identical(getOption("renv.config.synchronized.check"), FALSE)
    expect_identical(getOption("renv.config.repos.override"), character(0))
    expect_identical(getOption("repos"), c(CRAN = ""))
  })
})

test_that("W3.2 : ts_enforce_renv_offline_guard() repond aussi a BIOMENTAT_OFFLINE_STRICT=1", {
  withr::local_options(list(
    biomentat.offline_strict = NULL,
    renv.config.auto.snapshot = NULL,
    renv.config.synchronized.check = NULL,
    renv.config.repos.override = NULL,
    repos = getOption("repos")
  ))
  withr::local_envvar(c(BIOMENTAT_OFFLINE_STRICT = "1"))

  expect_true(isTRUE(.w32_enforce_guard()))
  expect_identical(getOption("renv.config.auto.snapshot"), FALSE)
  expect_identical(getOption("renv.config.synchronized.check"), FALSE)
  expect_identical(getOption("renv.config.repos.override"), character(0))
  expect_identical(getOption("repos"), c(CRAN = ""))
})

test_that("W3.2 : ts_renv_restore_guarded() bloque avec biomentat_offline_error en mode strict", {
  withr::local_options(list(biomentat.offline_strict = TRUE))
  withr::local_envvar(c(BIOMENTAT_OFFLINE_STRICT = ""))

  err <- expect_error(
    .w32_restore_guarded(prompt = FALSE),
    class = "biomentat_offline_error"
  )
  expect_identical(err$egress_id, "renv_cran")

  withr::with_options(list(biomentat.offline_strict = NULL), {
    withr::with_envvar(c(BIOMENTAT_OFFLINE_STRICT = "1"), {
      err_env <- expect_error(
        .w32_restore_guarded(prompt = FALSE),
        class = "biomentat_offline_error"
      )
      expect_identical(err_env$egress_id, "renv_cran")
    })
  })
})

test_that("W3.2 : global.R appelle ts_enforce_renv_offline_guard() au demarrage", {
  lines <- readLines(file.path(ts_project_root(), "global.R"), warn = FALSE)
  expect_true(any(grepl("^\\s*ts_enforce_renv_offline_guard\\(\\s*\\)\\s*$", lines)))
})
