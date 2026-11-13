# =============================================================================
# test-flaw-vfi-drive.R — P2.3 : FLAW V-FI (sans fileInput) & Drive sans navigateur
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")

test_that("P2.3 V-FI : ts_velocity_resolve_loom_path() resout et confine les chemins", {
  root <- withr::local_tempdir()
  loom_f <- file.path(root, "sample.loom")
  h5_f   <- file.path(root, "sample.h5ad")
  txt_f  <- file.path(root, "bad.exe")
  writeLines("dummy", loom_f)
  writeLines("dummy", h5_f)
  writeLines("dummy", txt_f)

  outside_dir <- withr::local_tempdir()
  outside_f <- file.path(outside_dir, "outside.loom")
  writeLines("dummy", outside_f)

  expect_identical(
    normalizePath(ts_velocity_resolve_loom_path(loom_f, roots = root), winslash = "/"),
    normalizePath(loom_f, winslash = "/")
  )
  withr::with_envvar(c(TS_VELOCITY_LOOM_PATH = h5_f), {
    expect_identical(
      normalizePath(ts_velocity_resolve_loom_path(roots = root), winslash = "/"),
      normalizePath(h5_f, winslash = "/")
    )
  })

  expect_error(
    ts_velocity_resolve_loom_path(file.path(root, "..", basename(outside_f)), roots = root),
    class = "velocity_validation_error"
  )
  expect_error(
    ts_velocity_resolve_loom_path(outside_f, roots = root),
    class = "velocity_validation_error"
  )
  expect_error(
    ts_velocity_resolve_loom_path(txt_f, roots = root),
    class = "velocity_validation_error"
  )
})

test_that("P2.3 Drive : TRANSCRIPTO_DRIVE_ROOT isole ts_drive_root()/ts_drive_path()", {
  iso_root <- withr::local_tempdir()
  withr::local_envvar(c(TRANSCRIPTO_DRIVE_ROOT = iso_root))
  .ts_drive_state$root <- NULL
  expect_identical(
    normalizePath(ts_drive_root(), winslash = "/"),
    normalizePath(iso_root, winslash = "/")
  )
  expect_true(startsWith(
    normalizePath(ts_drive_path("state.json"), winslash = "/", mustWork = FALSE),
    normalizePath(iso_root, winslash = "/")
  ))
})

test_that("P2.3 Option B : tools/ws_drive_client.mjs existe et passe node --check", {
  ws_script <- file.path(ts_project_root(), "tools", "ws_drive_client.mjs")
  expect_true(file.exists(ws_script))
  node_bin <- Sys.which("node")
  skip_if(!nzchar(node_bin), "node absent du PATH")
  if (file.exists(ws_script)) {
    out <- system2(node_bin, args = c("--check", shQuote(ws_script)),
                   stdout = TRUE, stderr = TRUE)
    code <- attr(out, "status") %||% 0L
    expect_identical(as.integer(code), 0L, info = paste(out, collapse = "\n"))
  }
})
