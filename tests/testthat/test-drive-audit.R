# =============================================================================
# test-drive-audit.R — P1.3a : journal JSONL append-only agent_audit.jsonl
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")

test_that("P1.3a : ts_drive_append_audit ecrit des lignes JSONL valides sans ecraser", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile("audit_")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  p <- file.path(tmp, "agent_audit.jsonl")

  ok1 <- ts_drive_append_audit(list(
    seq = 1L, status = "done", module = "bulk_de",
    action = "run_de", elapsed_s = 0.42
  ), audit_path = p)
  ok2 <- ts_drive_append_audit(list(
    seq = 2L, status = "invalid", module = NA_character_,
    action = NA_character_, error = "bad_seq"
  ), audit_path = p)

  expect_true(isTRUE(ok1))
  expect_true(isTRUE(ok2))
  lines <- readLines(p, warn = FALSE, encoding = "UTF-8")
  expect_length(lines, 2L)

  r1 <- jsonlite::fromJSON(lines[[1]], simplifyVector = TRUE)
  r2 <- jsonlite::fromJSON(lines[[2]], simplifyVector = TRUE)
  expect_identical(as.integer(r1$seq), 1L)
  expect_identical(r1$status, "done")
  expect_identical(r1$module, "bulk_de")
  expect_true(!is.null(r1$timestamp) && nzchar(r1$timestamp))
  expect_identical(as.integer(r2$seq), 2L)
  expect_identical(r2$status, "invalid")
})

test_that("P1.3a : ts_drive_append_audit ne leve jamais d'erreur si le chemin est non inscriptible", {
  tmp_dir <- tempfile("audit_dir_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  # Ecrire dans un chemin qui pointe vers un repertoire existant (echoue a l'ouverture)
  res <- expect_no_error(
    ts_drive_append_audit(list(seq = 1L, status = "done"), audit_path = tmp_dir)
  )
  expect_false(isTRUE(res))
})

test_that("P1.3a : ts_drive_tick appelle ts_drive_append_audit sur les 3 sites de verdict", {
  src <- paste(readLines(file.path(ts_project_root(), "tools/pkg_skeleton/R/drive_watcher.R"),
                         warn = FALSE), collapse = "\n")
  hits <- lengths(regmatches(src, gregexpr("ts_drive_append_audit(", src, fixed = TRUE)))[1]
  # 3 sites d'appel dans ts_drive_tick()
  expect_gte(hits, 3L)
})
