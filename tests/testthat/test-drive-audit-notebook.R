# =============================================================================
# test-drive-audit-notebook.R — BK-12 : lecture JSONL & export Quarto/Rmd
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")

test_that("BK-12 : ts_drive_read_audit gere fichier manquant et fichier vide", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile("audit_nb_")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  expected_cols <- c("timestamp", "seq", "tool", "status", "summary", "raw_json")

  missing_path <- file.path(tmp, "does_not_exist.jsonl")
  df_missing <- ts_drive_read_audit(missing_path)
  expect_s3_class(df_missing, "data.frame")
  expect_identical(names(df_missing), expected_cols)
  expect_identical(nrow(df_missing), 0L)

  empty_path <- file.path(tmp, "empty.jsonl")
  file.create(empty_path)
  df_empty <- ts_drive_read_audit(empty_path)
  expect_s3_class(df_empty, "data.frame")
  expect_identical(names(df_empty), expected_cols)
  expect_identical(nrow(df_empty), 0L)

  ws_path <- file.path(tmp, "whitespace.jsonl")
  writeLines(c("", "   ", "\t"), ws_path)
  df_ws <- ts_drive_read_audit(ws_path)
  expect_identical(nrow(df_ws), 0L)
})

test_that("BK-12 : ts_drive_read_audit parse JSONL valide et tolere les lignes malformees", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile("audit_nb_")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  p <- file.path(tmp, "agent_audit.jsonl")
  lines <- c(
    '{"timestamp":"2026-11-21T10:00:00Z","seq":1,"tool":"run_de","status":"done","summary":"DE ok | 12 genes"}',
    "{not valid json",
    '{"timestamp":"2026-11-21T10:01:00Z","seq":2,"action":"load_data","status":"done","module":"import","elapsed_s":0.5}',
    '{"timestamp":"2026-11-21T10:02:00Z","seq":3,"command":"export_tsv","status":"error","error":"disk full"}'
  )
  writeLines(lines, p, useBytes = TRUE)

  expect_warning(
    df <- ts_drive_read_audit(p),
    "Skipping malformed audit JSONL line 2"
  )
  expect_identical(nrow(df), 3L)
  expect_identical(df$seq, c(1L, 2L, 3L))
  expect_identical(df$tool, c("run_de", "load_data", "export_tsv"))
  expect_identical(df$status, c("done", "done", "error"))
  expect_identical(df$summary[1], "DE ok | 12 genes")
  expect_match(df$summary[2], "module=import")
  expect_match(df$summary[3], "error=disk full")
  expect_identical(df$raw_json[1], lines[1])
})

test_that("BK-12 : ts_drive_read_audit respecte limit et rejette limit invalide", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile("audit_nb_")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  p <- file.path(tmp, "agent_audit.jsonl")
  for (i in 1:5) {
    ts_drive_append_audit(
      list(seq = i, tool = paste0("tool_", i), status = "done",
           timestamp = sprintf("2026-11-21T10:0%d:00Z", i)),
      audit_path = p
    )
  }

  df2 <- ts_drive_read_audit(p, limit = 2)
  expect_identical(nrow(df2), 2L)
  expect_identical(df2$seq, c(4L, 5L))
  expect_identical(df2$tool, c("tool_4", "tool_5"))

  expect_error(ts_drive_read_audit(p, limit = 0), "positive integer")
  expect_error(ts_drive_read_audit(p, limit = -1), "positive integer")
})

test_that("BK-12 : ts_drive_export_audit_qmd genere .qmd et .Rmd en memoire et sur disque", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile("audit_nb_")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  p <- file.path(tmp, "agent_audit.jsonl")
  ts_drive_append_audit(
    list(seq = 1L, tool = "qc_pca", status = "done",
         summary = "pipe|escaped", timestamp = "2026-11-21T09:00:00Z"),
    audit_path = p
  )
  ts_drive_append_audit(
    list(seq = 2L, action = "de_limma", status = "done",
         summary = "ok", timestamp = "2026-11-21T09:05:00Z"),
    audit_path = p
  )

  res_qmd <- ts_drive_export_audit_qmd(p, out_path = NULL, format = "qmd")
  expect_null(res_qmd$path)
  expect_identical(res_qmd$format, "qmd")
  expect_identical(res_qmd$n_events, 2L)
  qmd_txt <- paste(res_qmd$content, collapse = "\n")
  expect_match(qmd_txt, "format: html", fixed = TRUE)
  expect_match(qmd_txt, "Total events**: 2", fixed = TRUE)
  expect_match(qmd_txt, "2026-11-21T09:00:00Z .. 2026-11-21T09:05:00Z", fixed = TRUE)
  expect_match(qmd_txt, "pipe\\|escaped", fixed = TRUE)
  expect_match(qmd_txt, "library(jsonlite)", fixed = TRUE)

  rmd_out <- file.path(tmp, "nested", "audit_report.Rmd")
  res_rmd <- ts_drive_export_audit_qmd(
    p,
    out_path = rmd_out,
    title = "Custom Audit",
    format = "rmd"
  )
  expect_identical(res_rmd$path, rmd_out)
  expect_identical(res_rmd$format, "rmd")
  expect_true(file.exists(rmd_out))
  disk_lines <- readLines(rmd_out, warn = FALSE, encoding = "UTF-8")
  expect_identical(disk_lines, res_rmd$content)
  rmd_txt <- paste(disk_lines, collapse = "\n")
  expect_match(rmd_txt, 'title: "Custom Audit"', fixed = TRUE)
  expect_match(rmd_txt, "output: html_document", fixed = TRUE)

  empty_qmd <- ts_drive_export_audit_qmd(file.path(tmp, "missing.jsonl"), format = "qmd")
  expect_identical(empty_qmd$n_events, 0L)
  expect_match(paste(empty_qmd$content, collapse = "\n"), "_No audit events recorded._", fixed = TRUE)
})
