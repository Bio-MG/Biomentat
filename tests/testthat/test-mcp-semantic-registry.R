# =============================================================================
# test-mcp-semantic-registry.R — BK-12 Formal MCP Semantic Tool Schema Tests
# =============================================================================

.mcp_sem_env <- function() {
  e <- new.env(parent = globalenv())
  old <- commandArgs
  # Load scripts/mcp_server.R without triggering --check or .ts_serve()
  src_lines <- readLines(file.path(ts_project_root(), "scripts", "mcp_server.R"),
                         encoding = "UTF-8", warn = FALSE)
  # Strip final .ts_serve() call so sourcing only loads definitions
  serve_idx <- which(grepl("^\\.ts_serve\\(\\)", src_lines))
  if (length(serve_idx)) src_lines <- src_lines[-serve_idx]
  tmp <- tempfile(fileext = ".R")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(src_lines, tmp, useBytes = TRUE)

  # Provide --file= via mocked commandArgs in local env
  e$.cmd_args_override <- c("Rscript", paste0("--file=", file.path(ts_project_root(), "scripts", "mcp_server.R")))
  e$commandArgs <- function(trailingOnly = FALSE) {
    if (isTRUE(trailingOnly)) character(0L) else e$.cmd_args_override
  }
  sys.source(tmp, envir = e)
  e
}

test_that("BK-12 semantic registry registers all 12 tools cleanly", {
  e <- .mcp_sem_env()
  expect_identical(e$.ts_mcp_semantic_problems(), character(0L))
  tools <- e$.ts_tools()
  expect_length(tools, 12L)
  for (tl in tools) {
    expect_true("annotations" %in% names(tl))
    expect_true("x_semantic" %in% names(tl))
    expect_identical(tl$x_semantic$schema_version, "1.0.0")
    expect_true(tl$x_semantic$category %in% e$TS_MCP_SEMANTIC_CATEGORIES)
    expect_true(length(tl$x_semantic$intents) >= 1L)
  }
})

test_that("BK-12 .ts_register_tool enforces schema, handler, and readOnlyHint invariants", {
  e <- .mcp_sem_env()
  valid_ann <- list(title = "T", readOnlyHint = TRUE, destructiveHint = FALSE,
                    idempotentHint = TRUE, openWorldHint = FALSE)
  valid_sem <- list(category = "observation", modalities = "session",
                    intents = "test_intent", side_effects = "none",
                    requires_armed = FALSE, requires_session = FALSE)
  valid_schema <- list(type = "object", properties = list(), additionalProperties = FALSE)

  expect_error(
    e$.ts_register_tool("", "desc", valid_schema, function() {}, valid_ann, valid_sem),
    "non-empty scalar string"
  )
  expect_error(
    e$.ts_register_tool("x", "desc", list(type = "array"), function() {}, valid_ann, valid_sem),
    "object schema"
  )
  expect_error(
    e$.ts_register_tool("x", "desc", valid_schema, "not_a_fn", valid_ann, valid_sem),
    "must be a function"
  )
  bad_sem <- valid_sem
  bad_sem$category <- "unknown_cat"
  expect_error(
    e$.ts_register_tool("x", "desc", valid_schema, function() {}, valid_ann, bad_sem),
    "invalid semantic category"
  )
  # Contradiction: readOnlyHint = TRUE with side_effects != "none"
  bad_side <- valid_sem
  bad_side$side_effects <- "writes_scenario_json"
  expect_error(
    e$.ts_register_tool("x", "desc", valid_schema, function() {}, valid_ann, bad_side),
    "contradicts side_effects"
  )
})

test_that("BK-12 .ts_resolve_semantic_tools filters by read_only, category, modality, and intent", {
  e <- .mcp_sem_env()
  ro <- e$.ts_resolve_semantic_tools(read_only = TRUE)
  expect_length(ro, 5L)
  ro_names <- vapply(ro, function(x) x$name, character(1L))
  expect_setequal(ro_names, c(
    "transcripto_drive_status", "transcripto_drive_read_result",
    "transcripto_drive_snapshot", "transcripto_drive_state",
    "transcripto_drive_datasets"
  ))

  exec_tools <- e$.ts_resolve_semantic_tools(category = "execution")
  expect_length(exec_tools, 1L)
  expect_identical(exec_tools[[1L]]$name, "transcripto_drive_run")

  spatial_tools <- e$.ts_resolve_semantic_tools(modality = "spatial")
  expect_true(length(spatial_tools) >= 7L)

  by_intent <- e$.ts_resolve_semantic_tools(intent = "import_dataset")
  expect_length(by_intent, 1L)
  expect_identical(by_intent[[1L]]$name, "transcripto_drive_import")
})

test_that("scripts/mcp_server.R --check passes with semantic check OK", {
  rscript <- file.path(R.home("bin"), "Rscript.exe")
  if (!file.exists(rscript)) rscript <- "Rscript"
  srv <- file.path(ts_project_root(), "scripts", "mcp_server.R")
  out <- suppressWarnings(system2(rscript, c("--no-init-file", srv, "--check"),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  expect_true(is.null(status) || identical(status, 0L))
  expect_true(any(grepl("semantic check: OK", out, fixed = TRUE)))
})
