# =============================================================================
# test-mod-sc-markers-presto.R — W2.3 wiring of sc_fast_find_all_markers()
#   into modules/sc/mod_sc_markers.R
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_helpers.R")
source_project_file("modules/sc/mod_sc_markers.R")

.SM_PRESTO_FILE <- "modules/sc/mod_sc_markers.R"

.sm_presto_mock_obj <- function() {
  if (!methods::isClass("TSMockSeuratMarkersPresto")) {
    methods::setClass("TSMockSeuratMarkersPresto", slots = c(meta.data = "data.frame"))
  }
  methods::new(
    "TSMockSeuratMarkersPresto",
    meta.data = data.frame(
      seurat_clusters = c("0", "0", "1", "1"),
      stringsAsFactors = FALSE
    )
  )
}

.sm_presto_env <- function(parent_env = globalenv(), include_fast = TRUE) {
  e <- new.env(parent = parent_env)
  e$.tr    <- function(key, ...) key
  e$.t_fmt <- ts_ast_assignment("global.R", ".t_fmt")
  e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
  e$req    <- shiny::req

  e$calls <- list()
  e$marker_log_rv <- function(x) invisible(NULL)
  e$markers_rv    <- function(x) {
    e$last_markers <- x
    invisible(NULL)
  }
  e$showNotification <- function(msg, ...) invisible(NULL)
  e$`Idents<-` <- function(obj, value) obj
  e$subsample_seurat_for_analysis <- function(obj, max_per_group, group_col) {
    list(was_subsampled = FALSE, object = obj, n_before = 4L, n_after = 4L)
  }
  e$optimize_bpcells_for_markers <- function(obj) NULL

  sample_df <- data.frame(
    gene = "G1", cluster = "0", avg_log2FC = 1.8,
    p_val = 1e-4, p_val_adj = 0.002, pct.1 = 0.9, pct.2 = 0.1,
    stringsAsFactors = FALSE
  )

  e$FindAllMarkers <- function(object, test.use = "wilcox", min.pct = 0.1,
                               logfc.threshold = 0.25, only.pos = TRUE,
                               verbose = FALSE, ...) {
    e$calls[[length(e$calls) + 1L]] <- list(
      fn = "FindAllMarkers",
      test.use = test.use,
      min.pct = min.pct,
      logfc.threshold = logfc.threshold,
      only.pos = only.pos,
      verbose = verbose
    )
    sample_df
  }

  if (isTRUE(include_fast)) {
    e$sc_fast_find_all_markers <- function(object, test.use = "wilcox",
                                           min.pct = 0.1,
                                           logfc.threshold = 0.25,
                                           only.pos = TRUE,
                                           verbose = FALSE, ...) {
      e$calls[[length(e$calls) + 1L]] <- list(
        fn = "sc_fast_find_all_markers",
        test.use = test.use,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        only.pos = only.pos,
        verbose = verbose
      )
      sample_df
    }
  }

  e$global_data <- list(i18n = NULL, sc_obj = .sm_presto_mock_obj())
  e$shared_rv   <- new.env(parent = emptyenv())
  e$shared_rv$max_cells_heavy <- Inf
  e$input <- list(
    marker_test = "wilcox",
    marker_min_pct = 0.15,
    marker_logfc = 0.30,
    sort_by = "logfc"
  )
  e$session <- list(onSessionEnded = function(f) invisible(NULL))
  e
}

test_that("mod_sc_markers.R wires sc_fast_find_all_markers with FindAllMarkers fallback", {
  p <- if (file.exists(.SM_PRESTO_FILE)) .SM_PRESTO_FILE else file.path("../..", .SM_PRESTO_FILE)
  src <- paste(readLines(p, warn = FALSE), collapse = "\n")

  expect_match(src, 'exists\\("sc_fast_find_all_markers",\\s*mode\\s*=\\s*"function"\\)')
  expect_match(src, "sc_fast_find_all_markers\\s*\\(")
  expect_match(src, "FindAllMarkers\\s*\\(")
})

test_that("run_markers dispatches to sc_fast_find_all_markers when available", {
  skip_if_not(
    exists("MockShinySession", envir = asNamespace("shiny")),
    "shiny::MockShinySession indisponible"
  )

  env <- .sm_presto_env(parent_env = globalenv(), include_fast = TRUE)
  shiny::withReactiveDomain(
    shiny::MockShinySession$new(),
    eval(ts_ast_observe_block(.SM_PRESTO_FILE, "run_markers"), envir = env)
  )

  expect_length(env$calls, 1L)
  expect_identical(env$calls[[1L]]$fn, "sc_fast_find_all_markers")
  expect_identical(env$calls[[1L]]$test.use, "wilcox")
  expect_equal(env$calls[[1L]]$min.pct, 0.15)
  expect_equal(env$calls[[1L]]$logfc.threshold, 0.30)
  expect_true(isTRUE(env$calls[[1L]]$only.pos))
  expect_false(isTRUE(env$calls[[1L]]$verbose))
  expect_true(is.data.frame(env$shared_rv$markers_data))
  expect_equal(nrow(env$shared_rv$markers_data), 1L)
})

test_that("run_markers falls back to FindAllMarkers when sc_fast_find_all_markers is absent", {
  skip_if_not(
    exists("MockShinySession", envir = asNamespace("shiny")),
    "shiny::MockShinySession indisponible"
  )

  # Parent environment without sc_fast_find_all_markers (inherits from baseenv)
  iso_parent <- new.env(parent = baseenv())
  env <- .sm_presto_env(parent_env = iso_parent, include_fast = FALSE)
  env$.normalize_marker_cols <- get(".normalize_marker_cols", envir = globalenv())

  shiny::withReactiveDomain(
    shiny::MockShinySession$new(),
    eval(ts_ast_observe_block(.SM_PRESTO_FILE, "run_markers"), envir = env)
  )

  expect_length(env$calls, 1L)
  expect_identical(env$calls[[1L]]$fn, "FindAllMarkers")
  expect_identical(env$calls[[1L]]$test.use, "wilcox")
  expect_equal(env$calls[[1L]]$min.pct, 0.15)
  expect_equal(env$calls[[1L]]$logfc.threshold, 0.30)
  expect_true(isTRUE(env$calls[[1L]]$only.pos))
  expect_false(isTRUE(env$calls[[1L]]$verbose))
})

test_that("real sc_fast_find_all_markers respects mocked FindAllMarkers in run_markers", {
  skip_if_not(
    exists("MockShinySession", envir = asNamespace("shiny")),
    "shiny::MockShinySession indisponible"
  )

  env <- .sm_presto_env(parent_env = globalenv(), include_fast = FALSE)
  env$sc_fast_find_all_markers <- get("sc_fast_find_all_markers", envir = globalenv())

  shiny::withReactiveDomain(
    shiny::MockShinySession$new(),
    eval(ts_ast_observe_block(.SM_PRESTO_FILE, "run_markers"), envir = env)
  )

  expect_length(env$calls, 1L)
  expect_identical(env$calls[[1L]]$fn, "FindAllMarkers")
  expect_true(is.data.frame(env$shared_rv$markers_data))
})
