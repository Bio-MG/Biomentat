# =============================================================================
# test-mod-bulk-de-cancor.R — W2.3 #3a check_design_cancor() in mod_bulk_de_run
# =============================================================================

.MBC_FILE   <- "modules/bulk_de/mod_bulk_de_run.R"
.MBC_BUTTON <- "bulk-de-run_de"

.mbc_env <- function(metadata, covariates = character(0)) {
  e <- new.env(parent = globalenv())
  sys.source(file.path(ts_project_root(), "tools/pkg_skeleton/R/drive_allowlist.R"), envir = e)
  sys.source(file.path(ts_project_root(), "tools/pkg_skeleton/R/validation.R"), envir = e)

  e$`%||%` <- function(a, b) if (is.null(a)) b else a
  e$calls <- list()
  e$notifs <- list()

  e$ts_drive_job_busy   <- function() TRUE
  e$ts_drive_job_state  <- function() list(button = .MBC_BUTTON)
  e$ts_drive_job_finish <- function(button, status = "done", error = NULL) {
    e$calls[[length(e$calls) + 1L]] <- list(button = button, status = status, error = error)
    invisible(TRUE)
  }
  e$req <- function(...) invisible(TRUE)
  e$drive_counter <- function() 0L

  e$input <- list(
    run_de = 1L,
    condition_col = "cond",
    group_ref = "a",
    group_target = "b",
    de_engine = "deseq2",
    covariates = covariates,
    no_rep_enable = FALSE,
    no_rep_attest = FALSE,
    no_rep_bcv = 0.4,
    padj_method = "BH",
    shrink_lfc = FALSE,
    contrast_name = "",
    padj_thresh = 0.05,
    lfc_thresh = 1
  )
  e$session <- list()
  e$global_data <- list(
    bulk_obj = list(metadata = metadata),
    i18n = NULL
  )
  e$shared_rv <- list2env(list(
    filtered_counts = matrix(10L, nrow = 4L, ncol = nrow(metadata),
                             dimnames = list(paste0("g", 1:4), rownames(metadata))),
    contrasts = list(),
    active_contrast = NULL,
    active_condition_col = NULL,
    lfc_thresh = 1,
    padj_thresh = 0.05,
    dds_full = NULL
  ), parent = emptyenv())
  e$helpers <- list(
    design_str = function() "~ cond",
    register_contrast = function(name, res, engine = NULL) {
      e$shared_rv$contrasts[[name]] <- res
    },
    remember_padj_ctx = function(...) invisible(NULL)
  )
  e$TS_PADJ_METHOD_DEFAULT <- "BH"
  e$build_dds <- function(...) structure(list(fake = TRUE), class = "dds")
  e$.normalize_de_cols <- function(res, ...) res
  e$updateSelectInput <- function(...) invisible(NULL)
  e$showNotification <- function(ui, type = "message", duration = NULL, ...) {
    e$notifs[[length(e$notifs) + 1L]] <- list(ui = as.character(ui), type = type)
    invisible(NULL)
  }
  e$removeNotification <- function(...) invisible(NULL)
  e$.tr <- function(key) key
  e$.t_fmt <- function(fmt, ...) {
    args <- list(...)
    out <- fmt
    for (nm in names(args)) {
      out <- gsub(paste0("{", nm, "}"), as.character(args[[nm]]), out, fixed = TRUE)
    }
    out
  }
  e$run_bulk_de_dispatch <- function(...) {
    data.frame(
      log2FoldChange = c(2, -2, 0.1, 3),
      padj           = c(0.001, 0.002, 0.9, 0.003),
      baseMean       = rep(10, 4L)
    )
  }
  e$ensure_app_daemons <- function() FALSE
  e$de_drive_job <- new.env(parent = emptyenv())
  e$de_drive_job$cell <- NULL
  e$de_drive_task <- list(invoke = function(...) invisible(NULL), status = function() NULL, result = function() NULL)
  e$run_de_core <- ts_ast_assignment("tools/pkg_skeleton/R/bulk_helpers.R", "run_de_core", e)
  e$.de_drive_commit <- tryCatch(ts_ast_assignment(.MBC_FILE, ".de_drive_commit", e), error = function(err) NULL)
  e$.de_drive_fail   <- tryCatch(ts_ast_assignment(.MBC_FILE, ".de_drive_fail", e), error = function(err) NULL)
  e
}

.mbc_run <- function(envir) {
  f <- function() NULL
  body(f) <- ts_ast_observe_block(.MBC_FILE, "drive_trigger")
  environment(f) <- envir
  if (requireNamespace("shiny", quietly = TRUE)) {
    shiny::withReactiveDomain(shiny::MockShinySession$new(), f())
  } else {
    f()
  }
}

test_that("mod_bulk_de_run warns (non-blocking) when numeric covariate is collinear via check_design_cancor", {
  meta <- data.frame(
    cond = factor(c("a", "a", "a", "b", "b", "b")),
    age  = c(20, 21, 22, 60, 61, 62),
    row.names = paste0("s", 1:6),
    stringsAsFactors = FALSE
  )
  e <- .mbc_env(meta, covariates = "age")
  .mbc_run(e)

  warn_msgs <- vapply(
    Filter(function(n) identical(n$type, "warning"), e$notifs),
    function(n) n$ui,
    character(1)
  )
  expect_true(any(grepl("age", warn_msgs, fixed = TRUE)))
  expect_true(any(grepl("cor", warn_msgs, fixed = TRUE)))
  expect_identical(e$calls[[1]]$status, "done")
})

test_that("mod_bulk_de_run emits no cancor warning when numeric covariate is uncorrelated", {
  meta <- data.frame(
    cond = factor(c("a", "a", "a", "b", "b", "b")),
    rnd  = c(10, 50, 20, 12, 48, 22),
    row.names = paste0("s", 1:6),
    stringsAsFactors = FALSE
  )
  e <- .mbc_env(meta, covariates = "rnd")
  .mbc_run(e)

  warn_msgs <- vapply(
    Filter(function(n) identical(n$type, "warning"), e$notifs),
    function(n) n$ui,
    character(1)
  )
  expect_false(any(grepl("rnd", warn_msgs, fixed = TRUE)))
  expect_identical(e$calls[[1]]$status, "done")
})
