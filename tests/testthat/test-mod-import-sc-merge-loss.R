# =============================================================================
# test-mod-import-sc-merge-loss.R — W2.4 #2: wire summarize_merge_gene_loss()
# into modules/import/mod_import_sc.R (.sc_log_merge_gene_report)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/sc_helpers.R")
source_project_file("tools/pkg_skeleton/R/sc_pipeline.R")
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(bsicons)
  library(shinyFiles)
  library(Seurat)
})
source_project_file("tools/pkg_skeleton/R/rdata_io.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("modules/import/mod_rdata_picker.R")
source_project_file("modules/import/mod_import_sc.R")

if (!exists(".tr", envir = globalenv(), inherits = FALSE)) {
  assign(".tr", function(key) key, envir = globalenv())
}
if (!exists(".t_fmt", envir = globalenv(), inherits = FALSE)) {
  assign(".t_fmt", function(template, ...) {
    vals <- list(...)
    for (nm in names(vals)) {
      template <- gsub(paste0("{", nm, "}"), format(vals[[nm]]), template, fixed = TRUE)
    }
    template
  }, envir = globalenv())
}

.make_tiny_sc <- function(genes, cells) {
  mat <- Matrix::Matrix(1, nrow = length(genes), ncol = length(cells), sparse = TRUE,
                        dimnames = list(genes, cells))
  SeuratObject::CreateSeuratObject(counts = mat)
}

test_that("mod_import_sc wires summarize_merge_gene_loss with exists + tryCatch guards", {
  root <- ts_project_root()
  src <- paste(readLines(file.path(root, "modules/import/mod_import_sc.R"), warn = FALSE),
               collapse = "\n")

  expect_match(src, 'exists\\("summarize_merge_gene_loss",\\s*mode\\s*=\\s*"function"\\)')
  expect_match(src, "summarize_merge_gene_loss\\(")
  expect_match(src, "tryCatch\\(")
})

test_that(".sc_log_merge_gene_report reports union retention and gene loss via summarize_merge_gene_loss", {
  o1 <- .make_tiny_sc(c("G1", "G2", "G3", "G4"), c("c1", "c2"))
  o2 <- .make_tiny_sc(c("G2", "G3", "G5"), c("c3", "c4"))
  merged_union <- merge(o1, y = o2, add.cell.ids = c("S1", "S2"))
  merged_common <- subset(merged_union, features = c("G2", "G3"))

  gd <- new.env(parent = emptyenv())
  gd$i18n <- NULL

  shiny::testServer(mod_import_sc_server, args = list(global_data = gd), expr = {
    notified <- list()
    mod_env <- environment(.sc_log_merge_gene_report)
    orig_notify <- get(".notify_sc", envir = mod_env)
    on.exit(assign(".notify_sc", orig_notify, envir = mod_env), add = TRUE)
    assign(".notify_sc", function(msg, type = "message", duration = 5) {
      notified[[length(notified) + 1L]] <<- list(msg = msg, type = type, duration = duration)
    }, envir = mod_env)

    # Case 1: union merge (0 genes lost)
    res_union <- .sc_log_merge_gene_report(list(S1 = o1, S2 = o2), merged_union, log_per_sample = TRUE)
    expect_false(is.null(res_union$gene_loss))
    expect_identical(res_union$gene_loss$total_union, 5L)
    expect_identical(res_union$gene_loss$total_retained, 5L)
    expect_identical(res_union$gene_loss$total_lost_vs_union, 0L)
    expect_match(logs(), "5/5 retained vs union \\(0 lost\\)")
    expect_identical(length(notified), 0L)

    # Case 2: subset/intersection merge (3 genes lost vs union, per-sample breakdown + warning notification)
    res_loss <- .sc_log_merge_gene_report(list(S1 = o1, S2 = o2), merged_common)
    expect_false(is.null(res_loss$gene_loss))
    expect_identical(res_loss$gene_loss$total_union, 5L)
    expect_identical(res_loss$gene_loss$total_retained, 2L)
    expect_identical(res_loss$gene_loss$total_lost_vs_union, 3L)
    expect_match(logs(), "2/5 retained vs union \\(3 lost\\)")
    expect_match(logs(), "S1=-2")
    expect_match(logs(), "S2=-1")
    expect_identical(length(notified), 1L)
    expect_identical(notified[[1L]]$type, "warning")
    expect_match(notified[[1L]]$msg, "2/5 retained vs union \\(3 lost\\)")
  })
})

test_that(".sc_log_merge_gene_report degrades gracefully if summarize_merge_gene_loss errors", {
  o1 <- .make_tiny_sc(c("G1", "G2"), c("c1", "c2"))
  o2 <- .make_tiny_sc(c("G2", "G3"), c("c3", "c4"))
  merged_union <- merge(o1, y = o2, add.cell.ids = c("S1", "S2"))

  gd <- new.env(parent = emptyenv())
  gd$i18n <- NULL

  orig_fn <- get("summarize_merge_gene_loss", envir = globalenv())
  on.exit(assign("summarize_merge_gene_loss", orig_fn, envir = globalenv()), add = TRUE)
  assign("summarize_merge_gene_loss", function(...) stop("boom"), envir = globalenv())

  shiny::testServer(mod_import_sc_server, args = list(global_data = gd), expr = {
    res <- .sc_log_merge_gene_report(list(S1 = o1, S2 = o2), merged_union)
    expect_identical(res$union, 3L)
    expect_identical(res$common, 1L)
    expect_null(res$gene_loss)
  })
})
