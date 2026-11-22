# =============================================================================
# test-mod-bulk-multi-meta.R — BK-07 UI & server wiring in mod_bulk_multi.R
# =============================================================================

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(DT))

source_project_file("tools/pkg_skeleton/R/bulk_helpers.R")
source_project_file("tools/pkg_skeleton/R/bulk_multi.R")
source_project_file("tools/pkg_skeleton/R/bulk_multi_compare.R")
source_project_file("tools/pkg_skeleton/R/bulk_meta_analysis.R")
source_project_file("modules/bulk/mod_bulk_multi.R")

test_that("mod_bulk_multi_ui exposes BK-07 meta_method and tab_meta outputs", {
  local_i18n <- list(t = function(x) x)
  env <- new.env(parent = globalenv())
  env$i18n <- local_i18n
  env$layout_sidebar <- bslib::layout_sidebar
  env$sidebar <- bslib::sidebar
  env$navset_card_underline <- bslib::navset_card_underline
  env$nav_panel <- bslib::nav_panel
  env$DTOutput <- DT::DTOutput
  sys.source(file.path(ts_project_root(), "modules/bulk/mod_bulk_multi.R"), envir = env)

  html <- as.character(env$mod_bulk_multi_ui("multi"))
  expect_match(html, "multi-meta_method", fixed = TRUE)
  expect_match(html, "tab_meta", fixed = TRUE)
  expect_match(html, "multi-meta_summary", fixed = TRUE)
  expect_match(html, "multi-meta_table", fixed = TRUE)
})

test_that("mod_bulk_multi_server computes BK-07 horizontal meta-analysis across saved studies", {
  mk_de <- function(pvals, lfcs) {
    df <- data.frame(
      gene = paste0("G", seq_along(pvals)),
      log2FoldChange = lfcs,
      pvalue = pvals,
      padj = stats::p.adjust(pvals, method = "BH"),
      stringsAsFactors = FALSE
    )
    rownames(df) <- df$gene
    df
  }

  mk_entry <- function(label, pvals, lfcs) {
    df <- mk_de(pvals, lfcs)
    list(
      label = label, producer = "pipeline_save",
      registered_at = Sys.time(), updated_at = Sys.time(),
      obj = list(
        counts = matrix(1:12, nrow = 3L, ncol = 4L,
                        dimnames = list(df$gene, paste0("S", 1:4))),
        metadata = data.frame(row.names = paste0("S", 1:4),
                              condition = c("A", "A", "B", "B")),
        project = label, type = "bulk", import_mode = "merged_matrix"
      ),
      pipeline = list(lfc_thresh = 1, padj_thresh = 0.05,
                      contrasts = list(KO_vs_WT = df))
    )
  }

  ds <- list(
    StudyA = mk_entry("StudyA", c(1e-5, 0.4, 0.8), c(2.1, 0.3, -0.1)),
    StudyB = mk_entry("StudyB", c(2e-4, 0.5, 0.7), c(1.8, 0.2, 0.1))
  )
  gd <- shiny::reactiveValues(bulk_datasets = ds, i18n = NULL, language = "fr")

  shiny::testServer(mod_bulk_multi_server, args = list(id = "multi", global_data = gd), {
    session$setInputs(
      datasets = c("StudyA", "StudyB"),
      contrast = "KO_vs_WT",
      lfc_thresh = 1,
      padj_thresh = 0.05,
      meta_method = "stouffer"
    )
    session$setInputs(run = 1L)
    session$flushReact()

    m <- meta_result()
    expect_identical(m$method, "stouffer")
    expect_identical(m$n_studies, 2L)
    expect_identical(m$n_genes, 3L)
    expect_identical(m$table$gene[1L], "G1")
    expect_true(m$table$FDR[1L] < 0.01)

    summary_txt <- output$meta_summary
    expect_match(summary_txt, "stouffer", fixed = TRUE)
    expect_match(summary_txt, "2", fixed = TRUE)

    # Switching to Fisher updates reactively
    session$setInputs(meta_method = "fisher")
    session$flushReact()
    expect_identical(meta_result()$method, "fisher")
  })
})
