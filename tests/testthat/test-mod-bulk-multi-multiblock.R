# ==============================================================================
# test-mod-bulk-multi-multiblock.R — BK-06 DIABLO/NEMO UI & server in mod_bulk_multi
# ==============================================================================

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
suppressPackageStartupMessages(library(DT))

.mb_test_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("mod_bulk_multi_ui exposes tab_multiblock and BK-06 controls", {
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)
  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }

  source(.mb_test_root("modules/bulk/mod_bulk_multi.R"), local = globalenv())
  ui_html <- as.character(mod_bulk_multi_ui("bm_mb"))
  expect_true(grepl("tab_multiblock", ui_html, fixed = TRUE))
  expect_true(grepl("bm_mb-multiblock_mode", ui_html, fixed = TRUE))
  expect_true(grepl("bm_mb-run_multiblock", ui_html, fixed = TRUE))
  expect_true(grepl("bm_mb-multiblock_table", ui_html, fixed = TRUE))
})

test_that("mod_bulk_multi_server runs NEMO and DIABLO via tab_multiblock", {
  source(.mb_test_root("tools/pkg_skeleton/R/bulk_multi_compare.R"), local = globalenv())
  source(.mb_test_root("tools/pkg_skeleton/R/bulk_multiblock.R"), local = globalenv())
  source(.mb_test_root("modules/bulk/mod_bulk_multi.R"), local = globalenv())

  if (!exists("ts_datatable", envir = globalenv())) {
    assign("ts_datatable", function(df, ...) DT::datatable(df), envir = globalenv())
  }
  if (!exists(".strip_i18n_html", envir = globalenv())) {
    assign(".strip_i18n_html", function(x) x, envir = globalenv())
  }

  set.seed(123)
  samples <- paste0("S", seq_len(8))
  m1 <- matrix(rnorm(80), nrow = 10, ncol = 8,
               dimnames = list(paste0("G", 1:10), samples))
  m1[1:3, 1:4] <- m1[1:3, 1:4] + 3
  m2 <- matrix(rnorm(64), nrow = 8, ncol = 8,
               dimnames = list(paste0("P", 1:8), samples))
  m2[1:3, 1:4] <- m2[1:3, 1:4] + 3

  meta <- data.frame(
    condition = rep(c("ctrl", "case"), each = 4),
    row.names = samples,
    stringsAsFactors = FALSE
  )
  fake_de <- data.frame(
    log2FoldChange = c(2, -2, 0.1),
    pvalue = c(1e-4, 2e-4, 0.8),
    padj = c(1e-3, 2e-3, 0.9),
    row.names = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )

  ds_list <- list(
    RNA = list(
      obj = list(counts = m1, metadata = meta),
      pipeline = list(
        vst_mat = m1,
        active_condition_col = "condition",
        contrasts = list(case_vs_ctrl = fake_de)
      )
    ),
    Proteomics = list(
      obj = list(counts = m2, metadata = meta),
      pipeline = list(
        vst_mat = m2,
        active_condition_col = "condition",
        contrasts = list(case_vs_ctrl = fake_de)
      )
    )
  )

  global_data <- shiny::reactiveValues(
    language = "fr",
    i18n = list(t = function(k) k),
    bulk_datasets = ds_list,
    bulk_multi_comparison = NULL
  )

  shiny::testServer(
    mod_bulk_multi_server,
    args = list(global_data = global_data),
    {
      session$setInputs(
        datasets = c("RNA", "Proteomics"),
        multiblock_mode = "nemo",
        multiblock_y_col = "condition"
      )
      session$setInputs(run_multiblock = 1L)
      sum_nemo <- output$multiblock_summary
      expect_true(grepl("NEMO", sum_nemo, fixed = TRUE))
      expect_true(grepl("clusters", sum_nemo, fixed = TRUE))
      expect_false(is.null(output$multiblock_table))

      session$setInputs(multiblock_mode = "diablo", run_multiblock = 2L)
      sum_diablo <- output$multiblock_summary
      expect_true(grepl("DIABLO", sum_diablo, fixed = TRUE))
      expect_true(grepl("Exactitude CV", sum_diablo, fixed = TRUE))
      expect_false(is.null(output$multiblock_table))
    }
  )
})
