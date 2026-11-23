# ==============================================================================
# test-bulk-tx2gene.R — In-app tx2gene generator (TxDb + GTF) and DTU UI
# ==============================================================================

.tx2g_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("build_tx2gene maps transcripts from GTF and reports unmapped IDs", {
  source(.tx2g_root("tools/pkg_skeleton/R/bulk_dtu.R"), local = globalenv())

  gtf <- tempfile(fileext = ".gtf")
  on.exit(unlink(gtf), add = TRUE)
  writeLines(c(
    'chr1\thavana\texon\t100\t200\t.\t+\t.\tgene_id "GENE1"; transcript_id "ENST0001";',
    'chr1\thavana\texon\t300\t400\t.\t+\t.\tgene_id "GENE1"; transcript_id "ENST0002";',
    'chr2\thavana\texon\t500\t600\t.\t-\t.\tgene_id "GENE2"; transcript_id "ENST0003";'
  ), gtf)

  cnt <- matrix(10, nrow = 4, ncol = 2,
                dimnames = list(c("ENST0001", "ENST0002", "ENST0003", "ENST_UNMAPPED"), c("S1", "S2")))
  res <- build_tx2gene(cnt, source = "gtf", gtf_path = gtf)

  expect_equal(res$n_input, 4L)
  expect_equal(res$n_mapped, 3L)
  expect_equal(res$n_unmapped, 1L)
  expect_equal(res$unmapped_ids, "ENST_UNMAPPED")
  expect_equal(res$table$gene_id[match("ENST0002", res$table$transcript_id)], "GENE1")
  expect_equal(res$source, "gtf")
})

test_that("build_tx2gene throws classed bulk_dtu_error on missing GTF or unknown TxDb", {
  source(.tx2g_root("tools/pkg_skeleton/R/bulk_dtu.R"), local = globalenv())

  expect_error(
    build_tx2gene(c("TX1", "TX2"), source = "gtf", gtf_path = tempfile()),
    class = "bulk_dtu_error"
  )
  expect_error(
    build_tx2gene(c("TX1", "TX2"), source = "txdb", txdb_pkg = "TxDb.NonExistent.Pkg"),
    class = "bulk_dtu_error"
  )
  expect_error(
    build_tx2gene(character(0L), source = "txdb"),
    class = "bulk_dtu_error"
  )
})

test_that("mod_bulk_de_dtu_ui exposes btn_build_tx2gene and dtu_tx2gene_source", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  suppressPackageStartupMessages({
    library(shiny)
    library(bslib)
    library(DT)
  })
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)
  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }

  source(.tx2g_root("modules/bulk_de/mod_bulk_de_ui.R"), local = globalenv())
  ui_html <- as.character(mod_bulk_de_dtu_ui("dtu_test"))
  expect_true(grepl("dtu_test-btn_build_tx2gene", ui_html, fixed = TRUE))
  expect_true(grepl("dtu_test-dtu_tx2gene_source", ui_html, fixed = TRUE))
  expect_true(grepl("dtu_test-dtu_tx2gene_preview", ui_html, fixed = TRUE))
})
