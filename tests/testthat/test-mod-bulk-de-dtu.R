# ==============================================================================
# test-mod-bulk-de-dtu.R — DTU (edgeR::diffSpliceDGE) UI & server wiring
# ==============================================================================

skip_if_not_installed("shiny")
skip_if_not_installed("bslib")
skip_if_not_installed("DT")
skip_if_not_installed("edgeR")
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(bslib))
suppressPackageStartupMessages(library(DT))

.dtu_test_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("mod_bulk_de_dtu_ui exposes DTU controls and tables", {
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)
  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }

  source(.dtu_test_root("modules/bulk_de/mod_bulk_de_ui.R"), local = globalenv())
  ui_html <- as.character(mod_bulk_de_dtu_ui("de_dtu"))
  expect_true(grepl("de_dtu-dtu_tx2gene_file", ui_html, fixed = TRUE))
  expect_true(grepl("de_dtu-dtu_test_method", ui_html, fixed = TRUE))
  expect_true(grepl("de_dtu-run_dtu", ui_html, fixed = TRUE))
  expect_true(grepl("de_dtu-dtu_gene_table", ui_html, fixed = TRUE))
  expect_true(grepl("de_dtu-dtu_tx_table", ui_html, fixed = TRUE))
})

test_that(".de_dtu_server runs DTU on transcript counts + tx2gene", {
  source(.dtu_test_root("tools/pkg_skeleton/R/bulk_dtu.R"), local = globalenv())
  source(.dtu_test_root("modules/bulk_de/mod_bulk_de_summary.R"), local = globalenv())

  if (!exists("ts_datatable", envir = globalenv())) {
    assign("ts_datatable", function(df, ...) DT::datatable(df), envir = globalenv())
  }
  if (!exists(".strip_i18n_html", envir = globalenv())) {
    assign(".strip_i18n_html", function(x) x, envir = globalenv())
  }

  set.seed(7)
  n_tx <- 12L
  tx_ids <- paste0("TX", seq_len(n_tx))
  gene_ids <- rep(paste0("G", seq_len(4L)), each = 3L)
  samples <- paste0("S", seq_len(6L))
  counts <- matrix(rpois(n_tx * 6L, lambda = 80), nrow = n_tx, ncol = 6L,
                   dimnames = list(tx_ids, samples))
  # Inject strong isoform switch in G1 between ctrl (1:3) and case (4:6)
  counts[1, 1:3] <- c(250, 260, 240)
  counts[1, 4:6] <- c(20, 25, 30)
  counts[2, 1:3] <- c(25, 30, 20)
  counts[2, 4:6] <- c(260, 245, 255)

  tx2gene <- data.frame(transcript_id = tx_ids, gene_id = gene_ids, stringsAsFactors = FALSE)
  meta <- data.frame(
    condition = rep(c("ctrl", "case"), each = 3L),
    row.names = samples,
    stringsAsFactors = FALSE
  )

  shared_rv <- shiny::reactiveValues(
    filtered_counts = counts,
    active_condition_col = "condition",
    contrasts = list(),
    lfc_thresh = 1,
    padj_thresh = 0.05
  )
  global_data <- shiny::reactiveValues(
    language = "fr",
    i18n = list(t = function(k) k),
    bulk_obj = list(counts = counts, metadata = meta, tx2gene = tx2gene)
  )

  dtu_wrapper_server <- function(id, global_data, shared_rv) {
    shiny::moduleServer(id, function(input, output, session) {
      .de_dtu_server(input, output, session, session$ns, global_data, shared_rv)
    })
  }

  shiny::testServer(
    dtu_wrapper_server,
    args = list(global_data = global_data, shared_rv = shared_rv),
    {
      session$setInputs(condition_col = "condition", dtu_test_method = "Simes")
      session$setInputs(run_dtu = 1L)
      st <- output$dtu_status
      expect_true(grepl("DTU (Simes)", st, fixed = TRUE))
      expect_true(grepl("4 gènes testés", st, fixed = TRUE))
      expect_false(is.null(output$dtu_gene_table))
      expect_false(is.null(output$dtu_tx_table))
    }
  )
})
