# =============================================================================
# test-bulk-import-tpm-guard.R — W2.4 Bulk import TPM/FPKM & normalized matrix
# detection + de_allowed guard in modules/import/mod_import_bulk.R
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/rdata_io.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("tools/pkg_skeleton/R/provenance.R")
source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/batch_correction.R")
source_project_file("tools/pkg_skeleton/R/bulk_import_engine.R")
source_project_file("modules/import/mod_rdata_picker.R")
source_project_file("modules/import/mod_import_bulk.R")

suppressPackageStartupMessages(library(shiny))

.tpm_guard_tmp <- tempfile("tpm_guard_")
dir.create(.tpm_guard_tmp, showWarnings = FALSE, recursive = TRUE)

.write_matrix_csv <- function(mat, filename) {
  path <- file.path(.tpm_guard_tmp, filename)
  df <- data.frame(gene = rownames(mat), as.data.frame(mat, check.names = FALSE),
                   stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.csv(df, path, row.names = FALSE)
  path
}

# ── 1. Pure helper .bulk_detect_matrix_scale() ──────────────────────────────

test_that(".bulk_detect_matrix_scale detects raw integer counts (de_allowed = TRUE)", {
  expect_true(exists(".bulk_detect_matrix_scale", mode = "function"),
              info = ".bulk_detect_matrix_scale helper must be defined in mod_import_bulk.R")
  if (!exists(".bulk_detect_matrix_scale", mode = "function")) return(invisible(NULL))

  mat <- matrix(c(0L, 15L, 42L, 120L, 8L, 300L, 19L, 64L), nrow = 4, ncol = 2,
                dimnames = list(paste0("G", 1:4), c("S1", "S2")))
  res <- .bulk_detect_matrix_scale(mat)

  expect_identical(res$verdict, "raw")
  expect_true(isTRUE(res$de_allowed))
  expect_true(isTRUE(res$is_raw_counts))
  expect_false(isTRUE(res$has_fractional))
  expect_false(isTRUE(res$has_negative))
  expect_false(isTRUE(res$is_tpm_like))
  expect_identical(res$n_non_integer, 0L)
})

test_that(".bulk_detect_matrix_scale detects fractional TPM/FPKM matrices (de_allowed = FALSE + warning)", {
  expect_true(exists(".bulk_detect_matrix_scale", mode = "function"))
  if (!exists(".bulk_detect_matrix_scale", mode = "function")) return(invisible(NULL))

  # TPM matrix with column sums ~ 1e6 and values in (0, 1)
  c1 <- c(0, 0.25, 0.85, 14.6, 880000.3, 119984.0)
  c2 <- c(0, 0.40, 0.75, 22.1, 895000.25, 104976.5)
  mat_tpm <- cbind(S1 = c1, S2 = c2)
  rownames(mat_tpm) <- paste0("G", seq_len(nrow(mat_tpm)))

  res_tpm <- .bulk_detect_matrix_scale(mat_tpm)
  expect_identical(res_tpm$verdict, "tpm_fpkm_likely")
  expect_false(isTRUE(res_tpm$de_allowed))
  expect_false(isTRUE(res_tpm$is_raw_counts))
  expect_true(isTRUE(res_tpm$has_fractional))
  expect_true(isTRUE(res_tpm$is_tpm_like))
  expect_match(res_tpm$message, "TPM/FPKM", fixed = TRUE)

  # Fractional FPKM/normalized subset without 1e6 column sum
  mat_frac <- matrix(c(12.3, 45.6, 0.75, 88.1, 14.2, 39.8, 0.42, 91.7), nrow = 4, ncol = 2,
                     dimnames = list(paste0("G", 1:4), c("S1", "S2")))
  res_frac <- .bulk_detect_matrix_scale(mat_frac)
  expect_identical(res_frac$verdict, "transformed_generic")
  expect_false(isTRUE(res_frac$de_allowed))
  expect_true(isTRUE(res_frac$has_fractional))
  expect_match(res_frac$message, "counts bruts", fixed = TRUE)
})

test_that(".bulk_detect_matrix_scale detects negative log-normalized values (de_allowed = FALSE)", {
  expect_true(exists(".bulk_detect_matrix_scale", mode = "function"))
  if (!exists(".bulk_detect_matrix_scale", mode = "function")) return(invisible(NULL))

  mat_neg <- matrix(c(-1.5, 0.8, 2.3, -0.4, 1.2, -2.1, 0.5, 1.9), nrow = 4, ncol = 2,
                    dimnames = list(paste0("G", 1:4), c("S1", "S2")))
  res <- .bulk_detect_matrix_scale(mat_neg)

  expect_identical(res$verdict, "transformed_generic")
  expect_false(isTRUE(res$de_allowed))
  expect_false(isTRUE(res$is_raw_counts))
  expect_true(isTRUE(res$has_negative))
  expect_match(res$message, "n\u00e9gatives", fixed = TRUE)
})

test_that(".bulk_detect_matrix_scale handles edge cases (all-zero, small matrix, empty/NULL)", {
  expect_true(exists(".bulk_detect_matrix_scale", mode = "function"))
  if (!exists(".bulk_detect_matrix_scale", mode = "function")) return(invisible(NULL))

  # All-zero matrix: no NaN/crash, flagged all_zero and de_allowed = FALSE
  mat_zero <- matrix(0, nrow = 4, ncol = 3,
                     dimnames = list(paste0("G", 1:4), paste0("S", 1:3)))
  res_zero <- .bulk_detect_matrix_scale(mat_zero)
  expect_true(isTRUE(res_zero$all_zero))
  expect_false(isTRUE(res_zero$de_allowed))
  expect_false(isTRUE(res_zero$has_fractional))
  expect_false(isTRUE(res_zero$has_negative))

  # Small 2x2 raw integer count matrix
  mat_small <- matrix(c(10L, 25L, 12L, 30L), nrow = 2, ncol = 2,
                      dimnames = list(c("G1", "G2"), c("S1", "S2")))
  res_small <- .bulk_detect_matrix_scale(mat_small)
  expect_identical(res_small$verdict, "raw")
  expect_true(isTRUE(res_small$de_allowed))
  expect_false(isTRUE(res_small$all_zero))

  # Small 1x1 fractional matrix
  res_1x1_frac <- .bulk_detect_matrix_scale(matrix(3.14, nrow = 1, ncol = 1))
  expect_false(isTRUE(res_1x1_frac$de_allowed))
  expect_true(isTRUE(res_1x1_frac$has_fractional))

  # Invalid / empty inputs never throw
  expect_no_error(res_null <- .bulk_detect_matrix_scale(NULL))
  expect_identical(res_null$verdict, "unknown")
  expect_false(isTRUE(res_null$de_allowed))

  expect_no_error(res_empty <- .bulk_detect_matrix_scale(matrix(numeric(0), nrow = 0, ncol = 0)))
  expect_identical(res_empty$verdict, "unknown")
  expect_false(isTRUE(res_empty$de_allowed))
})

# ── 2. Server integration: merged_matrix, generic_omics, and per_sample ─────

test_that("mod_import_bulk_server sets de_allowed = TRUE on raw integer counts and FALSE + warning on TPM/negative", {
  # 2a. Raw integer counts -> de_allowed = TRUE
  mat_raw <- matrix(c(50L, 120L, 80L, 60L, 140L, 95L, 45L, 110L, 70L, 55L, 130L, 85L),
                    nrow = 3, ncol = 4,
                    dimnames = list(c("G1", "G2", "G3"), c("S1", "S2", "S3", "S4")))
  csv_raw <- .write_matrix_csv(mat_raw, "raw_counts.csv")

  gd_raw <- new.env(parent = emptyenv())
  gd_raw$i18n <- NULL
  gd_raw$drive_registry <- new.env(parent = emptyenv())
  logs_raw <- NULL

  shiny::testServer(mod_import_bulk_server, args = list(global_data = gd_raw), expr = {
    session$setInputs(
      bulk_import_mode = "merged_matrix",
      counts_format = "rows",
      counts_has_header = TRUE,
      counts_has_rownames = TRUE,
      bulk_data_type = "rna_counts",
      min_counts = 10,
      project_name = "RawProj",
      multi_label = "",
      counts_file = list(name = "raw_counts.csv", datapath = csv_raw,
                         size = file.info(csv_raw)$size, type = "text/csv"),
      btn_load = 1L
    )
    session$flushReact()
    logs_raw <<- paste(logs(), collapse = "\n")
  })

  expect_false(is.null(gd_raw$bulk_obj))
  expect_true(isTRUE(gd_raw$bulk_obj$de_allowed))
  expect_false(grepl("TPM/FPKM|valeurs n\u00e9gatives|Matrice continue", logs_raw))

  # 2b. Fractional TPM matrix -> de_allowed = FALSE + user-facing warning in logs
  mat_tpm <- cbind(
    S1 = c(0, 0.35, 15.2, 880000.45, 119984.0),
    S2 = c(0, 0.50, 18.8, 890000.20, 109980.5),
    S3 = c(0, 0.25, 12.4, 870000.35, 129987.0),
    S4 = c(0, 0.60, 20.1, 895000.30, 104979.0)
  )
  rownames(mat_tpm) <- paste0("G", 1:5)
  csv_tpm <- .write_matrix_csv(mat_tpm, "tpm_matrix.csv")

  gd_tpm <- new.env(parent = emptyenv())
  gd_tpm$i18n <- NULL
  gd_tpm$drive_registry <- new.env(parent = emptyenv())
  logs_tpm <- NULL

  shiny::testServer(mod_import_bulk_server, args = list(global_data = gd_tpm), expr = {
    session$setInputs(
      bulk_import_mode = "merged_matrix",
      counts_format = "rows",
      counts_has_header = TRUE,
      counts_has_rownames = TRUE,
      bulk_data_type = "rna_counts",
      min_counts = 10,
      project_name = "TPMProj",
      multi_label = "",
      counts_file = list(name = "tpm_matrix.csv", datapath = csv_tpm,
                         size = file.info(csv_tpm)$size, type = "text/csv"),
      btn_load = 1L
    )
    session$flushReact()
    logs_tpm <<- paste(logs(), collapse = "\n")
  })

  expect_false(is.null(gd_tpm$bulk_obj))
  expect_false(isTRUE(gd_tpm$bulk_obj$de_allowed))
  expect_match(logs_tpm, "TPM/FPKM", fixed = TRUE)

  # 2c. Negative log-normalized matrix -> de_allowed = FALSE + warning (and rows preserved)
  mat_neg <- matrix(c(-1.8, 0.4, 1.2, -0.7, 0.9, -1.5, 2.1, -0.3,
                      -0.9, 1.4, -1.1, 0.6, 1.7, -0.8, -0.2, 1.0),
                    nrow = 4, ncol = 4,
                    dimnames = list(paste0("G", 1:4), paste0("S", 1:4)))
  csv_neg <- .write_matrix_csv(mat_neg, "neg_log.csv")

  gd_neg <- new.env(parent = emptyenv())
  gd_neg$i18n <- NULL
  gd_neg$drive_registry <- new.env(parent = emptyenv())
  logs_neg <- NULL

  shiny::testServer(mod_import_bulk_server, args = list(global_data = gd_neg), expr = {
    session$setInputs(
      bulk_import_mode = "merged_matrix",
      counts_format = "rows",
      counts_has_header = TRUE,
      counts_has_rownames = TRUE,
      bulk_data_type = "rna_counts",
      min_counts = 10,
      project_name = "NegProj",
      multi_label = "",
      counts_file = list(name = "neg_log.csv", datapath = csv_neg,
                         size = file.info(csv_neg)$size, type = "text/csv"),
      btn_load = 1L
    )
    session$flushReact()
    logs_neg <<- paste(logs(), collapse = "\n")
  })

  expect_false(is.null(gd_neg$bulk_obj))
  expect_false(isTRUE(gd_neg$bulk_obj$de_allowed))
  expect_equal(nrow(gd_neg$bulk_obj$counts), 4L)
  expect_match(logs_neg, "n\u00e9gatives", fixed = TRUE)
})

test_that("mod_import_bulk_server sets de_allowed = FALSE in per_sample mode when files contain fractional TPMs", {
  ps_dir <- file.path(.tpm_guard_tmp, "ps_tpm")
  dir.create(ps_dir, showWarnings = FALSE, recursive = TRUE)
  samples <- c("Ctrl_1", "Ctrl_2", "Treat_1", "Treat_2")
  paths <- character(length(samples))
  for (i in seq_along(samples)) {
    p <- file.path(ps_dir, paste0(samples[i], ".tsv"))
    writeLines(c(
      "gene_id\tcount",
      sprintf("G1\t%.2f", 0.35 + i * 0.1),
      sprintf("G2\t%.2f", 880000.25 - i * 100),
      sprintf("G3\t%.2f", 119999.40 + i * 100)
    ), p)
    paths[i] <- p
  }

  gd_ps <- new.env(parent = emptyenv())
  gd_ps$i18n <- NULL
  gd_ps$drive_registry <- new.env(parent = emptyenv())
  logs_ps <- NULL

  shiny::testServer(mod_import_bulk_server, args = list(global_data = gd_ps), expr = {
    session$setInputs(
      bulk_import_mode = "per_sample",
      project_name = "PS_TPM",
      ps_dup_threshold = 0.05,
      ps_fill_zero = TRUE,
      btn_load = 0L
    )
    session$setInputs(ps_files = data.frame(
      name = basename(paths), size = file.size(paths), type = "text/tsv",
      datapath = paths, stringsAsFactors = FALSE
    ))
    session$flushReact()
    session$setInputs(btn_load = 1L)
    session$flushReact()
    logs_ps <<- paste(logs(), collapse = "\n")
  })

  expect_false(is.null(gd_ps$bulk_obj))
  expect_false(isTRUE(gd_ps$bulk_obj$de_allowed))
  expect_match(logs_ps, "non entiers|TPM/FPKM")
})
