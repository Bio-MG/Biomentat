# =============================================================================
# test-mod-sc-pipeline-aberrations.R — Multi-replicate aberration wiring (W2.4 §2.1)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/sc_helpers.R")
source_project_file("modules/sc/mod_sc_pipeline.R")

.make_aberration_test_meta <- function(rep_sizes = c(S1 = 120L, S2 = 115L, S3 = 125L, S4 = 118L),
                                       nFeature_means = c(S1 = 1500, S2 = 1500, S3 = 1500, S4 = 1500),
                                       nCount_means = c(S1 = 4000, S2 = 4000, S3 = 4000, S4 = 4000),
                                       mt_means = c(S1 = 5, S2 = 5, S3 = 5, S4 = 5)) {
  dfs <- lapply(names(rep_sizes), function(nm) {
    n <- rep_sizes[[nm]]
    data.frame(
      orig.ident = rep(nm, n),
      nFeature_RNA = stats::rnorm(n, mean = nFeature_means[[nm]], sd = 40),
      nCount_RNA = stats::rnorm(n, mean = nCount_means[[nm]], sd = 80),
      percent.mt = pmax(0, stats::rnorm(n, mean = mt_means[[nm]], sd = 0.4)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, dfs)
}

.make_aberration_mock_obj <- function(meta_df) {
  if (requireNamespace("SeuratObject", quietly = TRUE)) {
    n_cells <- nrow(meta_df)
    counts <- Matrix::Matrix(
      stats::rpois(20L * n_cells, lambda = 2),
      nrow = 20L, ncol = n_cells, sparse = TRUE
    )
    rownames(counts) <- paste0("Gene", seq_len(20L))
    colnames(counts) <- paste0("Cell", seq_len(n_cells))
    rownames(meta_df) <- colnames(counts)
    SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta_df)
  } else {
    list(meta.data = meta_df)
  }
}

test_that("QC summary emits concise warning when multi-replicate QC detects an aberrant replicate", {
  set.seed(123)
  meta <- .make_aberration_test_meta(
    rep_sizes = c(S1 = 120L, S2 = 115L, S3 = 25L, S4 = 118L)
  )
  obj <- .make_aberration_mock_obj(meta)
  snap <- data.frame(
    Echantillon = c("S1", "S2", "S3", "S4"),
    Cellules_avant = c(140L, 135L, 130L, 135L),
    Cellules_apres = c(120L, 115L, 25L, 118L),
    stringsAsFactors = FALSE
  )
  params <- list(qc_min_gene = 200, qc_max_gene = 6000, qc_mt = 20)

  out <- .sc_pipeline_qc_summary(snap, params = params, obj = obj)
  expect_identical(out$status, "applied")
  expect_identical(out$flagged_replicates, "S3")
  expect_true(is.character(out$aberrations_text) && length(out$aberrations_text) == 1L && nzchar(out$aberrations_text))
  expect_match(out$aberrations_text, "S3", fixed = TRUE)
  expect_match(out$aberrations_text, "low_cell_count", fixed = TRUE)
  expect_match(out$text, "S3", fixed = TRUE)
})

test_that("QC summary emits no aberration warning on balanced multi-replicate runs", {
  set.seed(42)
  meta <- .make_aberration_test_meta(
    rep_sizes = c(S1 = 200L, S2 = 210L, S3 = 190L, S4 = 205L),
    nFeature_means = c(S1 = 1500, S2 = 1550, S3 = 1480, S4 = 1520),
    nCount_means = c(S1 = 4000, S2 = 4100, S3 = 3950, S4 = 4050),
    mt_means = c(S1 = 5, S2 = 5.2, S3 = 4.8, S4 = 5.1)
  )
  obj <- .make_aberration_mock_obj(meta)
  snap <- data.frame(
    Echantillon = c("S1", "S2", "S3", "S4"),
    Cellules_avant = c(210L, 220L, 200L, 215L),
    Cellules_apres = c(200L, 210L, 190L, 205L),
    stringsAsFactors = FALSE
  )

  out <- .sc_pipeline_qc_summary(snap, obj = obj)
  expect_identical(out$status, "applied")
  expect_null(out$aberrations_text)
  expect_identical(out$flagged_replicates, character(0))
})

test_that("QC summary skips aberration check when fewer than 2 replicates exist", {
  set.seed(123)
  meta <- .make_aberration_test_meta(
    rep_sizes = c(S1 = 30L),
    nFeature_means = c(S1 = 1500),
    nCount_means = c(S1 = 4000),
    mt_means = c(S1 = 5)
  )
  obj <- .make_aberration_mock_obj(meta)
  snap <- data.frame(
    Echantillon = "S1",
    Cellules_avant = 50L,
    Cellules_apres = 30L,
    stringsAsFactors = FALSE
  )

  out <- .sc_pipeline_qc_summary(snap, obj = obj)
  expect_identical(out$status, "applied")
  expect_null(out$aberrations_text)
  expect_identical(out$flagged_replicates, character(0))
})

test_that("QC summary handles objects without meta.data or orig.ident gracefully", {
  snap <- data.frame(
    Echantillon = c("S1", "S2"),
    Cellules_avant = c(100L, 100L),
    Cellules_apres = c(90L, 90L),
    stringsAsFactors = FALSE
  )
  out_mat <- .sc_pipeline_qc_summary(snap, obj = matrix(0, nrow = 50L, ncol = 180L))
  expect_null(out_mat$aberrations_text)
  expect_identical(out_mat$flagged_replicates, character(0))

  out_no_col <- .sc_pipeline_qc_summary(
    snap,
    obj = list(meta.data = data.frame(nFeature_RNA = 1:10))
  )
  expect_null(out_no_col$aberrations_text)
  expect_identical(out_no_col$flagged_replicates, character(0))
})

test_that("mod_sc_pipeline.R guards detect_replicate_aberrations and wires UI + notification", {
  src <- paste(
    readLines(file.path(ts_project_root(), "modules/sc/mod_sc_pipeline.R"),
              warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )
  expect_match(src, 'exists("detect_replicate_aberrations", mode = "function")', fixed = TRUE)
  expect_match(src, "detect_replicate_aberrations(", fixed = TRUE)
  expect_match(src, "tryCatch(", fixed = TRUE)
  expect_match(src, "aberrations_text", fixed = TRUE)
  expect_match(src, 'type = "warning"', fixed = TRUE)
})
