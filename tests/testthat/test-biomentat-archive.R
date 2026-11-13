# =============================================================================
# test-biomentat-archive.R — P3.3: Portable .biomentat archive + QC fields
# =============================================================================

source_project_file("tools/pkg_skeleton/R/biomentat_archive.R")

test_that("create_biomentat_archive validates arguments and extension", {
  err1 <- tryCatch(create_biomentat_archive(NULL, tempfile(fileext = ".biomentat")), error = function(e) e)
  expect_s3_class(err1, "biomentat_archive_error")

  err2 <- tryCatch(create_biomentat_archive(list(a = 1), tempfile(fileext = ".zip")), error = function(e) e)
  expect_s3_class(err2, "biomentat_archive_error")
})

test_that("create_biomentat_archive and read_biomentat_archive round-trip with QC fields & SHA-256 verification", {
  arch_path <- tempfile(fileext = ".biomentat")
  on.exit(unlink(arch_path), add = TRUE)

  payload <- list(counts = matrix(1:12, nrow = 3L), notes = "GSE_test")
  qc <- list(n_genes = 3L, n_samples = 4L, min_lib_size = 10L, mito_pct_median = 2.5)
  meta <- list(modality = "bulk", organism = "hsapiens")

  out <- create_biomentat_archive(payload, arch_path, qc_metrics = qc, metadata = meta)
  expect_true(file.exists(arch_path))
  expect_identical(out$schema, "biomentat-archive/v1")
  expect_identical(out$exchange_schema, "biomentat-exchange/v1")
  expect_true(nzchar(out$payload_sha256))

  rd <- read_biomentat_archive(arch_path, verify_hash = TRUE)
  expect_identical(rd$manifest$qc_metrics$n_genes, 3L)
  expect_identical(rd$manifest$qc_metrics$mito_pct_median, 2.5)
  expect_identical(rd$manifest$metadata$modality, "bulk")
  expect_identical(rd$payload, payload)
})

test_that("read_biomentat_archive rejects missing archive", {
  err <- tryCatch(read_biomentat_archive(tempfile(fileext = ".biomentat")), error = function(e) e)
  expect_s3_class(err, "biomentat_archive_error")
})
