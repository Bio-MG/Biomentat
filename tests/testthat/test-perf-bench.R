# =============================================================================
# test-perf-bench.R — P1.3b : harnais de benchmark offline + sortie JSON
# =============================================================================

test_that("P1.3b : benchmarks/run_perf_bench.R produit un JSON valide avec metriques finies", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")

  bench_script <- file.path(ts_project_root(), "benchmarks", "run_perf_bench.R")
  expect_true(file.exists(bench_script), info = "benchmarks/run_perf_bench.R absent")
  if (!file.exists(bench_script)) return(invisible(NULL))

  out_json <- tempfile(fileext = ".json")
  on.exit(unlink(out_json), add = TRUE)

  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  status <- system2(rscript, args = c(shQuote(bench_script), shQuote(paste0("--out=", out_json))),
                    stdout = TRUE, stderr = TRUE)
  exit_code <- attr(status, "status") %||% 0L
  expect_identical(as.integer(exit_code), 0L,
                   info = paste("Echec run_perf_bench.R:", paste(status, collapse = "\n")))
  expect_true(file.exists(out_json))

  parsed <- jsonlite::fromJSON(out_json, simplifyVector = FALSE)
  expect_true(is.list(parsed))
  expect_true(!is.null(parsed$timestamp) && nzchar(parsed$timestamp))
  expect_true(!is.null(parsed$r_version) && nzchar(parsed$r_version))
  expect_true(is.list(parsed$tasks) && length(parsed$tasks) >= 1L)

  first_task <- parsed$tasks[[1]]
  expect_true(is.character(first_task$name) && nzchar(first_task$name))
  expect_true(is.numeric(first_task$elapsed_s) && first_task$elapsed_s >= 0)
  expect_true(is.numeric(first_task$rss_before_mb) && is.finite(first_task$rss_before_mb))
  expect_true(is.numeric(first_task$rss_after_mb) && is.finite(first_task$rss_after_mb))
})
