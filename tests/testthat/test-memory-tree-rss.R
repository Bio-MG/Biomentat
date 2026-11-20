# =============================================================================
# test-memory-tree-rss.R — W1.4 process-tree RSS & W1.2 .rds size estimator
# =============================================================================
source_project_file("tools/pkg_skeleton/R/memory.R")

test_that("ts_process_rss_mb returns a finite scalar MB or NA_real_", {
  v <- ts_process_rss_mb()
  expect_true(is.numeric(v) && length(v) == 1L)
  expect_true(is.na(v) || (is.finite(v) && v > 0))
})

test_that("ts_process_tree_rss_mb works for self-only and include_children = TRUE", {
  skip_if_not(requireNamespace("ps", quietly = TRUE))

  self_only <- ts_process_tree_rss_mb(include_children = FALSE)
  expect_named(self_only, c("total_mb", "self_mb", "children_mb", "n_children"))
  expect_true(is.finite(self_only$self_mb) && self_only$self_mb > 0)
  expect_equal(self_only$children_mb, 0)
  expect_equal(self_only$n_children, 0L)
  expect_equal(self_only$total_mb, self_only$self_mb)

  with_kids <- ts_process_tree_rss_mb(include_children = TRUE)
  expect_named(with_kids, c("total_mb", "self_mb", "children_mb", "n_children"))
  expect_true(is.finite(with_kids$total_mb) && is.finite(with_kids$self_mb))
  expect_gte(with_kids$total_mb, with_kids$self_mb)
  expect_gte(with_kids$children_mb, 0)
  expect_gte(with_kids$n_children, 0L)
})

test_that("ts_process_tree_rss_mb returns NA fallback on invalid PID without throwing", {
  bad_pids <- list(-1, 0, NA_real_, Inf, "not-a-pid", 999999999L, integer(0))
  for (p in bad_pids) {
    res <- ts_process_tree_rss_mb(pid = p)
    expect_named(res, c("total_mb", "self_mb", "children_mb", "n_children"))
    expect_true(is.na(res$total_mb))
    expect_true(is.na(res$self_mb))
    expect_true(is.na(res$children_mb))
    expect_true(is.na(res$n_children))
  }
})

test_that("ts_system_ram_mb returns total_mb and available_mb", {
  ram <- ts_system_ram_mb()
  expect_named(ram, c("total_mb", "available_mb"))
  expect_true(is.numeric(ram$total_mb) && length(ram$total_mb) == 1L)
  expect_true(is.numeric(ram$available_mb) && length(ram$available_mb) == 1L)
  if (requireNamespace("ps", quietly = TRUE)) {
    expect_true(is.finite(ram$total_mb) && ram$total_mb > 0)
    expect_true(is.finite(ram$available_mb) && ram$available_mb >= 0)
    expect_lte(ram$available_mb, ram$total_mb)
  }
})

test_that("ts_ram_budget_check classifies none, warn, and block purely", {
  mb <- 1024^2
  res_na <- ts_ram_budget_check(1000 * mb, factor = 3,
                                available_mb = NA_real_, total_mb = NA_real_)
  expect_equal(res_na$level, "none")
  expect_equal(res_na$message, "")

  res_ok <- ts_ram_budget_check(100 * mb, factor = 3,
                                available_mb = 2000, total_mb = 8000)
  expect_equal(res_ok$level, "none")
  expect_equal(res_ok$projected_mb, 300)

  res_warn <- ts_ram_budget_check(1000 * mb, factor = 3,
                                  available_mb = 2000, total_mb = 8000)
  expect_equal(res_warn$level, "warn")
  expect_true(nzchar(res_warn$message))

  res_block <- ts_ram_budget_check(3000 * mb, factor = 3,
                                   available_mb = 2000, total_mb = 8000)
  expect_equal(res_block$level, "block")
  expect_true(nzchar(res_block$message))
})

test_that("ts_estimate_rds_size_mb calibrates compressed .rds size vs ~5x overestimate", {
  raw_bytes <- 100 * 1024^2
  est <- ts_estimate_rds_size_mb(raw_bytes)
  expect_equal(est, 22)
  expect_lt(est, 100 / 4)

  # object.size() class compatibility
  obj_sz <- object.size(numeric(100000))
  est_obj <- ts_estimate_rds_size_mb(obj_sz, compression_ratio = 0.22)
  expect_true(is.finite(est_obj) && est_obj > 0)

  # Zero and custom compression ratio
  expect_equal(ts_estimate_rds_size_mb(0), 0)
  expect_equal(ts_estimate_rds_size_mb(100 * 1024^2, compression_ratio = 0.5), 50)
})

test_that("ts_estimate_rds_size_mb handles NA and invalid inputs without throwing", {
  bad_inputs <- list(NA_real_, NaN, Inf, -Inf, -1, "100", NULL, c(1, 2))
  for (x in bad_inputs) {
    expect_identical(ts_estimate_rds_size_mb(x), NA_real_)
  }
  expect_identical(ts_estimate_rds_size_mb(1024^2, compression_ratio = -0.1), NA_real_)
  expect_identical(ts_estimate_rds_size_mb(1024^2, compression_ratio = NA_real_), NA_real_)
  expect_identical(ts_estimate_rds_size_mb(1024^2, compression_ratio = Inf), NA_real_)
})

test_that("ts_estimate_rds_size_mb matches realistic saveRDS compressed size on synthetic object", {
  set.seed(42)
  syn <- list(
    expr = matrix(round(rnorm(4e5), 3), nrow = 800, ncol = 500),
    meta = data.frame(cell = paste0("c", seq_len(500)), grp = rep(letters[1:5], 100))
  )
  raw_bytes <- as.numeric(object.size(syn))
  raw_mb <- raw_bytes / 1024^2

  tmp_gz <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_gz), add = TRUE)
  saveRDS(syn, tmp_gz, compress = TRUE)
  actual_gz_mb <- file.info(tmp_gz)[["size"]] / 1024^2

  est_mb <- ts_estimate_rds_size_mb(object.size(syn))
  expect_true(is.finite(est_mb) && est_mb > 0)
  # Must avoid the ~5x overestimate of raw object.size()
  expect_lt(est_mb, raw_mb * 0.30)
  # Realistic relative proximity to actual compressed .rds size
  expect_gt(est_mb, actual_gz_mb * 0.4)
  expect_lt(est_mb, actual_gz_mb * 1.5)
})

test_that("ts_estimate_rds_size_mb supports optional compress mode presets", {
  raw_bytes <- 100 * 1024^2
  expect_equal(ts_estimate_rds_size_mb(raw_bytes, compress = TRUE), 22)
  expect_equal(ts_estimate_rds_size_mb(raw_bytes, compress = "gzip"), 22)
  expect_equal(ts_estimate_rds_size_mb(raw_bytes, compress = "xz"), 20)
  expect_equal(ts_estimate_rds_size_mb(raw_bytes, compress = "bzip2"), 21)
  expect_equal(ts_estimate_rds_size_mb(raw_bytes, compress = FALSE), 100)
  expect_equal(ts_estimate_rds_size_mb(raw_bytes, compress = "none"), 100)

  # Explicit compression_ratio takes precedence over compress preset
  expect_equal(ts_estimate_rds_size_mb(raw_bytes, compression_ratio = 0.5, compress = "xz"), 50)

  # Invalid compress values return NA_real_ without throwing
  for (bad in list(NA, "lz4", 123, c("gzip", "xz"))) {
    expect_identical(ts_estimate_rds_size_mb(raw_bytes, compress = bad), NA_real_)
  }
})
