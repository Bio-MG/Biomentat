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
