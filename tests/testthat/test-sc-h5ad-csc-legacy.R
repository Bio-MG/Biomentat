# =============================================================================
# test-sc-h5ad-csc-legacy.R — tests for sc_normalize_h5ad_sparse_layout()
# (W2.4 §2.7, tools/pkg_skeleton/R/sc_bpcells.R)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_bpcells.R")

test_that("sc_normalize_h5ad_sparse_layout handles modern csr_matrix and csc_matrix", {
  res_csr <- sc_normalize_h5ad_sparse_layout(
    encoding_type = "csr_matrix",
    dims = c(500L, 2000L),
    n_obs = 500L,
    n_vars = 2000L
  )
  expect_type(res_csr, "list")
  expect_named(res_csr, c("layout", "needs_transpose", "legacy_fallback", "note"))
  expect_identical(res_csr$layout, "csr")
  expect_true(res_csr$needs_transpose)
  expect_false(res_csr$legacy_fallback)
  expect_true(nzchar(res_csr$note))

  res_csc <- sc_normalize_h5ad_sparse_layout(
    encoding_type = "csc_matrix",
    dims = c(500L, 2000L),
    n_obs = 500L,
    n_vars = 2000L
  )
  expect_identical(res_csc$layout, "csc")
  expect_true(res_csc$needs_transpose)
  expect_false(res_csc$legacy_fallback)
})

test_that("sc_normalize_h5ad_sparse_layout detects already genes-x-cells orientation", {
  res <- sc_normalize_h5ad_sparse_layout(
    encoding_type = "csc_matrix",
    dims = c(2000L, 500L),
    n_obs = 500L,
    n_vars = 2000L
  )
  expect_identical(res$layout, "csc")
  expect_false(res$needs_transpose)
  expect_false(res$legacy_fallback)
})

test_that("sc_normalize_h5ad_sparse_layout falls back cleanly when encoding-type is missing (legacy h5ad)", {
  for (enc in list(NULL, "", NA_character_, character(0))) {
    res <- sc_normalize_h5ad_sparse_layout(
      encoding_type = enc,
      dims = c(300L, 1200L),
      n_obs = 300L,
      n_vars = 1200L
    )
    expect_identical(res$layout, "csr")
    expect_true(res$needs_transpose)
    expect_true(res$legacy_fallback)
    expect_match(res$note, "legacy", ignore.case = TRUE)
  }
})

test_that("sc_normalize_h5ad_sparse_layout handles NULL dims and square matrices", {
  res_no_dims <- sc_normalize_h5ad_sparse_layout(encoding_type = "csr_matrix")
  expect_identical(res_no_dims$layout, "csr")
  expect_true(res_no_dims$needs_transpose)
  expect_false(res_no_dims$legacy_fallback)

  res_sq <- sc_normalize_h5ad_sparse_layout(
    encoding_type = "csc_matrix",
    dims = c(100L, 100L),
    n_obs = 100L,
    n_vars = 100L
  )
  expect_identical(res_sq$layout, "csc")
  expect_true(res_sq$needs_transpose)
})

test_that("sc_normalize_h5ad_sparse_layout handles dense array encoding and rejects invalid inputs", {
  res_dense <- sc_normalize_h5ad_sparse_layout(
    encoding_type = "array",
    dims = c(100L, 50L),
    n_obs = 100L,
    n_vars = 50L
  )
  expect_identical(res_dense$layout, "dense")
  expect_true(res_dense$needs_transpose)
  expect_false(res_dense$legacy_fallback)

  expect_error(
    sc_normalize_h5ad_sparse_layout(encoding_type = "coo_matrix"),
    class = "sc_bpcells_error"
  )
  expect_error(
    sc_normalize_h5ad_sparse_layout(
      encoding_type = "csr_matrix",
      dims = c(10L, 20L),
      n_obs = 30L,
      n_vars = 40L
    ),
    class = "sc_bpcells_error"
  )
  expect_error(
    sc_normalize_h5ad_sparse_layout(encoding_type = "csr_matrix", dims = c(10L)),
    class = "sc_bpcells_error"
  )
})
