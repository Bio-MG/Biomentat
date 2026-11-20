# =============================================================================
# test-sc-umap-render-budget.R — Adaptive UMAP point-budget & stratified subsampling
# =============================================================================
# Tests unitaires pour sc_umap_render_budget() et sc_subsample_umap_indices()
# dans tools/pkg_skeleton/R/sc_plotting.R (W1.1 #13).
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
load_project_pkg()
source_project_file("tools/pkg_skeleton/R/sc_plotting.R")

test_that("sc_umap_render_budget classifies small, 28k-50k, and >50k datasets", {
  # Small dataset: no subsampling, no raster
  b_small <- sc_umap_render_budget(5000L)
  expect_false(b_small$should_subsample)
  expect_false(b_small$use_raster)
  expect_identical(b_small$target_points, 5000L)
  expect_identical(b_small$mode, "full")
  expect_equal(b_small$subsample_fraction, 1)

  # Boundary at max_interactive (25k)
  b_25k <- sc_umap_render_budget(25000L)
  expect_false(b_25k$should_subsample)
  expect_false(b_25k$use_raster)
  expect_identical(b_25k$target_points, 25000L)
  expect_identical(b_25k$mode, "full")

  # 28k cells: subsample to 25k, still interactive vector (<= 50k)
  b_28k <- sc_umap_render_budget(28000L)
  expect_true(b_28k$should_subsample)
  expect_false(b_28k$use_raster)
  expect_identical(b_28k$target_points, 25000L)
  expect_identical(b_28k$mode, "subsampled")
  expect_equal(b_28k$subsample_fraction, 25000 / 28000)

  # Boundary at raster_threshold (50k): still <= 50k so use_raster = FALSE
  b_50k <- sc_umap_render_budget(50000L)
  expect_true(b_50k$should_subsample)
  expect_false(b_50k$use_raster)
  expect_identical(b_50k$target_points, 25000L)
  expect_identical(b_50k$mode, "subsampled")

  # 500k cells: both subsample and rasterize
  b_500k <- sc_umap_render_budget(500000L)
  expect_true(b_500k$should_subsample)
  expect_true(b_500k$use_raster)
  expect_identical(b_500k$target_points, 25000L)
  expect_identical(b_500k$mode, "raster")
  expect_equal(b_500k$subsample_fraction, 25000 / 500000)
})

test_that("sc_umap_render_budget handles zero cells and custom thresholds", {
  b_zero <- sc_umap_render_budget(0L)
  expect_false(b_zero$should_subsample)
  expect_false(b_zero$use_raster)
  expect_identical(b_zero$target_points, 0L)
  expect_identical(b_zero$mode, "full")

  b_custom <- sc_umap_render_budget(15000L, max_interactive = 10000L, raster_threshold = 12000L)
  expect_true(b_custom$should_subsample)
  expect_true(b_custom$use_raster)
  expect_identical(b_custom$target_points, 10000L)
  expect_identical(b_custom$mode, "raster")
})

test_that("sc_umap_render_budget validates arguments with sc_plotting_error", {
  expect_error(sc_umap_render_budget(-1L), class = "sc_plotting_error")
  expect_error(sc_umap_render_budget(NA_integer_), class = "sc_plotting_error")
  expect_error(sc_umap_render_budget(c(10L, 20L)), class = "sc_plotting_error")
  expect_error(sc_umap_render_budget(100L, max_interactive = 0L), class = "sc_plotting_error")
  expect_error(sc_umap_render_budget(100L, raster_threshold = 0L), class = "sc_plotting_error")
})

test_that("sc_subsample_umap_indices returns all indices when n_cells <= max_points", {
  expect_identical(sc_subsample_umap_indices(0L), integer(0))
  expect_identical(sc_subsample_umap_indices(100L, max_points = 25000L), seq_len(100L))
  expect_identical(sc_subsample_umap_indices(25000L, max_points = 25000L), seq_len(25000L))
})

test_that("sc_subsample_umap_indices is deterministic and bounded for 28k-500k cells", {
  idx1 <- sc_subsample_umap_indices(28000L, max_points = 25000L, seed = 42L)
  idx2 <- sc_subsample_umap_indices(28000L, max_points = 25000L, seed = 42L)
  idx3 <- sc_subsample_umap_indices(28000L, max_points = 25000L, seed = 99L)

  expect_type(idx1, "integer")
  expect_length(idx1, 25000L)
  expect_identical(idx1, idx2)
  expect_false(identical(idx1, idx3))
  expect_true(all(idx1 >= 1L & idx1 <= 28000L))
  expect_identical(idx1, sort(unique(idx1)))

  idx_500k <- sc_subsample_umap_indices(500000L, max_points = 25000L, seed = 42L)
  expect_length(idx_500k, 25000L)
  expect_true(all(idx_500k >= 1L & idx_500k <= 500000L))
  expect_identical(idx_500k, sort(unique(idx_500k)))
})

test_that("sc_subsample_umap_indices preserves rare clusters via stratified sampling", {
  # 28,000 cells: 27,980 in C1, 17 in C2, 3 in C3 (ultra-rare)
  groups <- factor(c(rep("C1", 27980L), rep("C2", 17L), rep("C3", 3L)))
  idx1 <- sc_subsample_umap_indices(28000L, max_points = 25000L, groups = groups, seed = 42L)
  idx2 <- sc_subsample_umap_indices(28000L, max_points = 25000L, groups = groups, seed = 42L)

  expect_identical(idx1, idx2)
  expect_length(idx1, 25000L)
  expect_identical(idx1, sort(unique(idx1)))

  sampled_groups <- groups[idx1]
  expect_true(sum(sampled_groups == "C3") >= 1L)
  expect_true(sum(sampled_groups == "C2") >= 1L)
  expect_true(sum(sampled_groups == "C1") >= 1L)

  # Even with aggressive subsampling (500k -> 1000), a 1-cell rare cluster is kept
  groups_500k <- c(rep("MajorA", 300000L), rep("MajorB", 199999L), "Singleton")
  idx_aggr <- sc_subsample_umap_indices(500000L, max_points = 1000L, groups = groups_500k, seed = 7L)
  expect_length(idx_aggr, 1000L)
  expect_true(500000L %in% idx_aggr)
})

test_that("sc_subsample_umap_indices handles NA groups and K > max_points", {
  groups_na <- c(rep("A", 40L), rep(NA_character_, 10L))
  idx_na <- sc_subsample_umap_indices(50L, max_points = 20L, groups = groups_na, seed = 1L)
  expect_length(idx_na, 20L)
  expect_true(any(is.na(groups_na[idx_na])))

  # More groups (30) than max_points (10): still returns exactly 10 unique indices
  groups_many <- paste0("G", seq_len(30L))
  idx_many <- sc_subsample_umap_indices(30L, max_points = 10L, groups = groups_many, seed = 1L)
  expect_length(idx_many, 10L)
  expect_identical(idx_many, sort(unique(idx_many)))
})

test_that("sc_subsample_umap_indices validates arguments with sc_plotting_error", {
  expect_error(sc_subsample_umap_indices(-5L), class = "sc_plotting_error")
  expect_error(sc_subsample_umap_indices(100L, max_points = 0L), class = "sc_plotting_error")
  expect_error(sc_subsample_umap_indices(100L, max_points = 10L, groups = c("A", "B")),
               class = "sc_plotting_error")
  expect_error(sc_subsample_umap_indices(100L, max_points = 10L, seed = NA_integer_),
               class = "sc_plotting_error")
})
