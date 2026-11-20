# =============================================================================
# test-sc-replicate-aberrations.R — Multi-replicate aberration report (W2.4 §2.1)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/sc_helpers.R")

.make_replicate_meta <- function(rep_sizes = c(s1 = 200L, s2 = 210L, s3 = 190L, s4 = 205L),
                                 nFeature_means = c(s1 = 1500, s2 = 1550, s3 = 1480, s4 = 1520),
                                 nCount_means = c(s1 = 4000, s2 = 4100, s3 = 3950, s4 = 4050),
                                 mt_means = c(s1 = 5, s2 = 5.2, s3 = 4.8, s4 = 5.1)) {
  dfs <- lapply(names(rep_sizes), function(nm) {
    n <- rep_sizes[[nm]]
    data.frame(
      orig.ident = rep(nm, n),
      nFeature_RNA = stats::rnorm(n, mean = nFeature_means[[nm]], sd = 50),
      nCount_RNA = stats::rnorm(n, mean = nCount_means[[nm]], sd = 100),
      percent.mt = pmax(0, stats::rnorm(n, mean = mt_means[[nm]], sd = 0.5)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, dfs)
}

test_that("detect_replicate_aberrations returns clean report on balanced replicates", {
  set.seed(42)
  meta <- .make_replicate_meta()
  res <- detect_replicate_aberrations(meta)

  expect_type(res, "list")
  expect_true(all(c("summary_df", "flagged_replicates", "has_aberrations") %in% names(res)))
  expect_false(res$has_aberrations)
  expect_identical(res$flagged_replicates, character(0))
  expect_s3_class(res$summary_df, "data.frame")
  expect_equal(nrow(res$summary_df), 4L)
  expect_true(all(c("replicate", "n_cells", "flagged", "reasons") %in% names(res$summary_df)))
  expect_true(all(c("nFeature_RNA_median", "nCount_RNA_median", "percent.mt_median") %in% names(res$summary_df)))
  expect_false(any(res$summary_df$flagged))
})

test_that("detect_replicate_aberrations flags low absolute cell count (< 50)", {
  set.seed(42)
  meta <- .make_replicate_meta(
    rep_sizes = c(s1 = 120L, s2 = 110L, s3 = 30L, s4 = 115L)
  )
  res <- detect_replicate_aberrations(meta)

  expect_true(res$has_aberrations)
  expect_identical(res$flagged_replicates, "s3")
  row_s3 <- res$summary_df[res$summary_df$replicate == "s3", , drop = FALSE]
  expect_true(row_s3$flagged)
  expect_match(row_s3$reasons, "low_cell_count")
})

test_that("detect_replicate_aberrations flags low relative cell count (< 0.25 * median)", {
  set.seed(42)
  meta <- .make_replicate_meta(
    rep_sizes = c(s1 = 1000L, s2 = 1050L, s3 = 80L, s4 = 980L)
  )
  res <- detect_replicate_aberrations(meta)

  expect_true(res$has_aberrations)
  expect_identical(res$flagged_replicates, "s3")
  row_s3 <- res$summary_df[res$summary_df$replicate == "s3", , drop = FALSE]
  expect_true(row_s3$flagged)
  expect_match(row_s3$reasons, "low_cell_count")
})

test_that("detect_replicate_aberrations flags extreme metric median via MAD (>= 3 replicates)", {
  set.seed(42)
  meta <- .make_replicate_meta(
    mt_means = c(s1 = 5.0, s2 = 5.2, s3 = 4.8, s4 = 35.0)
  )
  res <- detect_replicate_aberrations(meta, mad_k = 3)

  expect_true(res$has_aberrations)
  expect_identical(res$flagged_replicates, "s4")
  row_s4 <- res$summary_df[res$summary_df$replicate == "s4", , drop = FALSE]
  expect_true(row_s4$flagged)
  expect_match(row_s4$reasons, "extreme_percent\\.mt")
})

test_that("detect_replicate_aberrations skips MAD check when fewer than 3 replicates or MAD == 0", {
  set.seed(42)
  meta_2 <- .make_replicate_meta(
    rep_sizes = c(s1 = 200L, s2 = 200L),
    nFeature_means = c(s1 = 1500, s2 = 5000),
    nCount_means = c(s1 = 4000, s2 = 4000),
    mt_means = c(s1 = 5, s2 = 40)
  )
  res_2 <- detect_replicate_aberrations(meta_2)
  expect_false(res_2$has_aberrations)

  # Identical medians across 3 replicates -> MAD == 0, no division/zero-MAD false positives
  meta_const <- data.frame(
    orig.ident = rep(c("a", "b", "c"), each = 60L),
    nFeature_RNA = rep(1000, 180L),
    stringsAsFactors = FALSE
  )
  res_const <- detect_replicate_aberrations(meta_const)
  expect_false(res_const$has_aberrations)
})

test_that("detect_replicate_aberrations handles missing metrics, custom sample_col, and empty input", {
  set.seed(42)
  meta <- .make_replicate_meta()
  names(meta)[names(meta) == "orig.ident"] <- "sample_id"
  meta$percent.mt <- NULL

  res <- detect_replicate_aberrations(
    meta,
    sample_col = "sample_id",
    metrics = c("nFeature_RNA", "nonexistent_metric")
  )
  expect_false(res$has_aberrations)
  expect_true("nFeature_RNA_median" %in% names(res$summary_df))
  expect_false("nonexistent_metric_median" %in% names(res$summary_df))

  empty_meta <- meta[0, , drop = FALSE]
  res_empty <- detect_replicate_aberrations(empty_meta, sample_col = "sample_id")
  expect_false(res_empty$has_aberrations)
  expect_identical(res_empty$flagged_replicates, character(0))
  expect_equal(nrow(res_empty$summary_df), 0L)
})

test_that("detect_replicate_aberrations validates arguments with errorCondition", {
  meta <- data.frame(orig.ident = c("a", "b"), nFeature_RNA = c(100, 200))
  expect_error(detect_replicate_aberrations("not_a_df"), class = "sc_replicate_aberration_error")
  expect_error(detect_replicate_aberrations(meta, sample_col = "missing_col"), class = "sc_replicate_aberration_error")
  expect_error(detect_replicate_aberrations(meta, mad_k = -1), class = "sc_replicate_aberration_error")
})
