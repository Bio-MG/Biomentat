# ==============================================================================
# test-bulk-de-mnar-dispatch.R — BK-05 MNAR/MAR DE dispatch & daemon preload
# ==============================================================================

.mnar_test_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("run_de_core dispatches engine='mnar' on non-integer matrix with NAs", {
  source(.mnar_test_root("tools/pkg_skeleton/R/batch_correction.R"), local = globalenv())
  source(.mnar_test_root("tools/pkg_skeleton/R/bulk_mnar_de.R"), local = globalenv())
  source(.mnar_test_root("tools/pkg_skeleton/R/bulk_helpers.R"), local = globalenv())

  set.seed(42)
  mat <- matrix(rnorm(60, mean = 10, sd = 0.4), nrow = 10, ncol = 6)
  rownames(mat) <- paste0("PROT", seq_len(10))
  colnames(mat) <- paste0("S", seq_len(6))
  # Inject MNAR pattern in feature 1 (all NA in ctrl, high in case)
  mat[1, 1:3] <- NA_real_
  mat[1, 4:6] <- c(14.1, 14.3, 13.9)

  meta <- data.frame(
    sample = colnames(mat),
    condition = c("ctrl", "ctrl", "ctrl", "case", "case", "case"),
    stringsAsFactors = FALSE
  )

  out <- run_de_core(
    counts_matrix = mat,
    metadata = meta,
    engine = "mnar",
    condition_col = "condition",
    group_target = "case",
    group_ref = "ctrl",
    p_adjust_method = "BH"
  )

  expect_true(is.list(out))
  expect_s3_class(out$res, "data.frame")
  expect_true(all(c("log2FoldChange", "pvalue", "padj", "baseMean", "regime") %in% colnames(out$res)))
  expect_equal(nrow(out$res), 10L)
  expect_equal(out$res["PROT1", "regime"], "MNAR")
  expect_equal(bulk_de_engine_label("mnar"), "PolySTest (MNAR/MAR)")
})

test_that("APP_DAEMON_SOURCE_FILES preloads bulk_mnar_de.R before bulk_helpers.R", {
  source(.mnar_test_root("tools/pkg_skeleton/R/spatial_async.R"), local = globalenv())
  idx_mnar <- match("tools/pkg_skeleton/R/bulk_mnar_de.R", APP_DAEMON_SOURCE_FILES)
  idx_help <- match("tools/pkg_skeleton/R/bulk_helpers.R", APP_DAEMON_SOURCE_FILES)
  expect_false(is.na(idx_mnar))
  expect_false(is.na(idx_help))
  expect_lt(idx_mnar, idx_help)
})
