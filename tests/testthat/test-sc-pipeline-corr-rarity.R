# =============================================================================
# test-sc-pipeline-corr-rarity.R — Regressions from live SC auto-pipeline run
# =============================================================================
# 1. find_correlated_genes() must handle S4 matrix backends (e.g. BPCells
#    IterableMatrix, dgCMatrix) without throwing:
#    "cannot coerce type 'S4' to vector of type 'double'".
# 2. APP_DAEMON_SOURCE_FILES must preload sc_population_rarity.R so
#    compute_population_rarity() is available inside mirai workers.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_helpers.R")
source_project_file("tools/pkg_skeleton/R/spatial_async.R")

test_that("APP_DAEMON_SOURCE_FILES preloads sc_population_rarity.R before sc_pipeline.R", {
  idx_rarity <- match("tools/pkg_skeleton/R/sc_population_rarity.R", APP_DAEMON_SOURCE_FILES)
  idx_pipe   <- match("tools/pkg_skeleton/R/sc_pipeline.R", APP_DAEMON_SOURCE_FILES)
  expect_false(is.na(idx_rarity))
  expect_false(is.na(idx_pipe))
  expect_true(idx_rarity < idx_pipe)
})

test_that("find_correlated_genes coerces non-dgCMatrix S4 matrices (BPCells-style) cleanly", {
  # Minimal S4 matrix class mimicking BPCells IterableMatrix behavior:
  # - inherits(x, "dgCMatrix") is FALSE
  # - x[i, ] returns another S4 object (so as.numeric(x[i, ]) fails without as.matrix())
  if (!methods::isClass("MockBPCellsMat")) {
    methods::setClass("MockBPCellsMat", slots = c(data = "matrix"))
  }
  methods::setMethod("rownames", "MockBPCellsMat", function(x) rownames(x@data))
  methods::setMethod("[", c("MockBPCellsMat", "ANY", "ANY", "ANY"),
    function(x, i, j, ..., drop = TRUE) {
      methods::new("MockBPCellsMat", data = x@data[i, j, drop = FALSE])
    }
  )
  methods::setMethod("as.matrix", "MockBPCellsMat", function(x, ...) x@data)

  mat <- matrix(
    c(1, 2, 3, 4, 5, 6,
      2, 4, 6, 8, 10, 12,
      6, 5, 4, 3, 2, 1),
    nrow = 3L, byrow = TRUE,
    dimnames = list(c("S100A8", "S100A9", "ACTB"), paste0("C", 1:6))
  )
  s4_mat <- methods::new("MockBPCellsMat", data = mat)

  # Stub .get_norm_matrix locally to return the S4 matrix
  old_fn <- .get_norm_matrix
  on.exit(assign(".get_norm_matrix", old_fn, envir = globalenv()), add = TRUE)
  assign(".get_norm_matrix", function(obj) s4_mat, envir = globalenv())

  res <- find_correlated_genes(s4_mat, target_gene = "S100A8", method = "pearson", threshold = 0.5)
  expect_s3_class(res, "data.frame")
  expect_true("S100A9" %in% res$gene)
  expect_equal(res$correlation[res$gene == "S100A9"], 1, tolerance = 1e-8)
})
