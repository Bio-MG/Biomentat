# =============================================================================
# test-sc-spatial-regulon.R — P3.3: decoupleR TF/pathway regulons (SC & Spatial)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_spatial_regulon.R")

.mk_regulon_fixture <- function(n_genes = 20L, n_cells = 16L, seed = 42L) {
  set.seed(seed)
  genes <- paste0("G", seq_len(n_genes))
  cells <- paste0("C", seq_len(n_cells))
  mat <- matrix(stats::rnorm(n_genes * n_cells, mean = 0, sd = 0.3), nrow = n_genes, ncol = n_cells,
                dimnames = list(genes, cells))
  # TF1 activates G1..G6 in GroupA (cells 1..8)
  mat[1:6, 1:8] <- mat[1:6, 1:8] + 3.0
  # TF2 activates G7..G12 in GroupB (cells 9..16)
  mat[7:12, 9:16] <- mat[7:12, 9:16] + 3.0

  net <- data.frame(
    source = c(rep("TF1", 6L), rep("TF2", 6L)),
    target = paste0("G", 1:12),
    mor = rep(1, 12L),
    stringsAsFactors = FALSE
  )
  groups <- stats::setNames(c(rep("GroupA", 8L), rep("GroupB", 8L)), cells)
  list(mat = mat, net = net, groups = groups)
}

test_that("run_sc_spatial_regulons validates inputs and network schema", {
  fx <- .mk_regulon_fixture()

  err1 <- tryCatch(run_sc_spatial_regulons(NULL, fx$net), error = function(e) e)
  expect_s3_class(err1, "sc_spatial_regulon_error")

  err2 <- tryCatch(run_sc_spatial_regulons(fx$mat, fx$net[, c("source", "target")]), error = function(e) e)
  expect_s3_class(err2, "sc_spatial_regulon_error")

  err3 <- tryCatch(run_sc_spatial_regulons(fx$mat, fx$net, method = "invalid"), error = function(e) e)
  expect_s3_class(err3, "sc_spatial_regulon_error")
})

test_that("run_sc_spatial_regulons computes ULM and MLM regulon scores + group summary", {
  skip_if_not_installed("decoupleR")
  fx <- .mk_regulon_fixture()

  res_ulm <- run_sc_spatial_regulons(fx$mat, fx$net, method = "ulm", minsize = 5L, group_by = fx$groups)
  expect_identical(res_ulm$method, "ulm")
  expect_identical(res_ulm$n_regulons, 2L)
  expect_true(is.matrix(res_ulm$score_matrix))
  expect_identical(dim(res_ulm$score_matrix), c(2L, 16L))
  expect_true(is.data.frame(res_ulm$group_summary))

  ga_tf1 <- res_ulm$group_summary$mean_score[res_ulm$group_summary$source == "TF1" &
                                             res_ulm$group_summary$group == "GroupA"]
  gb_tf1 <- res_ulm$group_summary$mean_score[res_ulm$group_summary$source == "TF1" &
                                             res_ulm$group_summary$group == "GroupB"]
  expect_gt(ga_tf1, gb_tf1 + 1.0)

  res_mlm <- run_sc_spatial_regulons(fx$mat, fx$net, method = "mlm", minsize = 5L)
  expect_identical(res_mlm$method, "mlm")
  expect_identical(dim(res_mlm$score_matrix), c(2L, 16L))
})
