# =============================================================================
# test-bulk-multiblock.R — BK-06 Multi-Block DIABLO & NEMO Partial-Block Tests
# =============================================================================

source_project_file("tools/pkg_skeleton/R/bulk_multiblock.R")

.make_multiblock_fixture <- function(n_per_class = 12L, n_classes = 2L, seed = 101L) {
  set.seed(seed)
  n <- n_per_class * n_classes
  snames <- paste0("S", seq_len(n))
  y <- factor(rep(paste0("Grp", seq_len(n_classes)), each = n_per_class))
  names(y) <- snames

  # Shared latent factor separating classes
  latent <- as.numeric(y) - mean(as.numeric(y))

  mk_block <- function(p, prefix, signal_idx = 1:5, signal_strength = 2.5) {
    mat <- matrix(stats::rnorm(p * n, sd = 0.6), nrow = p, ncol = n,
                  dimnames = list(paste0(prefix, seq_len(p)), snames))
    for (g in signal_idx) {
      mat[g, ] <- mat[g, ] + signal_strength * latent
    }
    mat
  }

  list(
    blocks = list(
      rna  = mk_block(30L, "G_", signal_idx = 1:5, signal_strength = 2.5),
      prot = mk_block(20L, "P_", signal_idx = 1:4, signal_strength = 2.2),
      met  = mk_block(15L, "M_", signal_idx = 1:3, signal_strength = 2.0)
    ),
    y = y
  )
}

test_that("BK-06 validation rejects malformed multi-block inputs with typed errors", {
  fx <- .make_multiblock_fixture()

  expect_error(
    run_multiblock_diablo(list(rna = fx$blocks$rna), fx$y),
    class = "bulk_multiblock_error"
  )
  expect_error(
    run_multiblock_diablo(list(fx$blocks$rna, fx$blocks$prot), fx$y),
    class = "bulk_multiblock_error"
  )
  bad_inf <- fx$blocks
  bad_inf$rna[1L, 1L] <- Inf
  expect_error(
    run_multiblock_diablo(bad_inf, fx$y),
    class = "bulk_multiblock_error"
  )

  # Single-class outcome rejected
  y_one <- factor(rep("OnlyOne", length(fx$y)))
  names(y_one) <- names(fx$y)
  expect_error(
    run_multiblock_diablo(fx$blocks, y_one),
    class = "bulk_multiblock_error"
  )

  # Missing sample in named y rejected
  y_short <- fx$y[-1L]
  expect_error(
    run_multiblock_diablo(fx$blocks, y_short),
    class = "bulk_multiblock_error"
  )

  # Disjoint samples rejected for DIABLO (requires common samples)
  b_disjoint <- list(
    a = fx$blocks$rna[, 1:6],
    b = fx$blocks$prot[, 7:12]
  )
  expect_error(
    run_multiblock_diablo(b_disjoint, fx$y),
    class = "bulk_multiblock_error"
  )

  expect_error(
    run_multiblock_diablo(fx$blocks, fx$y, design_weight = 1.5),
    class = "bulk_multiblock_error"
  )
  expect_error(
    run_multiblock_nemo(fx$blocks, k = 1L),
    class = "bulk_multiblock_error"
  )
})

test_that("run_multiblock_diablo integrates 3 blocks, selects signal features, and achieves high CV accuracy", {
  fx <- .make_multiblock_fixture(n_per_class = 12L, n_classes = 2L, seed = 202L)

  res <- run_multiblock_diablo(
    blocks = fx$blocks,
    y = fx$y,
    ncomp = 2L,
    keepX = c(rna = 5L, prot = 4L, met = 3L),
    design_weight = 0.2,
    cv_folds = 4L,
    seed = 42L
  )

  expect_s3_class(res, "biomentat_diablo")
  expect_identical(res$method_used, "pure_r_sgcca_diablo")
  expect_identical(res$ncomp, 2L)
  expect_identical(names(res$variates), c("rna", "prot", "met"))
  expect_identical(dim(res$consensus_variates), c(24L, 2L))

  # Inter-block correlation on comp1 is strongly positive
  expect_true(res$block_correlations["rna", "prot"] > 0.6)
  expect_true(res$block_correlations["rna", "met"] > 0.6)

  # Sparsity respected and planted signal features recovered
  expect_true(length(res$selected_features$rna) <= 10L)
  expect_true(length(res$selected_features$prot) <= 8L)
  expect_true(any(paste0("G_", 1:5) %in% res$selected_features$rna))
  expect_true(any(paste0("P_", 1:4) %in% res$selected_features$prot))

  # Cross-validation separates planted classes cleanly
  expect_true(res$cv_accuracy >= 0.85)
  expect_true(res$cv_balanced_accuracy >= 0.85)
})

test_that("run_multiblock_diablo supports 3-class outcomes, NAs, zero-variance features, and determinism", {
  fx <- .make_multiblock_fixture(n_per_class = 10L, n_classes = 3L, seed = 303L)
  fx$blocks$rna[30L, ] <- 5.0 # zero-variance row
  fx$blocks$prot[10L, 1L] <- NA_real_ # isolated NA

  res1 <- run_multiblock_diablo(fx$blocks, fx$y, ncomp = 2L, cv_folds = 5L, seed = 77L)
  res2 <- run_multiblock_diablo(fx$blocks, fx$y, ncomp = 2L, cv_folds = 5L, seed = 77L)

  expect_identical(res1$predictions, res2$predictions)
  expect_equal(res1$consensus_variates, res2$consensus_variates)
  expect_false("G_30" %in% rownames(res1$loadings$rna))
  expect_true(res1$cv_accuracy >= 0.80)
})

test_that("run_multiblock_nemo fuses partial-block samples and recovers true clusters", {
  fx <- .make_multiblock_fixture(n_per_class = 8L, n_classes = 3L, seed = 404L)
  # Introduce partial sample coverage: drop 4 samples from prot and 4 different from met
  partial_blocks <- list(
    rna  = fx$blocks$rna,
    prot = fx$blocks$prot[, -c(2L, 10L, 18L, 22L), drop = FALSE],
    met  = fx$blocks$met[,  -c(5L, 13L, 20L, 24L), drop = FALSE]
  )

  res <- run_multiblock_nemo(partial_blocks, k = 3L, seed = 42L)

  expect_s3_class(res, "biomentat_nemo")
  expect_identical(res$method_used, "pure_r_nemo_partial")
  expect_identical(res$n_clusters, 3L)
  expect_length(res$clusters, 24L)
  expect_identical(dim(res$block_coverage), c(24L, 3L))
  expect_false(res$block_coverage["S2", "prot"])
  expect_false(res$block_coverage["S5", "met"])
  expect_identical(res$n_shared_blocks["S2", "S5"], 1L) # only rna measures both S2 and S5
  expect_true(isSymmetric(res$fused_similarity))
  expect_equal(unname(diag(res$fused_similarity)), rep(1, 24L))
  expect_true(res$silhouette_mean > 0.15)

  # Cluster purity against true 3 classes should be high despite missing blocks
  tab <- table(res$clusters, fx$y[names(res$clusters)])
  purity <- sum(apply(tab, 1L, max)) / length(res$clusters)
  expect_true(purity >= 0.90)

  # Automatic eigengap mode also runs cleanly
  res_auto <- run_multiblock_nemo(partial_blocks, k = NULL, max_k = 5L, seed = 42L)
  expect_true(res_auto$n_clusters >= 2L && res_auto$n_clusters <= 5L)
})
