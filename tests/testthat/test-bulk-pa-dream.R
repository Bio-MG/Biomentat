# =============================================================================
# test-bulk-pa-dream.R — GEMS PA v2a : moteur mixte variancePartition::dream()
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/validation.R")
source_project_file("tools/pkg_skeleton/R/batch_correction.R")
source_project_file("tools/pkg_skeleton/R/memory.R")
source_project_file("tools/pkg_skeleton/R/bulk_helpers.R")

.dream_fixture <- function(n_pairs = 4L, subject_sd = 2.5, n_de = 40L, n_genes = 400L) {
  set.seed(20261113)
  patient <- factor(rep(sprintf("P%d", seq_len(n_pairs)), each = 2L))
  cond <- factor(rep(c("ref", "target"), n_pairs), levels = c("ref", "target"))
  subj_offset <- stats::rnorm(n_pairs, 0, subject_sd)
  base <- stats::rgamma(n_genes, shape = 2, scale = 50) + 5
  lfc_true <- rep(0, n_genes)
  lfc_true[seq_len(n_de)] <- 1.5
  samples <- sprintf("S%d", seq_len(2 * n_pairs))
  counts <- vapply(seq_len(2 * n_pairs), function(i) {
    p <- ceiling(i / 2)
    is_target <- (i %% 2L == 0L)
    lambda <- base * 2^(lfc_true * is_target) * exp(subj_offset[p])
    stats::rpois(n_genes, lambda)
  }, numeric(n_genes))
  rownames(counts) <- sprintf("G%04d", seq_len(n_genes))
  colnames(counts) <- samples
  meta <- data.frame(condition = cond, patient = patient, row.names = samples)
  list(counts = counts, meta = meta, de_idx = seq_len(n_de), n_genes = n_genes)
}

.dream_cols_ok <- function(res) {
  all(c("log2FoldChange", "pvalue", "padj", "gene") %in% colnames(res))
}

.dream_lme4_incompat <- function() {
  requireNamespace("lme4", quietly = TRUE) &&
    requireNamespace("variancePartition", quietly = TRUE) &&
    utils::packageVersion("lme4") >= "2.0.0" &&
    utils::packageVersion("variancePartition") < "1.38.0"
}

test_that("PA-v2a : run_dream_de retourne le schema canonique et gagne de la puissance vs limma non-apparie", {
  skip_if_not_installed("variancePartition")
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")
  skip_if(.dream_lme4_incompat(),
          "variancePartition < 1.38 incompatible avec lme4 >= 2.0 (as_lmerModLmerTest2 attend formal 'theta', lme4 2.0.6 expose 'par')")

  d <- .dream_fixture()
  rss_before <- ts_process_rss_mb()
  unpaired <- run_limma_voom_de(d$counts, d$meta, "condition", "target", "ref")
  dream_res <- run_dream_de(d$counts, d$meta, "condition", "target", "ref",
                            pair_col = "patient")
  rss_after <- ts_process_rss_mb()

  expect_true(.dream_cols_ok(dream_res))
  expect_identical(nrow(dream_res), d$n_genes)
  expect_true(all(is.finite(dream_res$log2FoldChange)))
  expect_true(all(is.finite(dream_res$pvalue)))
  expect_true(is.numeric(rss_before) && is.numeric(rss_after))

  p_un <- unpaired$pvalue[match(dream_res$gene, unpaired$gene)]
  p_dr <- dream_res$pvalue
  de_genes <- sprintf("G%04d", d$de_idx)
  is_de <- dream_res$gene %in% de_genes
  expect_lt(median(p_dr[is_de], na.rm = TRUE),
            median(p_un[is_de], na.rm = TRUE))
})

test_that("PA-v2a : run_bulk_de_dispatch(engine = 'dream') et bulk_de_engine_label('dream')", {
  expect_identical(bulk_de_engine_label("dream"), "limma-dream")
  skip_if_not_installed("variancePartition")
  skip_if(.dream_lme4_incompat(),
          "variancePartition < 1.38 incompatible avec lme4 >= 2.0 (formal 'par' vs 'theta')")
  d <- .dream_fixture()
  res <- run_bulk_de_dispatch(
    engine = "dream",
    counts_matrix = d$counts,
    metadata = d$meta,
    condition_col = "condition",
    group_target = "target",
    group_ref = "ref",
    pair_col = "patient"
  )
  expect_true(.dream_cols_ok(res))
})

test_that("PA-v2a : echec dream() encapsule en bulk_de_error (jamais d'erreur brute)", {
  skip_if_not_installed("variancePartition")
  if (.dream_lme4_incompat()) {
    d <- .dream_fixture()
    expect_error(
      run_bulk_de_dispatch(
        engine = "dream",
        counts_matrix = d$counts,
        metadata = d$meta,
        condition_col = "condition",
        group_target = "target",
        group_ref = "ref",
        pair_col = "patient"
      ),
      class = "bulk_de_error"
    )
  } else {
    succeed("lme4/variancePartition compatibles sur cet hote")
  }
})

test_that("PA-v2a : refus honnete bulk_de_error si pair_col absent, singleton ou treat_lfc > 0", {
  skip_if_not_installed("variancePartition")
  d <- .dream_fixture()
  expect_error(
    run_dream_de(d$counts, d$meta, "condition", "target", "ref", pair_col = NULL),
    class = "bulk_de_error"
  )
  expect_error(
    run_dream_de(d$counts, d$meta, "condition", "target", "ref", pair_col = "missing_col"),
    class = "bulk_de_error"
  )
  meta_single <- d$meta
  meta_single$patient <- as.character(meta_single$patient)
  meta_single$patient[1] <- "P_SINGLE"
  expect_error(
    run_dream_de(d$counts, meta_single, "condition", "target", "ref", pair_col = "patient"),
    class = "bulk_de_error"
  )
  expect_error(
    run_dream_de(d$counts, d$meta, "condition", "target", "ref", pair_col = "patient", treat_lfc = 1),
    class = "bulk_de_error"
  )
})
