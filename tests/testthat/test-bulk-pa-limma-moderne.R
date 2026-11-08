# =============================================================================
# test-bulk-pa-limma-moderne.R — GEMS PA v1 (fiche
# PROPOSAL_GEMS_IMPROVEMENT_PACKAGES_2026-11-04.md §PA, décisions D1/D2) :
# « limma moderne apparié & pondéré ».
# =============================================================================
# Portée du jalon PA v1 (fiche §PA + D1/D2) :
#   1. designs appariés/bloqués : bloc `duplicateCorrelation` côté moteur —
#      les designs appariés (structure GSE52778 4×2) cessent d'être analysés
#      comme indépendants ;
#   2. `voomWithQualityWeights` (pondération qualité) ;
#   3. option `treat_lfc` (limma::treat) — test contre un seuil de FC
#      (kill la classe « significatif mais trivial ») ;
#   4. garde `check_design_cancor` (limma::canCorPairs) — complète la garde
#      catégorielle `check_design_confounding` (validation.R) qui ne voit
#      que le catégoriel.
# D1 : `dream()` est HORS périmètre v1 (seconde vague, derrière re-mesure).
#
# FIXTURE : appariée SYNTHÉTIQUE 3 patients × 2 conditions (structure
# GSE52778 documentée dans la fiche — le subset n'est pas sur disque ; le
# moteur est validé sur une appariée synthétique DÉTERMINISTE, seed posée).
# Effet sujet FORT (offset log par patient) : l'analyse appariée doit
# GAGNER de la puissance (p-values plus petites) sur les gènes réellement DE.
# La fumée du chemin non-apparié reste couverte par test-bulk-helpers.R.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/validation.R")
# bulk_assert_raw_counts vit dans le domaine correction de batch (même motif
# que test-bulk-helpers.R:33) — run_bulk_de_dispatch l'appelle systématiquement.
source_project_file("tools/pkg_skeleton/R/batch_correction.R")
source_project_file("tools/pkg_skeleton/R/bulk_helpers.R")

.pa_fixture <- function(n_pairs = 3L, subject_sd = 2.5, n_de = 40L, n_genes = 400L) {
  set.seed(20261108)
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
  colnames(counts) <- samples
  meta <- data.frame(condition = cond, patient = patient, row.names = samples)
  list(counts = counts, meta = meta, de_idx = seq_len(n_de), n_genes = n_genes)
}

.pa_cols_ok <- function(res) {
  all(c("log2FoldChange", "pvalue", "padj", "gene") %in% colnames(res))
}

test_that("PA-1 : run_limma_voom_de accepte pair_col et applique duplicateCorrelation (gagne de la puissance)", {
  d <- .pa_fixture()
  unpaired <- run_limma_voom_de(d$counts, d$meta, "condition", "target", "ref")
  paired   <- run_limma_voom_de(d$counts, d$meta, "condition", "target", "ref",
                                pair_col = "patient")
  expect_true(.pa_cols_ok(paired))
  expect_identical(nrow(paired), d$n_genes)
  # Puissance appariée : sur les gènes VRAIMENT DE, les p-values appariées
  # (effet sujet modélisé par duplicateCorrelation) sont plus petites que
  # celles de l'analyse qui traite les sujets comme indépendants.
  p_un <- unpaired$pvalue[match(rownames(paired), rownames(unpaired))]
  # NB : expect_lt (testthat 3.3.2) n'a PAS de paramètre info= — le message
  # vit dans le titre du bloc et le commentaire ci-dessus.
  # Assertion : le bloc apparié doit gagner de la puissance sur les gènes DE.
  expect_lt(median(paired$pvalue[d$de_idx], na.rm = TRUE),
            median(p_un[d$de_idx], na.rm = TRUE))
})

test_that("PA-2 : treat_lfc teste contre un seuil de FC (kill « significatif mais trivial »)", {
  d <- .pa_fixture()
  r0  <- run_limma_voom_de(d$counts, d$meta, "condition", "target", "ref", treat_lfc = 0)
  r25 <- run_limma_voom_de(d$counts, d$meta, "condition", "target", "ref", treat_lfc = 2.5)
  n0  <- sum(r0$padj < 0.05, na.rm = TRUE)
  n25 <- sum(r25$padj < 0.05, na.rm = TRUE)
  expect_true(n0 >= 1L, info = "le fixture doit produire au moins un gène sig. sans seuil")
  # Assertion : treat_lfc élevé doit réduire les « significatifs » (test
  # contre un seuil de FC — expect_lt n'a pas de info=, testthat 3.3.2).
  expect_lt(n25, n0)
  expect_true(.pa_cols_ok(r25))
})

test_that("PA-3 : refus honnêtes — pair_col inconnue / blocs singleton / moteur non-limma", {
  d <- .pa_fixture()
  expect_error(run_limma_voom_de(d$counts, d$meta, "condition", "target", "ref",
                                 pair_col = "introuvable"), class = "bulk_de_error")
  # Chaque bloc = 1 échantillon : duplicateCorrelation n'a aucun sens.
  d$meta$singleton <- factor(sprintf("b%d", seq_len(nrow(d$meta))))
  expect_error(run_limma_voom_de(d$counts, d$meta, "condition", "target", "ref",
                                 pair_col = "singleton"), class = "bulk_de_error")
  # Le bloc apparié n'est supporté que par le moteur limma en v1 (D1).
  expect_error(run_bulk_de_dispatch("edger", d$counts, d$meta, "condition", "target", "ref",
                                    pair_col = "patient"), class = "bulk_de_error")
})

test_that("PA-4 : dispatch route pair_col/treat_lfc vers le moteur limma", {
  d <- .pa_fixture()
  out <- run_bulk_de_dispatch("limma", d$counts, d$meta, "condition", "target", "ref",
                              pair_col = "patient", treat_lfc = 0)
  expect_true(.pa_cols_ok(out))
})

test_that("PA-5 : check_design_cancor détecte la confusion numérique (canCorPairs)", {
  d <- .pa_fixture()
  # Covariable numérique quasi-copiée sur la condition : canCorPairs doit tirer.
  d$meta$age <- ifelse(d$meta$condition == "target", 60, 40) + stats::rnorm(nrow(d$meta), 0, 0.5)
  # Covariable numérique indépendante : ne doit PAS tirer.
  d$meta$rnd <- stats::rnorm(nrow(d$meta))
  hit  <- check_design_cancor(d$meta, "condition", "age")
  miss <- check_design_cancor(d$meta, "condition", "rnd")
  expect_true(isTRUE(hit$confounded), info = "la covariable numérique collée à la condition est une confusion")
  expect_true(hit$correlation > 0.9, info = paste("correlation mesurée :", hit$correlation))
  # NB : la gemme canCorPairs n existe pas dans limma 3.62.2 (mesuré) — la
  expect_false(isTRUE(miss$confounded))
  # Covariable catégorielle : hors périmètre (garde catégorielle existante).
  expect_false(isTRUE(check_design_cancor(d$meta, "condition", "patient")$confounded))
})
