# =============================================================================
# test-golden-bulk-de-parity.R — GOLDEN / parité numérique (P1.4, PRE_ROADMAP
# 2026-11-12 §P1.4 ; audit externe #5 exigence S-tier « parité numérique vs
# scripts de référence dans une tolérance explicite »)
# =============================================================================
# Principe : le moteur d'un onglet et la recette de référence (le package
# upstream appelé DIRECTEMENT, recette miroir documentée file:line) doivent
# produire les MÊMES nombres sur une fixture déterministe. Toute dérive
# silencieuse de recette (changements d'arguments, réordonnancement, defaults)
# casse la parité et rend ces tests rouges AVANT qu'un utilisateur ne voie un
# résultat différent de la littérature.
#
# Fixture : counts ENTIERS déterministes (graine 20261113), 200 gènes × 12
# échantillons (2 groupes × 6), 30 gènes différentiels à effet connu — aucun
# jeu externe requis (reproductible sur tout hôte).
#
# Tolérances EXPLICITES : les moteurs appellent les mêmes packages avec les
# mêmes entrées → parité attendue quasi exacte (1e-8 sur logFC, 1e-6 relatif
# sur p-values). Un écart au-delà = dérive de recette.
#
# Transformations spatiales : propriétés affines de apply_coord_orientation()
# (spatial_io.R) — involution des miroirs, préservation de l'emprise bbox en
# repli sans image (Xenium/CosMx), miroir autour de l'étendue image sinon.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/bulk_helpers.R")
source_project_file("tools/pkg_skeleton/R/spatial_io.R")

.golden_fixture <- function(seed = 20261113L) {
  set.seed(seed)
  ng <- 200L
  n_per_group <- 6L
  ns <- 2L * n_per_group
  base_expr <- 2^stats::rnorm(ng, mean = 5, sd = 1)   # counts moyens ~30..3000
  effect <- rep(0, ng)
  effect[1:30] <- 1.2                                  # log2FC vrai (up dans target)
  cond <- factor(rep(c("ref", "target"), each = n_per_group), levels = c("ref", "target"))
  sf <- stats::runif(ns, 0.8, 1.2)                     # tailles de librairie relatives
  mu <- matrix(NA_real_, nrow = ng, ncol = ns)
  for (i in seq_len(ns)) {
    mu[, i] <- base_expr * (if (cond[i] == "target") 2^effect else 1) * sf[i]
  }
  counts <- matrix(stats::rpois(ng * ns, lambda = as.vector(mu)), nrow = ng, ncol = ns)
  rownames(counts) <- paste0("gene_", seq_len(ng))
  colnames(counts) <- paste0("sample_", seq_len(ns))
  metadata <- data.frame(row.names = colnames(counts), condition = cond)
  list(counts = counts, metadata = metadata)
}

.golden_align <- function(eng, ref) {
  eng <- eng[order(eng$gene), ]
  ref <- ref[order(rownames(ref)), ]
  stopifnot(identical(eng$gene, rownames(ref)))
  list(lfc = max(abs(eng$log2FoldChange - ref$logFC)),
       p = max(abs(eng$pvalue - ref$pvalue) / pmax(ref$pvalue, 1e-300)),
       padj = max(abs(eng$padj - ref$adj.P.Val) / pmax(ref$adj.P.Val, 1e-300)))
}

test_that("GOLDEN limma-voom : le moteur est bit-parité avec la recette limma miroir", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")
  fx <- .golden_fixture()
  eng <- run_limma_voom_de(fx$counts, fx$metadata, "condition", "target", "ref")
  # Recette miroir EXACTE (bulk_helpers.R:702-728, branche sans pair_col, treat_lfc = 0)
  grp <- factor(fx$metadata$condition, levels = c("ref", "target"))
  design <- stats::model.matrix(~grp)
  y <- edgeR::DGEList(counts = round(fx$counts))
  y <- edgeR::calcNormFactors(y)
  v <- limma::voomWithQualityWeights(y, design)
  fit <- limma::eBayes(limma::lmFit(v, design))
  ref <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P", adjust.method = "BH")
  d <- .golden_align(eng, ref)
  expect_lt(d$lfc, 1e-8)
  expect_lt(d$p, 1e-6)
  expect_lt(d$padj, 1e-6)
})

test_that("GOLDEN edgeR : parité QL (estimateDisp → glmQLFit → glmQLFTest coef=2)", {
  skip_if_not_installed("edgeR")
  fx <- .golden_fixture()
  eng <- run_edger_de(fx$counts, fx$metadata, "condition", "target", "ref")
  # Recette miroir EXACTE (bulk_helpers.R:558-585, branche dispersion estimée)
  grp <- factor(fx$metadata$condition, levels = c("ref", "target"))
  design <- stats::model.matrix(~grp)
  y <- edgeR::DGEList(counts = round(fx$counts), group = grp)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)
  ref <- edgeR::topTags(qlf, n = Inf, adjust.method = "BH")$table
  refdf <- data.frame(gene = rownames(ref), logFC = ref$logFC, pvalue = ref$PValue,
                      adj.P.Val = ref$FDR, stringsAsFactors = FALSE)
  rownames(refdf) <- refdf$gene   # .golden_align trie par rownames côté référence
  d <- .golden_align(eng, refdf)
  expect_lt(d$lfc, 1e-8)
  expect_lt(d$p, 1e-6)
  expect_lt(d$padj, 1e-6)
})

test_that("GOLDEN DESeq2 : parité Wald (build_dds → extract_deseq2_contrast shrink=FALSE)", {
  skip_if_not_installed("DESeq2")
  # build_dds dépend de bulk_assert_raw_counts (batch_correction.R — fichier compagnon)
  source_project_file("tools/pkg_skeleton/R/batch_correction.R")
  fx <- .golden_fixture()
  dds <- build_dds(fx$counts, fx$metadata, design_formula = "~ condition", run_deseq = TRUE)
  eng <- extract_deseq2_contrast(dds, "condition", "target", "ref", shrink = FALSE)
  # Recette miroir EXACTE (bulk_helpers.R:149-159 + :228-231) — DESeq2 est déterministe
  dds_ref <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(fx$counts), colData = fx$metadata, design = ~condition)
  dds_ref <- DESeq2::DESeq(dds_ref, quiet = TRUE)
  res <- DESeq2::results(dds_ref, contrast = c("condition", "target", "ref"),
                         alpha = 0.05, pAdjustMethod = "BH")
  refdf <- data.frame(gene = rownames(res), log2FoldChange = res$log2FoldChange,
                      pvalue = res$pvalue, adj.P.Val = res$padj, stringsAsFactors = FALSE)
  rownames(refdf) <- refdf$gene   # .golden_align trie par rownames côté référence
  eng2 <- data.frame(gene = eng$gene, log2FoldChange = eng$log2FoldChange,
                     pvalue = eng$pvalue, adj.P.Val = eng$padj, stringsAsFactors = FALSE)
  d <- .golden_align(eng2, refdf)
  expect_lt(d$lfc, 1e-8)
  expect_lt(d$p, 1e-6)
  expect_lt(d$padj, 1e-6)
})

test_that("GOLDEN spatial : apply_coord_orientation — identité et involutions", {
  coords <- data.frame(id = paste0("spot_", 1:5),
                       x = c(10, 20, 30, 40, 50),
                       y = c(5, 15, 25, 35, 45), stringsAsFactors = FALSE)
  # Tous drapeaux FALSE → objet rendu tel quel
  expect_identical(apply_coord_orientation(coords), coords)
  # swap_xy est une involution
  sw <- apply_coord_orientation(coords, swap_xy = TRUE)
  expect_identical(apply_coord_orientation(sw, swap_xy = TRUE), coords)
  expect_identical(sw$x, coords$y)
  expect_identical(sw$y, coords$x)
  # flip_x avec étendue image : miroir autour de [0, w] — appliqué deux fois = identité
  fl <- apply_coord_orientation(coords, flip_x = TRUE, img_width = 100, img_height = 100)
  expect_equal(fl$x, 100 - coords$x)
  expect_identical(fl$y, coords$y)
  expect_equal(apply_coord_orientation(fl, flip_x = TRUE, img_width = 100, img_height = 100)$x,
               coords$x)
  # swap + les DEUX miroirs avec étendue image (les dims suivent le swap)
  combo <- apply_coord_orientation(coords, swap_xy = TRUE, flip_x = TRUE, flip_y = TRUE,
                                   img_width = 60, img_height = 80)
  # après swap, x=ancien y ∈ [5,45], miroir w=80 → 80-y ; y=ancien x, miroir h=60 → 60-x
  expect_equal(combo$x, 80 - coords$y)
  expect_equal(combo$y, 60 - coords$x)
})

test_that("GOLDEN spatial : miroir SANS image préserve l'emprise bbox (repli Xenium/CosMx)", {
  # Nuage décalé du bord (cas réel : coords pixels avec offset) — le miroir
  # bbox doit permuter min et max, PAS déplacer le nuage hors de son emprise.
  coords <- data.frame(id = paste0("spot_", 1:5),
                       x = c(110, 120, 130, 140, 150),
                       y = c(210, 220, 230, 240, 250), stringsAsFactors = FALSE)
  fx <- apply_coord_orientation(coords, flip_x = TRUE)
  expect_equal(range(fx$x), range(coords$x))          # emprise x PRESERVÉE
  expect_equal(sort(fx$x), sort(coords$x))            # même valeurs, miroir
  fy <- apply_coord_orientation(coords, flip_y = TRUE)
  expect_equal(range(fy$y), range(coords$y))          # emprise y PRESERVÉE
  # Combiné swap+flip sans image : les deux emprises restent celles du nuage
  both <- apply_coord_orientation(coords, swap_xy = TRUE, flip_x = TRUE, flip_y = TRUE)
  expect_equal(range(both$x), range(coords$y))
  expect_equal(range(both$y), range(coords$x))
})
