# =============================================================================
# test-sc-cell-cycle.R — SF (2026-11-09) : cell cycle scoring + régression
# (fiche PROPOSAL_SEURAT_V56_GEMS §3-SF, moteur sc_cell_cycle_scores)
# =============================================================================
# Contrats : scoring en metadata (S.Score/G2M.Score/Phase + CC.Difference),
# régression ScaleData optionnelle via smart_scale_data(vars.to.regress=),
# refus classés (mouse sans mapping orthologue, IDs non-remappés), défaut
# pipeline OFF (zéro changement).
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_helpers.R")
source_project_file("tools/pkg_skeleton/R/sc_bpcells.R")
source_project_file("tools/pkg_skeleton/R/sc_pipeline.R")

skip_if_not_installed("Seurat")
skip_if_not_installed("SeuratObject")
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Seurat)
})

# fixture aux VRAIS symboles de cycle : 15 gènes S + 12 gènes G2M réels +
# remplissage. cc.genes vit dans le namespace Seurat (lazy data).
.cc_fixture <- function(n_cells = 80) {
  cc <- get("cc.genes", envir = asNamespace("Seurat"))
  genes <- c(head(cc$s.genes, 15), head(cc$g2m.genes, 12),
             paste0("filler", 1:600))
  # expression hétérogène des fillers (MESURE : des fillers à lambda uniforme
  # laissent les gènes S/G2M SEULS dans leur bin d'expression et le ctrl
  # sampling d'AddModuleScore échoue — « sample larger than population »)
  lambdas <- c(rep(2, 27), 0.5 + (seq_len(600) %% 40) * 0.25)
  counts <- matrix(stats::rpois(length(genes) * n_cells,
                                lambda = rep(lambdas, times = n_cells)),
                   length(genes),
                   dimnames = list(genes, paste0("cell", 1:n_cells)))
  # cycle réel : 1/3 des cellules en S (gènes S hauts), 1/3 en G2M
  s_cells <- 1:(n_cells %/% 3)
  g2m_cells <- (n_cells %/% 3 + 1):(2 * n_cells %/% 3)
  counts[1:15, s_cells] <- counts[1:15, s_cells] + 8
  counts[16:27, g2m_cells] <- counts[16:27, g2m_cells] + 8
  obj <- CreateSeuratObject(counts)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj
}

test_that("SF : scoring en metadata — S.Score, G2M.Score, Phase, CC.Difference", {
  obj <- .cc_fixture()
  res <- sc_cell_cycle_scores(obj, organism = "human", ctrl = 10)
  expect_s4_class(res$obj, "Seurat")
  expect_true(all(c("S.Score", "G2M.Score", "Phase", "CC.Difference") %in%
                  colnames(res$obj@meta.data)))
  expect_true(all(is.numeric(res$obj$S.Score), is.numeric(res$obj$G2M.Score)))
  expect_identical(as.numeric(res$obj$CC.Difference),
                   as.numeric(res$obj$S.Score - res$obj$G2M.Score))
  expect_gt(res$n_s, 10)
  expect_gt(res$n_g2m, 10)
  # les cellules "S" du fixture portent des scores S majoritairement hauts :
  expect_gt(mean(res$obj$S.Score[1:(ncol(obj) %/% 3)]),
            mean(res$obj$S.Score[((2 * ncol(obj) %/% 3) + 1):ncol(obj)]))
})

test_that("SF : refus classés — mouse (pas d'orthologue) et IDs non-remappés", {
  obj <- .cc_fixture()
  e1 <- tryCatch(sc_cell_cycle_scores(obj, organism = "mouse"),
                 error = function(e) e)
  expect_s3_class(e1, "sc_helpers_error")
  expect_match(conditionMessage(e1), "orthologue")
  # objet aux gènes renommés (IDs non-humains) : overlap insuffisant
  obj_bad <- obj
  counts <- GetAssayData(obj_bad, layer = "counts")
  rownames(counts) <- paste0("ENSG", seq_len(nrow(counts)))
  obj_bad <- CreateSeuratObject(counts)
  obj_bad <- NormalizeData(obj_bad, verbose = FALSE)
  e2 <- tryCatch(sc_cell_cycle_scores(obj_bad, organism = "human"),
                 error = function(e) e)
  expect_s3_class(e2, "sc_helpers_error")
  expect_match(conditionMessage(e2), "remapp")
})

test_that("SF : régression optionnelle — smart_scale_data(vars.to.regress) pose la régression", {
  obj <- .cc_fixture(n_cells = 60)
  res <- sc_cell_cycle_scores(obj, organism = "human", ctrl = 10)
  obj <- res$obj
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  # voie différence (défaut fiche) : CC.Difference régressé, pas les scores bruts
  out <- smart_scale_data(obj, vars.to.regress = "CC.Difference")
  expect_true(nrow(GetAssayData(out, layer = "scale.data")) > 0)
  # voie complète
  out2 <- smart_scale_data(obj, vars.to.regress = c("S.Score", "G2M.Score"))
  expect_true(nrow(GetAssayData(out2, layer = "scale.data")) > 0)
  # défaut inchangé : sans vars, le comportement historique est intact
  out3 <- smart_scale_data(obj)
  expect_true(nrow(GetAssayData(out3, layer = "scale.data")) > 0)
})

test_that("SF : défaut pipeline OFF — sc_ap_normalize_params ferme le jeu de clés", {
  d <- sc_ap_normalize_params()
  expect_false(isTRUE(d$sc_ap_cellcycle))
  expect_identical(d$sc_ap_cellcycle_regress, "difference")
  # override possible :
  d2 <- sc_ap_normalize_params(list(sc_ap_cellcycle = TRUE,
                                    sc_ap_cellcycle_regress = "complete"))
  expect_true(isTRUE(d2$sc_ap_cellcycle))
  expect_identical(d2$sc_ap_cellcycle_regress, "complete")
})
