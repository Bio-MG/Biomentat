# =============================================================================
# test-sc-integration-native.R — voie native v5 IntegrateLayers (fiche
# PROPOSAL_SEURAT_V56_GEMS_2026-11-05 §3-SC, jalon SC, 2026-11-09).
# =============================================================================
# Le sélecteur d'intégration offre désormais rpca/cca (voie native v5 :
# split par batch → normalisation/PCA par couche → IntegrateLayers →
# JoinLayers) à côté d'Harmony (comportement historique inchangé) et de
# `none`. INVARIANTS testés ici sur le cœur réel (run_sc_auto_pipeline_core,
# fixture déterministe 2 groupes = 2 batchs) :
#   1. la voie rpca produit une réduction `integrated.rpca` ET une réduction
#      `pca` intacte (garde de comparabilité) ;
#   2. les clusters non intégrés survivent sous `unintegrated_clusters` et les
#      intégrés sous `<methode>_clusters` — côte à côte (garde fiche §3-SC) ;
#   3. `seurat_clusters` = clusters INTÉGRÉS (ordre mesuré : FindClusters
#      écrase seurat_clusters au dernier appel — l'intégration doit être la
#      dernière) ;
#   4. les couches sont re-jointes (JoinLayers) — l'aval (marqueurs) lit un
#      seul layer `data` ;
#   5. les skips sont HONNÊTES et par message : SCT non activée au jalon 1,
#      batch absent, et `mnn` MORT SUR MESURE (FastMNNIntegration absent de
#      Seurat 5.6.0 — getNamespaceExports mesuré) ;
#   6. la voie `none` est une régression stricte : aucun artefact d'intégration.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/sc_helpers.R")   # resolve_sketch_preset, robust_find_clusters
source_project_file("tools/pkg_skeleton/R/sc_bpcells.R")   # smart_scale_data
source_project_file("tools/pkg_skeleton/R/sc_trajectory.R")  # calculate_pseudotime
source_project_file("tools/pkg_skeleton/R/sc_pipeline.R")

# Fixture 2 groupes × 2 batchs (motif test-sc-auto-pipeline.R : blocs
# marqueurs séparables + 3 gènes MT), orig.ident = la variable de batch.
.make_native_fixture <- function(n_per_group = 80L, n_genes = 120L, seed = 42) {
  set.seed(seed)
  n_cells <- 2L * n_per_group
  mat <- matrix(stats::rpois(n_genes * n_cells, lambda = 1),
                nrow = n_genes, dimnames = list(paste0("G", seq_len(n_genes)), NULL))
  mat[1:15, seq_len(n_per_group)] <-
    mat[1:15, seq_len(n_per_group)] + stats::rpois(15L * n_per_group, 8)
  mat[16:30, (n_per_group + 1L):n_cells] <-
    mat[16:30, (n_per_group + 1L):n_cells] + stats::rpois(15L * n_per_group, 8)
  mt <- matrix(stats::rpois(3L * n_cells, lambda = 1), nrow = 3L,
               dimnames = list(c("MT-A", "MT-B", "MT-C"), NULL))
  mat <- rbind(mat, mt)
  obj <- SeuratObject::CreateSeuratObject(counts = mat)
  obj$orig.ident <- factor(rep(c("s1", "s2"), each = n_per_group))
  obj
}

.ap_base_params <- function() {
  list(
    run_qc             = FALSE,
    sc_ap_norm         = "log",
    sc_ap_pca_dim      = 10,
    sc_ap_res          = 0.5,
    sc_ap_cluster_algo = "1",
    sc_ap_markers      = FALSE,
    sc_ap_compute_umap = TRUE,
    sc_ap_batch_var    = "orig.ident"
  )
}

test_that("voie native rpca : réduction intégrée + garde side-by-side + JoinLayers", {
  testthat::skip_if_not_installed("Seurat")
  library(Seurat)   # génériques non qualifiés du cœur (NormalizeData, ScaleData, ...)
  testthat::skip_if(utils::packageVersion("Seurat") < "5.0.0")

  res <- run_sc_auto_pipeline_core(
    .make_native_fixture(n_per_group = 150L),  # batchs > k.weight rpca par défaut (100)
    params = modifyList(.ap_base_params(), list(sc_ap_integration = "rpca"))
  )

  # Pas de marqueur d'échec ; la voie native s'est appliquée.
  # NB : res$log est un VECTEUR de lignes (cœur) — any() sinon expect_false
  # compare une longueur (mesuré 2026-11-09 : échec de longueur, pas de valeur).
  testthat::expect_false(any(grepl("\u274c", res$log)), info = paste(res$log, collapse = "\n"))
  testthat::expect_match(paste(res$log, collapse = "\n"),
                         "Intégration native rpca appliquée", fixed = TRUE)

  obj <- res$obj
  # 1. réduction intégrée + PCA brute intacte
  testthat::expect_true("integrated.rpca" %in% names(obj@reductions))
  testthat::expect_true("pca" %in% names(obj@reductions))

  # 2. garde side-by-side (fiche §3-SC)
  testthat::expect_true("unintegrated_clusters" %in% colnames(obj@meta.data))
  testthat::expect_true("rpca_clusters" %in% colnames(obj@meta.data))

  # 3. seurat_clusters = clusters INTÉGRÉS (miroir identique à rpca_clusters)
  testthat::expect_identical(as.integer(obj$seurat_clusters),
                             as.integer(obj$rpca_clusters))

  # 4. couches re-jointes : un seul layer `data` (pas de data.s1/data.s2)
  layers <- SeuratObject::Layers(obj[["RNA"]])
  testthat::expect_false(any(grepl("^data\\.", layers)),
               info = paste("layers après JoinLayers :", paste(layers, collapse = ",")))

  # 5. le clustering suit l'espace intégré : UMAP sur integrated.rpca
  testthat::expect_true("umap" %in% names(obj@reductions))
  testthat::expect_equal(res$metrics$integration, "rpca")
})

test_that("voie native cca : s'applique et produit integrated.cca", {
  testthat::skip_if_not_installed("Seurat")
  library(Seurat)   # génériques non qualifiés du cœur (NormalizeData, ScaleData, ...)
  testthat::skip_if(utils::packageVersion("Seurat") < "5.0.0")

  res <- run_sc_auto_pipeline_core(
    .make_native_fixture(n_per_group = 150L, seed = 43),
    params = modifyList(.ap_base_params(), list(sc_ap_integration = "cca"))
  )
  testthat::expect_false(any(grepl("\u274c", res$log)), info = paste(res$log, collapse = "\n"))
  testthat::expect_match(paste(res$log, collapse = "\n"),
                         "Intégration native cca appliquée", fixed = TRUE)
  testthat::expect_true("integrated.cca" %in% names(res$obj@reductions))
  testthat::expect_true("cca_clusters" %in% colnames(res$obj@meta.data))
  testthat::expect_true("unintegrated_clusters" %in% colnames(res$obj@meta.data))
})

test_that("skips honnêtes : SCT non activée, batch absent, mnn mort sur mesure", {
  testthat::skip_if_not_installed("Seurat")
  library(Seurat)   # génériques non qualifiés du cœur (NormalizeData, ScaleData, ...)
  testthat::skip_if(utils::packageVersion("Seurat") < "5.0.0")

  # SCT : documentée fiche §3-SC mais NON activée au jalon 1 — skip affiché,
  # jamais silencieux, et AUCUN artefact d'intégration.
  res_sct <- run_sc_auto_pipeline_core(
    .make_native_fixture(n_per_group = 150L),  # batchs > k.weight rpca par défaut (100)
    params = modifyList(.ap_base_params(),
                        list(sc_ap_integration = "rpca", sc_ap_norm = "sct"))
  )
  testthat::expect_match(paste(res_sct$log, collapse = "\n"),
                         "Intégration native ignorée : compat SCT", fixed = TRUE)
  testthat::expect_false("integrated.rpca" %in% names(res_sct$obj@reductions))
  testthat::expect_false("unintegrated_clusters" %in% colnames(res_sct$obj@meta.data))

  # Variable de batch absente : skip affiché.
  res_nobatch <- run_sc_auto_pipeline_core(
    .make_native_fixture(n_per_group = 150L),  # batchs > k.weight rpca par défaut (100)
    params = modifyList(.ap_base_params(),
                        list(sc_ap_integration = "rpca", sc_ap_batch_var = "n_existe_pas"))
  )
  testthat::expect_match(paste(res_nobatch$log, collapse = "\n"),
                         "variable de batch 'n_existe_pas' absente", fixed = TRUE)
  testthat::expect_false("integrated.rpca" %in% names(res_nobatch$obj@reductions))

  # mnn : mort sur mesure (FastMNNIntegration absent de Seurat 5.6.0) — le
  # refus est AFFICHÉ, jamais un repli silencieux.
  res_mnn <- run_sc_auto_pipeline_core(
    .make_native_fixture(n_per_group = 150L),  # batchs > k.weight rpca par défaut (100)
    params = modifyList(.ap_base_params(), list(sc_ap_integration = "mnn"))
  )
  testthat::expect_match(paste(res_mnn$log, collapse = "\n"),
                         "'mnn' indisponible", fixed = TRUE)
  testthat::expect_false(any(grepl("^integrated\\.", names(res_mnn$obj@reductions))))
})

test_that("régression voie none : aucun artefact d'intégration", {
  testthat::skip_if_not_installed("Seurat")
  library(Seurat)   # génériques non qualifiés du cœur (NormalizeData, ScaleData, ...)

  res <- run_sc_auto_pipeline_core(
    .make_native_fixture(n_per_group = 150L),  # batchs > k.weight rpca par défaut (100)
    params = .ap_base_params()
  )
  testthat::expect_false(any(grepl("^integrated\\.", names(res$obj@reductions))))
  testthat::expect_false("unintegrated_clusters" %in% colnames(res$obj@meta.data))
  testthat::expect_true("seurat_clusters" %in% colnames(res$obj@meta.data))
  testthat::expect_equal(res$metrics$integration, "none")
})
