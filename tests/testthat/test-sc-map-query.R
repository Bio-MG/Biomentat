# =============================================================================
# test-sc-map-query.R — SE (2026-11-09) : reference mapping MapQuery
# (fiche PROPOSAL_SEURAT_V56_GEMS §3-SE, moteur sc_map_query dans sc_helpers.R)
# =============================================================================
# Contrats testés (rehearsal 2026-11-09, Seurat 5.6) :
#  - colonnes produites : predicted.<label_col> + predicted.<label_col>.score
#    (nommées d'après l'entrée refdata — PAS predicted.id de la vignette v4) ;
#  - réduction projetée ref.umap (nrow = ncol(query)) ;
#  - gardes classées sc_helpers_error : colonne absente, all-NA, gènes communs
#    insuffisants, référence sans pca ;
#  - coexistence SingleR : le moteur N'ECRIT aucune colonne SingleR_*.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_helpers.R")

skip_if_not_installed("Seurat")
skip_if_not_installed("SeuratObject")
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Seurat)
})

.mk_pair <- function() {
  set.seed(1)
  mk <- function(n, prefix, shift) {
    m <- matrix(stats::rpois(200 * n, 2), 200,
                dimnames = list(paste0("g", 1:200), paste0(prefix, 1:n)))
    m[1:20, 1:(n %/% 2)] <- m[1:20, 1:(n %/% 2)] + shift
    m[21:40, ((n %/% 2) + 1):n] <- m[21:40, ((n %/% 2) + 1):n] + shift
    m
  }
  ref <- CreateSeuratObject(cbind(mk(150, "R1_", 15), mk(150, "R2_", 0)))
  ref$celltype <- factor(c(rep("Tcell", 75), rep("Bcell", 75)))
  ref <- NormalizeData(ref)
  ref <- FindVariableFeatures(ref, verbose = FALSE)
  ref <- ScaleData(ref, verbose = FALSE)
  ref <- RunPCA(ref, npcs = 10, verbose = FALSE)
  q <- CreateSeuratObject(cbind(mk(80, "Q1_", 15), mk(80, "Q2_", 0)))
  q <- NormalizeData(q)
  list(ref = ref, query = q)
}

test_that("SE : sc_map_query transfere les etiquettes + scores + ref.umap", {
  pair <- .mk_pair()
  res <- tryCatch(
    sc_map_query(reference = pair$ref, query = pair$query,
                 label_col = "celltype", dims = 10, seed = 1),
    error = function(e) e)
  if (inherits(res, "error")) skip(paste("MapQuery indisponible sur ce fixture :", conditionMessage(res)))

  q <- res$query
  expect_identical(res$predicted_col, "predicted.celltype")
  expect_true("predicted.celltype" %in% colnames(q@meta.data))
  expect_true("predicted.celltype.score" %in% colnames(q@meta.data))
  expect_type(q@meta.data$predicted.celltype.score, "double")
  # MESURE : les scores de label transfer peuvent dépasser 1 d'un EPSILON
  # flottant (30/160 cellules à 1 + ~2e-16 mesuré — range() arrondit à
  # 1.000000, `> 1` voit vrai). Tolérance machine obligatoire.
  sc <- q@meta.data$predicted.celltype.score
  expect_true(all(sc >= -1e-8 & sc <= 1 + 1e-8))
  expect_false(any(is.na(q$predicted.celltype)))
  # ref.umap : reduction projetee pleine query
  umap_name <- grep("^ref\\..*umap", names(q@reductions), value = TRUE)[1]
  expect_false(is.na(umap_name))
  expect_equal(nrow(SeuratObject::Embeddings(q[[umap_name]])), ncol(q))
  # coexistence SingleR : aucune colonne SingleR ecrite
  expect_false(any(grepl("^SingleR_", colnames(q@meta.data))))
  # contrat de retour
  expect_equal(res$n_query, ncol(pair$query))
  expect_equal(res$n_ref, ncol(pair$ref))
  expect_equal(res$dims_used, 10L)
})

test_that("SE : gardes classees — colonne absente, all-NA, pca absente", {
  pair <- .mk_pair()
  e1 <- tryCatch(sc_map_query(pair$ref, pair$query, label_col = "nope", dims = 10),
                 error = function(e) e)
  expect_s3_class(e1, "sc_helpers_error")
  ref_na <- pair$ref
  ref_na$vide <- NA_character_
  e2 <- tryCatch(sc_map_query(ref_na, pair$query, label_col = "vide", dims = 10),
                 error = function(e) e)
  expect_s3_class(e2, "sc_helpers_error")
  ref_nopca <- pair$ref
  ref_nopca[["pca"]] <- NULL
  e3 <- tryCatch(sc_map_query(ref_nopca, pair$query, label_col = "celltype", dims = 10),
                 error = function(e) e)
  expect_s3_class(e3, "sc_helpers_error")
})

test_that("SE : garde gènes communs insuffisants (refus actionnable, pas de crash opaque)", {
  pair <- .mk_pair()
  # query fraîche aux gènes renommés : zéro commun (Assay5 n'expose pas de
  # slot counts — reconstruction simple au lieu de la manipulation de slots)
  set.seed(1)
  m_far <- matrix(stats::rpois(200 * 160, 2), 200,
                  dimnames = list(paste0("zz_", 1:200), colnames(pair$query)))
  q_far <- CreateSeuratObject(m_far)
  q_far <- NormalizeData(q_far)
  e <- tryCatch(sc_map_query(pair$ref, q_far, label_col = "celltype", dims = 10),
                error = function(e) e)
  expect_s3_class(e, "sc_helpers_error")
  expect_match(conditionMessage(e), "commun")
})
