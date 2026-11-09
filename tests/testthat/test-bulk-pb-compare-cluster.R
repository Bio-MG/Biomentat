# =============================================================================
# test-bulk-pb-compare-cluster.R — PB v2 (GEMS D5) : compareCluster
# =============================================================================
# Jalon PB v2 (« comparaison d'enrichissement inter-conditions » — fiche §PB
# D5 : « compareCluster après »). Même discipline RED→GREEN que PB v1 :
# fond OBLIGATOIRE (M12), règle composites B1-C, class pathway_error,
# map_fun injecté. GO testé EN RÉEL (org.Hs.eg.db au lock, offline) sur des
# listes T-cell vs B-cell — les directions biologiques doivent se séparer ;
# KEGG est structurellement testé (requiert le réseau à l'exécution — la
# suite standard reste offline, test 429→retry mis à part).
# =============================================================================

source_project_file("tools/pkg_skeleton/R/pathway_helpers.R")

# listes réalistes (marqueurs immuno classiques) + fond étendu
.t_lists <- function() {
  tcell <- c("CD3D", "CD3E", "CD8A", "CD4", "LCK", "TRAC", "CD2", "IL7R",
             "CCR7", "SELL", "TCF7", "LEF1", "MAL", "NKG7", "GNLY", "PRF1")
  bcell <- c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "CD19", "BLK",
             "BANK1", "PAX5", "EBF1", "TCL1A", "FCRL5", "XBP1", "SDC1",
             "DERL3", "IGHG1")
  myeloid <- c("LYZ", "S100A8", "S100A9", "CD14", "FCGR3A", "MS4A7",
               "CST3", "CTSS", "FCN1", "VCAN", "S100A12", "LST1", "AIF1",
               "CD68", "MRC1", "CSF1R")
  filler <- paste0("GENEFILL", 1:60)
  list(Tcell = tcell, Bcell = bcell, Myeloid = myeloid)
}
.t_universe <- function() {
  c(unlist(.t_lists()), paste0("GENEFILL", 1:60),
    c("ACTB", "GAPDH", "B2M", "HLA-A", "HLA-B", "PTPRC", "VIM", "TUBB",
      "RPL13A", "EEF1A1"))
}

test_that("PB v2 : compareCluster GO (réel org.Hs, offline) sépare les listes", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  res <- run_pathway_compare_cluster(
    gene_lists = .t_lists(), universe_ids = .t_universe(),
    db = "GO", organism = "human", ont = "BP",
    pvalue_cutoff = 0.1, qvalue_cutoff = 0.5)
  expect_false(is.null(res$go))
  expect_true(nrow(res$go) > 0)
  expect_true(all(c("Cluster", "ID", "Description", "p.adjust") %in% colnames(res$go)))
  expect_true(all(res$go$Cluster %in% names(.t_lists())))
  expect_true(all(res$go$db == "GO"))
  # séparation biologique : les T cell genes présentent au moins un terme
  # « T cell » (falsifiable — un moteur muet ne le trouverait pas)
  t_terms <- res$go$Description[res$go$Cluster == "Tcell"]
  expect_true(any(grepl("T cell", t_terms, ignore.case = TRUE)))
  # meta honnête : 3 listes gardées, comptes > seuil, fond > listes
  expect_equal(sum(vapply(res$meta$per_list, function(x) isTRUE(x$kept), logical(1))), 3L)
  # MESURE : les 60 GENEFILL synthétiques ne mappent pas en Entrez — le fond
  # effectif = les 48 marqueurs + 10 housekeeping ≈ 58 (honnêteté M12)
  expect_gte(res$meta$universe_n_entrez, 50L)
})

test_that("PB v2 : liste trop petite DÉPOSÉE avec note, jamais testée muette", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  lists <- .t_lists()
  lists$Tiny <- c("CD3D", "CD3E", "CD8A")   # 3 symboles réels -> < 5 Entrez
  res <- run_pathway_compare_cluster(
    gene_lists = lists, universe_ids = .t_universe(),
    db = "GO", organism = "human", pvalue_cutoff = 0.1, qvalue_cutoff = 0.5)
  expect_true("Tiny" %in% res$meta$dropped)
  expect_false(is.null(res$go))
  expect_false(any(res$go$Cluster == "Tiny"))
  expect_true(any(grepl("déposée", res$meta$notes)))
})

test_that("PB v2 : refus classés — <2 listes, liste non nommée, fond vide, seuils invalides", {
  skip_if_not_installed("clusterProfiler")
  e1 <- tryCatch(run_pathway_compare_cluster(
    gene_lists = list(A = c("CD3D", "CD8A")), universe_ids = .t_universe()),
    error = function(e) e)
  expect_s3_class(e1, "pathway_error")
  e2 <- tryCatch(run_pathway_compare_cluster(
    gene_lists = c("CD3D", "CD8A"), universe_ids = .t_universe()),
    error = function(e) e)
  expect_s3_class(e2, "pathway_error")
  e3 <- tryCatch(run_pathway_compare_cluster(
    gene_lists = .t_lists(), universe_ids = character(0)),
    error = function(e) e)
  expect_s3_class(e3, "pathway_error")
  e4 <- tryCatch(run_pathway_compare_cluster(
    gene_lists = .t_lists(), universe_ids = .t_universe(), pvalue_cutoff = 1.5),
    error = function(e) e)
  expect_s3_class(e4, "pathway_error")
  # toutes les listes déposées -> < 2 valides -> refus avec comptes
  e5 <- tryCatch(run_pathway_compare_cluster(
    gene_lists = list(A = c("CD3D"), B = c("CD8A"), C = c("LYZ")),
    universe_ids = .t_universe()),
    error = function(e) e)
  expect_s3_class(e5, "pathway_error")
  expect_match(conditionMessage(e5), "Moins de 2 listes valides")
})

test_that("PB v2 : règle composites B1-C — les ENSG..._SYMBOL sortent du compte, jamais du silence", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  lists <- .t_lists()
  # les 3 premiers gènes de chaque liste passent en composite : le strip
  # B1-C garde la partie ENSG (sub("[_.].*$","",...)) qui NE mappe PAS en
  # SYMBOL→Entrez — ils sortent HONNÊTEMENT des comptes (déclarés 16 ->
  # Entrez 13), la liste reste valide (13 >= 5) et le calcul continue.
  lists_mix <- lapply(lists, function(x) {
    c(paste0("ENSG00000", seq_along(head(x, 3)), "_", head(x, 3)), tail(x, -3))
  })
  # seuils assouplis : 13 gènes/liste (vs 16) → les termes perdent en
  # significativité, le couple 0.1/0.5 du test principal ne passe plus
  res <- run_pathway_compare_cluster(lists_mix, .t_universe(), db = "GO",
                                     organism = "human",
                                     pvalue_cutoff = 0.5, qvalue_cutoff = 1)
  expect_false(is.null(res$go))
  for (nm in names(lists_mix)) {
    pl <- res$meta$per_list[[nm]]
    expect_equal(pl$declared, 16L)
    expect_equal(pl$entrez, 13L)
    expect_true(isTRUE(pl$kept))
  }
})
