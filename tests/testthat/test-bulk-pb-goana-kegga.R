# =============================================================================
# test-bulk-pb-goana-kegga.R — PB (GEMS D5) : ORA directionnelle limma
# =============================================================================
# Jalon 1 du paquet PB (« enrichissement — le 3ᵉ idiome ») : goana/kegga +
# barcodeplot. L'ORA clusterProfiler (run_pathway_enrichment) travaille sur la
# LISTE des gènes significatifs ; goana/kegga travaillent sur la TABLE DE
# complète et testent la DIRECTION (Up/Down) par voie — le 3ᵉ idiome, d'après
# la fiche GEMS §PB (GO opérateur A1, séquence A→C→B→D ; D5 : goana/kegga
# d'abord, compareCluster après ; pilier 4 fgsea/fry/camera = PACKS, jamais
# double-owner).
#
# Discipline RED→GREEN : ce fichier a été écrit AVANT le moteur (run #B :
# « could not find function run_pathway_direction » = RED confirmé), puis le
# moteur pathway_helpers.R a été écrit pour le faire passer.
#
# Structure testée :
#   - refus classés pathway_error (colonnes manquantes, table vide, organisme
#     inconnu, db inconnu, univers vide / zéro-mappé, zéro DE au seuil) ;
#   - nettoyage des identifiants composites (règle B1-C, même règle que
#     run_pathway_enrichment) ;
#   - déduplication SYMBOLE→ENTREZ (deux symboles → un même Entrez : le plus
#     significatif garde la ligne — sinon data.frame rownames dupliqués) ;
#   - intégration org.Hs.eg.db réelle : tables goana/kegga non vides, colonnes
#     de direction présentes, meta honnête (comptes, note de fond) ;
#   - plot_pathway_barcode : refus classé si aucun gène de la voie dans le
#     classement, tracé OK sinon (graphics off).
# =============================================================================

source_project_file("tools/pkg_skeleton/R/pathway_helpers.R")

.pbd_fixture_de <- function() {
  data.frame(
    gene = c("TP53", "IL6", "TNF", "CD4", "CD8A", "FOXP3", "GZMB", "IFNG",
             "STAT1", "CXCL8"),
    log2FoldChange = c(-2.1, 3.4, 2.8, 1.2, -1.5, 2.2, 3.0, 2.5, 1.9, 3.6),
    pvalue = c(1e-6, 1e-8, 1e-5, 0.03, 0.04, 2e-6, 3e-7, 5e-6, 0.01, 4e-9),
    padj = c(0.001, 0.0002, 0.002, 0.2, 0.25, 0.001, 0.0004, 0.003, 0.08, 0.0001),
    stringsAsFactors = FALSE
  )
}

# Fond synthétique : le fixture + un grand pool de symboles factices qui ne
# mappent PAS (vérifie que la note de fond honnête reflète le non-mappé).
.pbd_fixture_universe <- function() {
  c(.pbd_fixture_de()$gene, paste0("ZZ_FAKE_", seq_len(200)))
}

# map_fun injecté : PURE, sans OrgDb — TP53/IL6/… → Entrez ; le reste ne
# mappe pas (0 ligne retournée par bitr simulé).
.pbd_map_fun <- function(keys) {
  known <- c(TP53 = "7157", IL6 = "3569", TNF = "7124", CD4 = "920",
             CD8A = "925", FOXP3 = "50951", GZMB = "1104", IFNG = "3458",
             STAT1 = "6772", CXCL8 = "3576")
  hit <- intersect(keys, names(known))
  if (!length(hit)) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0),
                      stringsAsFactors = FALSE))
  }
  data.frame(SYMBOL = hit, ENTREZID = unname(known[hit]),
             stringsAsFactors = FALSE)
}

# ── Refus classés (PURES — aucun paquet Bioconductor requis) ────────────────

test_that("PB : refus classés pathway_error sur entrées invalides", {
  de <- .pbd_fixture_de()
  uni <- .pbd_fixture_universe()

  expect_error(
    run_pathway_direction(de[, c("gene", "log2FoldChange")], uni),
    class = "pathway_error")           # colonne pvalue manquante
  expect_error(
    run_pathway_direction(de[0, ], uni),
    class = "pathway_error")           # table DE vide
  expect_error(
    run_pathway_direction(de, character(0)),
    class = "pathway_error")           # univers vide
  expect_error(
    run_pathway_direction(de, uni, db = "REACTOME"),
    class = "pathway_error")           # db inconnu (v1 : GO/KEGG seulement)
  expect_error(
    run_pathway_direction(de, uni, organism = "rat"),
    class = "pathway_error")           # organisme hors portée v1 (human/mouse)
})

test_that("PB : zéro gène DE au seuil = refus honnête, jamais une table vide muette", {
  de <- .pbd_fixture_de()
  de$pvalue <- 0.9   # rien de significatif au seuil par défaut (FDR 0.05)
  expect_error(
    run_pathway_direction(de, .pbd_fixture_universe(),
                          map_fun = .pbd_map_fun),
    class = "pathway_error")
})

test_that("PB : identifiants composites B1-C acceptés et nettoyés", {
  de <- .pbd_fixture_de()
  de$gene[1] <- "ENSG00000141510_TP53"   # composite ENSG…_SYMBOL
  # Ne doit PAS refuser : le nettoyage composite précède le mapping. Le strip
  # B1-C garde la partie ENSG (sub("[_.].*$","",...)) qui ne mappe PAS en
  # SYMBOL→Entrez : ce gène sort HONNÊTEMENT du compte (9/10), jamais au
  # silence dans un entre-deux.
  expect_no_error(
    res <- run_pathway_direction(de, .pbd_fixture_universe(),
                                 db = "GO", map_fun = .pbd_map_fun))
  expect_identical(res$meta$de_n_entrez, 9L)
})

# ── Intégration org.Hs.eg.db réelle (le moteur complet) ─────────────────────

test_that("PB : goana rend une table directionnelle avec meta honnête", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("limma")

  de <- .pbd_fixture_de()
  uni <- c(de$gene, paste0("GENE_NOMAP_", seq_len(80)))
  res <- run_pathway_direction(de, uni, db = "GO", map_fun = .pbd_map_fun)

  expect_false(is.null(res$go))
  expect_true(is.data.frame(res$go))
  expect_true(nrow(res$go) > 0L)
  # Colonnes de direction — le cœur du 3ᵉ idiome : effectifs Up/Down par voie
  # (API limma moderne = liste par direction) + p-valeurs + FDR appliqué côté
  # moteur (BH par direction, même politique que le module).
  expect_true(all(c("Term", "N", "Up", "Down", "P.Up", "P.Down",
                    "FDR.Up", "FDR.Down") %in% names(res$go)))
  # Meta honnête : comptes réels, note de fond conservée (M12).
  expect_identical(res$meta$de_n_entrez, 10L)
  expect_identical(res$meta$universe_n_entrez, 10L)  # les ZZ_FAKE_ ne mappent pas
  expect_match(res$meta$organism, "human")
})

test_that("PB : kegga rend une table directionnelle (db=KEGG)", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("limma")

  res <- run_pathway_direction(.pbd_fixture_de(), .pbd_fixture_universe(),
                               db = "KEGG", map_fun = .pbd_map_fun)
  expect_false(is.null(res$kegg))
  expect_true(is.data.frame(res$kegg))
  expect_true(nrow(res$kegg) > 0L)
  expect_true(all(c("Up", "Down", "P.Up", "P.Down", "FDR.Up", "FDR.Down")
                  %in% names(res$kegg)))
  expect_null(res$go)   # db=KEGG seul : go absent (pas de calcul fantôme)
})

test_that("PB : note de fond M12 estampillée quand des gènes ne mappent pas", {
  skip_if_not_installed("org.Hs.eg.db")
  res <- run_pathway_direction(.pbd_fixture_de(), .pbd_fixture_universe(),
                               db = "GO", map_fun = .pbd_map_fun)
  expect_false(is.null(res$meta$universe_note))
  expect_match(res$meta$universe_note, "fond|universe|repli", ignore.case = TRUE)
})

# ── barcodeplot ─────────────────────────────────────────────────────────────

test_that("PB : barcodeplot refuse honnêtement sans gène de la voie, trace sinon", {
  skip_if_not_installed("limma")

  stats <- setNames(.pbd_fixture_de()$log2FoldChange, .pbd_fixture_de()$gene)
  expect_error(
    plot_pathway_barcode(stats, character(0)),
    class = "pathway_error")                     # index vide
  expect_error(
    plot_pathway_barcode(stats, c("ABSENTE_1", "ABSENTE_2")),
    class = "pathway_error")                     # aucun gène de la voie classé

  png(tempfile(fileext = ".png"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(
    plot_pathway_barcode(stats, c("TP53", "IL6", "TNF")))
})
