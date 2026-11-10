# =============================================================================
# test-packs-convergence.R — PK-3 (pilier 3 « couche de convergence », fiche
# PROPOSAL_PACKS_ECOSYSTEM_2026-11-09.md §3 jalon PK-3 ; note de conception
# PACKS_PK3_CONVERGENCE_DESIGN_2026-11-10.md — D5 GO opérateur 2026-11-10 =
# OPTION A : BH globale sur l'union des tests + effondrement de redondance
# enrichplot::pairwise_termsim, 1 cluster = 1 vote).
#
# Créé AVANT l'implémentation (règle maison : test d'abord, rouge d'abord) —
# ce fichier doit être ROUGE tant que la section « PACKS (PK-3) » n'existe pas
# dans pathway_helpers.R (packs_split_libraries / packs_converge /
# build_convergence_dt).
#
# Contrat PK-3 (note de conception §3) : PK-3 consomme les sorties PK-2 SANS
# re-définir ni re-runner les moteurs (pilier 4 = GEMS PB, jamais double-owner,
# fiche §4/D1) — la BH globale est recalculée à partir des p BRUTES portées par
# le contrat PK-2 (colonne `pvalue`) et l'effondrement consomme l'attribut
# `enrich_obj` (contrat STAT-S2). Zéro dépendance nouvelle : enrichplot /
# clusterProfiler / igraph sont déjà au lock (igraph déjà consommé par
# bulk_network.R et mod_sc_corr.R).
#
# Preuves exigées par l'ask (vérificateur) :
#   - la p.adjust GLOBALE de packs_converge est identique (bit à bit) à une
#     recomputation manuelle stats::p.adjust(, "BH") sur l'union des tests ;
#   - le cluster termsim : 1 cluster = 1 vote — le score compte les
#     LIBRAIRIES DISTINCTES significatives sous le seuil global.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")   # strip_composite_gene_ids (dépendance run_packs_enrichment — motif helper-sccoda-fixtures.R:31)
source_project_file("tools/pkg_skeleton/R/datatable.R")    # ts_datatable (dépendance build_convergence_dt)
source_project_file("tools/pkg_skeleton/R/pathway_helpers.R")

# -- fixtures -----------------------------------------------------------------
# 2 librairies à chevauchement CONTRÔLÉ :
#   libA : « Alpha pathway » (G1-30) + « Alpha pathway bis » (MÊMES gènes —
#          cluster 3 termes, 2 librairies seulement) + « Beta pathway »
#          (G31-60) + « Gamma pathway » (G61-90, aucun hit → absent du run)
#   libB : « Alpha pathway » (G1-30, MÊMES gènes que libA — JC = 1, ET
#          description IDENTIQUE à libA — la matrice termsim enrichplot est
#          étiquetée par Description : le cas dupliqué doit rester exact) +
#          « Delta pathway » (G91-120, disjoint)

.tmp_conv_dir <- function() {
  d <- file.path(tempdir(), sprintf("packs_pk3_%d", as.integer(Sys.time()) + sample(1e6, 1)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.conv_gmt_lines <- function(term, desc, genes) paste(c(term, desc, genes), collapse = "\t")

.conv_pack_registry <- function(d) {
  writeLines(c(
    .conv_gmt_lines("Alpha pathway", "Signalisation alpha", paste0("G", 1:30)),
    .conv_gmt_lines("Alpha pathway bis", "Signalisation alpha bis", paste0("G", 1:30)),
    .conv_gmt_lines("Beta pathway", "Signalisation beta", paste0("G", 31:60)),
    .conv_gmt_lines("Gamma pathway", "Signalisation gamma", paste0("G", 61:90))
  ), file.path(d, "libA.gmt"))
  writeLines(c(
    .conv_gmt_lines("Alpha pathway", "Signalisation alpha", paste0("G", 1:30)),
    .conv_gmt_lines("Delta pathway", "Signalisation delta", paste0("G", 91:120))
  ), file.path(d, "libB.gmt"))
  structure(
    list(id = "pk_conv", type = "libraries", name = "Conv fixture",
         version = "1.0", attribution = "test", source = NULL,
         downloaded = NULL, license = NULL,
         files = c("libA.gmt", "libB.gmt"), palettes = NULL,
         dir = d, manifest = file.path(d, "manifest.json")),
    class = "biomentat_pack"
  )
}

# Genes à tester : Alpha plein (30/30), Beta 5/30, Delta 5/30, Gamma 0/30.
# pval/qval cutoffs à 1 : la table PK-2 porte TOUS les termes testés (l'union
# est complète — cf. test d'équivalence BH PK-2 / BH globale ci-dessous).
.conv_genes <- c(paste0("G", 1:30), paste0("G", 31:35), paste0("G", 91:95))
.conv_universe <- paste0("G", 1:150)
.conv_ids_expected <- c("pk_conv::libA::Alpha pathway", "pk_conv::libA::Alpha pathway bis",
                        "pk_conv::libA::Beta pathway", "pk_conv::libB::Alpha pathway",
                        "pk_conv::libB::Delta pathway")

#' Fixture complète : registre + run PK-2 combiné (exactement la voie du
#' module : packs_term2gene multi-librairies → UN SEUL run_packs_enrichment —
#' aucun moteur re-run, PK-3 consomme ce que PK-2 a produit).
.conv_fixture <- function() {
  d <- .tmp_conv_dir()
  reg <- .conv_pack_registry(d)
  t2g <- packs_term2gene(list(reg), c("pk_conv||libA.gmt", "pk_conv||libB.gmt"))
  combined <- run_packs_enrichment(genes = .conv_genes,
                                   term2gene = t2g$term2gene,
                                   term2name = t2g$term2name,
                                   universe = .conv_universe,
                                   pval_cutoff = 1, qval_cutoff = 1)
  list(registry = reg, t2g = t2g, combined = combined)
}

# -- packs_split_libraries : découpe de la table combinée PK-2 ----------------

test_that("packs_split_libraries découpe le run combiné par préfixe pack::fichier et reporte enrich_obj", {
  fx <- .conv_fixture()
  split <- packs_split_libraries(fx$combined)
  expect_type(split, "list")
  expect_setequal(names(split), c("pk_conv::libA", "pk_conv::libB"))
  expect_identical(sum(vapply(split, nrow, integer(1))), nrow(fx$combined))
  # chaque part porte l'objet enrichResult brut du run combiné (contrat STAT-S2)
  for (p in split) {
    expect_s3_class(p, "data.frame")
    expect_identical(attr(p, "enrich_obj"), attr(fx$combined, "enrich_obj"))
  }
  # les IDs d'une part partagent le même préfixe librairie
  expect_true(all(grepl("^pk_conv::libA::", split[["pk_conv::libA"]]$ID)))
  expect_true(all(grepl("^pk_conv::libB::", split[["pk_conv::libB"]]$ID)))
  # table 0 ligne (contrat PK-2 « 0 hit ») → liste vide, jamais d'erreur
  empty <- fx$combined[0, , drop = FALSE]
  expect_identical(packs_split_libraries(empty), list())
})

test_that("packs_split_libraries refuse honnêtement : table sans ID, IDs sans préfixe pack::fichier::", {
  expect_error(packs_split_libraries(data.frame(x = 1)), class = "pathway_error")
  # un run NON-Packs (KEGG : IDs sans préfixe) n'est PAS un run multi-librairies
  kegg_like <- data.frame(ID = c("hsa04110", "hsa04150"),
                          Description = c("a", "b"), pvalue = c(0.01, 0.02))
  expect_error(packs_split_libraries(kegg_like), class = "pathway_error")
})

# -- packs_converge : BH GLOBALE (option A, D5 GO 2026-11-10) -----------------

test_that("packs_converge : p.adjust GLOBALE identique bit à bit à la recomputation manuelle p.adjust(, 'BH') sur l'union", {
  fx <- .conv_fixture()
  split <- packs_split_libraries(fx$combined)
  conv <- packs_converge(split)
  expect_s3_class(conv, "data.frame")
  # les 5 termes testés (Gamma sans hit n'existe pas)
  expect_identical(sort(conv$ID), sort(.conv_ids_expected))
  # colonnes de la note de conception §3 (ID préfixé, Description, p brut,
  # p.adjust GLOBAL, librairie(s), cluster, score) + audit PK-2
  expect_true(all(c("ID", "Description", "library", "pvalue", "p.adjust",
                    "p.adjust_pk2", "cluster", "n_libraries",
                    "cluster_libraries") %in% names(conv)))
  # PREUVE EXIGÉE : BH globale = p.adjust() manuel sur l'union des tests
  expected_bh <- stats::p.adjust(fx$combined$pvalue, method = "BH")
  names(expected_bh) <- fx$combined$ID
  expect_identical(conv$p.adjust, unname(expected_bh[conv$ID]))
  # indépendance d'ordre : l'union dans l'autre ordre de librairies donne la
  # même BH par ID (BH est invariante à l'ordre des tests)
  conv_rev <- packs_converge(split[rev(names(split))])
  a <- setNames(conv$p.adjust, conv$ID)
  b <- setNames(conv_rev$p.adjust, conv_rev$ID)
  expect_identical(a[sort(names(a))], b[sort(names(b))])
})

test_that("packs_converge : p.adjust PK-2 porté en colonne dédiée ; union complète → même dénominateur BH (équivalence mesurée), union partielle → divergence honnête", {
  fx <- .conv_fixture()
  conv <- packs_converge(packs_split_libraries(fx$combined))
  # cutoffs à 1 : l'union des lignes EST l'ensemble des tests du run PK-2 →
  # la BH PK-2 et la BH globale partagent le même dénominateur
  expect_equal(conv$p.adjust, conv$p.adjust_pk2)
  # avec les cutoffs par défaut, des tests sont éliminés du contrat PK-2 →
  # la BH globale (recalculée sur l'union PORTÉE) diverge de la BH PK-2, dans
  # la direction documentée : MOINS de tests → BH moins conservatrice
  # (all.equal est INUTILISABLE ici : les valeurs < ~1.5e-8 en valeur absolue
  # sont ignorées par sa tolérance — comparaison exacte exigée)
  strict <- run_packs_enrichment(genes = .conv_genes,
                                 term2gene = fx$t2g$term2gene,
                                 term2name = fx$t2g$term2name,
                                 universe = .conv_universe, pval_cutoff = 0.05)
  if (nrow(strict) >= 1L && nrow(strict) < nrow(fx$combined)) {
    conv_strict <- packs_converge(packs_split_libraries(strict))
    expect_false(identical(conv_strict$p.adjust, conv_strict$p.adjust_pk2))
    expect_true(all(conv_strict$p.adjust <= conv_strict$p.adjust_pk2),
                info = "BH sur l'union portée (moins de tests) <= BH PK-2 (tous tests)")
  }
})

# -- packs_converge : effondrement pairwise_termsim, 1 cluster = 1 vote -------

test_that("packs_converge : termes à jeux de gènes identiques → 1 cluster inter-librairies ; score = librairies DISTINCTES significatives (pas le nombre de lignes)", {
  fx <- .conv_fixture()
  conv <- packs_converge(packs_split_libraries(fx$combined))
  row_of <- function(id) conv[match(id, conv$ID), , drop = FALSE]
  a1 <- row_of("pk_conv::libA::Alpha pathway")
  a2 <- row_of("pk_conv::libA::Alpha pathway bis")
  a3 <- row_of("pk_conv::libB::Alpha pathway")
  # JC(libA Alpha, libB Alpha) = 1 (mêmes gènes) → un SEUL cluster
  expect_identical(a1$cluster, a3$cluster)
  expect_identical(a2$cluster, a3$cluster)
  # Beta / Delta disjoints → clusters propres
  expect_false(a1$cluster == row_of("pk_conv::libA::Beta pathway")$cluster)
  expect_false(a1$cluster == row_of("pk_conv::libB::Delta pathway")$cluster)
  expect_identical(length(unique(conv$cluster)), 3L)
  # « 1 cluster = 1 vote » : le score compte les LIBRAIRIES DISTINCTES (2),
  # pas les lignes du cluster (3)
  expect_identical(a1$n_libraries, 2L)
  expect_identical(a2$n_libraries, 2L)
  expect_identical(a3$n_libraries, 2L)
  expect_identical(a1$cluster_libraries, "pk_conv::libA, pk_conv::libB")
  # Beta/Delta : pas significatifs sous le seuil global 0.05 → score 0
  beta <- row_of("pk_conv::libA::Beta pathway")
  delta <- row_of("pk_conv::libB::Delta pathway")
  expect_true(beta$p.adjust > 0.05)
  expect_true(delta$p.adjust > 0.05)
  expect_identical(beta$n_libraries, 0L)
  expect_identical(delta$n_libraries, 0L)
  # le seuil global GATE le vote : à p_threshold = 1, chaque terme isolé
  # redevient significatif dans SA librairie (score 1), Alpha reste 2
  conv1 <- packs_converge(packs_split_libraries(fx$combined), p_threshold = 1)
  expect_identical(conv1$n_libraries[match("pk_conv::libA::Beta pathway", conv1$ID)], 1L)
  expect_identical(conv1$n_libraries[match("pk_conv::libB::Delta pathway", conv1$ID)], 1L)
  expect_identical(conv1$n_libraries[match("pk_conv::libB::Alpha pathway", conv1$ID)], 2L)
})

test_that("packs_converge : déterminisme strict (même entrée → même table) + paramètres d'audit portés", {
  fx <- .conv_fixture()
  split <- packs_split_libraries(fx$combined)
  expect_identical(packs_converge(split), packs_converge(split))
  conv <- packs_converge(split, p_threshold = 0.1, sim_cutoff = 0.2)
  prm <- attr(conv, "converge_params")
  expect_identical(prm$p_threshold, 0.1)
  expect_identical(prm$p_adjust_method, "BH")
  expect_identical(prm$method, "JC")
})

test_that("packs_converge : flux à objets enrichResult PAR librairie — aucune arête inter-librairies (conservateur, documenté), BH globale inchangée", {
  fx <- .conv_fixture()
  # DEUX runs PK-2 distincts (un par librairie) : chacun porte SON enrich_obj
  t2gA <- packs_term2gene(list(fx$registry), "pk_conv||libA.gmt")
  t2gB <- packs_term2gene(list(fx$registry), "pk_conv||libB.gmt")
  runA <- run_packs_enrichment(.conv_genes, t2gA$term2gene, t2gA$term2name,
                               universe = .conv_universe,
                               pval_cutoff = 1, qval_cutoff = 1)
  runB <- run_packs_enrichment(.conv_genes, t2gB$term2gene, t2gB$term2name,
                               universe = .conv_universe,
                               pval_cutoff = 1, qval_cutoff = 1)
  conv <- packs_converge(list("pk_conv::libA" = runA, "pk_conv::libB" = runB))
  expect_identical(sort(conv$ID), sort(.conv_ids_expected))
  # objets distincts → la similarité INTER-OBJETS n'est pas calculable par
  # pairwise_termsim → 0-fill conservateur : pas d'arête inter-librairies ;
  # les arêtes INTRA-objet subsistent (AlphaA ≡ AlphaBisA de libA → 1 cluster
  # commun) → 4 clusters, pas 5
  expect_identical(length(unique(conv$cluster)), 4L)
  expect_identical(conv$n_libraries[match("pk_conv::libA::Alpha pathway", conv$ID)], 1L)
  expect_identical(conv$n_libraries[match("pk_conv::libB::Alpha pathway", conv$ID)], 1L)
  # MESURE HONNÊTE (probe 2026-11-10) : les p brutes NE SONT PAS transversales
  # aux deux flux — enricher() restreint la liste ET l'univers aux gènes
  # annotés par TERM2GENE (un run par librairie teste un univers différent) :
  # p(Alpha) combiné 4.99e-20 vs ~7e-19 (libA seul : univers 90) vs ~2.7e-12
  # (libB seul : liste rognée à 35). C'est exactement pourquoi la voie
  # nominale multi-librairies PK-2 est le run COMBINÉ (union des annotations)
  # — PK-3 consomme ce run combiné ; le flux par-librairie reste supporté
  # (plomberie prouvée ci-dessus) mais ses p brutes sont les siennes.
  conv_cmb <- packs_converge(packs_split_libraries(fx$combined))
  expect_false(identical(sort(conv$p.adjust), sort(conv_cmb$p.adjust)))
})

# -- packs_converge : refus honnêtes (classe pathway_error, style PK-2) -------

test_that("packs_converge refuse honnêtement : entrée vide, moins de 2 librairies, clés absentes", {
  expect_error(packs_converge(list()), class = "pathway_error")
  expect_error(packs_converge(NULL), class = "pathway_error")
  fx <- .conv_fixture()
  split <- packs_split_libraries(fx$combined)
  # une seule librairie → pas de convergence possible
  one <- split[1]
  expect_error(packs_converge(one), class = "pathway_error")
  expect_error(packs_converge(one), "2 librairies")
  # liste non nommée : les clés de librairie sont le contrat
  expect_error(packs_converge(unname(split)), class = "pathway_error")
})

test_that("packs_converge refuse honnêtement : colonne pvalue absente, lignes malformées, doublons (librairie, terme)", {
  fx <- .conv_fixture()
  split <- packs_split_libraries(fx$combined)
  k1 <- names(split)[1]
  .rebuild <- function(part1) {
    l <- c(list(part1), split[-1]); names(l) <- names(split); l
  }
  # colonne pvalue absente (contrat PK-2 : pvalue est portée)
  bad <- split[[k1]][, setdiff(names(split[[k1]]), "pvalue"), drop = FALSE]
  expect_error(packs_converge(.rebuild(bad)), "pvalue", class = "pathway_error")
  # pvalue NA / hors [0;1] / ID vide → lignes malformées
  m1 <- split[[k1]]; m1$pvalue[1] <- NA_real_
  expect_error(packs_converge(.rebuild(m1)), class = "pathway_error")
  m2 <- split[[k1]]; m2$pvalue[1] <- 1.5
  expect_error(packs_converge(.rebuild(m2)), class = "pathway_error")
  m3 <- split[[k1]]; m3$ID[1] <- "  "
  expect_error(packs_converge(.rebuild(m3)), class = "pathway_error")
  # doublon (librairie, terme) dans la même librairie
  m4 <- rbind(split[[k1]], split[[k1]][1, , drop = FALSE])
  expect_error(packs_converge(.rebuild(m4)), "dupliqu", class = "pathway_error")
})

test_that("packs_converge refuse honnêtement : enrich_obj absent, couverture incomplète, méthode non supportée", {
  fx <- .conv_fixture()
  split <- packs_split_libraries(fx$combined)
  # aucun attribut enrich_obj (contrat STAT-S2 violé) → la similarité
  # enrichplot est incalculable sans re-run (interdit) → refus
  bare <- lapply(split, function(p) { attr(p, "enrich_obj") <- NULL; p })
  expect_error(packs_converge(bare), class = "pathway_error")
  # couverture incomplète : une ligne dont l'ID n'existe dans AUCUN enrich_obj
  g <- split[[1]][1, , drop = FALSE]
  g$ID <- "pk_conv::libA::Terme fantôme"
  split_ghost <- split
  split_ghost[[1]] <- rbind(split[[1]], g)
  expect_error(packs_converge(split_ghost), "Terme fantôme", class = "pathway_error")
  # méthode de similarité non supportée en v1 (JC seul : sans semData)
  expect_error(packs_converge(split, method = "Wang"), class = "pathway_error")
})

# -- build_convergence_dt : rendu DT au style maison ---------------------------

test_that("build_convergence_dt rend un DT-widget au style maison (ts_datatable)", {
  fx <- .conv_fixture()
  conv <- packs_converge(packs_split_libraries(fx$combined))
  w <- build_convergence_dt(conv, tr = function(x) x)
  expect_s3_class(w, "htmlwidget")
  expect_true("datatables" %in% class(w))
  expect_false(is.null(w$x))
})

# -- câblage modules (assertions source-level, style freeze tests) ------------

test_that("les deux modules pathway portent la vue « Convergence » PK-3 (garde source-level)", {
  .mod_src <- function(p) paste(readLines(file.path(ts_project_root(), p), warn = FALSE), collapse = "\n")
  for (m in c("modules/bulk/mod_bulk_pathways.R", "modules/sc/mod_sc_pathways.R")) {
    src <- .mod_src(m)
    expect_match(src, 'i18n$t("Convergence")', fixed = TRUE,
                 info = paste(m, "— nav_panel Convergence absent"))
    expect_match(src, "packs_split_libraries(", fixed = TRUE,
                 info = paste(m, "— appel packs_split_libraries absent"))
    expect_match(src, "packs_converge(", fixed = TRUE,
                 info = paste(m, "— appel packs_converge absent"))
    expect_match(src, "build_convergence_dt(", fixed = TRUE,
                 info = paste(m, "— appel build_convergence_dt absent"))
    expect_match(src, "PK-3", fixed = TRUE,
                 info = paste(m, "— marqueur jalon PK-3 absent"))
  }
  # la section PK-3 vit dans le MÊME fichier moteur (section PACKS), pas ailleurs
  ph <- paste(readLines(file.path(ts_project_root(), "tools/pkg_skeleton/R/pathway_helpers.R"),
                        warn = FALSE), collapse = "\n")
  expect_match(ph, "packs_converge <- function", fixed = TRUE)
  expect_match(ph, "packs_split_libraries <- function", fixed = TRUE)
  expect_match(ph, "build_convergence_dt <- function", fixed = TRUE)
  expect_match(ph, "PK-3", fixed = TRUE)
})

test_that("i18n : les clés PK-3 sont présentes avec fr ET en non vides (appends comptés)", {
  tr <- jsonlite::fromJSON(file.path(ts_project_root(), "i18n", "translation.json"),
                           simplifyVector = FALSE)$translation
  fr_keys <- vapply(tr, function(e) e$fr %||% NA_character_, character(1))
  en_vals <- setNames(vapply(tr, function(e) e$en %||% "", character(1)), fr_keys)
  for (k in c("Convergence",
              "La convergence s'applique aux runs Packs multi-librairies (≥ 2 librairies sélectionnées).",
              "BH globale sur l'union des tests (toutes librairies confondues) + effondrement de redondance (similarité de Jaccard enrichplot) — 1 cluster = 1 vote.",
              "Librairie", "p (brut)", "p.adjust GLOBAL", "p.adjust PK-2 (run)",
              "Cluster", "Score (librairies)", "Librairies convergentes")) {
    expect_true(k %in% fr_keys, info = paste("clé i18n (fr) absente :", k))
    # && paresseux : l'indexation [[k]] ne doit pas masquer l'absence de la clé
    expect_true(k %in% fr_keys && nzchar(en_vals[[k]]),
                info = paste("traduction en absente :", k))
  }
})
