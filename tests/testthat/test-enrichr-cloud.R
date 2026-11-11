# =============================================================================
# test-enrichr-cloud.R — PK-5 (pilier 2bis, fiche PROPOSAL_PACKS_ECOSYSTEM
# 2026-11-09.md §3 jalon PK-5, idéation §4.2 ; GO vie privée opérateur
# 2026-11-11 = B-1) : mode cloud Enrichr OPT-IN EXPLICITE.
#
# Créé AVANT l'implémentation (règle maison : test d'abord, rouge d'abord).
#
# Contrat vie privée (D3) : l'argument `consent` DOIT être littéralement TRUE
# pour qu'UN SEUL octet quitte la machine — le refus est levé AVANT toute
# construction de requête, donc il est testable HORS réseau. Les tests de
# cette suite ne contactent JAMAIS maayanlab.cloud : parsing sur fixtures,
# construction d'URL pure, gardes d'arguments avant vol ; le seul test qui
# touche une URL réseau vise localhost:1 (refus instantané, aucun paquet vers
# l'extérieur) pour prouver l'enveloppement des erreurs réseau en
# pathway_error. Défaut = local strict (PK-2) ; échec réseau = pathway_error
# VISIBLE, jamais un repli silencieux vers le local.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/pathway_helpers.R")

# -- fixtures ------------------------------------------------------------------

# TSV /export natif Enrichr (colonnes mesurées via le client officiel enrichR) :
# Term, Overlap, P.value, Adjusted.P.value, Old.P.value, Old.Adjusted.P.value,
# Odds.Ratio, Combined.Score, Genes. Deux blocs (2 librairies) concaténés
# comme le ferait enrichr_cloud_enrich().
.enrichr_tsv <- function() paste0(
  "Term\tOverlap\tP.value\tAdjusted.P.value\tOld.P.value\tOld.Adjusted.P.value\tOdds.Ratio\tCombined.Score\tGenes\n",
  "Breast cancer\t3/70\t0.001\t0.01\t0.001\t0.01\t12.5\t800\tTP53;BRCA1;EGFR\n",
  "Cell cycle\t2/120\t0.008\t0.04\t0.008\t0.04\t5.2\t300\tTP53;EGFR\n",
  "DNA repair\t2/90\t0.03\t0.06\t0.03\t0.06\t4.1\t150\tBRCA1;TP53\n")

# -- garde de consentement (LE contrat vie privée — testé hors réseau) ---------

test_that("consent autre que TRUE = pathway_error AVANT tout vol (contrat D3)", {
  # consent FALSE / NA / NULL : le message porte le refus de vie privée et le
  # fait que RIEN n'a été envoyé — si le gate était après la construction de
  # requête, un appel avec base_url irroutable renverrait une erreur curl,
  # pas ce message (l'ordre gate→validation→réseau est testé au bloc suivant).
  for (bad in list(FALSE, NA, NULL)) {
    expect_error(
      enrichr_cloud_submit(c("TP53", "BRCA1"), consent = bad,
                           base_url = "https://localhost:1"),
      class = "pathway_error", info = paste("consent =", deparse(bad)))
    expect_error(
      enrichr_cloud_run(c("TP53", "BRCA1"), c("KEGG_2021_Human"), consent = bad,
                        base_url = "https://localhost:1"),
      class = "pathway_error", info = paste("run, consent =", deparse(bad)))
  }
  # le message du gate dit explicitement que rien n'a quitté la machine
  err <- tryCatch(enrichr_cloud_submit("TP53", consent = FALSE,
                                       base_url = "https://localhost:1"),
                  error = function(e) e)
  expect_match(conditionMessage(err), "RIEN", ignore.case = FALSE)
})

test_that("les gardes d'arguments passent AVANT le vol (base https + gènes non vides)", {
  # base_url non https refusée sans réseau (validation avant construction)
  expect_error(
    enrichr_cloud_url("http://insecure.example", "enrich"),
    class = "pathway_error")
  # gènes vides refusés avec consent=TRUE (validation avant la requête —
  # sinon ce test tenterait le réseau)
  expect_error(
    enrichr_cloud_submit(character(0), consent = TRUE,
                         base_url = "https://localhost:1"),
    class = "pathway_error")
  expect_error(
    enrichr_cloud_run(c("   ", NA), c("KEGG_2021_Human"), consent = TRUE,
                      base_url = "https://localhost:1"),
    class = "pathway_error")
  # librairies vides refusées côté enrich et run (avant vol)
  expect_error(
    enrichr_cloud_enrich(character(0), base_url = "https://localhost:1"),
    class = "pathway_error")
  expect_error(
    enrichr_cloud_run("TP53", character(0), consent = TRUE,
                      base_url = "https://localhost:1"),
    class = "pathway_error")
  # session absente = refus honnête (sans session l'API répond « expired »,
  # mesuré 2026-11-11 — on refuse avant le vol plutôt que de rendre un
  # résultat vide qui ressemblerait à « rien d'enrichi »)
  expect_error(
    enrichr_cloud_enrich("KEGG_2021_Human", base_url = "https://localhost:1",
                         cookie_jar = NULL),
    class = "pathway_error", info = "session requise pour /export")
})

test_that("échec réseau = pathway_error visible (localhost:1, aucun paquet sortant)", {
  skip_if_not(requireNamespace("httr2", quietly = TRUE))
  expect_error(
    enrichr_cloud_submit("TP53", consent = TRUE, base_url = "https://localhost:1"),
    class = "pathway_error")
  expect_error(
    enrichr_cloud_libraries(base_url = "https://localhost:1"),
    class = "pathway_error")
})

# -- constantes + URL builder (purs) -------------------------------------------

test_that("TS_ENRICHR_BASE_URLS : 5 hôtes https uniques chez maayanlab.cloud", {
  expect_length(TS_ENRICHR_BASE_URLS, 5L)
  expect_true(all(startsWith(TS_ENRICHR_BASE_URLS, "https://")))
  expect_true(all(grepl("maayanlab.cloud", TS_ENRICHR_BASE_URLS)))
  expect_length(unique(TS_ENRICHR_BASE_URLS), 5L)
  expect_true("Enrichr" %in% names(TS_ENRICHR_BASE_URLS))
})

test_that("enrichr_cloud_url construit l'URL finale (slash normalisé + query)", {
  expect_identical(
    enrichr_cloud_url("https://maayanlab.cloud/Enrichr", "enrich"),
    "https://maayanlab.cloud/Enrichr/enrich")
  expect_identical(
    enrichr_cloud_url("https://maayanlab.cloud/Enrichr/", "export",
                      list(file = "API", backgroundType = "KEGG_2021_Human")),
    paste0("https://maayanlab.cloud/Enrichr/export",
           "?file=API&backgroundType=KEGG_2021_Human"))
})

# -- parsing des résultats (purs, fixtures) ------------------------------------

test_that("enrichr_cloud_parse_export : TSV /enrich→export natif en data.frame contractuel", {
  df <- enrichr_cloud_parse_export(.enrichr_tsv(), "KEGG_2021_Human")
  expect_identical(nrow(df), 3L)
  expect_identical(df$Term, c("Breast cancer", "Cell cycle", "DNA repair"))
  expect_identical(df$Overlap, c("3/70", "2/120", "2/90"))
  expect_identical(df$P.value, c(0.001, 0.008, 0.03))
  expect_identical(df$library, c("KEGG_2021_Human", "KEGG_2021_Human", "KEGG_2021_Human"))
  # entrées dégénérées → contrat vide, jamais d'erreur
  expect_identical(nrow(enrichr_cloud_parse_export("", "")), 0L)
  expect_identical(nrow(enrichr_cloud_parse_export("<!DOCTYPE html>garbage", "X")), 0L)
})

test_that("enrichr_cloud_results_table : colonnes du contrat PK-2 + p.adjust BH GLOBAL recalculé", {
  tbl <- enrichr_cloud_results_table(
    enrichr_cloud_parse_export(.enrichr_tsv(), "KEGG_2021_Human"),
    pval_cutoff = 1, p_adjust_method = "BH")
  expect_setequal(c("ID", "Description", "GeneRatio", "Count", "pvalue",
                    "p.adjust", "library", "Genes"), colnames(tbl))
  expect_identical(nrow(tbl), 3L)
  # BH recalculée localement sur l'UNION des termes renvoyés (motif PK-3),
  # la p-ajustée native Enrichr est ignorée — vérification terme à terme :
  expect_identical(tbl$p.adjust,
                   unname(stats::p.adjust(c(0.001, 0.008, 0.03), "BH")))
  expect_identical(tbl$Description, c("Breast cancer", "Cell cycle", "DNA repair"))
  expect_identical(tbl$Count, c(3L, 2L, 2L))
  expect_identical(tbl$GeneRatio, c("3/70", "2/120", "2/90"))
  expect_identical(tbl$Genes, c("TP53;BRCA1;EGFR", "TP53;EGFR", "BRCA1;TP53"))
})

test_that("enrichr_cloud_results_table : cutoff p.adjust filtre + zéro résultat = contrat vide", {
  # p.adjust BH de c(0.001, 0.008, 0.03) = (0.003, 0.012, 0.03) :
  # cutoff 0.005 ne garde que le 1er terme.
  tbl0 <- enrichr_cloud_results_table(
    enrichr_cloud_parse_export(.enrichr_tsv(), "KEGG_2021_Human"),
    pval_cutoff = 0.005, p_adjust_method = "BH")
  expect_identical(nrow(tbl0), 1L)
  expect_identical(tbl0$Description, "Breast cancer")
  # entrée vide → contrat vide, jamais NULL
  tbl1 <- enrichr_cloud_results_table(
    enrichr_cloud_parse_export(character(0), ""), pval_cutoff = 1,
    p_adjust_method = "BH")
  expect_identical(nrow(tbl1), 0L)
  expect_s3_class(tbl1, "data.frame")
})

test_that("enrichr_cloud_parse_libraries : noms triés + dédupliqués, nid {statistics:[…]} mesuré live", {
  # forme mesurée 2026-11-11 : top-level {statistics: [...], categories: [...]}
  raw <- list(
    statistics = list(
      list(libraryName = "KEGG_2021_Human", geneCoverage = 27884L),
      list(libraryName = "Reactome_2022", geneCoverage = 25000L),
      list(libraryName = "KEGG_2021_Human", geneCoverage = 27884L)
    ),
    categories = list(list(categoryId = 1L))
  )
  expect_identical(enrichr_cloud_parse_libraries(raw),
                   c("KEGG_2021_Human", "Reactome_2022"))
  # forme plate (robustesse) et vide
  expect_identical(enrichr_cloud_parse_libraries(list(list(libraryName = "A"))), "A")
  expect_identical(enrichr_cloud_parse_libraries(list()), character(0))
})

test_that("enrichr_cloud_libs_output : corps renderUI pur, dépendances injectées", {
  seen <- list()
  notify <- function(msg, ...) { seen$msg <<- msg; seen$type <<- list(...)$type }
  tr <- identity; t_fmt <- function(tpl, ...) {
    args <- list(...)
    out <- tpl
    for (k in names(args)) out <- gsub(paste0("{", k, "}"), args[[k]], out, fixed = TRUE)
    out
  }
  # mode inactif → NULL (pas de rendu)
  expect_null(enrichr_cloud_libs_output(FALSE, NULL, function(x) NULL,
                                        function() character(0), tr, t_fmt,
                                        notify, function(x) x))
  # fetch en échec → notification honnête + slot discret, cache posé (NA)
  cache <- NULL
  store <- function(x) cache <<- x
  out <- enrichr_cloud_libs_output(TRUE, cache, store,
                                   function() stop("boom"), tr, t_fmt,
                                   notify, function(x) x)
  expect_match(as.character(out), "indisponibles", fixed = TRUE)
  expect_match(seen$msg, "boom", fixed = TRUE)
  expect_identical(seen$type, "warning")
  expect_true(is.na(cache))
  # 2ᵉ rendu : le cache NA est re-posé mais PAS re-fetché (fetch compté)
  fetches <- 0L
  fetch <- function() { fetches <<- fetches + 1L; c("KEGG_2021_Human", "Reactome_2022") }
  out2 <- enrichr_cloud_libs_output(TRUE, NULL, store, fetch, tr, t_fmt,
                                    notify, function(x) x)
  expect_s3_class(out2, "shiny.tag")
  expect_match(as.character(out2), "KEGG_2021_Human", fixed = TRUE)
  expect_identical(fetches, 1L)
  # 3ᵉ rendu avec le cache plein → zéro fetch supplémentaire
  out3 <- enrichr_cloud_libs_output(TRUE, cache, store, fetch, tr, t_fmt,
                                    notify, function(x) x)
  expect_identical(fetches, 1L)
  expect_match(as.character(out3), "Reactome_2022", fixed = TRUE)
})

# -- câblage UI (assertions source-level, style freeze tests) ------------------

test_that("les 2 modules pathway câblent le mode cloud opt-in (choix + moteur + panneau partagé)", {
  .src <- function(p) paste(readLines(file.path(ts_project_root(), p), warn = FALSE), collapse = "\n")
  for (m in c("modules/bulk/mod_bulk_pathways.R", "modules/sc/mod_sc_pathways.R")) {
    src <- .src(m)
    expect_match(src, "EnrichrCloud", fixed = TRUE,
                 info = paste(m, "— base EnrichrCloud absente du sélecteur"))
    expect_match(src, "enrichr_cloud_optin_ui(", fixed = TRUE,
                 info = paste(m, "— panneau opt-in moteur non branché"))
    expect_match(src, "enrichr_cloud_libs_output(", fixed = TRUE,
                 info = paste(m, "— corps renderUI libs non branché"))
    expect_match(src, "enrichr_cloud_run(", fixed = TRUE,
                 info = paste(m, "— moteur cloud non appelé au run"))
  }
  # le MOTEUR porte le contrat vie privée littéral (bannière destination
  # réelle + refus consentement avant tout vol)
  eng <- .src("tools/pkg_skeleton/R/pathway_helpers.R")
  expect_match(eng, "maayanlab.cloud", fixed = TRUE,
               info = "la bannière ne nomme pas la destination réelle")
  expect_match(eng, "pathway_enrichr_consent", fixed = TRUE,
               info = "l'id de la case de consentement est absent du moteur")
  expect_match(eng, "RIEN n'a été envoyé", fixed = TRUE,
               info = "le message de refus vie privée attendu est absent du moteur")
})
