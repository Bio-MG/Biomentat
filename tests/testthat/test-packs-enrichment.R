# =============================================================================
# test-packs-enrichment.R — PK-2 (pilier 2 local, fiche PROPOSAL_PACKS_ECOSYSTEM
# §3 jalon PK-2 ; détail idéation §4.1/§4.3) : ORA locale sur les GMT de packs.
#
# Créé AVANT l'implémentation (règle maison : test d'abord, rouge d'abord) —
# ce fichier doit être ROUGE tant que pathway_helpers.R ne porte pas la
# section « PACKS (PK-2) » (packs_read_gmt / packs_library_choices /
# packs_term2gene / run_packs_enrichment).
#
# Contrat (fiche §4.1) : parsing GMT → TERM2GENE → clusterProfiler::enricher()
# — même objet enrichResult, même contrat de table, mêmes exports/plots.
# Espace SYMBOLES (les librairies Enrichr sont en symboles humains — aucun bitr,
# zéro perte) ; univers = gènes détectés de l'expérience (fiche §5 étape 1).
# Orthologie souris : refus honnête côté modules (D4 = décision JOINT GEMS D5,
# les packs souris attendent) — le moteur reste organisme-free.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/pathway_helpers.R")

# -- fixtures -----------------------------------------------------------------

.tmp_gmt_dir <- function() {
  d <- file.path(tempdir(), sprintf("packs_pk2_%d", as.integer(Sys.time()) + sample(1e6, 1)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

.write_gmt <- function(path, lines) writeLines(lines, path)

# 2 librairies : libA (3 termes, descriptions), libB (1 terme, description vide
# — le format texte Enrichr vérifié en fiche §2.7 : `term \t \t gene1 …`).
.libA_lines <- c(
  paste(c("Term A", "Activation pathway", paste0("GENE", 1:20)), collapse = "\t"),
  paste(c("Term B", "Metabolism pathway", paste0("GENE", 21:40)), collapse = "\t"),
  paste(c("Term C", "Immune pathway", paste0("GENE", 41:60)), collapse = "\t")
)
.libB_lines <- paste(c("Term D", "", paste0("GENE", 61:80)), collapse = "\t")

.pack_entry <- function(id, dir, files, type = "libraries", name = id) {
  structure(
    list(id = id, type = type, name = name, version = "1.0",
         attribution = "test", source = NULL, downloaded = NULL,
         license = NULL, files = files, palettes = NULL,
         dir = dir, manifest = file.path(dir, "manifest.json")),
    class = "biomentat_pack"
  )
}

# -- packs_read_gmt : parsing pur ---------------------------------------------

test_that("packs_read_gmt parse un GMT standard (termes, descriptions, gènes)", {
  d <- .tmp_gmt_dir()
  p <- file.path(d, "libA.gmt"); .write_gmt(p, .libA_lines)
  out <- packs_read_gmt(p)
  expect_type(out, "list")
  expect_setequal(names(out), c("term2gene", "term2name"))
  expect_setequal(out$term2gene$term, c("Term A", "Term B", "Term C"))
  expect_setequal(out$term2gene$gene[out$term2gene$term == "Term A"],
                  paste0("GENE", 1:20))
  expect_setequal(out$term2name$term, c("Term A", "Term B", "Term C"))
  expect_identical(out$term2name$name[out$term2name$term == "Term A"],
                   "Activation pathway")
  # gènes uniques et nettoyés (le même gène répété sur la ligne ne duplique pas)
  d2 <- .tmp_gmt_dir()
  p2 <- file.path(d2, "dup.gmt")
  .write_gmt(p2, paste(c("T", "desc", "G1", "G1", " G2 "), collapse = "\t"))
  out2 <- packs_read_gmt(p2)
  expect_identical(sort(out2$term2gene$gene), c("G1", "G2"))
})

test_that("packs_read_gmt accepte la description vide (format texte Enrichr)", {
  d <- .tmp_gmt_dir()
  p <- file.path(d, "libB.gmt"); .write_gmt(p, .libB_lines)
  out <- packs_read_gmt(p)
  expect_setequal(out$term2gene$term, "Term D")
  expect_identical(out$term2name$name[out$term2name$term == "Term D"], "Term D")
})

test_that("packs_read_gmt refuse honnêtement : fichier absent, ligne malformée, terme dupliqué, fichier vide, terme sans gène", {
  # fichier absent
  expect_error(packs_read_gmt(file.path(tempdir(), "inexistant_pk2.gmt")),
               class = "pathway_error")
  # ligne malformée (< 3 champs) — le message NOMME la ligne
  d <- .tmp_gmt_dir(); p <- file.path(d, "bad.gmt")
  .write_gmt(p, c(.libA_lines[1], "SeulementUnChamp"))
  expect_error(packs_read_gmt(p), "ligne 2", class = "pathway_error")
  # terme dupliqué dans le même fichier
  d2 <- .tmp_gmt_dir(); p2 <- file.path(d2, "dup_term.gmt")
  .write_gmt(p2, c(.libA_lines[1], paste(c("Term A", "x", "G9"), collapse = "\t")))
  expect_error(packs_read_gmt(p2), "dupliqu", class = "pathway_error")
  # fichier vide
  d3 <- .tmp_gmt_dir(); p3 <- file.path(d3, "empty.gmt")
  .write_gmt(p3, character(0))
  expect_error(packs_read_gmt(p3), class = "pathway_error")
  # terme sans aucun gène valide
  d4 <- .tmp_gmt_dir(); p4 <- file.path(d4, "nogene.gmt")
  .write_gmt(p4, paste(c("T", "desc", "", "  "), collapse = "\t"))
  expect_error(packs_read_gmt(p4), class = "pathway_error")
})

# -- packs_library_choices / packs_term2gene : registre → sélection ----------

test_that("packs_library_choices liste les librairies des packs « libraries » uniquement", {
  d <- .tmp_gmt_dir(); .write_gmt(file.path(d, "libA.gmt"), .libA_lines)
  .write_gmt(file.path(d, "libB.gmt"), .libB_lines)
  reg <- list(
    .pack_entry("chem_packs", d, c("libA.gmt", "libB.gmt"), name = "Chem & Gems"),
    .pack_entry("palettes_wes", d, character(0), type = "palette", name = "Palettes")
  )
  ch <- packs_library_choices(reg)
  expect_setequal(unname(ch), ch)  # vecteur nommé
  expect_setequal(names(ch),
                  c("chem_packs||libA.gmt", "chem_packs||libB.gmt"))
  expect_match(ch[["chem_packs||libA.gmt"]], "Chem & Gems", fixed = TRUE)
  expect_match(ch[["chem_packs||libA.gmt"]], "libA", fixed = TRUE)
  # registre vide → vecteur vide (pas d'erreur)
  expect_identical(packs_library_choices(list()), character(0))
  expect_identical(packs_library_choices(NULL), character(0))
})

test_that("packs_term2gene combine les librairies sélectionnées avec préfixe pack::fichier::terme", {
  d <- .tmp_gmt_dir()
  .write_gmt(file.path(d, "libA.gmt"), .libA_lines)
  .write_gmt(file.path(d, "libB.gmt"), .libB_lines)
  reg <- list(.pack_entry("chem_packs", d, c("libA.gmt", "libB.gmt"), name = "Chem"))
  out <- packs_term2gene(reg, c("chem_packs||libA.gmt", "chem_packs||libB.gmt"))
  expect_true("chem_packs::libA::Term A" %in% out$term2gene$term)
  expect_true("chem_packs::libB::Term D" %in% out$term2gene$term)
  expect_setequal(out$term2gene$gene[out$term2gene$term == "chem_packs::libA::Term A"],
                  paste0("GENE", 1:20))
  # la description GMT survive dans term2name (brute — la provenance vit dans l'ID)
  expect_identical(out$term2name$name[out$term2name$term == "chem_packs::libA::Term A"],
                   "Activation pathway")
  # même nom de terme dans deux packs → PAS de collision (préfixe)
  d2 <- .tmp_gmt_dir(); .write_gmt(file.path(d2, "libA.gmt"), .libA_lines)
  reg2 <- list(reg[[1]], .pack_entry("other_packs", d2, "libA.gmt", name = "Other"))
  out2 <- packs_term2gene(reg2, c("chem_packs||libA.gmt", "other_packs||libA.gmt"))
  expect_length(unique(out2$term2gene$term), 6L)
})

test_that("packs_term2gene refuse honnêtement : sélection vide, clé inconnue (liste les clés valides), GMT absent du disque", {
  d <- .tmp_gmt_dir(); .write_gmt(file.path(d, "libA.gmt"), .libA_lines)
  reg <- list(.pack_entry("chem_packs", d, "libA.gmt", name = "Chem"))
  # sélection vide
  expect_error(packs_term2gene(reg, character(0)),
               class = "pathway_error")
  expect_error(packs_term2gene(reg, NULL), class = "pathway_error")
  # clé inconnue — le message liste les clés valides
  expect_error(packs_term2gene(reg, "ghost_pack||x.gmt"), "chem_packs\\|\\|libA\\.gmt",
               class = "pathway_error")
  # GMT déclaré présent au scan mais disparu du disque — refus honnête
  d2 <- .tmp_gmt_dir()
  reg2 <- list(.pack_entry("vanish", d2, "libX.gmt", name = "Vanish"))
  expect_error(packs_term2gene(reg2, "vanish||libX.gmt"), class = "pathway_error")
})

# -- run_packs_enrichment : enricher() en espace symboles ---------------------

test_that("run_packs_enrichment rend le contrat enrichResult (ID/Description/p.adjust/Count/GeneRatio + attr enrich_obj)", {
  term2gene <- data.frame(
    term = rep(c("T1", "T2", "T3"), each = 20),
    gene = c(paste0("GENE", 1:20), paste0("GENE", 21:40), paste0("GENE", 41:60))
  )
  term2name <- data.frame(term = c("T1", "T2", "T3"),
                          name = c("Activation", "Metabolism", "Immune"))
  genes <- c(paste0("GENE", 1:15), paste0("GENE", 61:65))  # 15/20 du T1
  universe <- paste0("GENE", 1:100)
  res <- run_packs_enrichment(genes, term2gene, term2name,
                              universe = universe, pval_cutoff = 0.05)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("ID", "Description", "p.adjust", "Count", "GeneRatio")
                  %in% names(res)))
  expect_true(nrow(res) >= 1L)
  expect_match(res$ID[1], "T1", fixed = TRUE)
  expect_true(res$p.adjust[which(res$ID == "T1")] < 0.05)
  expect_s4_class(attr(res, "enrich_obj"), "enrichResult")
  # la description GMT passe au rendu (TERM2NAME consommé)
  expect_true("Activation" %in% res$Description)
  # l'univers est honoré : le dénominateur BgRatio ne dépasse pas l'univers
  expect_true(all(unlist(strsplit(res$BgRatio[1], "/")) |> as.integer() <= c(100, 100)))
})

test_that("run_packs_enrichment : 0 hit → data.frame vide au même contrat (jamais NULL)", {
  term2gene <- data.frame(term = "T1", gene = paste0("GENE", 1:20))
  genes <- paste0("GENE", 61:70)  # hors des termes, dans l'univers
  res <- run_packs_enrichment(genes, term2gene, universe = paste0("GENE", 1:100),
                              pval_cutoff = 0.05)
  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 0L)
  expect_true(all(c("ID", "Description", "p.adjust", "Count", "GeneRatio")
                  %in% names(res)))
})

test_that("run_packs_enrichment refuse honnêtement : gènes vides, qval hors borne, TERM2GENE malformé", {
  term2gene <- data.frame(term = "T1", gene = paste0("GENE", 1:20))
  expect_error(run_packs_enrichment(character(0), term2gene),
               class = "pathway_error")
  expect_error(run_packs_enrichment(c("A", "B"), term2gene, qval_cutoff = 1.5),
               class = "pathway_error")
  expect_error(run_packs_enrichment(c("A", "B"), data.frame(x = 1, y = 2)),
               class = "pathway_error")
})

# -- câblage modules (assertions source-level, style freeze tests) ------------

test_that("les deux modules pathway portent la branche Packs (garde souris D4 + sélection + moteur packs)", {
  .mod_src <- function(p) paste(readLines(file.path(ts_project_root(), p), warn = FALSE), collapse = "\n")
  for (m in c("modules/bulk/mod_bulk_pathways.R", "modules/sc/mod_sc_pathways.R")) {
    src <- .mod_src(m)
    # le choix sélecteur + l'UI multi-sélection + la branche moteur
    expect_match(src, '"Packs"', fixed = TRUE,
                 info = paste(m, "— choix sélecteur Packs absent"))
    expect_match(src, "packs_libs_ui", fixed = TRUE,
                 info = paste(m, "— uiOutput multi-sélection absent"))
    expect_match(src, "packs_term2gene(", fixed = TRUE,
                 info = paste(m, "— appel packs_term2gene absent"))
    expect_match(src, "run_packs_enrichment(", fixed = TRUE,
                 info = paste(m, "— appel run_packs_enrichment absent"))
    # refus honnête souris (D4 = décision JOINT GEMS D5) — côté module, pas moteur
    expect_match(src, 'identical(input$pathway_org, "human")', fixed = TRUE,
                 info = paste(m, "— garde organisme absente"))
    expect_match(src, "arbitrage D4 (GEMS D5)", fixed = TRUE,
                 info = paste(m, "— message de refus D4 absent"))
    expect_match(src, "PK-2", fixed = TRUE,
                 info = paste(m, "— marqueur jalon PK-2 absent"))
  }
  # i18n : les 5 clés nouvelles sont présentes avec fr ET en non vides
  # (modèle i18n : la chaîne FR est la clé ; translation = array {fr, en})
  tr <- jsonlite::fromJSON(file.path(ts_project_root(), "i18n", "translation.json"),
                           simplifyVector = FALSE)$translation
  fr_keys <- vapply(tr, function(e) e$fr %||% NA_character_, character(1))
  en_vals <- setNames(vapply(tr, function(e) e$en %||% "", character(1)), fr_keys)
  for (k in c("Librairies de packs",
              "Sélectionnez une ou plusieurs librairies",
              "Aucun pack de librairies installé — déposez un dossier dans data_sc/packs/ (effectif à la session suivante).",
              "Packs : sélectionnez au moins une librairie.",
              "Packs : orthologie souris en attente de l'arbitrage D4 (GEMS D5) — librairies symboles humains uniquement.")) {
    expect_true(k %in% fr_keys, info = paste("clé i18n (fr) absente :", k))
    expect_true(nzchar(en_vals[[k]]), info = paste("traduction en absente :", k))
  }
})
