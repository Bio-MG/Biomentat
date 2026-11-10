# =============================================================================
# helpers_pathway.R — Pathway enrichment (ORA + GSEA), shared sc + bulk
# =============================================================================
# Extracted from global.R (refactor, session post-v1.0 Bulk) — single source
# of truth for pathway analysis, called by BOTH mod_sc_pathways.R and
# mod_bulk_pathways.R so the two domains never diverge in methodology.
#
# Contents:
#   - run_pathway_enrichment() : ORA (GO/KEGG/Reactome) via clusterProfiler
#   - run_gsea_enrichment()    : pre-ranked GSEA (GO/KEGG/Reactome)
#   - plot_pathway_barplot(), plot_pathway_dotplot(), build_pathway_dt()
#
# Depends on: clusterProfiler, org.Hs.eg.db / org.Mm.eg.db, ReactomePA
# (KEGG/Reactome only), KEGGREST (optional), ggplot2, viridis, DT.
#
# D2 M5 lot 5 (2026-10-25) : le fichier vit desormais dans l'arbre du paquet
# (tools/pkg_skeleton/R/). Les 4 sites de lecture des constantes config
# TS_PATHWAY_* (thresholds.R:94-96) sont passes en modele « carry + pin »
# (precedents .ts_dt_const / .ts_report_const) : sous globalenv (app.R Bloc 03
# et les tests), la constante config est PORTEE ; sous namespace, le repli
# EPINGLE — EGAL a la valeur config (thresholds.R:94-96) — prend le relais.
# =============================================================================

# Modele « carry + pin » (precedents .ts_dt_const datatable.R:33,
# .ts_report_const report_render.R:38) : porte la constante definie par
# config/thresholds.R quand elle est visible depuis globalenv (app.R et les
# tests qui la sourcent), sinon retombe sur le repli EPINGLE ci-dessous.
.ts_pathway_const <- function(name, fallback) {
  if (exists(name, envir = globalenv(), inherits = TRUE)) {
    get(name, envir = globalenv(), inherits = TRUE)
  } else {
    fallback
  }
}




# M12 (2026-10-17) — resolution of the ORA background, extracted so that it can
# be DECIDED without clusterProfiler.
#
# Two things were wrong when this lived inline inside run_pathway_enrichment():
#   1. `universe = NULL` — or a universe that cleaned down to nothing, or one
#      that mapped to zero Entrez ids — left `universe_entrez` at NULL and was
#      forwarded as-is to enrichGO/enrichKEGG/enrichPathway. clusterProfiler
#      then silently used its DEFAULT background (every annotated gene), which
#      is a different analysis wearing the same name (see the measured gap in
#      test-stats-quickwins.R, "the enrichment BACKGROUND"). Only the
#      "universe supplied but nothing converted" branch warned; the two other
#      fallbacks were mute.
#   2. The note was a bare `warning()`, so a downstream consumer (export,
#      report, drive verdict) could not tell a default-background run from a
#      declared-background one without scraping the console.
#
# `map_fun` is injected (it is `bitr` in production) so the resolution is a
# pure function: no OrgDb, no network, no package needed to test it.
#
# Returns list(entrez = <character or NULL>, note = <character or NULL>).
# `note` is non-NULL exactly when the background FELL BACK; `entrez` is then
# NULL, which is what the enrichment functions receive as "use the default".
.resolve_ora_universe <- function(universe, ids, map_fun) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids)]

  if (is.null(universe)) {
    return(list(entrez = NULL, note = paste0(
      "ORA : aucun univers fourni — le fond par défaut de clusterProfiler ",
      "(tous les gènes annotés) est utilisé. Déclarez l'univers (les gènes ",
      "détectables) pour une ORA comparable.")))
  }

  universe_clean <- unique(trimws(universe[!is.na(universe) &
                                             nchar(trimws(universe)) > 0]))
  if (length(universe_clean) == 0L) {
    return(list(entrez = NULL, note = paste0(
      "ORA : l'univers fourni est VIDE après nettoyage — repli sur le fond ",
      "par défaut de clusterProfiler (tous les gènes annotés).")))
  }

  universe_map <- tryCatch(map_fun(universe_clean), error = function(e) NULL)
  if (is.null(universe_map) || nrow(universe_map) == 0L) {
    return(list(entrez = NULL, note = paste0(
      "ORA : aucun ID de l'univers n'a pu être converti en Entrez — repli sur ",
      "le fond par défaut de clusterProfiler (tous les gènes annotés).")))
  }

  list(entrez = union(unique(universe_map$ENTREZID), ids), note = NULL)
}


# One place that stamps the fallback note on the returned table, so the two
# return paths of run_pathway_enrichment() cannot drift apart (the M8 lesson:
# a duplicated expression is two truths waiting to diverge).
.ora_attach_universe_note <- function(df, note) {
  if (!is.null(note)) attr(df, "universe_fallback") <- note
  df
}


#' Pathway Enrichment Analysis (Robust v3)

#' @param genes Vecteur de gènes (Symbols)

#' @param organism Organisme ("human", "mouse")

#' @param database Base de données ("GOBP", "KEGG", "Reactome")

#' @param pval_cutoff Seuil p-value

#' @param qval_cutoff Seuil q-value, dans ]0, 1]. C'était une CONSTANTE CACHÉE
#'   (0.2, écrite en dur sur le seul site GOBP) alors que `pval_cutoff` était
#'   exposé : un second filtre que l'appelant ne pouvait ni voir ni régler.
#'   Défaut 0.2 = la valeur que les trois bases appliquaient déjà, donc le
#'   comportement des appelants existants est inchangé. Validé avant tout
#'   chargement de paquet (M12).

#' @param p_adjust_method Méthode de correction pour tests multiples ("BH",
#'   "BY", "bonferroni", "holm", ... — voir stats::p.adjust.methods).
#'   Défaut "BH" : comportement historique strictement inchangé.

#' @return data.frame avec pathways enrichis

run_pathway_enrichment <- function(genes, organism = "human",

                                   database = "GOBP",

                                   pval_cutoff = 0.05, universe = NULL,
                                   qval_cutoff = 0.2,
                                   p_adjust_method = "BH") {

  

  # STAT-Q1 : méthode de correction exposée à l'utilisateur. Validation contre
  # stats::p.adjust.methods (liste de référence de R) plutôt que contre
  # TS_PADJ_METHODS, qui n'est que le sous-ensemble offert dans l'UI — un
  # appelant programmatique reste libre d'utiliser "none" ou "hommel".
  p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)

  # B1-C (O21, 2026-10-28) : identifiants composites ENSG…_SYMBOL / ENSG….N
  # acceptés — même règle que le chemin de mapping, sinon bitr rend un df
  # vide en silence et l'univers désaligne les gènes (dossier round 5).
  genes <- strip_composite_gene_ids(genes)
  if (!is.null(universe)) universe <- strip_composite_gene_ids(universe)

  # M12 : le q-filtre était une CONSTANTE CACHÉE (0.2, écrite en dur sur le seul
  # site GOBP) alors que `pvalueCutoff` est exposé. Il devient un argument, validé
  # AVANT tout chargement de paquet : le refus est immédiat et testable sans
  # Bioconductor. Défaut 0.2 = la valeur que les trois bases appliquaient déjà.
  qval_cutoff <- suppressWarnings(as.numeric(qval_cutoff))
  if (length(qval_cutoff) != 1L || is.na(qval_cutoff) ||
      qval_cutoff <= 0 || qval_cutoff > 1) {
    stop(errorCondition("qval_cutoff doit être un nombre dans ]0, 1].",
                        class = "pathway_error"))
  }

  # 1. Check Core Dependencies

  if (!requireNamespace("clusterProfiler", quietly = TRUE))

    stop(errorCondition("Package 'clusterProfiler' requis. Installez-le via BiocManager.", class = "pathway_error"))

  

  library(clusterProfiler)

  

  # 2. Prepare Organism Database & Convert Genes

  gene_entrez <- NULL

  

  if (organism == "human") {

    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))

      stop(errorCondition("Package 'org.Hs.eg.db' requis pour l'analyse humaine.", class = "pathway_error"))

    library(org.Hs.eg.db)

    orgdb <- org.Hs.eg.db

    

  } else if (organism == "mouse") {

    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))

      stop(errorCondition("Package 'org.Mm.eg.db' requis pour l'analyse souris.", class = "pathway_error"))

    library(org.Mm.eg.db)

    orgdb <- org.Mm.eg.db

    

  } else {

    stop(errorCondition("Organisme non supporté (choisir 'human' ou 'mouse')", class = "pathway_error"))

  }

  

  # Clean input genes: remove NAs, empty strings, trim whitespace

  genes_clean <- unique(trimws(genes[!is.na(genes) & nchar(trimws(genes)) > 0]))

  

  if (length(genes_clean) == 0) {

    stop(errorCondition("Aucun gène valide fourni après nettoyage.", class = "pathway_error"))

  }

  

  # Attempt conversion Symbol -> Entrez ID

  # Note: bitr can fail if keys are invalid. We wrap it in tryCatch.

  tryCatch({

    gene_entrez <- bitr(

      genes_clean, 

      fromType = "SYMBOL", 

      toType   = "ENTREZID", 

      OrgDb    = orgdb

    )

  }, error = function(e) {

    # If bitr fails completely, return empty

    warning(paste("Erreur lors de la conversion des gènes:", e$message))

    gene_entrez <<- data.frame(SYMBOL=character(0), ENTREZID=character(0))

  })

  

  # Check if any genes were successfully converted

  if (is.null(gene_entrez) || nrow(gene_entrez) == 0) {

    stop(errorCondition("Aucun gène n'a pu être converti en Entrez ID. Vérifiez que les noms de gènes sont des symboles officiels (ex: 'TP53', 'Actb') et correspondent à l'organisme sélectionné.", class = "pathway_error"))

  }

  

  ids <- gene_entrez$ENTREZID

  # M12 : résolution PURE du fond (cf. .resolve_ora_universe() en tête de
  # fichier). La note de repli est conservée : elle est à la fois avertie ici
  # et estampillée sur la table rendue, pour qu'un consommateur en aval (export,
  # rapport, verdict drive) puisse distinguer un fond DÉCLARÉ d'un fond par
  # défaut sans lire la console.
  universe_res <- .resolve_ora_universe(
    universe, ids,
    function(keys) bitr(keys, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb))
  universe_entrez <- universe_res$entrez
  universe_note <- universe_res$note
  if (!is.null(universe_note)) warning(universe_note)

  # 3. Run Enrichment based on Database

  enrich_result <- NULL

  

  if (database == "GOBP") {

    enrich_result <- enrichGO(

      gene          = ids,

      OrgDb         = orgdb,

      ont           = "BP",

      pAdjustMethod = p_adjust_method,

      pvalueCutoff  = pval_cutoff,

      qvalueCutoff  = qval_cutoff,

      readable      = FALSE, # Keep as Entrez for consistency, or TRUE for Symbols

      universe      = universe_entrez

    )

    

  } else if (database == "KEGG") {

    # KEGG requires organism code (hsa/mmu)

    org_code <- if (organism == "human") "hsa" else "mmu"

    

    # Check if KEGGREST is available for newer clusterProfiler versions

    if (!requireNamespace("KEGGREST", quietly = TRUE)) {

      warning("KEGGREST recommandé pour KEGG. Installation suggérée.")

    }

    

    enrich_result <- enrichKEGG(

      gene          = ids,

      organism      = org_code,

      pAdjustMethod = p_adjust_method,

      pvalueCutoff  = pval_cutoff,

      qvalueCutoff  = qval_cutoff,

      universe      = universe_entrez

    )

    

  } else if (database == "Reactome") {

    if (!requireNamespace("ReactomePA", quietly = TRUE))

      stop(errorCondition("Package 'ReactomePA' requis pour l'analyse Reactome.", class = "pathway_error"))

    

    library(ReactomePA)

    

    # ReactomePA::enrichPathway expects organism name exactly as per its DB

    # Usually "human" or "mouse" works if the DB is loaded.

    enrich_result <- enrichPathway(

      gene          = ids,

      organism      = organism, 

      pAdjustMethod = p_adjust_method,

      pvalueCutoff  = pval_cutoff,

      qvalueCutoff  = qval_cutoff,

      universe      = universe_entrez

    )

  } else {

    stop(errorCondition("Base de données non supportée (GOBP/KEGG/Reactome)", class = "pathway_error"))

  }

  

  # 4. Format Output

  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {

    return(.ora_attach_universe_note(data.frame(

      ID          = character(0),

      Description = character(0),

      p.adjust    = numeric(0),

      Count       = integer(0),

      GeneRatio   = character(0)

    ), universe_note))

  }

  

  res_df <- as.data.frame(enrich_result)

  # STAT-S2 : attach the raw enrichResult S4 object as an attribute (additive —
  # mirrors the `gsea_obj` pattern of run_gsea_enrichment()). enrichplot network
  # views (emapplot/cnetplot + pairwise_termsim) need this raw object (gene
  # sets, GeneRatio, geneID); consumers that only need the table can ignore it —
  # the data.frame contract for existing callers is unchanged.
  attr(res_df, "enrich_obj") <- enrich_result

  # M12 : même estampille sur le chemin NOMINAL que sur le chemin vide — deux
  # retours, un seul helper, donc aucune dérive possible entre les deux.
  return(.ora_attach_universe_note(res_df, universe_note))

}



#' Gene Set Enrichment Analysis (GSEA, pre-ranked) — complements run_pathway_enrichment()

#'

#' Unlike ORA (run_pathway_enrichment), GSEA does not require an arbitrary

#' significance threshold: it ranks ALL tested genes by a continuous score

#' (here, signed -log10(p) * sign(log2FC), the standard ranking metric) and

#' tests for enrichment along that full ranking. More statistically robust

#' when few genes pass a hard significance cutoff.

#'

#' @param de_results data.frame with at least columns: gene, log2FoldChange, pvalue.

#' @param organism "human" or "mouse".

#' @param database "GOBP", "KEGG", or "Reactome".

#' @param pval_cutoff p-value cutoff passed to the GSEA call.

#' @param p_adjust_method Méthode de correction pour tests multiples ("BH",
#'   "BY", "bonferroni", "holm", ... — voir stats::p.adjust.methods).
#'   Défaut "BH" : comportement historique strictement inchangé.

#' @return data.frame: ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, ...

run_gsea_enrichment <- function(de_results, organism = "human",

                                 database = "GOBP", pval_cutoff = 0.05,
                                 p_adjust_method = "BH") {

  # STAT-Q1 : cf. run_pathway_enrichment() — même validation, même défaut.
  p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)

  # B1-C (O21, 2026-10-28) : même strip composite que le funnel ORA — les
  # composites ENSG…_SYMBOL de la colonne gene écourtaient le rangement GSEA.
  if (!is.null(de_results) && "gene" %in% names(de_results)) {
    de_results$gene <- strip_composite_gene_ids(de_results$gene)
  }



  if (!requireNamespace("clusterProfiler", quietly = TRUE))

    stop(errorCondition("Package 'clusterProfiler' requis. Installez-le via BiocManager.", class = "pathway_error"))

  library(clusterProfiler)



  if (organism == "human") {

    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))

      stop(errorCondition("Package 'org.Hs.eg.db' requis pour l'analyse humaine.", class = "pathway_error"))

    library(org.Hs.eg.db)

    orgdb <- org.Hs.eg.db

  } else if (organism == "mouse") {

    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))

      stop(errorCondition("Package 'org.Mm.eg.db' requis pour l'analyse souris.", class = "pathway_error"))

    library(org.Mm.eg.db)

    orgdb <- org.Mm.eg.db

  } else {

    stop(errorCondition("Organisme non supporté (choisir 'human' ou 'mouse')", class = "pathway_error"))

  }



  required_cols <- c("gene", "log2FoldChange", "pvalue")

  missing_cols <- setdiff(required_cols, colnames(de_results))

  if (length(missing_cols) > 0) {

    stop(errorCondition(paste0("Colonnes manquantes dans les résultats DE pour GSEA : ", paste(missing_cols, collapse = ", ")), class = "pathway_error"))

  }



  de_clean <- de_results[!is.na(de_results$pvalue) & !is.na(de_results$log2FoldChange), ]

  de_clean <- de_clean[!is.na(de_clean$gene) & nchar(trimws(de_clean$gene)) > 0, ]

  if (nrow(de_clean) < 10) {

    stop(errorCondition(paste0("Trop peu de gènes valides (", nrow(de_clean), ") après nettoyage pour GSEA."), class = "pathway_error"))

  }



  gene_entrez <- tryCatch({

    AnnotationDbi::select(orgdb, keys = unique(de_clean$gene), keytype = "SYMBOL", columns = "ENTREZID")

  }, error = function(e) {

    stop(errorCondition(paste0("Erreur lors de la conversion des gènes pour GSEA : ", conditionMessage(e)), class = "pathway_error"))

  })

  gene_entrez <- gene_entrez[!is.na(gene_entrez$ENTREZID), ]

  gene_entrez <- gene_entrez[!duplicated(gene_entrez$SYMBOL), ]



  de_merged <- merge(de_clean, gene_entrez, by.x = "gene", by.y = "SYMBOL")

  if (nrow(de_merged) < 10) {

    stop(errorCondition("Aucun gène n'a pu être converti en Entrez ID pour GSEA. Vérifiez l'organisme sélectionné.", class = "pathway_error"))

  }



  # Standard pre-ranking metric: signed -log10(p), tie-broken by averaging

  # duplicated Entrez IDs (can happen when multiple symbols map to one gene).

  de_merged$rank_metric <- -log10(pmax(de_merged$pvalue, 1e-300)) * sign(de_merged$log2FoldChange)

  de_merged <- stats::aggregate(rank_metric ~ ENTREZID, data = de_merged, FUN = mean)



  ranked <- sort(setNames(de_merged$rank_metric, de_merged$ENTREZID), decreasing = TRUE)



  gsea_result <- if (database == "GOBP") {

    clusterProfiler::gseGO(geneList = ranked, OrgDb = orgdb, ont = "BP",

                           pvalueCutoff = pval_cutoff, pAdjustMethod = p_adjust_method, verbose = FALSE)

  } else if (database == "KEGG") {

    org_code <- if (organism == "human") "hsa" else "mmu"

    clusterProfiler::gseKEGG(geneList = ranked, organism = org_code,

                             pvalueCutoff = pval_cutoff, pAdjustMethod = p_adjust_method, verbose = FALSE)

  } else if (database == "Reactome") {

    if (!requireNamespace("ReactomePA", quietly = TRUE))

      stop(errorCondition("Package 'ReactomePA' requis pour l'analyse Reactome.", class = "pathway_error"))

    library(ReactomePA)

    ReactomePA::gsePathway(geneList = ranked, organism = organism,

                          pvalueCutoff = pval_cutoff, pAdjustMethod = p_adjust_method, verbose = FALSE)

  } else {

    stop(errorCondition("Base de données non supportée pour GSEA (GOBP/KEGG/Reactome)", class = "pathway_error"))

  }



  res_df <- as.data.frame(gsea_result)

  if (nrow(res_df) == 0) {

    return(data.frame(ID = character(0), Description = character(0), setSize = integer(0),

                      enrichmentScore = numeric(0), NES = numeric(0),

                      pvalue = numeric(0), p.adjust = numeric(0)))

  }

  # Normalize column naming to stay compatible with build_pathway_dt() /

  # plot_pathway_barplot() / plot_pathway_dotplot(), which expect "Count"

  # and "GeneRatio" — GSEA uses "setSize" instead, so we alias it.

  res_df$Count     <- res_df$setSize

  res_df$GeneRatio <- paste0(res_df$setSize, "/", length(ranked))



  # Attach the raw gseaResult S4 object as an attribute (additive — does not

  # change the data.frame contract for existing callers). enrichplot::gseaplot2()

  # needs this raw object (with @geneList, @result, etc.) to draw the running

  # enrichment-score curve; consumers that only need the table can ignore it.

  attr(res_df, "gsea_obj") <- gsea_result

  res_df

}








#' Top-N pathway barplot (shared by mod_sc_pathways.R and mod_bulk.R)

plot_pathway_barplot <- function(df, db_label = "", top_n = 15, tr = NULL,
                                 palette = "default", manual_gradient = NULL) {

  tr <- tr %||% function(x) x

  df_top <- head(df, top_n)

  df_top$Description <- factor(df_top$Description, levels = rev(df_top$Description))

  ggplot(df_top, aes(x = Count, y = Description, fill = -log10(p.adjust))) +

    geom_bar(stat = "identity", width = 0.7) +

    expression_continuous_scale(palette, "fill", manual_gradient, base_option = "plasma", direction = -1) +

    labs(title = paste(tr("Top"), top_n, tr("Pathways"), "-", db_label),

         x = tr("Nombre de gènes"), y = NULL, fill = "-log10(P-adj)") +

    theme_minimal(base_size = 12) +

    theme(axis.text.y = element_text(size = 10), plot.title = element_text(face = "bold", size = 14),

          legend.position = "right", panel.grid.major.y = element_blank())

}



#' Pathway dotplot (shared)

plot_pathway_dotplot <- function(df, db_label = "", top_n = 20, tr = NULL,
                                 palette = "default", manual_gradient = NULL) {

  tr <- tr %||% function(x) x

  df_top <- head(df, top_n)

  df_top$Description <- factor(df_top$Description, levels = rev(df_top$Description))

  if ("GeneRatio" %in% colnames(df_top)) {

    df_top$GeneRatioNum <- sapply(df_top$GeneRatio, function(r) {

      parts <- strsplit(as.character(r), "/")[[1]]

      if (length(parts) == 2) as.numeric(parts[1]) / as.numeric(parts[2]) else NA

    })

  } else {

    df_top$GeneRatioNum <- df_top$Count / max(df_top$Count, na.rm = TRUE)

  }

  ggplot(df_top, aes(x = GeneRatioNum, y = Description, color = -log10(p.adjust), size = Count)) +

    geom_point(alpha = 0.85) +

    expression_continuous_scale(palette, "color", manual_gradient, base_option = "magma", direction = -1) +

    labs(title = paste(tr("Dotplot Pathways"), "-", db_label),

         x = tr("Ratio de gènes"), y = NULL, color = "-log10(P-adj)", size = tr("Nombre de gènes")) +

    theme_minimal(base_size = 12) +

    theme(axis.text.y = element_text(size = 10), plot.title = element_text(face = "bold", size = 14),

          legend.position = "right")

}



#' Pathway enrichment network (emapplot / cnetplot) — STAT-S2
#'
#' Pure consumer of the RAW enrichment object (enrichResult/gseaResult) stored
#' as attribute `enrich_obj` on the results data.frame by
#' `run_pathway_enrichment()` / `run_gsea_enrichment()` — the same additive
#' pattern as `gsea_obj`. The network is DESCRIPTIVE: edges encode gene
#' similarity between pathways (emap) or pathway↔gene membership (cnet),
#' never causality.
#'
#' @param df results data.frame carrying the `enrich_obj` attribute.
#' @param db_label label of the enrichment database, shown in the title.
#' @param top_n maximum number of pathways shown (selected on p.adjust).
#' @param mode "emap" (pathway–pathway similarity network) or "cnet"
#'   (pathway–gene bipartite network).
#' @param tr optional translation function.
#' @return a ggplot object.
#'
#' @export
plot_pathway_network <- function(df, db_label = "", top_n = 30,
                                 mode = c("emap", "cnet"), tr = NULL) {

  tr <- tr %||% function(x) x
  mode <- match.arg(mode)

  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 2) {
    stop("top_n doit être un nombre >= 2 (un réseau d'une seule voie n'a pas de sens).", call. = FALSE)
  }

  enrich_obj <- attr(df, "enrich_obj")

  if (is.null(enrich_obj)) {
    stop("Aucun objet d'enrichissement brut attaché à ce résultat (attribut 'enrich_obj' absent). Relancez l'enrichissement pour produire le réseau.", call. = FALSE)
  }

  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Package 'enrichplot' requis. Installez-le via BiocManager.", call. = FALSE)
  }

  if (nrow(as.data.frame(enrich_obj)) < 2) {
    stop("Au moins deux voies enrichies sont nécessaires pour tracer un réseau.", call. = FALSE)
  }

  if (mode == "emap") {
    obj <- enrichplot::pairwise_termsim(enrich_obj)
    p <- enrichplot::emapplot(obj, showCategory = top_n)
    subtitle <- tr("Voies reliées par similarité de gènes (descriptif)")
  } else {
    p <- enrichplot::cnetplot(enrich_obj, showCategory = top_n)
    subtitle <- tr("Voies reliées à leurs gènes (descriptif)")
  }

  p + ggplot2::labs(
    title = paste(tr("Réseau d'enrichissement"), "-", db_label),
    subtitle = subtitle
  )

}


# ── STAT-S2 V2 : réseau d'enrichissement INTERACTIF (plotly + igraph) ───────
# Cadrage §2bd.5 #3 (2026-09-16) : la V1 (`plot_pathway_network()`) est
# statique ; la V2 est un réseau INTERACTIF. Choix tranché le 2026-09-28
# (décision utilisateur) : plotly + igraph — ZÉRO dépendance nouvelle (toutes
# deux déjà au renv.lock et installées ; précédent PLOT-S6 pour le toggle
# plotly). visNetwork est ÉCARTÉ pour cette version (dépendance nouvelle —
# évaluation documentée dans `docs/mcp_propagation.md` §12 et le contrat).

#' Enrichment network data — STAT-S2 V2 (structure figée)
#'
#' Construit la structure nœuds/arêtes consommée par
#' `plot_pathway_network_interactive()`. R PUR, sans enrichplot : la
#' similarité terme–terme (mode "emap") est un Jaccard ordinaire sur les jeux
#' de gènes, ce qui rend le constructeur testable sans Bioconductor et
#' reproductible indépendamment du rendu.
#' Contrat : `docs/contracts/PATHWAY_NETWORK_CONTRACT.md`.
#'
#' @param df résultats portant l'attribut `enrich_obj` (même convention que
#'   `plot_pathway_network()`).
#' @param top_n nombre maximal de termes (triés sur `p.adjust`), plafonné à
#'   `TS_PATHWAY_NETWORK_MAX_TERMS`.
#' @param mode "emap" (termes ↔ termes) ou "cnet" (termes ↔ gènes).
#' @return une liste de classe `pathway_network_data` : `nodes`
#'   (id/label/kind/count/p.adjust), `edges` (from/to/weight), `meta`.
#'
#' @export
build_pathway_network_data <- function(df, top_n = 30, mode = c("emap", "cnet")) {

  mode <- match.arg(mode)

  if (!is.data.frame(df)) {
    stop(errorCondition(paste0("df doit être un data.frame de résultats de voies (reçu : ", class(df)[1], ")."),
                        class = "pathway_error"))
  }

  enrich_obj <- attr(df, "enrich_obj")
  if (is.null(enrich_obj)) {
    stop(errorCondition("Aucun objet d'enrichissement brut attaché à ce résultat (attribut 'enrich_obj' absent). Relancez l'enrichissement pour produire le réseau.",
                        class = "pathway_error"))
  }

  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n < 2) {
    stop(errorCondition("top_n doit être un nombre >= 2 (un réseau d'une seule voie n'a pas de sens).",
                        class = "pathway_error"))
  }
  top_n <- min(as.integer(top_n), .ts_pathway_const("TS_PATHWAY_NETWORK_MAX_TERMS", 300L))

  res <- as.data.frame(enrich_obj)
  if (nrow(res) < 2) {
    stop(errorCondition("Au moins deux voies enrichies sont nécessaires pour tracer un réseau.",
                        class = "pathway_error"))
  }

  needed <- c("ID", "Description", "Count", "p.adjust", "geneID")
  if (!all(needed %in% names(res))) {
    stop(errorCondition(paste0("Colonnes manquantes dans l'objet d'enrichissement : ",
                               paste(setdiff(needed, names(res)), collapse = ", "), "."),
                        class = "pathway_error"))
  }

  terms <- res[order(res$p.adjust), ][seq_len(min(top_n, nrow(res))), ]

  gene_sets <- strsplit(terms$geneID, "/", fixed = TRUE)
  term_nodes <- data.frame(
    id = as.character(terms$ID),
    label = ifelse(is.na(terms$Description) | !nzchar(terms$Description),
                   as.character(terms$ID), as.character(terms$Description)),
    kind = "term",
    count = as.integer(terms$Count),
    p.adjust = as.numeric(terms$p.adjust),
    stringsAsFactors = FALSE
  )

  if (mode == "emap") {
    edges <- .pathway_emap_edges(term_nodes$id, gene_sets)
    nodes <- term_nodes
    min_sim <- .ts_pathway_const("TS_PATHWAY_NET_MIN_SIM", 0.20)
  } else {
    edges <- .pathway_cnet_edges(term_nodes$id, gene_sets)
    gene_ids <- sort(unique(edges$to))
    nodes <- rbind(term_nodes,
                   data.frame(id = gene_ids, label = gene_ids, kind = "gene",
                              count = NA_integer_, p.adjust = NA_real_,
                              stringsAsFactors = FALSE))
    min_sim <- NA_real_
  }

  structure(
    list(
      nodes = nodes,
      edges = edges,
      meta = list(mode = mode, top_n = nrow(term_nodes),
                  n_edges = nrow(edges), min_similarity = min_sim,
                  gene_separator = "/")
    ),
    class = "pathway_network_data"
  )

}

#' Jaccard term–term edges (emap) — helper pur de `build_pathway_network_data`
.pathway_emap_edges <- function(term_ids, gene_sets) {
  n <- length(term_ids)
  if (n < 2L) {
    return(data.frame(from = character(0), to = character(0),
                      weight = numeric(0), stringsAsFactors = FALSE))
  }
  from <- character(0); to <- character(0); w <- numeric(0)
  for (i in seq_len(n - 1L)) {
    for (j in seq(i + 1L, n)) {
      inter <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
      if (inter == 0L) next
      sim <- inter / length(union(gene_sets[[i]], gene_sets[[j]]))
      if (sim >= .ts_pathway_const("TS_PATHWAY_NET_MIN_SIM", 0.20)) {
        pair <- sort(c(term_ids[[i]], term_ids[[j]]))
        from <- c(from, pair[[1]]); to <- c(to, pair[[2]]); w <- c(w, sim)
      }
    }
  }
  data.frame(from = from, to = to, weight = w, stringsAsFactors = FALSE)
}

#' Term–gene membership edges (cnet) — helper pur de `build_pathway_network_data`
.pathway_cnet_edges <- function(term_ids, gene_sets) {
  from <- unlist(Map(rep, term_ids, vapply(gene_sets, length, integer(1))), use.names = FALSE)
  to <- unlist(gene_sets, use.names = FALSE)
  data.frame(from = from, to = to, weight = 1,
             stringsAsFactors = FALSE)
}

#' Deterministic network layout (Fruchterman-Reingold, graine figée)
#'
#' Séparé du rendu pour être testable : deux appels avec la même entrée
#' rendent les MÊMES coordonnées (garantie par le contrat).
#'
#' @param net sortie de `build_pathway_network_data()`.
#' @return le data.frame `nodes` de `net`, augmenté de `x`/`y`.
#' @export
pathway_network_layout <- function(net) {

  if (!all(c("nodes", "edges") %in% names(net)) ||
      !is.data.frame(net$nodes) || !is.data.frame(net$edges)) {
    stop(errorCondition("net doit être la sortie de build_pathway_network_data().",
                        class = "pathway_error"))
  }

  g <- igraph::graph_from_data_frame(net$edges, directed = FALSE,
                                     vertices = net$nodes)
  coords <- withr::with_seed(
    .ts_pathway_const("TS_PATHWAY_NETWORK_LAYOUT_SEED", 42L),
    igraph::layout_with_fr(g)
  )
  net$nodes$x <- coords[, 1]
  net$nodes$y <- coords[, 2]
  net$nodes
}

#' Interactive enrichment network — STAT-S2 V2 (plotly)
#'
#' Rendu plotly de la structure `pathway_network_data` : survol d'un terme →
#' libellé, effectif et p.adjust ; survol d'un gène (mode cnet) →
#' identifiant. La barre de mode plotly fournit l'export PNG. Le CLIC
#' (filtrage du tableau par un terme cliqué) est l'item ouvert du contrat —
#' non câblé dans cette version.
#'
#' @param net sortie de `build_pathway_network_data()`.
#' @param title titre du tracé (peut être pré-traduit par l'appelant).
#' @param tr fonction de traduction optionnelle.
#' @return un widget plotly.
#' @export
plot_pathway_network_interactive <- function(net, title = "", tr = NULL) {

  tr <- tr %||% function(x) x

  if (!inherits(net, "pathway_network_data")) {
    stop(errorCondition("net doit être la sortie de build_pathway_network_data().",
                        class = "pathway_error"))
  }

  nodes <- pathway_network_layout(net)
  edges <- net$edges

  # arêtes : un seul trace de segments (x/y alternés, NA sépare les segments)
  edge_x <- as.vector(t(cbind(nodes$x[match(edges$from, nodes$id)],
                              nodes$x[match(edges$to, nodes$id)], NA)))
  edge_y <- as.vector(t(cbind(nodes$y[match(edges$from, nodes$id)],
                              nodes$y[match(edges$to, nodes$id)], NA)))

  terms <- nodes[nodes$kind == "term", ]
  genes <- nodes[nodes$kind == "gene", ]

  term_hover <- paste0(
    "<b>", terms$label, "</b><br>",
    terms$count, " ", tr("g\u00e8nes"), "<br>",
    "p.adjust = ", format.pval(terms$p.adjust, digits = 3)
  )

  p <- plotly::plot_ly(type = "scatter", mode = "lines")
  if (nrow(edges) > 0) {
    p <- plotly::add_trace(p, x = edge_x, y = edge_y, mode = "lines",
                           line = list(color = "#cccccc", width = 1),
                           hoverinfo = "none", showlegend = FALSE)
  }
  if (nrow(genes) > 0) {
    p <- plotly::add_trace(p, data = genes, x = ~x, y = ~y, type = "scatter",
                           mode = "markers", text = genes$label,
                           hoverinfo = "text",
                           marker = list(color = "#999999", size = 7),
                           showlegend = FALSE)
  }
  p <- plotly::add_trace(p, data = terms, x = ~x, y = ~y, type = "scatter",
                         mode = "markers", text = term_hover,
                         hoverinfo = "text",
                         marker = list(
                           color = "#2c7fb8",
                           size = pmax(8, 3.5 * sqrt(pmax(terms$count, 1)))
                         ),
                         showlegend = FALSE)
  p <- plotly::layout(
    p,
    title = title,
    showlegend = FALSE,
    xaxis = list(visible = FALSE),
    yaxis = list(visible = FALSE),
    hoverlabel = list(bgcolor = "white"),
    margin = list(t = 40)
  )
  p
}



#' Pathway results DT table (shared)

build_pathway_dt <- function(df, tr = NULL) {

  tr <- tr %||% function(x) x

  cols_available <- intersect(c("ID", "Description", "p.adjust", "Count", "GeneRatio"), colnames(df))

  df_display <- df[, cols_available, drop = FALSE]

  colnames(df_display) <- c("ID", "Description", "P-adj", tr("Nb Gènes"), "Ratio")[seq_along(cols_available)]

  ts_datatable(df_display, page_length = 15L, filename_base = "pathways") %>%

    DT::formatStyle("P-adj", color = DT::styleInterval(c(0.001, 0.01, 0.05),

                                                        c("darkgreen", "green", "orange", "red")))

}

# =============================================================================
# PB (GEMS D5, GO opérateur 2026-11-09) — ORA DIRECTIONNELLE limma :
# goana/kegga + barcodeplot, le « 3e idiome » d'enrichissement.
# =============================================================================
# L'ORA clusterProfiler (run_pathway_enrichment ci-dessus) travaille sur la
# LISTE des gènes significatifs ; goana/kegga travaillent sur la TABLE DE
# complète et testent la DIRECTION (Up/Down) PAR VOIE. Extension ADDITIVE :
# même discipline de fond que M12 (fond obligatoire, jamais le fond par
# défaut silencieux), même règle d'identifiants composites que B1-C, même
# classe d'erreur pathway_error, map_fun injecté pour la pureté testable.
# compareCluster est EXPLICITEMENT hors v1 (fiche §PB, D5 : « d'abord ») ;
# fgsea/fry/camera restent la propriété du pilier 4 PACKS (jamais
# double-owner). Aucune UX nouvelle : le branchement module arrive en jalon 2.
# =============================================================================

#' ORA directionnelle limma (goana/kegga) sur une table DE complète
#'
#' @param de_table data.frame DE du module bulk : colonnes requises `gene`
#'   (SYMBOL ou identifiants composites ENSG..._SYMBOL — règle B1-C),
#'   `log2FoldChange`, `pvalue` (p-valeur BRUTE — goana applique lui-même le
#'   BH interne ; passer `padj` comme p-valeur produirait une double
#'   correction, refusé). `padj` optionnel, non consommé par limma.
#' @param universe_ids vecteur caractère des identifiants du fond (les
#'   rownames de la matrice filtrée — SYMBOLs comme côté DE). OBLIGATOIRE :
#'   l'ORA directionnelle refuse le fond par défaut de limma (doctrine M12
#'   poussée au 3e idiome — un fond non déclaré est une analyse différente
#'   sous le même nom).
#' @param db "GO", "KEGG" ou c("GO","KEGG") — chaque base demandée est
#'   calculée ; une base non demandée reste NULL (pas de calcul fantôme).
#' @param organism "human" ou "mouse" (org.Hs.eg.db / org.Mm.eg.db).
#' @param fdr_cutoff seuil BH sur la p-valeur brute, appliqué par goana/kegga
#'   pour sélectionner l'ensemble DE dont les directions sont comptées.
#' @param map_fun injection de la fonction de mapping (production :
#'   `function(keys) bitr(keys, fromType="SYMBOL", toType="ENTREZID",
#'   OrgDb=orgdb)`). Injectée, la résolution reste pure et testable sans
#'   OrgDb (motif M12, §2cp).
#'
#' @return list(go =, kegg =, meta =) — `go`/`kegg` : tables limma natives
#'   (colonnes de direction DE.Up/DE.Down/P.Up/P.Down) augmentées d'une
#'   colonne `db` ; NULL pour une base non demandée. `meta` : comptes
#'   honnêtes (de_n_entrez, universe_n_entrez), note de repli de fond (M12),
#'   organisme, seuil.
run_pathway_direction <- function(de_table, universe_ids,
                                  db = c("GO", "KEGG"),
                                  organism = "human",
                                  fdr_cutoff = 0.05,
                                  map_fun = NULL) {
  db <- tryCatch(match.arg(db, c("GO", "KEGG"), several.ok = TRUE),
                 error = function(e) NULL)
  if (is.null(db) || !length(db)) {
    stop(errorCondition("db doit être 'GO', 'KEGG' ou les deux (v1 — compareCluster et les autres bases sont explicitement hors portée).",
                        class = "pathway_error"))
  }
  organism <- tryCatch(match.arg(organism, c("human", "mouse")),
                       error = function(e) NULL)
  if (is.null(organism)) {
    stop(errorCondition("organism doit être 'human' ou 'mouse' (v1 — l'orthologie souris élargie reste une décision JOINT GEMS D5/PACKS).",
                        class = "pathway_error"))
  }
  fdr_cutoff <- suppressWarnings(as.numeric(fdr_cutoff))
  if (length(fdr_cutoff) != 1L || is.na(fdr_cutoff) ||
      fdr_cutoff <= 0 || fdr_cutoff >= 1) {
    stop(errorCondition("fdr_cutoff doit être un nombre dans ]0, 1[.",
                        class = "pathway_error"))
  }

  # -- 1. Validation de la table DE (avant tout paquet) ------------------------
  if (!is.data.frame(de_table) || nrow(de_table) == 0L) {
    stop(errorCondition("de_table vide ou non data.frame — rien à tester directionnellement.",
                        class = "pathway_error"))
  }
  required <- c("gene", "log2FoldChange", "pvalue")
  missing_cols <- setdiff(required, names(de_table))
  if (length(missing_cols)) {
    stop(errorCondition(paste0(
      "Colonnes DE manquantes pour l'ORA directionnelle : ",
      paste(missing_cols, collapse = ", "),
      ". La p-valeur BRUTE est requise — passer padj produirait un double BH chez goana/kegga."),
      class = "pathway_error"))
  }

  # -- 2. Fond obligatoire (doctrine M12 au 3e idiome) -------------------------
  universe_ids <- strip_composite_gene_ids(as.character(universe_ids))
  universe_ids <- unique(universe_ids[!is.na(universe_ids) & nzchar(universe_ids)])
  if (!length(universe_ids)) {
    stop(errorCondition("Fond (universe_ids) vide — l'ORA directionnelle refuse le fond par défaut de limma : fournissez les rownames de la matrice filtrée.",
                        class = "pathway_error"))
  }

  # -- 3. Mapping SYMBOL -> Entrez (map_fun injecté, même idiome que M12) ------
  species <- switch(organism, human = "Hs", mouse = "Mm")
  orgpkg <- switch(organism, human = "org.Hs.eg.db", mouse = "org.Mm.eg.db")
  if (is.null(map_fun)) {
    if (!requireNamespace(orgpkg, quietly = TRUE)) {
      stop(errorCondition(paste0("Package '", orgpkg,
                                  "' requis pour l'analyse ", organism, "."),
                          class = "pathway_error"))
    }
    map_fun <- function(keys) clusterProfiler::bitr(
      keys, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgpkg)
  }

  de_ids <- strip_composite_gene_ids(as.character(de_table$gene))
  de_map <- tryCatch(map_fun(unique(de_ids)), error = function(e) NULL)
  if (is.null(de_map) || !nrow(de_map)) {
    stop(errorCondition("Aucun gène DE n'a pu être converti en Entrez — l'ORA directionnelle refuse de continuer sur un ensemble vide.",
                        class = "pathway_error"))
  }
  universe_res <- .resolve_ora_universe(universe_ids, de_map$ENTREZID, map_fun)
  universe_entrez <- universe_res$entrez
  universe_note <- universe_res$note
  if (is.null(universe_entrez) || !length(universe_entrez)) {
    stop(errorCondition("Le fond déclaré ne mappe vers aucun Entrez — refus (jamais de repli silencieux sur le fond par défaut).",
                        class = "pathway_error"))
  }
  # Note PARTIELLE (au-delà des 3 cas du helper M12) : si une partie du fond
  # déclaré ne convertit pas en Entrez, le consommateur doit le savoir sans
  # lire la console — le fond effectif est plus petit que le fond déclaré.
  uni_map_probe <- tryCatch(map_fun(universe_ids), error = function(e) NULL)
  n_uni_declared <- length(universe_ids)
  n_uni_mapped <- if (!is.null(uni_map_probe)) length(unique(uni_map_probe$ENTREZID)) else 0L
  if (n_uni_mapped < n_uni_declared) {
    partial_note <- paste0(
      "Fond ORA directionnelle : ", n_uni_mapped, "/", n_uni_declared,
      " identifiants declares convertis en Entrez — les non-convertibles sont exclus du fond effectif.")
    universe_note <- if (is.null(universe_note)) partial_note
                     else paste(universe_note, partial_note, sep = " | ")
  }

  # -- 4. Déduplication Entrez (data.frame rownames uniques obligatoires) ------
  # Deux symboles peuvent converger vers un même Entrez : on garde, par
  # Entrez, la ligne la plus significative (min pvalue) — politique honnête
  # et déterministe, estampillée dans meta.
  de_map <- unique(de_map[, c("SYMBOL", "ENTREZID")])
  idx <- match(de_map$SYMBOL, de_ids)
  de_sel <- de_table[idx, , drop = FALSE]
  de_sel$.entrez <- de_map$ENTREZID
  de_sel <- de_sel[order(de_sel$.entrez, de_sel$pvalue), ]
  keep <- !duplicated(de_sel$.entrez)
  de_sel <- de_sel[keep, ]
  de_n_entrez <- nrow(de_sel)

  # -- 5. Pré-compte honnête du sous-ensemble DE (avant tout appel limma) ------
  bh <- p.adjust(de_sel$pvalue, method = "BH")
  n_de <- sum(bh < fdr_cutoff)
  if (n_de == 0L) {
    stop(errorCondition(paste0(
      "Zéro gène DE sous le seuil fdr_cutoff=", fdr_cutoff,
      " (BH interne = le même padj que le module) — l'ORA directionnelle refuse de produire une table vide muette."),
      class = "pathway_error"))
  }

  # -- 6. goana / kegga (API limma moderne : LISTE de vecteurs Entrez par
  #      direction — le data.frame logFC/P.Value de l'ancien goana n'existe
  #      plus dans limma 3.62 au lock, mesuré) -------------------------------
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop(errorCondition("Package 'limma' requis (moteur DE bulk — absent de la bibliothèque ?).",
                        class = "pathway_error"))
  }
  sig <- bh < fdr_cutoff
  de_up <- as.character(de_sel$.entrez)[sig & de_sel$log2FoldChange > 0]
  de_down <- as.character(de_sel$.entrez)[sig & de_sel$log2FoldChange < 0]
  direction_list <- list()
  if (length(de_up)) direction_list[["Up"]] <- de_up
  if (length(de_down)) direction_list[["Down"]] <- de_down
  if (!length(direction_list)) {
    stop(errorCondition("Les gènes DE au seuil portent un logFC nul — aucune direction testable.",
                        class = "pathway_error"))
  }

  .norm_direction_table <- function(t, db_name) {
    # limma rend P.Up/P.Down (P.<name>) sans FDR : le BH par direction est
    # appliqué ICI, sur les colonnes de la table (même politique que le
    # module pour l'ORA classique — qval exposé, jamais muet).
    for (nm in names(direction_list)) {
      pcol <- paste0("P.", nm)
      fcol <- paste0("FDR.", nm)
      if (pcol %in% names(t)) t[[fcol]] <- p.adjust(t[[pcol]], method = "BH")
    }
    t$db <- db_name
    t
  }

  out <- list(go = NULL, kegg = NULL)
  if ("GO" %in% db) {
    t <- limma::goana(direction_list, universe = universe_entrez,
                      species = species)
    out$go <- .norm_direction_table(t, "GO")
  }
  if ("KEGG" %in% db) {
    t <- limma::kegga(direction_list, universe = universe_entrez,
                      species = species)
    out$kegg <- .norm_direction_table(t, "KEGG")
  }

  list(go = out$go, kegg = out$kegg,
       meta = list(de_n_entrez = de_n_entrez,
                   universe_n_entrez = length(universe_entrez),
                   universe_note = universe_note,
                   n_de_at_cutoff = n_de,
                   directions = names(direction_list),
                   n_up = length(de_up),
                   n_down = length(de_down),
                   fdr_cutoff = fdr_cutoff,
                   organism = organism,
                   db = db))
}

#' Barcodeplot limma — rang logFC avec la voie surlignée
#'
#' @param rank_stats vecteur numérique NOMMÉ (names = mêmes identifiants que
#'   `index_genes`, typiquement Entrez) — le classement complet (logFC).
#' @param index_genes identifiants des gènes DE la voie à surligner.
#' @param ... passé à limma::barcodeplot (quantiles, labels, worm...).
#' @return la sortie de limma::barcodeplot (tracé base R, invisible) — à
#'   rendre via renderPlot côté module.
plot_pathway_barcode <- function(rank_stats, index_genes, ...) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop(errorCondition("Package 'limma' requis (moteur DE bulk — absent de la bibliothèque ?).",
                        class = "pathway_error"))
  }
  rank_stats <- rank_stats[!is.na(rank_stats)]
  if (!length(rank_stats)) {
    stop(errorCondition("Classement vide — rien à tracer pour le barcodeplot.",
                        class = "pathway_error"))
  }
  idx <- which(names(rank_stats) %in% as.character(index_genes))
  if (!length(idx)) {
    stop(errorCondition("Aucun gène de la voie n'apparaît dans le classement — le barcodeplot refuse de tracer un surlignage vide.",
                        class = "pathway_error"))
  }
  limma::barcodeplot(rank_stats, index = idx, ...)
}


# =============================================================================
# PB v2 (GEMS D5, GO opérateur 2026-11-10) — COMPARAISON INTER-CONDITIONS :
# compareCluster, le « 3e idiome » complété.
# =============================================================================
# PB v1 (goana/kegga) compare les DIRECTIONS Up/Down PAR VOIE sur UNE table
# DE ; compareCluster compare les ENRICHISSEMENTS ENTRE LISTES DE GÈNES
# (conditions/tissus/contrastes) — fiche §PB : « compareCluster après »,
# GSE164073 3 tissus. MÊME DISCIPLINE que PB v1 (commit 933ec9b) : fond
# OBLIGATOIRE (doctrine M12), règle d'identifiants composites B1-C, classe
# pathway_error, map_fun injecté pour la pureté testable. Moteur d'abord —
# le branchement module suivra en jalon 2 (précédent PB v1 : moteur
# 933ec9b → branchement 68172f9).
# =============================================================================

#' Comparaison d'enrichissement inter-conditions (compareCluster, PB v2)
#'
#' Compare les enrichissements de PLUSIEURS listes de gènes (une liste par
#' condition/tissu/contraste) dans un MÊME fond —
#' `clusterProfiler::compareCluster`, le 3e idiome complété après l'ORA
#' directionnelle PB v1. Gardes : >= 2 listes nommées ; listes trop petites
#' DÉPOSÉES avec note (jamais testées muettement) ; fond OBLIGATOIRE (M12) ;
#' règle composites B1-C ; map_fun injecté (motif M12, pureté testable).
#'
#' @param gene_lists liste NOMMÉE de vecteurs caractère (SYMBOLs ou
#'   composites ENSG..._SYMBOL).
#' @param universe_ids fond obligatoire partagé par toutes les listes.
#' @param db "GO", "KEGG" ou les deux — bases calculées indépendamment ;
#'   une base vide reste NULL avec note.
#' @param organism "human" ou "mouse".
#' @param ont "BP" (défaut), "CC" ou "MF".
#' @param pvalue_cutoff,qvalue_cutoff seuils passés à compareCluster.
#' @param min_list_size taille mini d'une liste mappée en Entrez.
#' @param map_fun injection du mapping SYMBOL -> Entrez.
#'
#' @return list(go =, kegg =, meta =) — data.frame compareCluster + colonne
#'   `db` ; NULL si la base n'a rien produit ; meta (comptes par liste,
#'   notes M12, seuils).
run_pathway_compare_cluster <- function(gene_lists, universe_ids,
                                        db = c("GO", "KEGG"),
                                        organism = "human",
                                        ont = "BP",
                                        pvalue_cutoff = 0.05,
                                        qvalue_cutoff = 0.2,
                                        min_list_size = 5L,
                                        map_fun = NULL) {
  db <- tryCatch(match.arg(db, c("GO", "KEGG"), several.ok = TRUE),
                 error = function(e) NULL)
  if (is.null(db) || !length(db)) {
    stop(errorCondition("db doit être 'GO', 'KEGG' ou les deux.",
                        class = "pathway_error"))
  }
  organism <- tryCatch(match.arg(organism, c("human", "mouse")),
                       error = function(e) NULL)
  if (is.null(organism)) {
    stop(errorCondition("organism doit être 'human' ou 'mouse'.",
                        class = "pathway_error"))
  }
  if (!is.list(gene_lists) || is.null(names(gene_lists)) ||
      any(!nzchar(names(gene_lists)))) {
    stop(errorCondition("gene_lists doit être une liste NOMMÉE de vecteurs de gènes (une entrée par condition).",
                        class = "pathway_error"))
  }
  if (length(gene_lists) < 2L) {
    stop(errorCondition("compareCluster compare PLUSIEURS listes : fournissez >= 2 conditions (une seule liste = enrichGO/enrichKEGG simple, voie déjà livrée).",
                        class = "pathway_error"))
  }
  ont <- tryCatch(match.arg(ont, c("BP", "CC", "MF")),
                  error = function(e) NULL)
  if (is.null(ont)) stop(errorCondition("ont doit être 'BP', 'CC' ou 'MF'.", class = "pathway_error"))
  pvalue_cutoff <- suppressWarnings(as.numeric(pvalue_cutoff)[1])
  qvalue_cutoff <- suppressWarnings(as.numeric(qvalue_cutoff)[1])
  if (is.na(pvalue_cutoff) || pvalue_cutoff <= 0 || pvalue_cutoff >= 1 ||
      is.na(qvalue_cutoff) || qvalue_cutoff <= 0 || qvalue_cutoff > 1) {
    stop(errorCondition("pvalue_cutoff (]0,1[) et qvalue_cutoff (]0,1]) invalides.",
                        class = "pathway_error"))
  }

  # -- 1. Fond obligatoire (doctrine M12, identique PB v1) ---------------------
  universe_ids <- strip_composite_gene_ids(as.character(universe_ids))
  universe_ids <- unique(universe_ids[!is.na(universe_ids) & nzchar(universe_ids)])
  if (!length(universe_ids)) {
    stop(errorCondition("Fond (universe_ids) vide — compareCluster refuse le fond par défaut : fournissez les rownames de la matrice filtrée.",
                        class = "pathway_error"))
  }

  # -- 2. Mapping injecté (même idiome M12 que PB v1) ---------------------------
  orgpkg <- switch(organism, human = "org.Hs.eg.db", mouse = "org.Mm.eg.db")
  if (is.null(map_fun)) {
    if (!requireNamespace(orgpkg, quietly = TRUE)) {
      stop(errorCondition(paste0("Package '", orgpkg,
                                 "' requis pour l'analyse ", organism, "."),
                          class = "pathway_error"))
    }
    map_fun <- function(keys) clusterProfiler::bitr(
      keys, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgpkg)
  }
  orgdb <- getExportedValue(orgpkg, orgpkg)

  # -- 3. Mapping PAR LISTE + dépôt honnête des listes trop petites ------------
  per_list <- list()
  dropped <- character(0)
  lists_entrez <- list()
  for (nm in names(gene_lists)) {
    ids <- strip_composite_gene_ids(as.character(gene_lists[[nm]]))
    ids <- unique(ids[!is.na(ids) & nzchar(ids)])
    n_declared <- length(ids)
    mp <- if (n_declared) tryCatch(map_fun(ids), error = function(e) NULL) else NULL
    n_mapped <- if (!is.null(mp)) length(unique(mp$ENTREZID)) else 0L
    if (n_mapped < min_list_size) {
      dropped <- c(dropped, nm)
      per_list[[nm]] <- list(declared = n_declared, entrez = n_mapped, kept = FALSE)
      next
    }
    lists_entrez[[nm]] <- unique(mp$ENTREZID)
    per_list[[nm]] <- list(declared = n_declared, entrez = n_mapped, kept = TRUE)
  }
  notes <- character(0)
  if (length(dropped)) {
    notes <- c(notes, sprintf(
      "%d liste(s) déposée(s) sous les %d gènes Entrez min. (%s) — déposées, jamais testées muettement.",
      length(dropped), min_list_size, paste(dropped, collapse = ", ")))
  }
  if (length(lists_entrez) < 2L) {
    stop(errorCondition(paste0(
      "Moins de 2 listes valides après mapping (min_list_size=", min_list_size,
      ") — compareCluster refuse de tourner sur une comparaison vide. Comptes : ",
      paste(vapply(per_list, function(x) paste0(x$declared, "->", x$entrez), ""),
            collapse = ", ")),
      class = "pathway_error"))
  }

  # -- 4. Fond mappé + notes (réutilise le résolveur M12 + sonde partielle) -----
  universe_res <- .resolve_ora_universe(universe_ids, unlist(lists_entrez), map_fun)
  universe_entrez <- universe_res$entrez
  universe_note <- universe_res$note
  if (is.null(universe_entrez) || !length(universe_entrez)) {
    stop(errorCondition("Le fond déclaré ne mappe vers aucun Entrez — refus (jamais de repli silencieux sur le fond par défaut).",
                        class = "pathway_error"))
  }
  uni_probe <- tryCatch(map_fun(unique(universe_ids)), error = function(e) NULL)
  n_uni_mapped <- if (!is.null(uni_probe)) length(unique(uni_probe$ENTREZID)) else 0L
  if (n_uni_mapped < length(universe_ids)) {
    notes <- c(notes, sprintf(
      "Fond : %d/%d identifiants déclarés convertis en Entrez — les non-convertibles sont exclus du fond effectif.",
      n_uni_mapped, length(universe_ids)))
  }
  if (!is.null(universe_note)) notes <- c(notes, universe_note)

  # -- 5. compareCluster par base (GO offline via OrgDb ; KEGG via enrichKEGG) --
  .run_one <- function(fun, extra) {
    res <- tryCatch(do.call(clusterProfiler::compareCluster, c(
      list(geneClusters = lists_entrez, fun = fun,
           universe = universe_entrez,
           pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff),
      extra)), error = function(e) e)
    if (inherits(res, "error")) {
      notes <<- c(notes, sprintf("compareCluster %s : %s", fun, conditionMessage(res)))
      return(NULL)
    }
    df <- tryCatch(as.data.frame(res), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0L) {
      notes <<- c(notes, sprintf("compareCluster %s : 0 voie sous les seuils.", fun))
      return(NULL)
    }
    df
  }
  out <- list(go = NULL, kegg = NULL)
  if ("GO" %in% db) {
    df <- .run_one("enrichGO", list(OrgDb = orgdb, keyType = "ENTREZID", ont = ont))
    if (!is.null(df)) out$go <- df
  }
  if ("KEGG" %in% db) {
    kegg_code <- switch(organism, human = "hsa", mouse = "mmu")
    df <- .run_one("enrichKEGG", list(organism = kegg_code))
    if (!is.null(df)) out$kegg <- df
  }
  # colonne db = nom affiché (MESURE : la boucle interne nommait "go"/"kegg")
  for (nm in c("go", "kegg")) if (!is.null(out[[nm]])) {
    out[[nm]]$db <- if (nm == "go") "GO" else "KEGG"
  }

  list(go = out$go, kegg = out$kegg,
       meta = list(organism = organism, ont = ont,
                   pvalue_cutoff = pvalue_cutoff, qvalue_cutoff = qvalue_cutoff,
                   universe_n_entrez = length(universe_entrez),
                   per_list = per_list, dropped = dropped,
                   notes = notes))
}

# ═════════════════════════════════════════════════════════════════════════════
# PACKS — ORA locale sur les GMT de packs (PK-2, pilier 2 local)
# fiche docs/proposals/PROPOSAL_PACKS_ECOSYSTEM_2026-11-09.md §3 jalon PK-2 ·
# détail docs/ideation/PACKS_ECOSYSTEM_ENRICHR_PALETTES.md §4.1/§4.3.
#
# Contrat (fiche §4.1) : parsing GMT → TERM2GENE → clusterProfiler::enricher()
# — le moteur existant fait tout : même objet enrichResult, même contrat de
# table (ID/Description/pvalue/padj/GeneRatio/Count), mêmes exports/plots
# (cnetplot/emapplot consomment l'attribut `enrich_obj`, motif STAT-S2).
# Espace SYMBOLES : les librairies Enrichr sont des symboles humains — AUCUN
# bitr SYMBOL→Entrez (zéro perte de gènes), donc l'univers (gènes détectés de
# l'expérience, fiche §5 étape 1) est consommé tel quel. Zéro nouvelle
# dépendance (parsing GMT = base R ; clusterProfiler déjà au lock).
# Orthologie souris = D4, décision JOINT GEMS D5 — refus honnête côté modules,
# le moteur reste organisme-free.
# ═════════════════════════════════════════════════════════════════════════════

#' Lecture pure d'un fichier GMT (librairie de gènes de pack)
#'
#' Format standard vérifié (fiche packs §2.7) : une ligne par terme,
#' `terme \t description \t gene1 \t gene2 …` — la description peut être vide
#' (mode texte Enrichr). Gènes uniques et nettoyés par terme.
#'
#' @param gmt_path Chemin du fichier .gmt (chaîne unique).
#' @return list(term2gene = data.frame(term, gene),
#'              term2name = data.frame(term, name) — name = description ou
#'              terme si vide).
#' @refusals pathway_error : chemin malformé, fichier absent, fichier vide,
#'   ligne < 3 champs (la ligne est nommée), terme vide, terme sans gène,
#'   terme dupliqué dans le même fichier.
packs_read_gmt <- function(gmt_path) {
  if (missing(gmt_path) || !is.character(gmt_path) || length(gmt_path) != 1L ||
      is.na(gmt_path) || !nzchar(trimws(gmt_path))) {
    stop(errorCondition(
      "packs_read_gmt : chemin de fichier GMT attendu (chaîne unique non vide).",
      class = "pathway_error"))
  }
  if (!file.exists(gmt_path)) {
    stop(errorCondition(sprintf("Fichier GMT introuvable : « %s ».", gmt_path),
                        class = "pathway_error"))
  }
  lines <- readLines(gmt_path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) == 0L) {
    stop(errorCondition(sprintf("Fichier GMT vide : « %s ».", gmt_path),
                        class = "pathway_error"))
  }
  seen <- new.env(parent = emptyenv())
  t2g_parts <- vector("list", length(lines))
  t2n_parts <- vector("list", length(lines))
  for (i in seq_along(lines)) {
    fields <- strsplit(lines[[i]], "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3L) {
      stop(errorCondition(sprintf(
        "GMT malformé « %s » ligne %d : attendus >= 3 champs tab-séparés (terme, description, gènes…), reçus %d.",
        gmt_path, i, length(fields)), class = "pathway_error"))
    }
    term <- trimws(fields[[1]])
    name_raw <- trimws(fields[[2]])
    genes <- unique(trimws(fields[-c(1L, 2L)]))
    genes <- genes[nzchar(genes)]
    if (!nzchar(term)) {
      stop(errorCondition(sprintf(
        "GMT malformé « %s » ligne %d : terme vide.", gmt_path, i),
        class = "pathway_error"))
    }
    if (length(genes) == 0L) {
      stop(errorCondition(sprintf(
        "GMT malformé « %s » ligne %d : le terme « %s » n'a aucun gène.",
        gmt_path, i, term), class = "pathway_error"))
    }
    if (!is.null(seen[[term]])) {
      stop(errorCondition(sprintf(
        "GMT malformé « %s » ligne %d : terme dupliqué « %s » — une librairie GMT doit avoir des termes uniques.",
        gmt_path, i, term), class = "pathway_error"))
    }
    seen[[term]] <- TRUE
    t2g_parts[[i]] <- data.frame(term = term, gene = genes,
                                 stringsAsFactors = FALSE)
    t2n_parts[[i]] <- data.frame(term = term,
                                 name = if (nzchar(name_raw)) name_raw else term,
                                 stringsAsFactors = FALSE)
  }
  list(term2gene = do.call(rbind, t2g_parts),
       term2name = do.call(rbind, t2n_parts))
}

#' Librairies disponibles dans le registre des packs (choix de sélection)
#'
#' @param packs Registre session (named list de biomentat_pack — scan PK-1) ;
#'   NULL/vide tolérés. Les packs dont le type n'est pas « libraries » sont
#'   ignorés (les palettes ne sont pas des librairies de gènes).
#' @return named character : valeur = « <pack_id>||<basename> » (clé de
#'   sélection stable), label = « <pack_name> — <fichier sans extension> ».
#'   Vecteur vide quand aucun pack de librairies.
packs_library_choices <- function(packs) {
  if (is.null(packs) || !is.list(packs) || length(packs) == 0L) return(character(0))
  vals <- character(0)
  # itération sur les VALEURS : le registre session est nommé (PK-1), mais un
  # registre anonyme (fixture de test, appel programmatique) reste fonctionnel.
  for (e in packs) {
    if (!is.list(e) || !identical(e$type, "libraries")) next
    files <- e$files
    if (!is.character(files) || length(files) == 0L) next
    pid <- if (is.character(e$id) && length(e$id) == 1L && nzchar(e$id)) e$id else "pack"
    for (f in files) {
      key <- paste0(pid, "||", basename(f))
      stem <- sub("\\.gmt$", "", basename(f), ignore.case = TRUE)
      label <- sprintf("%s — %s", e$name %||% pid, stem)
      vals[[key]] <- label
    }
  }
  vals
}

#' TERM2GENE/TERM2NAME combinés pour les librairies sélectionnées
#'
#' Les identifiants de terme sont préfixés `<pack_id>::<fichier>::<terme>` :
#' un même nom de terme dans deux librairies ne collisionne jamais, et la
#' provenance (pack + librairie) reste visible dans la table rendue et les
#' exports — la Description brute du GMT passe par TERM2NAME.
#'
#' @param packs Registre session (cf. packs_library_choices).
#' @param selected Character : clés de sélection telles que émises par
#'   packs_library_choices (« <pack_id>||<basename> »).
#' @return list(term2gene = data.frame(term, gene),
#'              term2name = data.frame(term, name)).
#' @refusals pathway_error : sélection vide/NULL, clé inconnue (les clés
#'   valides sont listées), fichier GMT déclaré par le pack absent du disque.
packs_term2gene <- function(packs, selected) {
  if (missing(selected) || is.null(selected) || !is.character(selected) ||
      length(selected) == 0L || all(!nzchar(trimws(selected)))) {
    stop(errorCondition(
      "Packs : sélectionnez au moins une librairie (sélection vide).",
      class = "pathway_error"))
  }
  valid <- packs_library_choices(packs)
  bad <- setdiff(trimws(selected), names(valid))
  if (length(bad) > 0L) {
    stop(errorCondition(sprintf(
      "Packs : librairie(s) inconnue(s) : « %s » — clés valides : %s. Ré-ouvrez le sélecteur (le registre a peut-être changé).",
      paste(bad, collapse = ", "),
      if (length(valid) == 0L) "(aucune — aucun pack de librairies installé)" else paste(names(valid), collapse = ", ")),
      class = "pathway_error"))
  }
  t2g_parts <- vector("list", length(selected))
  t2n_parts <- vector("list", length(selected))
  # recherche du pack PAR ID (les ids PK-1 sont uniques au scan ; le registre
  # peut être nommé ou non — l'indexation par nom ne suffit pas)
  .find_pack <- function(pid) {
    for (e in packs) if (is.list(e) && identical(e$id, pid)) return(e)
    NULL
  }
  for (i in seq_along(selected)) {
    key <- trimws(selected[[i]])
    pack_id <- sub("\\|\\|.*$", "", key)
    fname <- sub("^.*\\|\\|", "", key)
    e <- .find_pack(pack_id)
    if (is.null(e)) {
      stop(errorCondition(sprintf(
        "Packs : le pack « %s » est introuvable dans le registre session — ré-ouvrez le sélecteur (le registre a peut-être changé).",
        pack_id), class = "pathway_error"))
    }
    gmt_path <- file.path(e$dir, fname)
    if (!file.exists(gmt_path)) {
      stop(errorCondition(sprintf(
        "Packs : le fichier GMT « %s » déclaré par le pack « %s » est introuvable sur le disque — le dossier a peut-être bougé depuis le scan.",
        gmt_path, pack_id), class = "pathway_error"))
    }
    lib <- packs_read_gmt(gmt_path)
    stem <- sub("\\.gmt$", "", fname, ignore.case = TRUE)
    prefix <- paste0(pack_id, "::", stem, "::")
    t2g <- lib$term2gene
    t2g$term <- paste0(prefix, t2g$term)
    t2n <- lib$term2name
    t2n$term <- paste0(prefix, t2n$term)
    t2g_parts[[i]] <- t2g
    t2n_parts[[i]] <- t2n
  }
  list(term2gene = do.call(rbind, t2g_parts),
       term2name = do.call(rbind, t2n_parts))
}

#' ORA locale sur les librairies de packs — enricher(), contrat inchangé
#'
#' Même famille de sorties que run_pathway_enrichment() : data.frame
#' (ID/Description/GeneRatio/BgRatio/pvalue/p.adjust/qvalue/geneID/Count) ou
#' le data.frame vide au même contrat si aucun résultat, et l'objet
#' enrichResult brut attaché en attribut `enrich_obj` (STAT-S2 — les vues
#' réseau enrichplot cnetplot/emapplot le consomment). Espace SYMBOLES :
#' pas de conversion Entrez — l'univers déclaré (gènes détectés de
#' l'expérience) est consommé tel quel, en éliminant le biais de background
#' d'Enrichr (fiche §4.1 « bonus rigueur »).
#'
#' @param genes Character : gènes à tester (symboles ; identifiants composites
#'   ENSG…_SYMBOL acceptés, même règle B1-C que le moteur).
#' @param term2gene data.frame(term, gene) — cf. packs_term2gene().
#' @param term2name data.frame(term, name) optionnel (descriptions GMT).
#' @param universe Character optionnel : gènes mesurés de l'expérience.
#' @param pval_cutoff,qval_cutoff,p_adjust_method Identiques au moteur
#'   (validations miroir, avant tout chargement de paquet).
#' @return data.frame enrichi (ou vide au même contrat) + attribut enrich_obj.
#' @refusals pathway_error : TERM2GENE malformé, gènes vides après nettoyage,
#'   qval hors ]0;1], méthode d'ajustement inconnue.
run_packs_enrichment <- function(genes, term2gene, term2name = NULL,
                                 universe = NULL, pval_cutoff = 0.05,
                                 qval_cutoff = 0.2, p_adjust_method = "BH") {
  p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)
  qval_cutoff <- suppressWarnings(as.numeric(qval_cutoff))
  if (length(qval_cutoff) != 1L || is.na(qval_cutoff) ||
      qval_cutoff <= 0 || qval_cutoff > 1) {
    stop(errorCondition("qval_cutoff doit être un nombre dans ]0, 1].",
                        class = "pathway_error"))
  }
  if (missing(term2gene) || !is.data.frame(term2gene) ||
      !all(c("term", "gene") %in% names(term2gene)) ||
      nrow(term2gene) == 0L) {
    stop(errorCondition(
      "Packs : TERM2GENE attendu (data.frame non vide, colonnes term/gene) — cf. packs_term2gene().",
      class = "pathway_error"))
  }
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop(errorCondition(
      "Package 'clusterProfiler' requis. Installez-le via BiocManager.",
      class = "pathway_error"))
  }
  library(clusterProfiler)
  genes <- strip_composite_gene_ids(genes)
  genes <- unique(trimws(genes[!is.na(genes) & nchar(trimws(genes)) > 0]))
  if (length(genes) == 0L) {
    stop(errorCondition(
      "Packs : aucun gène valide fourni après nettoyage.", class = "pathway_error"))
  }
  universe_used <- NULL
  if (!is.null(universe)) {
    universe_used <- strip_composite_gene_ids(universe)
    universe_used <- unique(trimws(universe_used[!is.na(universe_used) &
                                                   nchar(trimws(universe_used)) > 0]))
    if (length(universe_used) == 0L) universe_used <- NULL
  }
  t2g <- term2gene[, c("term", "gene"), drop = FALSE]
  t2g$term <- trimws(as.character(t2g$term))
  t2g$gene <- trimws(as.character(t2g$gene))
  t2g <- t2g[nzchar(t2g$term) & nzchar(t2g$gene), , drop = FALSE]
  if (nrow(t2g) == 0L) {
    stop(errorCondition(
      "Packs : TERM2GENE vide après nettoyage (termes/gènes vides).",
      class = "pathway_error"))
  }
  t2n <- NULL
  if (!is.null(term2name) && is.data.frame(term2name) &&
      all(c("term", "name") %in% names(term2name)) && nrow(term2name) > 0L) {
    t2n <- term2name[, c("term", "name"), drop = FALSE]
    t2n$term <- trimws(as.character(t2n$term))
    t2n$name <- trimws(as.character(t2n$name))
    t2n <- t2n[nzchar(t2n$term) & nzchar(t2n$name), , drop = FALSE]
    if (nrow(t2n) == 0L) t2n <- NULL
  }
  enrich_result <- enricher(
    gene          = genes,
    TERM2GENE     = t2g,
    TERM2NAME     = t2n,
    pvalueCutoff  = pval_cutoff,
    qvalueCutoff  = qval_cutoff,
    pAdjustMethod = p_adjust_method,
    universe      = universe_used
  )
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0L) {
    return(data.frame(
      ID          = character(0),
      Description = character(0),
      p.adjust    = numeric(0),
      Count       = integer(0),
      GeneRatio   = character(0)
    ))
  }
  res_df <- as.data.frame(enrich_result)
  # STAT-S2 (additif — motif gsea_obj) : l'objet brut reste porté pour les
  # vues réseau enrichplot (cnetplot/emapplot + pairwise_termsim).
  attr(res_df, "enrich_obj") <- enrich_result
  res_df
}

