# =============================================================================
# R/sc/sc_pipeline.R — SC auto-pipeline execution & R-script export
# (extracted from mod_sc.R, Block 8 refactor)
# =============================================================================
# Pure logic: called from mod_sc_server() inside observeEvent().
# Depends on: helpers_sc.R (resolve_sketch_preset, robust_find_clusters,
#             standardize_sketch_reductions, subsample_seurat_for_analysis,
#             remap_seurat_ids_to_symbol), helpers_pathway.R,
#             R/sc/sc_trajectory.R
# =============================================================================

# ── M-2 phase 1 (2026-09-30) : normalisation des paramètres + cœur QC pur ──

#' Paramètres normalisés de l'autopipeline — LE point de passage unique
#'
#' Prend une liste brute (UI convertie, scénario drive, appel direct) et rend
#' un jeu FERMÉ de 28 clés : défauts canoniques appliqués aux clés absentes,
#' clés inconnues perdues. La classe « une clé manquante change silencieusement
#' un seuil » (audit C-1) est fermée STRUCTURELLEMENT : le cœur ne voit jamais
#' de NULL. Les défauts vivent ICI (couche pure) — le canal drive délègue.
#'
#' @param p Named list (ou vide). Les reactiveValues doivent être convertis
#'   par l'appelant (`shiny::reactiveValuesToList(input)`).
#' @return Named list de 28 clés `sc_ap_*` (M4 2026-10-26 : +4 clés doublets/ambiant).
sc_ap_normalize_params <- function(p = list()) {
  d <- list(
    sc_ap_mapping              = FALSE,
    sc_ap_mapping_org          = "human",
    sc_ap_bpcells              = FALSE,
    sc_ap_min_gene             = TS_SC_QC_MIN_GENES,
    sc_ap_max_gene             = TS_SC_QC_MAX_GENES,
    sc_ap_mt                   = TS_SC_QC_MAX_PCT_MT,
    sc_ap_norm                 = "log",
    sc_ap_pca_dim              = TS_SC_QC_PCA_DIMS,
    sc_ap_res                  = 0.5,
    sc_ap_cluster_algo         = "1",
    sc_ap_integration          = "none",
    sc_ap_batch_var            = "orig.ident",
    sc_ap_compute_umap         = TRUE,
    sc_ap_sketch_preset        = "max",
    sc_ap_sketch_ncells_custom = NA,
    sc_ap_singler              = FALSE,
    sc_ap_singler_ref          = "hpca",
    sc_ap_singler_level        = "main",
    sc_ap_markers              = TRUE,
    sc_ap_pathway              = FALSE,
    sc_ap_pathway_db           = "GOBP",
    sc_ap_pathway_org          = "human",
    sc_ap_correlation          = FALSE,
    sc_ap_trajectory           = TRUE,
    # M4 (2026-10-26, jalon BIOMENTAT) : QC doublets & ambiant. Défauts = source
    # unique defaults.R (convention QW-2). sc_ap_doublets est OPT-IN (FALSE) :
    # l'appel moteur est PARESSEUX (Step 1c) et le cœur reste exécutable sans
    # scDblFinder ; le canal drive hérite du même jeu fermé (délégation).
    sc_ap_doublets             = FALSE,
    sc_ap_run_ambient          = FALSE,
    sc_ap_doublet_score_cutoff = TS_SC_QC_DOUBLET_SCORE_CUTOFF,
    sc_ap_max_contamination    = TS_SC_QC_MAX_CONTAMINATION
  )
  for (k in names(d)) d[[k]] <- p[[k]] %||% d[[k]]
  d
}

.sc_tr <- function(msg, ...) {
  tr_fn <- if (exists(".tr", mode = "function", inherits = TRUE)) {
    get(".tr", mode = "function", inherits = TRUE)
  } else {
    NULL
  }
  raw <- if (is.function(tr_fn)) {
    tryCatch(tr_fn(msg), error = function(e) msg)
  } else {
    msg
  }
  if (...length() > 0L) sprintf(raw, ...) else raw
}

#' Cœur PUR de l'étape 1 (QC) de l'autopipeline — zéro Shiny
#'
#' Détection du pattern mito (humain `MT-` / souris `mt-`), calcul de
#' `percent.mt`, filtrage par seuils, garde de survie. Utilisable hors app
#' (tests, CLI) — c'est le premier pas de l'extraction du cœur pur (M-2).
#'
#' @param obj Seurat object brut.
#' @param min_genes,max_genes,max_mt Seuils QC (paramètres normalisés).
#' @param min_cells Plancher de survie (défaut 10, comportement historique).
#' @return List `object` (Seurat filtré), `n_before`, `n_removed`.
#' @raises Erreur classée `sc_pipeline_error` (message français) si moins de
#'   `min_cells` cellules survivent.
sc_pipeline_qc <- function(obj, min_genes, max_genes, max_mt, min_cells = 10) {
  mt_pat <- if (any(grepl("^MT-", rownames(obj)))) "^MT-"
            else if (any(grepl("^mt-", rownames(obj)))) "^mt-" else NULL
  obj[["percent.mt"]] <- if (!is.null(mt_pat))
    PercentageFeatureSet(obj, pattern = mt_pat) else 0
  n_before <- ncol(obj)
  # ⚠️ subset() lève sa propre erreur brute (« No cells found ») quand AUCUNE
  # cellule ne survit — convertie ici en erreur classée française (sinon la
  # garde min_cells ci-dessous n'est jamais atteinte ; trouvé par le test).
  obj <- tryCatch(
    subset(obj,
           subset = nFeature_RNA > min_genes &
                    nFeature_RNA < max_genes &
                    percent.mt   < max_mt),
    error = function(e) stop(errorCondition(
      .sc_tr("Seulement %d cellule(s) après QC (départ: %d). Réduisez les seuils.", 0L, n_before),
      class = "sc_pipeline_error")))
  if (ncol(obj) < min_cells) stop(errorCondition(
    .sc_tr("Seulement %d cellule(s) après QC (départ: %d). Réduisez les seuils.", ncol(obj), n_before),
    class = "sc_pipeline_error"))
  list(object = obj, n_before = n_before, n_removed = n_before - ncol(obj))
}

# Graine déclarée des étapes stochastiques (audit 2026-09-27 §1.8) — withr
# restaure l'état du RNG après chaque appel : deux exécutions identiques
# produisent les mêmes embeddings, sans fuiter la graine vers le reste de l'app.
.SC_AP_SEED <- if (exists("TS_SC_PIPELINE_SEED")) TS_SC_PIPELINE_SEED else 989L

.seeded_ap <- function(expr) {
  withr::with_seed(.SC_AP_SEED, force(expr))
}

# ── Jalon filets (2026-10-27) : filet de classification Step 2-5 ─────────────
# Les erreurs de dépendance BRUTES (Seurat/irlba/umap-rgl…) ne franchissent
# plus la frontière du cœur : converties en `sc_pipeline_error`, message
# d'origine préservé in extenso. Une condition DÉJÀ classée par le dépôt
# (`<domaine>_error` / `<domaine>_missing` — doctrine C5, erreurs M4 comme
# sc_ambient_package_missing) est RELANCÉE telle quelle (F20 : jamais de
# re-classement silencieux). Défauts mesurés : T2-QC-5 (FindClusters à haute
# résolution, « invalid first argument ») et T2-QC-6 (irlba, npcs >= cellules).
.ap_is_domain_classed <- function(e) {
  any(grepl("_(error|missing)$", class(e)))
}

.ap_classed_step <- function(expr) {
  tryCatch(
    expr,
    error = function(e) {
      if (.ap_is_domain_classed(e)) stop(e)
      stop(errorCondition(
        .sc_tr("Échec de l'autopipeline SC (normalisation/PCA/clustering) : %s",
               conditionMessage(e)),
        class = "sc_pipeline_error"))
    }
  )
}

#' Pure headless Single-Cell auto-pipeline core
#'
#' @param obj Seurat object (validated via assert_seurat)
#' @param params Named list of parameters (flexible aliases supported)
#' @param log_file Optional path to log file for inter-process progress streaming
#' @return list(obj = obj, metrics = metrics, timings = timings, log = log)
run_sc_auto_pipeline_core <- function(obj, params = list(), log_file = NULL) {
  # 1. Input assertions
  assert_seurat(obj, context = "run_sc_auto_pipeline_core")

  t_start <- Sys.time()
  step_times <- list(
    mapping_s    = 0,
    qc_s         = 0,
    doublets_s   = 0,
    norm_s       = 0,
    pca_s        = 0,
    harmony_s    = 0,
    clustering_s = 0,
    umap_s       = 0,
    tsne_s       = 0,
    markers_s    = 0,
    pathway_s    = 0,
    correlation_s= 0,
    trajectory_s = 0,
    total_s      = 0
  )
  log_lines <- character(0)
  log_fn <- function(msg) {
    line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
    log_lines <<- c(log_lines, line)
    if (!is.null(log_file) && is.character(log_file) && nzchar(log_file)) {
      tryCatch(cat(line, "\n", file = log_file, append = TRUE), error = function(e) NULL)
    }
  }

  # 2. Extract parameters with alias resolution
  p_norm <- sc_ap_normalize_params(params)

  pca_dim      <- as.integer(params$n_dims %||% params$dims %||% params$pca_dim %||% p_norm$sc_ap_pca_dim %||% 30L)
  resolution   <- as.numeric(params$resolution %||% params$res %||% params$clust_res %||% p_norm$sc_ap_res %||% 0.5)
  cluster_algo <- as.character(params$cluster_algo %||% params$algorithm %||% p_norm$sc_ap_cluster_algo %||% "louvain")
  norm_method  <- as.character(params$norm_method %||% params$norm %||% p_norm$sc_ap_norm %||% "log")
  nfeatures    <- as.integer(params$nfeatures %||% params$n_features %||% 2000L)
  assay        <- as.character(params$assay %||% Seurat::DefaultAssay(obj) %||% "RNA")
  integration  <- as.character(params$integration %||% params$sc_ap_integration %||% p_norm$sc_ap_integration %||% "none")
  batch_var    <- as.character(params$batch_var %||% params$sc_ap_batch_var %||% p_norm$sc_ap_batch_var %||% "orig.ident")
  compute_umap <- isTRUE(params$compute_umap %||% p_norm$sc_ap_compute_umap)
  compute_tsne <- if (!is.null(params$compute_tsne)) isTRUE(params$compute_tsne) else compute_umap
  run_qc       <- if (!is.null(params$run_qc)) isTRUE(params$run_qc) else if (!is.null(params$skip_qc)) !isTRUE(params$skip_qc) else TRUE
  seed         <- as.integer(params$seed %||% .SC_AP_SEED %||% 989L)

  n_cells_initial <- ncol(obj)
  n_genes_initial <- nrow(obj)
  qc_removed <- 0L
  qc_snapshot <- NULL

  # ── Step 0: ID Mapping (optional) ──────────────────────────────────────────
  if (isTRUE(p_norm$sc_ap_mapping)) {
    t0 <- Sys.time()
    detected <- tryCatch(detect_gene_id_type(rownames(obj)), error = function(e) "unknown")
    if (detected %in% c("ensembl", "entrez")) {
      log_fn(sprintf("Mapping IDs (%s)...", detected))
      map_res <- tryCatch(
        remap_seurat_ids_to_symbol(obj, from_type = detected,
                                   organism = p_norm$sc_ap_mapping_org %||% "human",
                                   collapse_method = "sum"),
        error = function(e) { log_fn(paste("Mapping ignoré:", e$message)); NULL }
      )
      if (!is.null(map_res)) {
        obj <- map_res$object
        log_fn(sprintf("✓ Mapping : %d gènes finaux", nrow(obj)))
      }
    } else {
      log_fn("Mapping IDs : symboles déjà détectés (ou type inconnu) — ignoré.")
    }
    step_times$mapping_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }

  # ── Step 1: QC ─────────────────────────────────────────────────────────────
  if (isTRUE(run_qc)) {
    t_qc <- Sys.time()
    log_fn("QC...")
    qcr <- sc_pipeline_qc(obj,
                          min_genes = p_norm$sc_ap_min_gene,
                          max_genes = p_norm$sc_ap_max_gene,
                          max_mt    = p_norm$sc_ap_mt,
                          min_cells = as.integer(params$min_cells %||% 10L))
    obj <- qcr$object
    qc_removed <- qcr$n_removed
    qc_snapshot <- list(min_genes = p_norm$sc_ap_min_gene, max_genes = p_norm$sc_ap_max_gene,
                        max_mt = p_norm$sc_ap_mt, n_removed = qc_removed, n_remaining = ncol(obj))
    log_fn(sprintf("✓ QC : %d cellules (retirées: %d)", ncol(obj), qc_removed))
    step_times$qc_s <- as.numeric(difftime(Sys.time(), t_qc, units = "secs"))
  } else {
    log_fn("QC ignoré (matrice pré-filtrée fournie).")
    step_times$qc_s <- 0
  }

  # ── Step 1c (M4, 2026-10-26) : QC doublets & ambiant — appel PARESSEUX ──────
  # scDblFinder est OPTIONNEL : le cœur doit rester exécutable sans le package
  # (log + skip quand le moteur n'est pas résident). Quand le moteur tourne,
  # ses erreurs classées (ex. sc_ambient_package_missing si run_ambient sans
  # SoupX/celda) se propagent TELLES QUELLES — jamais de repli silencieux
  # (feature 20). L'appel est AVANT la conversion BPCells (Step 1b) : la
  # détection exige des counts en mémoire. Provenance + distributions
  # rejoignent la liste de retour (consommables par le rapport).
  m4_provenance <- NULL
  m4_distributions <- NULL
  if (isTRUE(p_norm$sc_ap_doublets)) {
    if (exists("sc_filter_doublets_and_ambient", mode = "function")) {
      t_m4 <- Sys.time()
      log_fn("QC doublets & ambiant (scDblFinder)...")
      m4_res <- sc_filter_doublets_and_ambient(obj, params = list(
        seed = seed,
        batch_var = batch_var,
        doublet_score_cutoff = p_norm$sc_ap_doublet_score_cutoff,
        max_contamination = p_norm$sc_ap_max_contamination,
        run_ambient = isTRUE(p_norm$sc_ap_run_ambient)
      ))
      obj <- m4_res$obj
      m4_provenance <- m4_res$provenance
      m4_distributions <- m4_res$qc_metrics$distributions
      log_fn(sprintf("✓ QC doublets : %d cellules conservées (exclues: %d)",
                     ncol(obj),
                     length(as.character(m4_res$provenance$cells_excluded))))
      step_times$doublets_s <- as.numeric(difftime(Sys.time(), t_m4, units = "secs"))
    } else {
      log_fn("QC doublets ignoré : moteur sc_filter_doublets_and_ambient indisponible (scDblFinder optionnel).")
      step_times$doublets_s <- 0
    }
  }

  # ── Step 1b: Backend disque (BPCells) ──────────────────────────────────────
  if (isTRUE(p_norm$sc_ap_bpcells) && ncol(obj) > 10000L &&
      exists("sc_backend_status", mode = "function") && sc_backend_status(obj) == "memory") {
    if (exists(".bpcells_available", mode = "function") && .bpcells_available() &&
        exists("convert_seurat_to_bpcells", mode = "function")) {
      conv <- tryCatch(convert_seurat_to_bpcells(obj), error = function(e) NULL)
      if (!is.null(conv) && !is.null(conv$object)) {
        obj <- conv$object
        log_fn(sprintf("✓ Backend disque (BPCells) activé — %s cellules", format(ncol(obj), big.mark = " ")))
      }
    }
  }

  # ── Step 2-5: Normalisation / PCA / Clustering / UMAP ───────────────────────
  t_norm <- Sys.time()
  n_total_cells <- ncol(obj)
  sketch_params <- resolve_sketch_preset(
    p_norm$sc_ap_sketch_preset %||% "standard", n_total_cells,
    p_norm$sc_ap_sketch_ncells_custom)
  use_sketch <- !identical(norm_method, "sct") &&
                sketch_params$ncells < n_total_cells &&
                identical(integration, "none")

  # ── Garde de dimension (jalon filets, 2026-10-27) ──────────────────────────
  # T2-QC-6 (M5) : 12 cellules < pca_dim (20) faisaient sortir l'erreur irlba
  # brute (« max(nu, nv) must be strictly less than min(nrow(A), ncol(A)) »).
  # Refus CLASSÉ et actionnable avant tout calcul : les deux knobs nommés.
  if (!isTRUE(use_sketch) && n_total_cells <= pca_dim) {
    stop(errorCondition(
      .sc_tr("Pas assez de cellules (%d) pour %d dimensions PCA (pca_dim). Fournissez plus de cellules ou réduisez pca_dim.",
             n_total_cells, pca_dim),
      class = "sc_pipeline_error"))
  }

  if (exists("sc_backend_status", mode = "function") && sc_backend_status(obj) == "disk") {
    .ap_old_plan <- future::plan()
    on.exit(future::plan(.ap_old_plan), add = TRUE)
    future::plan("sequential")
    log_fn("Backend disque : future séquentiel forcé.")
  }

  clust_dims <- pca_dim
  clust_red <- "pca"
  mixing_before <- NULL
  mixing_after <- NULL

  # Filet de classification autour de TOUTE la branche sketch/complete (le
  # block est évalué dans le frame de l'appelant — même mécanisme que
  # .seeded_ap) : voir .ap_classed_step ci-dessus.
  .ap_classed_step({
  if (isTRUE(use_sketch)) {
    # ── Sketch Branch ────────────────────────────────────────────────────────
    log_fn(sprintf("Sketch : %s / %s cellules...", format(sketch_params$ncells, big.mark=" "), format(n_total_cells, big.mark=" ")))
    Seurat::DefaultAssay(obj) <- assay
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
    obj <- withr::with_seed(seed, Seurat::SketchData(object = obj, ncells = sketch_params$ncells,
                                                    method = "LeverageScore", sketched.assay = "sketch"))
    Seurat::DefaultAssay(obj) <- "sketch"
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
    obj <- smart_scale_data(obj)
    step_times$norm_s <- as.numeric(difftime(Sys.time(), t_norm, units = "secs"))

    t_pca <- Sys.time()
    obj <- Seurat::RunPCA(obj, npcs = sketch_params$npcs, verbose = FALSE)
    step_times$pca_s <- as.numeric(difftime(Sys.time(), t_pca, units = "secs"))
    log_fn(sprintf("✓ PCA (%d dims, sketch)", sketch_params$npcs))

    t_clust <- Sys.time()
    obj <- Seurat::FindNeighbors(obj, dims = 1:sketch_params$npcs, verbose = FALSE)
    obj <- robust_find_clusters(obj, resolution = resolution, algo = cluster_algo, log_fn = log_fn)
    step_times$clustering_s <- as.numeric(difftime(Sys.time(), t_clust, units = "secs"))
    log_fn(sprintf("✓ Clustering sketch OK (res %.1f)", resolution))

    if (isTRUE(compute_umap)) {
      t_umap <- Sys.time()
      obj <- withr::with_seed(seed, Seurat::RunUMAP(obj, dims = 1:sketch_params$npcs, reduction = "pca",
                                                   return.model = TRUE, verbose = FALSE))
      step_times$umap_s <- as.numeric(difftime(Sys.time(), t_umap, units = "secs"))
      log_fn("✓ UMAP sketch OK")
    }

    t_proj <- Sys.time()
    project_args <- list(
      object = obj, assay = assay, sketched.assay = "sketch",
      sketched.reduction = "pca", full.reduction = "pca.full",
      dims = 1:sketch_params$npcs,
      refdata = list(seurat_clusters = "seurat_clusters")
    )
    if (isTRUE(compute_umap)) project_args$umap.model <- "umap"
    obj <- do.call(Seurat::ProjectData, project_args)
    obj <- standardize_sketch_reductions(obj, full_pca_name = "pca.full")
    Seurat::DefaultAssay(obj) <- assay
    clust_red <- "pca"
    clust_dims <- sketch_params$npcs
    log_fn(sprintf("✓ Projection OK — %s cellules", format(ncol(obj), big.mark = " ")))

  } else {
    # ── Full Dataset Branch ──────────────────────────────────────────────────
    log_fn("Normalisation...")
    if (identical(norm_method, "sct")) {
      obj <- Seurat::SCTransform(obj, verbose = FALSE, variable.features.n = nfeatures, vst.flavor = "v2")
    } else {
      Seurat::DefaultAssay(obj) <- assay
      obj <- Seurat::NormalizeData(obj, verbose = FALSE)
      obj <- Seurat::FindVariableFeatures(obj, nfeatures = nfeatures, verbose = FALSE)
      obj <- smart_scale_data(obj)
    }
    step_times$norm_s <- as.numeric(difftime(Sys.time(), t_norm, units = "secs"))
    log_fn("✓ Normalisation OK")

    t_pca <- Sys.time()
    obj <- Seurat::RunPCA(obj, verbose = FALSE, npcs = pca_dim)
    step_times$pca_s <- as.numeric(difftime(Sys.time(), t_pca, units = "secs"))
    log_fn(sprintf("✓ PCA (%d dims)", pca_dim))

    # Harmony Integration
    if (identical(integration, "harmony")) {
      t_harm <- Sys.time()
      if (requireNamespace("harmony", quietly = TRUE) &&
          batch_var %in% colnames(obj@meta.data) &&
          length(unique(obj@meta.data[[batch_var]])) >= 2) {
        log_fn(sprintf("Harmony (variable de batch : %s)...", batch_var))
        obj <- withr::with_seed(seed, harmony::RunHarmony(obj, group.by.vars = batch_var,
                                                          dims.use = 1:min(30, pca_dim), verbose = FALSE))
        clust_red <- "harmony"
        clust_dims <- min(pca_dim, ncol(Seurat::Embeddings(obj, "harmony")))
        mixing_before <- tryCatch(sc_batch_mixing_score(obj, reduction = "pca", batch_col = batch_var)$score, error = function(e) NULL)
        mixing_after  <- tryCatch(sc_batch_mixing_score(obj, reduction = "harmony", batch_col = batch_var)$score, error = function(e) NULL)
        log_fn(sprintf("✓ Harmony appliquée (%d dims)", clust_dims))
      } else {
        log_fn("Harmony non appliquée : conditions non remplies.")
      }
      step_times$harmony_s <- as.numeric(difftime(Sys.time(), t_harm, units = "secs"))
    }

    t_clust <- Sys.time()
    obj <- Seurat::FindNeighbors(obj, dims = 1:clust_dims, reduction = clust_red, verbose = FALSE)
    obj <- robust_find_clusters(obj, resolution = resolution, algo = cluster_algo, log_fn = log_fn)
    step_times$clustering_s <- as.numeric(difftime(Sys.time(), t_clust, units = "secs"))
    log_fn(sprintf("✓ %d clusters (res %.1f)", length(unique(obj$seurat_clusters)), resolution))

    if (isTRUE(compute_umap)) {
      t_umap <- Sys.time()
      obj <- withr::with_seed(seed, Seurat::RunUMAP(obj, reduction = clust_red, dims = 1:clust_dims, verbose = FALSE))
      step_times$umap_s <- as.numeric(difftime(Sys.time(), t_umap, units = "secs"))
      log_fn("✓ UMAP OK")
    }

    if (isTRUE(compute_tsne)) {
      t_tsne <- Sys.time()
      max_tsne_cells <- if (exists(".AUTO_TSNE_MAX_CELLS")) .AUTO_TSNE_MAX_CELLS else 30000L
      if (ncol(obj) <= max_tsne_cells) {
        obj <- tryCatch(withr::with_seed(seed, Seurat::RunTSNE(obj, dims = 1:pca_dim, verbose = FALSE)),
                        error = function(e) { log_fn(paste("t-SNE ignoré:", e$message)); obj })
        log_fn("✓ t-SNE secondaire OK")
      } else {
        log_fn(sprintf("t-SNE secondaire ignoré (%d > %d cellules).", ncol(obj), max_tsne_cells))
      }
      step_times$tsne_s <- as.numeric(difftime(Sys.time(), t_tsne, units = "secs"))
    }
  }
  })

  n_cl <- length(unique(obj$seurat_clusters))

  # ── Step 6: SingleR (optional) ─────────────────────────────────────────────
  if (isTRUE(p_norm$sc_ap_singler)) {
    if (requireNamespace("SingleR", quietly = TRUE) && requireNamespace("celldex", quietly = TRUE) &&
        exists(".run_singler_safe", mode = "function")) {
      t_singler <- Sys.time()
      log_fn("Annotation SingleR...")
      result <- tryCatch(
        .run_singler_safe(obj, p_norm$sc_ap_singler_ref, p_norm$sc_ap_singler_level),
        error = function(e) { log_fn(paste("SingleR:", e$message)); NULL }
      )
      if (!is.null(result)) {
        col_name <- paste0("SingleR_", p_norm$sc_ap_singler_ref, "_", p_norm$sc_ap_singler_level)
        obj[[col_name]] <- result$labels
        log_fn(sprintf("✓ Annoté [%s] — %d types", result$method, length(unique(result$labels))))
      }
    } else {
      log_fn("SingleR/celldex non installés ou non disponibles — annotation ignorée.")
    }
  }

  # ── Step 7: FindAllMarkers & Pathways (optional) ───────────────────────────
  markers <- NULL
  pw <- NULL
  if (isTRUE(p_norm$sc_ap_markers) || isTRUE(p_norm$sc_ap_correlation)) {
    t_mark <- Sys.time()
    log_fn("FindAllMarkers...")
    cap_m <- as.integer(params$max_cells_heavy %||% 500L)
    sub_res <- subsample_seurat_for_analysis(obj, max_per_group = cap_m, group_col = "seurat_clusters")
    markers <- tryCatch({
      Seurat::Idents(sub_res$object) <- sub_res$object$seurat_clusters
      Seurat::FindAllMarkers(sub_res$object, only.pos = TRUE, min.pct = 0.1,
                             logfc.threshold = 0.25, verbose = FALSE)
    }, error = function(e) { log_fn(paste("Markers:", e$message)); NULL })

    if (!is.null(markers) && nrow(markers) > 0) {
      markers <- as.data.frame(markers)
      rownames(markers) <- NULL
      if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
      if (!"avg_log2FC" %in% colnames(markers)) markers$avg_log2FC <- markers$avg_logFC %||% 0
      has_padj <- "p_val_adj" %in% colnames(markers)
      if (!has_padj) {
        markers$p_val_adj <- NA_real_
      }
      if (!"cluster" %in% colnames(markers)) markers$cluster <- "Unknown"
      if (!"pct.1" %in% colnames(markers)) markers$pct.1 <- NA_real_
      if (!"pct.2" %in% colnames(markers)) markers$pct.2 <- NA_real_
      log_fn(sprintf("✓ %d marqueurs", nrow(markers)))
      step_times$markers_s <- as.numeric(difftime(Sys.time(), t_mark, units = "secs"))

      # Pathway ORA
      if (isTRUE(p_norm$sc_ap_pathway)) {
        t_pw <- Sys.time()
        log_fn("Pathway ORA...")
        pathway_org <- p_norm$sc_ap_pathway_org %||% "human"
        top_g_raw <- head(markers$gene[order(markers$p_val_adj)], 100)
        top_g <- if (exists(".remap_if_ensg", mode = "function")) {
          .remap_if_ensg(top_g_raw, pathway_org, notify_fn = function(m, ...) log_fn(m))
        } else {
          top_g_raw
        }
        if (length(top_g) > 0) {
          pw <- tryCatch(
            run_pathway_enrichment(top_g, organism = pathway_org,
                                   database = p_norm$sc_ap_pathway_db %||% "GOBP",
                                   pval_cutoff = 0.05,
                                   universe = rownames(obj)),
            error = function(e) { log_fn(paste("Pathway:", e$message)); NULL }
          )
          if (!is.null(pw) && nrow(pw) > 0) {
            log_fn(sprintf("✓ %d pathways", nrow(pw)))
          }
        }
        step_times$pathway_s <- as.numeric(difftime(Sys.time(), t_pw, units = "secs"))
      }
    } else {
      log_fn("Aucun marqueur trouvé.")
    }
  }

  # ── Step 8: Gene Correlation (optional) ────────────────────────────────────
  corr_res <- NULL
  if (isTRUE(p_norm$sc_ap_correlation)) {
    t_corr <- Sys.time()
    target_gene <- NULL
    if (!is.null(markers) && nrow(markers) > 0 && any(!is.na(markers$p_val_adj))) {
      ranked <- markers[order(markers$p_val_adj), ]
      target_gene <- ranked$gene[1]
    }
    if (!is.null(target_gene)) {
      log_fn(sprintf("Corrélation avec %s...", target_gene))
      cap_c <- as.integer(params$max_cells_heavy %||% 500L)
      sub_res_c <- subsample_seurat_for_analysis(obj, max_per_group = cap_c, group_col = "orig.ident")
      corr_res <- tryCatch(
        find_correlated_genes(sub_res_c$object, target_gene = target_gene,
                              method = "pearson", threshold = 0.3, top_n = 50),
        error = function(e) { log_fn(paste("Corrélation:", e$message)); NULL }
      )
      if (!is.null(corr_res) && nrow(corr_res) > 0) {
        log_fn(sprintf("✓ %d gènes corrélés avec %s", nrow(corr_res), target_gene))
      }
    }
    step_times$correlation_s <- as.numeric(difftime(Sys.time(), t_corr, units = "secs"))
  }

  # ── Step 9: Trajectory (optional) ──────────────────────────────────────────
  if (isTRUE(p_norm$sc_ap_trajectory)) {
    t_traj <- Sys.time()
    max_traj_cells <- if (exists(".MAX_TRAJECTORY_CELLS")) .MAX_TRAJECTORY_CELLS else 50000L
    if (ncol(obj) <= max_traj_cells) {
      log_fn("Trajectoire...")
      traj_red_use <- if ("umap" %in% names(obj@reductions)) "umap" else "pca"
      traj_embedding <- Seurat::Embeddings(obj, reduction = traj_red_use)
      attr(traj_embedding, "reduction") <- traj_red_use
      traj_res <- tryCatch(
        calculate_pseudotime(embeddings = traj_embedding, k = 15,
                             root_cells = NULL, root_method = "diameter"),
        error = function(e) { log_fn(paste("Trajectoire:", e$message)); NULL }
      )
      if (!is.null(traj_res)) {
        obj@meta.data$pseudotime        <- traj_res$pseudotime
        obj@meta.data$traj_in_component <- traj_res$in_root_component
        obj@meta.data$traj_method                <- rep("exploratory_knn", ncol(obj))
        obj@meta.data$traj_computation_reduction <- rep(traj_red_use, ncol(obj))
        obj@meta.data$traj_root_method           <- rep(traj_res$root_method, ncol(obj))
        obj@meta.data$traj_root_cell             <- rep(traj_res$root_cell, ncol(obj))
        obj@meta.data$traj_root_cluster          <- rep(NA_character_, ncol(obj))
        obj@meta.data$traj_root_component_size   <- rep(traj_res$root_component_size, ncol(obj))
        log_fn(sprintf("✓ Pseudotemps calculé (exploratoire kNN, racine auto/diamètre, réduction: %s)", toupper(traj_red_use)))
      }
    } else {
      log_fn(sprintf("Trajectoire ignorée (%d > %d cellules).", ncol(obj), max_traj_cells))
    }
    step_times$trajectory_s <- as.numeric(difftime(Sys.time(), t_traj, units = "secs"))
  }

  t_end <- Sys.time()
  step_times$total_s <- as.numeric(difftime(t_end, t_start, units = "secs"))

  # 3. Structured return
  metrics <- list(
    n_cells_initial     = n_cells_initial,
    n_cells_final       = ncol(obj),
    n_genes_initial     = n_genes_initial,
    n_genes_final       = nrow(obj),
    n_clusters          = n_cl,
    pca_dims            = clust_dims,
    resolution          = resolution,
    reductions          = names(obj@reductions),
    integration         = integration,
    batch_mixing_before = mixing_before,
    batch_mixing_after  = mixing_after,
    qc_removed          = qc_removed
  )

  list(
    obj              = obj,
    metrics          = metrics,
    timings          = step_times,
    log              = log_lines,
    markers_data     = markers,
    pathway_results  = pw,
    correlated_genes = corr_res,
    qc_snapshot      = qc_snapshot,
    doublets_provenance = m4_provenance,   # M4 : entrée new_provenance_entry (NULL si step désactivé)
    qc_distributions    = m4_distributions # M4 : distributions nCount/nFeature avant/après (NULL si step désactivé)
  )
}

#' Execute the SC auto-pipeline (called inside observeEvent)
#' @param input Shiny input object
#' @param global_data Global reactiveValues
#' @param shared_rv SC shared state
#' @param session Shiny session (for showNotification, log)
#' @param sc_log_rv reactiveVal for log
# ── P1-6a (2026-11-07, plan post-E2E #6) : cœur PUR du pipeline manuel ───────
#
# Transposition du corps du handler mod_sc_pipeline.R (qui tournait
# SYNCHRONE dans la session Shiny : 15-60 min de thread bloqué, inputs en
# file, fausse « session morte »). Le module délègue désormais à ce cœur via
# ExtendedTask + mirai (pool unique ensure_app_daemons(), pattern de
# run_sc_auto_pipeline_core déjà en daemon), avec repli synchrone tracé.
#
# AUCUN appel Shiny ici : les notifications sont COLLECTÉES (liste de
# list(type=, msg=)), l'état partagé est RETOURNÉ (le module l'applique dans
# sa session), la progression part par log_fn (lignes [TL] step|status
# consommées par le tracker du module) + on_progress (repli synchrone).
#
# @param obj Seurat object (assert_seurat)
# @param params liste nommée : qc_min_gene, qc_max_gene, qc_mt,
#   doublets_enabled, doublet_score_cutoff, max_contamination, run_ambient,
#   bpcells_mode, bpcells_threshold, norm_method, pca_dim, clust_res,
#   cluster_algo, reduction_method, seed — défauts = défauts UI du module.
# @param log_file chemin de log inter-process (streaming tracker)
# @param on_progress function(step, status) — repli synchrone (timeline)
# @return list(ok, obj, qc_snapshot, qc_dist_before, qc_dist_after,
#   provenance, bpcells_conv, skipped_steps, notifications, step, timings,
#   harmony_applied, mix_msg)
run_sc_manual_pipeline_core <- function(obj, params = list(), log_file = NULL,
                                        on_progress = NULL) {
  assert_seurat(obj, context = "run_sc_manual_pipeline_core")

  # ── helpers locaux (aucune pollution globale, aucune collision module) ──
  .man_mt_pattern <- function(o) {
    genes <- rownames(o)
    if (any(grepl("^MT-",  genes, ignore.case = FALSE))) return("^MT-")
    if (any(grepl("^mt-",  genes, ignore.case = FALSE))) return("^mt-")
    if (any(grepl("^MT-",  genes, ignore.case = TRUE)))  return("^MT-")
    NULL
  }
  .man_sequential_future <- function(expr) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan("sequential")
    force(expr)
  }
  .man_seeded <- function(expr, seed) withr::with_seed(seed, force(expr))
  # Copie daemon-safe de .dm_from_pca (module mod_sc_pipeline.R, verrouillé
  # par test-sc-pipeline-diffusion-map.R — on ne le déplace pas) ; mêmes
  # arguments, même orientation (CELLULES en lignes).
  .man_dm_from_pca <- function(pca_embeddings, npcs = NULL) {
    emb <- as.matrix(pca_embeddings)
    if (!is.null(npcs) && length(npcs) == 1L && is.finite(npcs)) {
      npcs <- as.integer(npcs)
      if (npcs >= 1L && npcs < ncol(emb)) emb <- emb[, seq_len(npcs), drop = FALSE]
    }
    dm <- destiny::DiffusionMap(emb)
    n_out <- min(2L, ncol(dm@eigenvectors))
    dm@eigenvectors[, seq_len(n_out), drop = FALSE]
  }
  .man_tl <- function(step, status) {
    line <- sprintf("[TL] %s|%s", step, status)
    log_fn(line)
    if (is.function(on_progress)) {
      tryCatch(on_progress(step, status), error = function(e) NULL)
    }
    invisible()
  }

  t_start <- Sys.time()
  log_lines <- character(0)
  log_fn <- function(msg) {
    line <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)
    log_lines <<- c(log_lines, line)
    if (!is.null(log_file) && is.character(log_file) && nzchar(log_file)) {
      tryCatch(cat(line, "\n", file = log_file, append = TRUE), error = function(e) NULL)
    }
  }

  # ── résolution des paramètres (défauts = défauts UI du module) ────────────
  qc_min_gene  <- as.numeric(params$qc_min_gene %||% 100)
  qc_max_gene  <- as.numeric(params$qc_max_gene %||% 8000)
  qc_mt        <- as.numeric(params$qc_mt %||% 20)
  doublets_enabled <- isTRUE(params$doublets_enabled %||%
                               requireNamespace("scDblFinder", quietly = TRUE))
  doublet_cutoff <- as.numeric(params$doublet_score_cutoff %||%
    if (exists("TS_SC_QC_DOUBLET_SCORE_CUTOFF")) TS_SC_QC_DOUBLET_SCORE_CUTOFF else 0.25)
  max_contam   <- as.numeric(params$max_contamination %||%
    if (exists("TS_SC_QC_MAX_CONTAMINATION")) TS_SC_QC_MAX_CONTAMINATION else 0.20)
  run_ambient  <- isTRUE(params$run_ambient)
  bpcells_mode <- as.character(params$bpcells_mode %||% "auto")
  bpcells_threshold <- as.numeric(params$bpcells_threshold %||%
    if (exists("TS_BPCELLS_AUTO_THRESHOLD")) TS_BPCELLS_AUTO_THRESHOLD else 150000L)
  norm_method  <- as.character(params$norm_method %||% "log")
  pca_dim      <- as.integer(params$pca_dim %||% 20)
  clust_res    <- as.numeric(params$clust_res %||% 0.5)
  cluster_algo <- suppressWarnings(as.integer(params$cluster_algo %||% "1"))
  reduction_method <- as.character(params$reduction_method %||% "umap")
  seed         <- as.integer(params$seed %||%
                               if (exists("TS_SC_PIPELINE_SEED")) TS_SC_PIPELINE_SEED else 989L)

  notifications <- list()
  note <- function(type, msg) {
    notifications <<- c(notifications, list(list(type = type, msg = msg)))
    log_fn(paste0("[", toupper(type), "] ", msg))
    invisible()
  }

  timings <- list(qc_s = 0, doublets_s = 0, backend_s = 0, norm_s = 0, pca_s = 0,
                  harmony_s = 0, clustering_s = 0, reduction_s = 0, tsne_s = 0,
                  total_s = 0)
  qc_snapshot <- NULL
  qc_dist_before <- NULL
  qc_dist_after <- NULL
  prov_list <- list()
  bpcells_conv <- NULL
  harmony_applied <- FALSE
  mix_msg <- NULL
  skipped_steps <- character(0)
  step <- "qc"

  out_fail <- function(e) {
    list(ok = FALSE, obj = obj, qc_snapshot = qc_snapshot,
         qc_dist_before = qc_dist_before, qc_dist_after = qc_dist_after,
         provenance = prov_list, bpcells_conv = bpcells_conv,
         skipped_steps = skipped_steps, notifications = notifications,
         step = step, timings = timings, harmony_applied = harmony_applied,
         mix_msg = mix_msg, error = e)
  }

  tryCatch({
    .man_tl("qc", "running")

    # ── Step 1 : QC ────────────────────────────────────────────────────────
    t0 <- Sys.time()
    mt_pattern <- .man_mt_pattern(obj)
    if (!is.null(mt_pattern)) {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
    } else {
      obj[["percent.mt"]] <- 0
      note("warning", "Aucun gène mitochondrial détecté. Filtre MT ignoré.")
    }

    n_before <- ncol(obj)
    meta_pre <- obj@meta.data
    n_ok_min <- sum(meta_pre$nFeature_RNA >  qc_min_gene, na.rm = TRUE)
    n_ok_max <- sum(meta_pre$nFeature_RNA <  qc_max_gene, na.rm = TRUE)
    n_ok_mt  <- sum(meta_pre$percent.mt   <  qc_mt,       na.rm = TRUE)
    n_ok_all <- sum(meta_pre$nFeature_RNA >  qc_min_gene &
                    meta_pre$nFeature_RNA <  qc_max_gene &
                    meta_pre$percent.mt   <  qc_mt, na.rm = TRUE)

    obj <- subset(obj,
      subset = nFeature_RNA > qc_min_gene &
               nFeature_RNA < qc_max_gene &
               percent.mt   < qc_mt)

    if (ncol(obj) < 10) {
      stop(errorCondition(sprintf(
        paste0(
          "Seulement %d cellule(s) après QC (départ : %d). Diagnostic par filtre :\n",
          "  • Min gènes > %d  : %d cellules passent (%d%%)\n",
          "  • Max gènes < %d  : %d cellules passent (%d%%)\n",
          "  • %% Mito < %d%%  : %d cellules passent (%d%%)\n",
          "  • Les 3 combinés  : %d cellules passent\n\n",
          "Suggestions : réduire Min gènes (essayez %d), augmenter Max gènes (%d) ou %% Mito (%d%%)."
        ),
        ncol(obj), n_before,
        qc_min_gene, n_ok_min, round(100 * n_ok_min / n_before),
        qc_max_gene, n_ok_max, round(100 * n_ok_max / n_before),
        qc_mt,       n_ok_mt,  round(100 * n_ok_mt  / n_before),
        n_ok_all,
        max(10,  qc_min_gene - 100),
        min(50000, qc_max_gene + 2000),
        min(50,  qc_mt + 15)
      ), class = "mod_sc_pipeline_error"))
    }

    # snapshot QC par échantillon (méta pré-filtre, cf. module Step-3.7A.2)
    samples    <- union(names(table(meta_pre$orig.ident)), unique(as.character(obj$orig.ident)))
    before_tab <- table(meta_pre$orig.ident)
    after_tab  <- table(obj$orig.ident)
    mt_med     <- tapply(meta_pre$percent.mt, meta_pre$orig.ident, median, na.rm = TRUE)
    qc_snapshot <- data.frame(
      Echantillon     = samples,
      Cellules_avant  = as.integer(ifelse(is.na(before_tab[samples]), 0L, before_tab[samples])),
      Cellules_apres  = as.integer(ifelse(is.na(after_tab[samples]),  0L, after_tab[samples])),
      Pct_Mito_Median = round(as.numeric(mt_med[samples]), 1),
      stringsAsFactors = FALSE
    )
    qc_snapshot$Pct_Conserve <- round(100 * qc_snapshot$Cellules_apres /
                                        pmax(qc_snapshot$Cellules_avant, 1), 1)
    timings$qc_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    .man_tl("qc", "done")

    # ── Step 1c : doublets & ambiant (M4) ───────────────────────────────────
    if (isTRUE(doublets_enabled)) {
      if (requireNamespace("scDblFinder", quietly = TRUE) &&
          exists("sc_filter_doublets_and_ambient", mode = "function")) {
        step <- "doublets"; .man_tl("doublets", "running")
        t0 <- Sys.time()
        m4_res <- sc_filter_doublets_and_ambient(obj, params = list(
          seed = seed,
          batch_var = "orig.ident",
          doublet_score_cutoff = doublet_cutoff,
          max_contamination    = max_contam,
          run_ambient          = run_ambient
        ))
        obj <- m4_res$obj
        qc_dist_before <- list(nCount_RNA   = m4_res$qc_metrics$distributions$nCount_RNA$before,
                               nFeature_RNA = m4_res$qc_metrics$distributions$nFeature_RNA$before)
        qc_dist_after  <- list(nCount_RNA   = m4_res$qc_metrics$distributions$nCount_RNA$after,
                               nFeature_RNA = m4_res$qc_metrics$distributions$nFeature_RNA$after)
        excl_cells <- as.character(m4_res$provenance$cells_excluded)
        if ("orig.ident" %in% colnames(meta_pre)) {
          exc_by_sample <- table(factor(meta_pre[excl_cells, "orig.ident"], levels = samples))
          qc_snapshot$Doublets_exclus <- as.integer(
            ifelse(is.na(exc_by_sample[samples]), 0L, exc_by_sample[samples]))
        } else {
          qc_snapshot$Doublets_exclus <- 0L
        }
        qc_snapshot$Pct_Doublets <- round(100 * qc_snapshot$Doublets_exclus /
                                            pmax(qc_snapshot$Cellules_apres, 1), 1)
        prov_list <- c(prov_list, list(m4_res$provenance))
        if (length(excl_cells) > 0) {
          note("message", sprintf("✓ Doublets : %d cellule(s) exclue(s).", length(excl_cells)))
        }
        timings$doublets_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        .man_tl("doublets", "done")
      } else {
        note("warning", "scDblFinder non installé — étape doublets ignorée.")
        skipped_steps <- c(skipped_steps, "doublets")
        .man_tl("doublets", "skip")
      }
    } else {
      skipped_steps <- c(skipped_steps, "doublets")
      .man_tl("doublets", "skip")
    }

    # ── Step 1b : backend disque (BPCells) ──────────────────────────────────
    step <- "backend"; .man_tl("backend", "running")
    t0 <- Sys.time()
    should_convert <- switch(bpcells_mode,
      "always" = TRUE,
      "never"  = FALSE,
      ncol(obj) > bpcells_threshold
    )
    if (should_convert && sc_backend_status(obj) == "memory") {
      if (!.bpcells_available()) {
        note("warning", paste0("Package 'BPCells' non installé — pipeline exécuté en RAM comme d'habitude. ",
                               "Installez-le pour les très gros datasets (remotes::install_github('bnprks/BPCells/r'))."))
      } else {
        conv <- tryCatch(convert_seurat_to_bpcells(obj), error = function(e) {
          note("warning", paste("BPCells: conversion échouée, poursuite en RAM —", e$message))
          NULL
        })
        if (!is.null(conv)) {
          obj <- conv$object
          bpcells_conv <- conv
          if (!isTRUE(conv$already_disk)) {
            note("message", sprintf("💽 Backend disque activé (BPCells) — %s cellules, RAM préservée pour la suite.",
                                    format(conv$n_cells, big.mark = " ")))
          }
        }
      }
    }
    # future séquentiel sur objets disque (cf. module — export de closures
    # géantes par le dispatch parallèle Seurat, observé sur 1.3M cellules)
    if (sc_backend_status(obj) == "disk") .man_sequential_future(NULL)
    timings$backend_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    .man_tl("backend", "done")

    # ── Step 2 : Normalisation ──────────────────────────────────────────────
    step <- "norm"; .man_tl("norm", "running")
    t0 <- Sys.time()
    if (identical(norm_method, "sct")) {
      obj <- SCTransform(obj, verbose = FALSE, variable.features.n = 3000, vst.flavor = "v2")
    } else {
      DefaultAssay(obj) <- "RNA"
      obj <- NormalizeData(obj)
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
      obj <- smart_scale_data(obj)
    }
    timings$norm_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    .man_tl("norm", "done")

    # ── Step 3 : PCA (+ libération scale.data) ─────────────────────────────
    step <- "pca"; .man_tl("pca", "running")
    t0 <- Sys.time()
    obj <- RunPCA(obj, verbose = FALSE, npcs = pca_dim)
    obj <- tryCatch({
      da <- DefaultAssay(obj)
      if (isTRUE("scale.data" %in% SeuratObject::Layers(obj[[da]]))) obj[[da]]$scale.data <- NULL
      obj
    }, error = function(e) obj)
    # libération RAM : clean_mem (memory.R, sourcé session) si présent — le
    # cœur doit rester exécutable hors app (smoke headless)
    if (exists("clean_mem", mode = "function")) clean_mem() else invisible(gc())
    timings$pca_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    .man_tl("pca", "done")

    # ── Step 4 : Harmony (AVANT clustering) ─────────────────────────────────
    if (identical(reduction_method, "harmony")) {
      if (requireNamespace("harmony", quietly = TRUE)) {
        if ("orig.ident" %in% colnames(obj@meta.data) &&
            length(unique(obj$orig.ident)) >= 2) {
          step <- "harmony"; .man_tl("harmony", "running")
          t0 <- Sys.time()
          if (ncol(obj) > 100000) {
            note("info", "Large dataset: paramètres Harmony réduits.")
            obj <- .man_sequential_future(
              RunHarmony(obj, group.by.vars = "orig.ident",
                         dims.use = 1:min(30, pca_dim), max.iter.harmony = 10, verbose = FALSE))
          } else {
            obj <- .man_sequential_future(
              RunHarmony(obj, group.by.vars = "orig.ident",
                         dims.use = 1:pca_dim, verbose = FALSE))
          }
          harmony_applied <- TRUE
          mix_msg <- tryCatch({
            before <- sc_batch_mixing_score(obj, reduction = "pca", batch_col = "orig.ident")
            after  <- sc_batch_mixing_score(obj, reduction = "harmony", batch_col = "orig.ident")
            sprintf("Mélange des batchs (0 = séparés, 1 = mélangés) : PCA brute %.0f%% → Harmony %.0f%%.",
                    100 * before$score, 100 * after$score)
          }, error = function(e) NULL)
          if (!is.null(mix_msg)) note("message", mix_msg)
          timings$harmony_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
          .man_tl("harmony", "done")
        } else {
          note("warning", "Harmony non appliquée : 0 ou 1 échantillon (orig.ident) — clustering et UMAP sur PCA brute.")
          skipped_steps <- c(skipped_steps, "harmony")
          .man_tl("harmony", "skip")
        }
      } else {
        note("warning", "Package 'harmony' introuvable — clustering et UMAP sur PCA brute.")
        skipped_steps <- c(skipped_steps, "harmony")
        .man_tl("harmony", "skip")
      }
    } else {
      skipped_steps <- c(skipped_steps, "harmony")
      .man_tl("harmony", "skip")
    }

    # espace de clustering = Harmony si appliquée, sinon PCA brute
    clust_red  <- if (harmony_applied) "harmony" else "pca"
    clust_dims <- if (harmony_applied) min(pca_dim, ncol(Seurat::Embeddings(obj, "harmony")))
                  else pca_dim

    # ── Step 5 : Clustering ─────────────────────────────────────────────────
    step <- "cluster"; .man_tl("cluster", "running")
    t0 <- Sys.time()
    obj <- FindNeighbors(obj, reduction = clust_red, dims = 1:clust_dims)
    obj <- tryCatch(
      FindClusters(obj, resolution = clust_res, algorithm = cluster_algo),
      error = function(e) {
        if (cluster_algo != 1L) {
          algo_name <- c("2" = "Louvain multilevel", "3" = "SLM", "4" = "Leiden")[as.character(cluster_algo)]
          note("warning", paste0(algo_name %||% "Algorithme",
                                 " indisponible/a échoué sur cet objet — repli sur Louvain standard. ",
                                 e$message))
          FindClusters(obj, resolution = clust_res, algorithm = 1)
        } else stop(e)
      }
    )
    timings$clustering_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    .man_tl("cluster", "done")

    # ── Step 6 : réduction principale ───────────────────────────────────────
    step <- "reduction"; .man_tl("reduction", "running")
    t0 <- Sys.time()
    if (harmony_applied) {
      obj <- .man_seeded(RunUMAP(obj, reduction = "harmony", dims = 1:clust_dims,
                                 verbose = FALSE, reduction.name = "umap_harmony"), seed)
    } else if (identical(reduction_method, "harmony")) {
      obj <- .man_seeded(RunUMAP(obj, dims = 1:pca_dim, verbose = FALSE), seed)
    } else if (identical(reduction_method, "umap")) {
      obj <- .man_seeded(RunUMAP(obj, dims = 1:pca_dim, verbose = FALSE), seed)
    } else if (identical(reduction_method, "tsne")) {
      obj <- .man_seeded(RunTSNE(obj, dims = 1:pca_dim, verbose = FALSE), seed)
    } else if (identical(reduction_method, "dm")) {
      if (requireNamespace("destiny", quietly = TRUE)) {
        dm_emb <- .man_dm_from_pca(Seurat::Embeddings(obj, "pca"), npcs = pca_dim)
        obj[["dm"]] <- CreateDimReducObject(embeddings = dm_emb,
                                            key = "DM_", assay = DefaultAssay(obj))
      } else {
        note("warning", "Package 'destiny' introuvable. UMAP utilisé.")
        obj <- .man_seeded(RunUMAP(obj, dims = 1:pca_dim, verbose = FALSE), seed)
      }
    }
    timings$reduction_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    .man_tl("reduction", "done")

    # ── Step 5b : t-SNE secondaire ──────────────────────────────────────────
    tsne_cap <- if (exists("TS_AUTO_TSNE_MAX_CELLS")) TS_AUTO_TSNE_MAX_CELLS else 10000L
    if (!"tsne" %in% names(obj@reductions) && !identical(reduction_method, "tsne")) {
      if (ncol(obj) > tsne_cap) {
        skipped_steps <- c(skipped_steps, "tsne")
        .man_tl("tsne", "skip")
        note("info", sprintf("t-SNE secondaire ignoré (%s cellules > %s max). Choisissez 't-SNE' comme méthode principale si nécessaire.",
                             format(ncol(obj), big.mark = " "), format(tsne_cap, big.mark = " ")))
      } else {
        step <- "tsne"; .man_tl("tsne", "running")
        t0 <- Sys.time()
        obj <- tryCatch(.man_seeded(RunTSNE(obj, dims = 1:pca_dim, verbose = FALSE), seed),
                        error = function(e) {
                          note("warning", paste("t-SNE secondaire ignoré:", e$message))
                          obj
                        })
        timings$tsne_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        .man_tl("tsne", "done")
      }
    } else {
      skipped_steps <- c(skipped_steps, "tsne")
      .man_tl("tsne", "skip")
    }

    timings$total_s <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    .man_tl("END", "ok")
    list(ok = TRUE, obj = obj, qc_snapshot = qc_snapshot,
         qc_dist_before = qc_dist_before, qc_dist_after = qc_dist_after,
         provenance = prov_list, bpcells_conv = bpcells_conv,
         skipped_steps = skipped_steps, notifications = notifications,
         step = step, timings = timings, harmony_applied = harmony_applied,
         mix_msg = mix_msg, error = NULL)
  }, error = function(e) {
    .man_tl("END", paste0("error:", step))
    out_fail(e)
  })
}

run_sc_auto_pipeline <- function(input, global_data, shared_rv, session, sc_log_rv) {
      # ── M-2 phase 1 (2026-09-30) : les paramètres passent par le normalisateur PUR —
      # liste fermée, défauts canoniques, clés inconnues perdues (la classe « une clé
      # manquante change silencieusement un seuil » — audit C-1 — est fermée
      # structurellement).
      params <- sc_ap_normalize_params(
        if (is.list(input)) input else shiny::reactiveValuesToList(input))
      removeModal()
      req(state_get(global_data, "sc_obj"))

      ll <- character(0)
      log_sc <- function(msg) {
        ll <<- c(ll, paste0("[", format(Sys.time(),"%H:%M:%S"), "] ", msg))
        sc_log_rv(paste(ll, collapse="\n"))
      }

      p <- shiny::Progress$new(); on.exit(p$close())

      # ── Préflight RAM (M-4, 2026-09-30) : ALERTE dans le log, JAMAIS
      # bloquante — le refus risquerait de bloquer des analyses légitimes
      # (MODE C, direction 2). Facteur à calibrer par benchmark (§2dr).
      ram_sys <- ts_system_ram_mb()
      ram_budget <- ts_ram_budget_check(
        object_bytes = as.numeric(object.size(state_get(global_data, "sc_obj"))),
        factor = TS_RAM_PREFLIGHT_FACTOR,
        available_mb = ram_sys$available_mb,
        total_mb = ram_sys$total_mb)
      if (ram_budget$level != "none") log_sc(ram_budget$message)

      tryCatch({
        obj <- state_get(global_data, "sc_obj")
        p$set(0.1, .tr("Calcul en cours..."))

        res <- run_sc_auto_pipeline_core(obj, params)

        for (l in res$log) {
          log_sc(sub("^\\[[0-9:]+\\]\\s*", "", l))
        }

        obj <- res$obj
        state_set(global_data, "sc_obj", obj)
        if (!is.null(res$markers_data)) state_set(shared_rv, "markers_data", res$markers_data)
        if (!is.null(res$pathway_results)) {
          state_set(shared_rv, "pathway_results", res$pathway_results)
          state_set(shared_rv, "pathway_db", params$sc_ap_pathway_db %||% "GOBP")
        }
        if (!is.null(res$correlated_genes)) {
          state_set(shared_rv, "correlated_genes", res$correlated_genes)
          if (!is.null(res$markers_data) && nrow(res$markers_data) > 0) {
            state_set(shared_rv, "corr_target_gene", res$markers_data$gene[1])
          }
        }
        if ("pseudotime" %in% colnames(obj@meta.data)) {
          traj_red_use <- if ("umap" %in% names(obj@reductions)) "umap" else "pca"
          state_set(shared_rv, "traj_reduction", traj_red_use)
          state_set(shared_rv, "traj_method", "exploratory_knn")
        }

        state_set(shared_rv, "active_tab", "tab_viz")
        showNotification(
          sprintf(.tr("✓ Pipeline SC : %d cellules, %d clusters"), ncol(obj), res$metrics$n_clusters),
          type="message", duration=6)

      }, error=function(e) {
        log_sc(paste(.tr("❌ Erreur:"), e$message))
        showNotification(paste(.tr("Erreur pipeline SC:"), e$message), type="error", duration=10)
        # Roadmap 5.2 (audit 2026-09-27 §1.9) : l'erreur est journalisée ET
        # notifiée, puis RELANCÉE — la classe sc_pipeline_error devient
        # observable (le verrou source qui interdisait le stop est réécrit
        # dans test-sc-pipeline.R). Les appelants UI enveloppent dans
        # tryCatch (mod_sc.R) / .sc_ap_run_drive pour ne pas casser la session.
        stop(e)
      })
}

