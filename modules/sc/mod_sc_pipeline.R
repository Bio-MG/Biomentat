# mod_sc_pipeline.R  —  Child 1: QC → Norm → PCA → Harmony → Clusters → UMAP
# Step-3.6: relaxed QC defaults for test datasets; actionable QC error messages.
# Step-3.7:
#   - Secondary t-SNE always computed (capped) alongside PCA/UMAP so it's
#     available in "3. Visualisation -> Réduction à visualiser" without an
#     extra manual run, mirroring PCA/UMAP's "always on" availability.
#   - "Max cellules/cluster pour analyses lourdes" slider -> shared_rv$max_cells_heavy,
#     read by mod_sc_markers.R / mod_sc_corr.R (and the auto-pipeline in mod_sc.R)
#     to cap FindAllMarkers / find_correlated_genes on large datasets.
#
# P1-6a (2026-11-07, plan post-E2E #6) : le CORPS du pipeline vit désormais
# dans run_sc_manual_pipeline_core() (tools/pkg_skeleton/R/sc_pipeline.R,
# moteur PUR, préchargé daemons via APP_DAEMON_SOURCE_FILES) — le module
# délègue via ExtendedTask + mirai (pool unique ensure_app_daemons(), même
# pattern que l'auto-pipeline mod_sc.R), repli synchrone TRACÉ si le pool
# est indisponible (jamais muet, convention run_job). Les helpers transférés
# dans le cœur (mt-pattern, future séquentiel, seed, big-dataset) ne vivent
# plus ici — 15-60 min de thread Shiny bloqué : terminé.

# Carte de diffusion depuis les embeddings PCA (chemin manuel "dm").
#
# 🔴 ORIENTATION : `Seurat::Embeddings()` rend CELLULES x DIMENSIONS, et c'est
# exactement ce que `destiny::DiffusionMap()` attend — il prend les ECHANTILLONS
# EN LIGNES et les variables en COLONNES (aide : « vars: Variables (columns) of
# the data »). Le `t()` historique inversait les axes : la carte portait sur
# `pca_dim` « echantillons » au lieu des n cellules, puis
# `CreateDimReducObject()` levait (les noms de lignes ne correspondaient a
# aucune cellule) et le run entier etait perdu — `sc_obj` n'etant ecrit qu'en
# fin de pipeline.
#
# Le helper existe pour etre TESTABLE sans Shiny : la mesure d'orientation
# (entree 6x3 => 6 lignes) est faite dans test-sc-pipeline-diffusion-map.R.
.dm_from_pca <- function(pca_embeddings, npcs = NULL) {
  emb <- as.matrix(pca_embeddings)
  if (!is.null(npcs) && length(npcs) == 1L && is.finite(npcs)) {
    npcs <- as.integer(npcs)
    if (npcs >= 1L && npcs < ncol(emb)) emb <- emb[, seq_len(npcs), drop = FALSE]
  }
  dm <- destiny::DiffusionMap(emb)
  n_out <- min(2L, ncol(dm@eigenvectors))
  dm@eigenvectors[, seq_len(n_out), drop = FALSE]
}

# Audit 2026-09-27 §1.8 : graine DÉCLARÉE pour les étapes stochastiques
# (UMAP / t-SNE) — restaurée après chaque appel (withr::with_seed), donc deux
# exécutions du pipeline produisent les mêmes embeddings sans fuiter l'état
# du RNG global vers le reste de l'application.
.SC_PIPELINE_SEED <- if (exists("TS_SC_PIPELINE_SEED")) TS_SC_PIPELINE_SEED else 989L

# t-SNE is quadratic-ish in cost vs UMAP's approximate NN — cap the automatic
# "secondary" run so it never becomes the bottleneck of a routine pipeline run
# on a 32Go CPU-only machine. Users who explicitly pick "t-SNE" as their PRIMARY
# reduction method below are not capped (informed, deliberate choice).
# fix 2026-11-03 (mesure E2E 16 échantillons GSE96583) : 27 550 cellules =
# >60 min dans la phase t-SNE secondaire (~7% CPU, aucun log, Progress
# invisible) — le seuil 30k laissait les jeux ~30k dans la zone coûteuse.
# 10k garde le confort « toutes réductions dispo » sur les petits jeux sans
# bloquer les moyens ; t-SNE explicite reste possible en méthode principale.
.AUTO_TSNE_MAX_CELLS <- if (exists("TS_AUTO_TSNE_MAX_CELLS")) TS_AUTO_TSNE_MAX_CELLS else 10000L

# Step-3.7A: default cell-count threshold above which "Auto" mode switches the
# pipeline to a disk-backed (BPCells) counts layer instead of RAM. Below this,
# a sparse dgCMatrix comfortably fits in 32Go even with generous headroom for
# ScaleData()/RunPCA() intermediates -- BPCells' overhead (slower per-access,
# disk I/O) only pays off once RAM would otherwise be the bottleneck.
.BPCELLS_AUTO_THRESHOLD <- if (exists("TS_BPCELLS_AUTO_THRESHOLD")) TS_BPCELLS_AUTO_THRESHOLD else 150000L

# M4 (features 19-22, jalon BIOMENTAT 2026-10-26) : défauts des contrôles
# contrefactuels doublets/ambiant — SOURCE UNIQUE defaults.R (contrat
# PROJECT.md:146-147, convention QW-2 anti-fork UI/drive) ; retombe sur le
# littéral contractuel si defaults.R n'est pas résidente (test isolé).
.PIPELINE_DOUBLET_CUTOFF_DEFAULT <- if (exists("TS_SC_QC_DOUBLET_SCORE_CUTOFF")) TS_SC_QC_DOUBLET_SCORE_CUTOFF else 0.25
.PIPELINE_MAX_CONTAMINATION_DEFAULT <- if (exists("TS_SC_QC_MAX_CONTAMINATION")) TS_SC_QC_MAX_CONTAMINATION else 0.20

# ── UI ────────────────────────────────────────────────────────────────────────

mod_sc_pipeline_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # QC filters — relaxed defaults for broader compatibility
    layout_columns(
      numericInput(ns("qc_min_gene"), i18n$t("Min Gènes"),  100,  step=50),   # was 200
      numericInput(ns("qc_max_gene"), i18n$t("Max Gènes"),  8000, step=500)   # was 6000
    ),
    sliderInput(ns("qc_mt"), i18n$t("% Mito Max"), 0, 50, 20, step=1),        # was 5
    div(class="alert alert-light", style="font-size:0.82em;padding:6px;",
        bsicons::bs_icon("info-circle"),
        span(i18n$t("Valeurs conseillées pour jeux de test réduits :")),
        " ",
        span(i18n$t("Min gènes = 50, Max gènes = 10 000, % Mito = 40."))),
    # --- M4 (features 19-22) : QC doublets & ambiant — contrôles contrefactuels.
    # L'activation par défaut suit la présence du package (opt-in honnête) ;
    # les défauts lisent defaults.R (source unique QW-2).
    hr(),
    h6(i18n$t("QC Doublets & Ambiant (scDblFinder)"), style="font-weight:bold;"),
    checkboxInput(ns("doublets_enabled"),
                  i18n$t("Filtrer les doublets (scDblFinder)"),
                  value = requireNamespace("scDblFinder", quietly = TRUE)),
    layout_columns(
      sliderInput(ns("doublet_score_cutoff"), i18n$t("Seuil score doublet"),
                  0, 1, .PIPELINE_DOUBLET_CUTOFF_DEFAULT, step = 0.05),
      sliderInput(ns("max_contamination"), i18n$t("Contamination ambiante max"),
                  0, 1, .PIPELINE_MAX_CONTAMINATION_DEFAULT, step = 0.05)
    ),
    checkboxInput(ns("run_ambient"),
                  i18n$t("Décontamination ambiante (SoupX/celda) — opt-in"),
                  value = FALSE),
    div(class="small text-muted mb-2",
        i18n$t("Exclut les cellules appelées doublets dont le score dépasse strictement le seuil ; la détection tourne une fois, le seuil ne pilote que le filtrage."),
        " ",
        i18n$t("Sans moteur SoupX/celda installé, l'étape ambiante est une erreur classée (jamais de repli silencieux).")),
    # --- UI: sketch size preset selector ---
    selectInput(
      ns("sketch_preset"),
      i18n$t("Taille du sketch (vitesse vs precision)"),
      choices = setNames(
        c("fast", "light", "medium", "standard", "high", "max", "custom"),
        c(.tr_plain("Rapide (test, 5 000 cellules)"),
          .tr_plain("Léger (10 000 cellules)"),
          .tr_plain("Moyen (25 000 cellules)"),
          .tr_plain("Standard (50 000 cellules)"),
          .tr_plain("Élevé (100 000 cellules)"),
          .tr_plain("Max (dataset complet)"),
          .tr_plain("Personnalisé"))
      ),
      selected = "standard"
    ),
    
    conditionalPanel(
      condition = sprintf("input['%s'] == 'custom'", ns("sketch_preset")),
      numericInput(ns("sketch_ncells_custom"), i18n$t("Nombre de cellules (sketch)"),
                   value = 20000, min = 1000, max = 500000, step = 1000)
    ),
    hr(),
    h6(i18n$t("Normalisation")),
    radioButtons(ns("norm_method"), NULL,
      choices=setNames(c("log", "sct"),
                       c(.tr_plain("LogNormalize (Standard)"), .tr_plain("SCTransform (Avancé)")))),
    hr(),
    h6(i18n$t("Réduction Dimensionnelle")),
    selectInput(ns("reduction_method"), i18n$t("Méthode de Réduction (principale)"),
      choices=setNames(c("umap", "pca", "tsne", "dm", "harmony"),
                       c("UMAP", .tr_plain("PCA"), "t-SNE", .tr_plain("Diffusion Maps"), "Harmony")), selected="umap"),
    div(class="small text-muted mb-2",
        i18n$t("t-SNE est en plus toujours calculé automatiquement en secondaire (si < "),
        format(.AUTO_TSNE_MAX_CELLS, big.mark=" "),
        i18n$t(" cellules) pour être disponible dans l'onglet Visualisation.")),
    layout_columns(
      sliderInput(ns("pca_dim"),    i18n$t("Dims PCA"),   5, 50, 20),
      numericInput(ns("clust_res"), i18n$t("Résolution"), 0.5, step=0.1)
    ),
    selectInput(ns("cluster_algo"), i18n$t("Algorithme de clustering"),
                choices = setNames(c("1", "2", "3", "4"),
                                   c(.tr_plain("Louvain (standard)"),
                                     .tr_plain("Louvain (multilevel refinement)"),
                                     .tr_plain("SLM (Smart Local Moving)"),
                                     .tr_plain("Leiden (nécessite reticulate + leidenalg)"))),
                selected = "1"),
    div(class="small text-muted mb-2",
        i18n$t("Leiden est souvent plus rapide/qualitatif sur les gros datasets, mais"),
        i18n$t(" nécessite un environnement Python (reticulate + leidenalg/igraph)"),
        i18n$t(" installé séparément — repli automatique sur Louvain si indisponible.")),
    actionButton(ns("run_pipeline"), i18n$t("Lancer Pipeline"),
                 class="btn-danger w-100", icon=icon("play")),
    # P1-10 (fix 2026-11-07, plan post-E2E #10) : timeline PERSISTANTE des
    # étapes du pipeline — contrairement à la Progress (fermée en fin de
    # course) et aux notifications éphémères, elle reste lisible après le run
    # (mesuré E2E : « Progress$new() ne montre rien pendant la phase t-SNE/
    # finale »). Alimentée par le serveur ci-dessous (.tl_set/.tl_end).
    uiOutput(ns("pipeline_timeline")),
    hr(),
    h6(i18n$t("Gros Dataset — Backend Disque (BPCells)"), style="font-weight:bold;"),
    selectInput(ns("bpcells_mode"), i18n$t("Mode"),
                choices = setNames(c("auto", "always", "never"),
                                   c(.tr_plain("Auto (recommandé)"),
                                     .tr_plain("Toujours (forcer disque)"),
                                     .tr_plain("Jamais (toujours RAM)"))),
                selected = "auto"),
    conditionalPanel(
      condition = "input.bpcells_mode == 'auto'", ns = ns,
      numericInput(ns("bpcells_threshold"), i18n$t("Seuil auto (cellules)"),
                   value = .BPCELLS_AUTO_THRESHOLD, min = 10000, step = 10000)
    ),
    div(class="small text-muted mb-2",
        i18n$t("Au-delà du seuil, les counts bruts sont écrits sur disque (BPCells) avant"),
        i18n$t(" Normalisation — la RAM n'est plus le facteur limitant. Nécessite le"),
        i18n$t(" package 'BPCells' (optionnel) ; sans lui, le pipeline continue en RAM"),
        i18n$t(" normalement avec un avertissement.")),
    hr(),
    h6(i18n$t("Performance / RAM"), style="font-weight:bold;"),
    numericInput(ns("max_cells_heavy"),
                 i18n$t("Max cellules/cluster pour analyses lourdes (Marqueurs, Corrélation)"),
                 value = 5000, min = 0, step = 500),
    div(class="small text-muted",
        i18n$t("0 = pas de sous-échantillonnage (toutes les cellules). Réduit uniquement"),
        i18n$t(" les calculs de FindAllMarkers/Corrélation — l'objet complet reste intact"),
        i18n$t(" pour la visualisation et l'export."))
  )
}

# ── Server ────────────────────────────────────────────────────────────────────

mod_sc_pipeline_server <- function(id, global_data, shared_rv) {
  moduleServer(id, function(input, output, session) {
    # ── i18n proxy ──────────────────────────────────────────────────────────
    .tr <- function(key) {
      tr <- isolate(global_data$i18n)
      if (is.null(tr)) return(key)
      tryCatch(.strip_i18n_html(tr$t(key)), error = function(e) key)
    }

    observeEvent(global_data$language, {
      updateNumericInput(session, "qc_min_gene", label = .tr("Min Gènes"))
      updateNumericInput(session, "qc_max_gene", label = .tr("Max Gènes"))
      updateSliderInput(session, "qc_mt", label = .tr("% Mito Max"))
      updateCheckboxInput(session, "doublets_enabled", label = .tr("Filtrer les doublets (scDblFinder)"))
      updateSliderInput(session, "doublet_score_cutoff", label = .tr("Seuil score doublet"))
      updateSliderInput(session, "max_contamination", label = .tr("Contamination ambiante max"))
      updateCheckboxInput(session, "run_ambient", label = .tr("Décontamination ambiante (SoupX/celda) — opt-in"))
      updateSelectInput(session, "sketch_preset",
        label = .tr("Taille du sketch (vitesse vs precision)"),
        choices = setNames(
          c("fast", "light", "medium", "standard", "high", "max", "custom"),
          c(.tr("Rapide (test, 5 000 cellules)"),
            .tr("Léger (10 000 cellules)"),
            .tr("Moyen (25 000 cellules)"),
            .tr("Standard (50 000 cellules)"),
            .tr("Élevé (100 000 cellules)"),
            .tr("Max (dataset complet)"),
            .tr("Personnalisé"))
        ),
        selected = isolate(input$sketch_preset) %||% "standard"
      )
      updateSelectInput(session, "reduction_method",
        label = .tr("Méthode de Réduction (principale)"),
        choices = setNames(
          c("umap", "pca", "tsne", "dm", "harmony"),
          c("UMAP", .tr("PCA"), "t-SNE", .tr("Diffusion Maps"), "Harmony")
        ),
        selected = isolate(input$reduction_method) %||% "umap"
      )
      updateSliderInput(session, "pca_dim", label = .tr("Dims PCA"))
      updateNumericInput(session, "clust_res", label = .tr("Résolution"))
      updateSelectInput(session, "cluster_algo",
        label = .tr("Algorithme de clustering"),
        choices = setNames(
          c("1", "2", "3", "4"),
          c(.tr("Louvain (standard)"),
            .tr("Louvain (multilevel refinement)"),
            .tr("SLM (Smart Local Moving)"),
            .tr("Leiden (nécessite reticulate + leidenalg)"))
        ),
        selected = isolate(input$cluster_algo) %||% "1"
      )
      updateActionButton(session, "run_pipeline", label = .tr("Lancer Pipeline"))
      updateSelectInput(session, "bpcells_mode",
        label = .tr("Mode"),
        choices = setNames(
          c("auto", "always", "never"),
          c(.tr("Auto (recommandé)"),
            .tr("Toujours (forcer disque)"),
            .tr("Jamais (toujours RAM)"))
        ),
        selected = isolate(input$bpcells_mode) %||% "auto"
      )
      updateNumericInput(session, "max_cells_heavy", label = .tr("Max cellules/cluster pour analyses lourdes (Marqueurs, Corrélation)"))
    }, ignoreInit = TRUE)


    if (is.null(getOption("future.rng.onMisuse")))
      options(future.rng.onMisuse = "ignore")
    if (is.null(getOption("future.globals.maxSize")))
      options(future.globals.maxSize = 2 * 1024^3)

    # ── Mirror the RAM-safety cap to shared_rv so mod_sc_markers.R / mod_sc_corr.R
    #    (and the auto-pipeline in mod_sc.R) can read a single source of truth.
    #    Fires once at startup too (default ignoreInit=FALSE) so the cap is set
    #    even if the user never opens this accordion panel. ──────────────────────
    observeEvent(input$max_cells_heavy, {
      state_set(shared_rv, "max_cells_heavy",
                if (is.na(input$max_cells_heavy) || input$max_cells_heavy <= 0)
                  Inf else input$max_cells_heavy)
    })

    # ── P1-10 : timeline persistante des étapes (plan post-E2E #10) ──────────
    # Étapes figées du pipeline manuel ; .tl_set() marque l'étape courante
    # (running) et clôt les précédentes (done) ; .tl_end() clôt la course
    # (tout done, ou étape en erreur + restantes pending). L'état vit dans un
    # reactiveVal local (rendu) ET dans shared_rv (visible hors module —
    # checkpoints/rapport). Pas d'étape « harmonie/t-SNE » exécutée → statut
    # honnête « skip » en fin de course.
    .TL_STEPS <- c(qc = "QC — filtres", doublets = "QC — doublets / ambiant",
                   backend = "Backend stockage", norm = "Normalisation",
                   pca = "PCA", harmony = "Harmony", cluster = "Clustering",
                   reduction = "Réduction principale", tsne = "t-SNE secondaire")
    timeline_rv <- shiny::reactiveVal(NULL)
    .tl_set <- function(step, status = "running") {
      idx <- match(step, names(.TL_STEPS))
      if (is.na(idx)) return(invisible())
      st <- ifelse(seq_along(.TL_STEPS) < idx, "done",
            ifelse(seq_along(.TL_STEPS) == idx, status, "pending"))
      tl <- data.frame(step = unname(.TL_STEPS), status = st, stringsAsFactors = FALSE)
      timeline_rv(tl)
      state_set(shared_rv, "pipeline_timeline", tl)
    }
    .tl_end <- function(ok, failed_step = NULL) {
      cur_tl <- timeline_rv()
      if (isTRUE(ok)) {
        # les étapes marquées « skip » restent skip (honnêteté), le reste done
        st <- if (!is.null(cur_tl) && nrow(cur_tl) == length(.TL_STEPS))
          ifelse(cur_tl$status == "skip", "skip", "done")
          else rep("done", length(.TL_STEPS))
      } else {
        failed_idx <- match(failed_step, names(.TL_STEPS))
        if (is.null(failed_step) || is.na(failed_idx)) failed_idx <- 1L
        st <- ifelse(seq_along(.TL_STEPS) < failed_idx, "done",
              ifelse(seq_along(.TL_STEPS) == failed_idx, "error", "skip"))
      }
      tl <- data.frame(step = unname(.TL_STEPS), status = st, stringsAsFactors = FALSE)
      timeline_rv(tl)
      state_set(shared_rv, "pipeline_timeline", tl)
    }
    # Étape sautée (honnêteté : jamais « done » pour une étape non exécutée)
    .tl_skip <- function(step) {
      idx <- match(step, names(.TL_STEPS))
      if (is.na(idx)) return(invisible())
      tl <- timeline_rv()
      if (is.null(tl)) return(invisible())
      tl$status[idx] <- "skip"
      timeline_rv(tl)
      state_set(shared_rv, "pipeline_timeline", tl)
    }
    output$pipeline_timeline <- renderUI({
      tl <- timeline_rv()
      shiny::req(tl)
      rows <- lapply(seq_len(nrow(tl)), function(i) {
        mark <- switch(tl$status[i], done = "✓", running = "▶", error = "✗",
                       skip = "—", "○")
        cls <- switch(tl$status[i], done = "text-success", running = "text-primary",
                      error = "text-danger", skip = "text-muted", "text-muted")
        shiny::div(class = cls, style = "font-size:0.85em; line-height:1.5;",
                   shiny::HTML(mark), " ", tl$step[i])
      })
      shiny::div(class = "alert alert-light mt-1", style = "padding:8px 12px;",
                 shiny::tags$strong(style = "font-size:0.82em;",
                                    .tr("Pipeline — dernière exécution")),
                 rows)
    })

    # ── P1-6a : ExtendedTask + mirai (pattern auto-pipeline mod_sc.R) ───────
    # Pool UNIQUE de l'application (ensure_app_daemons → init_spatial_daemons,
    # preload APP_DAEMON_SOURCE_FILES dont sc_pipeline.R) ; bind_task_button
    # désactive le bouton pendant la course (anti double-clic) ; le tracker
    # streame le log du daemon (lignes [TL] step|status → timeline P1-10).
    ext_task_cls <- if (exists("ExtendedTask", where = asNamespace("shiny"), inherits = FALSE)) {
      get("ExtendedTask", envir = asNamespace("shiny"))
    } else if (exists("ExtendedTask", where = asNamespace("bslib"), inherits = FALSE)) {
      get("ExtendedTask", envir = asNamespace("bslib"))
    } else {
      get0("ExtendedTask", envir = globalenv())
    }
    manual_task <- if (!is.null(ext_task_cls) && requireNamespace("mirai", quietly = TRUE)) {
      # Déps du cœur NON présentes dans le preload daemons (APP_DAEMON_SOURCE_
      # FILES — spatial_async.R, périmètre Lane A intouché) : memory.R (clean_mem)
      # et sc_qc_doublets.R (étape M4) partent avec la tâche, chemins absolus
      # calculés côté session ; source() est idempotent et les garde-fous
      # exists() du cœur restent honnêtes.
      core_deps <- c(
        file.path(ts_project_root(), "tools", "pkg_skeleton", "R", "memory.R"),
        file.path(ts_project_root(), "tools", "pkg_skeleton", "R", "sc_qc_doublets.R")
      )
      t <- ext_task_cls$new(function(obj, params, log_file, deps = core_deps) {
        ensure_app_daemons()
        mirai::mirai(
          {
            for (f in deps) if (file.exists(f)) source(f, local = FALSE)
            run_sc_manual_pipeline_core(obj, params, log_file = log_file)
          },
          obj = obj, params = params, log_file = log_file, deps = core_deps,
          .timeout = TS_SC_AUTO_PIPELINE_TIMEOUT_S * 1000L
        )
      })
      tryCatch(bslib::bind_task_button(t, "run_pipeline"), error = function(e) NULL)
      t
    } else {
      NULL
    }

    # Applique le résultat du cœur dans la session (async ET repli synchrone) :
    .apply_manual_result <- function(res) {
      if (!isTRUE(res$ok)) {
        for (n in res$notifications) {
          showNotification(n$msg, type = n$type, duration = 6)
        }
        err <- res$error
        .tl_end(ok = FALSE, failed_step = res$step)
        showNotification(paste(.tr("Erreur pipeline:"),
                               if (inherits(err, "error")) err$message else as.character(err)),
                         type = "error", duration = 10)
        return(invisible())
      }
      # P1-7 : point de pose ASR du pipeline MANUEL — le moteur
      # session_checkpoint (sourcé app.R:106) n'avait aucun publish depuis ce
      # module (imports, auto-pipeline mod_sc.R:1772, spatial en avaient) :
      # un reload page après 30+ min de pipeline perdait tout. Verbatim du
      # point de pose mod_sc.R:1772, APRÈS l'écriture de sc_obj.
      global_data$sc_obj <- res$obj
      if (!is.null(res$qc_snapshot)) state_set(shared_rv, "qc_snapshot", res$qc_snapshot)
      if (!is.null(res$qc_dist_before)) {
        state_set(shared_rv, "qc_dist_before", res$qc_dist_before)
        state_set(shared_rv, "qc_dist_after", res$qc_dist_after)
      }
      for (pr in res$provenance) provenance_append(shared_rv, pr)
      if (!is.null(res$bpcells_conv) && !isTRUE(res$bpcells_conv$already_disk)) {
        # One folder per conversion accumulates in tempdir() — cleanup at
        # session close (comportement historique du module).
        session$onSessionEnded(function() unlink(res$bpcells_conv$dir, recursive = TRUE))
      }
      for (n in res$notifications) {
        showNotification(n$msg, type = n$type, duration = if (identical(n$type, "warning")) 8 else 6)
      }
      if (is.function(global_data$session_checkpoint_publish)) global_data$session_checkpoint_publish("pipeline_step_complete")
      .tl_end(ok = TRUE)
      showNotification(.tr("✓ Pipeline terminé"), type = "message")
    }

    observeEvent(input$run_pipeline, {
      req(global_data$sc_obj)
      # garde re-entrance : une course déjà active ne relance pas
      if (!is.null(manual_task) && identical(manual_task$status(), "running")) return(invisible())
      obj <- shiny::isolate(global_data$sc_obj)

      params <- list(
        qc_min_gene = input$qc_min_gene, qc_max_gene = input$qc_max_gene,
        qc_mt = input$qc_mt,
        doublets_enabled = isTRUE(input$doublets_enabled),
        doublet_score_cutoff = input$doublet_score_cutoff %||% .PIPELINE_DOUBLET_CUTOFF_DEFAULT,
        max_contamination = input$max_contamination %||% .PIPELINE_MAX_CONTAMINATION_DEFAULT,
        run_ambient = isTRUE(input$run_ambient),
        bpcells_mode = input$bpcells_mode %||% "auto",
        bpcells_threshold = input$bpcells_threshold %||% .BPCELLS_AUTO_THRESHOLD,
        norm_method = input$norm_method %||% "log",
        pca_dim = input$pca_dim, clust_res = input$clust_res,
        cluster_algo = input$cluster_algo,
        reduction_method = input$reduction_method %||% "umap",
        seed = .SC_PIPELINE_SEED
      )

      pool_ok <- isTRUE(tryCatch(ensure_app_daemons(), error = function(e) FALSE))
      if (!is.null(manual_task) && pool_ok) {
        .tl_set("qc")
        # log de daemon → tracker → timeline (lignes [TL] step|status)
        manual_log_file <- tempfile(pattern = "sc_manual_log_", fileext = ".txt")
        manual_tracker <- create_reactive_tracker(session, manual_log_file, interval_ms = 800)
        observe({
          lines <- manual_tracker()
          if (length(lines)) for (ln in lines) {
            if (!grepl("^[[:space:]]*\\[TL\\] ", ln)) next
            body <- sub("^[[:space:]]*\\[TL\\] ", "", ln)
            stp <- sub("\\|.*$", "", body)
            sta <- sub("^[^|]*\\|", "", body)
            if (identical(stp, "END")) {
              if (grepl("^error", sta)) .tl_end(FALSE, failed_step = sub("^error:", "", sta))
              else .tl_end(TRUE)
            } else if (identical(sta, "running")) .tl_set(stp)
            else if (identical(sta, "skip")) .tl_skip(stp)
          }
        })
        manual_task$invoke(obj, params, manual_log_file)
      } else {
        # Repli synchrone TRACÉ (jamais muet — convention run_job) : le pool
        # mirai est indisponible ; l'interface se re-bloque pour la course,
        # la timeline est alimentée en direct par on_progress.
        showNotification(.tr("Worker mirai indisponible — exécution synchrone : l'interface restera bloquée pendant la course."),
                         type = "warning", duration = 8)
        .tl_set("qc")
        p <- shiny::Progress$new(); on.exit(p$close())
        p$set(message = .tr("Pipeline..."), value = 0)
        res <- tryCatch(
          run_sc_manual_pipeline_core(obj, params, on_progress = function(stp, sta) {
            if (identical(stp, "END")) {
              if (grepl("^error", sta)) .tl_end(FALSE, failed_step = sub("^error:", "", sta))
              else .tl_end(TRUE)
            } else if (identical(sta, "running")) .tl_set(stp)
            else if (identical(sta, "skip")) .tl_skip(stp)
          }),
          error = function(e) list(ok = FALSE, error = e, notifications = list(),
                                   step = "qc", obj = NULL)
        )
        .apply_manual_result(res)
      }
    })

    # Achèvement de la tâche async : application du résultat dans la session.
    if (!is.null(manual_task)) {
      observeEvent(manual_task$status(), {
        st <- manual_task$status()
        if (identical(st, "success")) {
          .apply_manual_result(manual_task$result())
        } else if (identical(st, "error")) {
          # échec hors cœur (daemon mort, timeout) : pas de ligne [TL] END —
          # l'étape faillible = la dernière « running » de la timeline.
          tl <- timeline_rv()
          failed <- if (!is.null(tl)) {
            idx <- which(tl$status == "running")
            if (length(idx)) names(.TL_STEPS)[idx[1]] else NULL
          } else NULL
          .tl_end(ok = FALSE, failed_step = failed)
          err <- tryCatch(manual_task$result(), error = function(e) e)
          showNotification(paste(.tr("Erreur pipeline:"),
                                 if (inherits(err, "error")) err$message else as.character(err)),
                           type = "error", duration = 10)
        }
      })
    }
  })
}
