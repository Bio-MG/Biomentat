# =============================================================================
# mod_sc_velocity.R — RNA Velocity Phase 3A Hardening + Stage 8 (3B-1) + 9 (3B-2)
# =============================================================================
# Strict RDS/MTX import, matrix validation, orientation, cell/gene alignment,
# precomputed UMAP vector validation, phase portraits, provenance export,
# reset when Seurat object changes. No inference, no loom/h5ad, no projection.
# Stage 8 : le module ORCHESTRE uniquement — la validation, le statut canonique
# (velocity_status_labels), la provenance et les exports residuent dans
# R/sc/sc_velocity.R. Validite technique affichee explicitement ; aucune
# validite biologique n'est suggeree.
# Stage 9 : les visualisations sont des CONSOMMATRICES pures du resultat
# canonique (assert_velocity_result) — fleches uniquement depuis des vecteurs
# valides, portee de donnees (preview/full) affichee, exports traces par
# analysis_id. Une figure ne cree, n'infere et ne repare jamais un resultat.
# =============================================================================

# Pure internal snapshot helper for drive observability (W1.1 #1).
# Returns a bounded descriptor with logical scalar `has_vectors` and
# character scalar `optin_status` in {idle, running, done, error, unavailable}.
.velocity_drive_probe_fields <- function(velocity_state, optin_status = NULL) {
  valid_statuses <- c("idle", "running", "done", "error", "unavailable")
  has_vectors <- if (is.null(velocity_state)) {
    FALSE
  } else {
    isTRUE(!is.null(velocity_state$umap_vectors))
  }
  resolved_status <- NULL
  if (is.character(optin_status) && length(optin_status) == 1L &&
      !is.na(optin_status) && optin_status %in% valid_statuses) {
    resolved_status <- optin_status
  } else if (!is.null(velocity_state)) {
    st_val <- velocity_state$optin_status
    if (is.character(st_val) && length(st_val) == 1L &&
        !is.na(st_val) && st_val %in% valid_statuses) {
      resolved_status <- st_val
    } else if (!is.null(velocity_state$optin_error)) {
      resolved_status <- "error"
    } else if (identical(velocity_state$vector_validation$field, "scvelo_optin") &&
               isTRUE(velocity_state$vector_validation$ok)) {
      resolved_status <- "done"
    }
  }
  if (is.null(resolved_status)) {
    resolved_status <- if (is.null(velocity_state)) "unavailable" else "idle"
  }
  list(
    has_vectors = has_vectors,
    optin_status = resolved_status
  )
}

mod_sc_velocity_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "alert alert-light", style = "font-size:0.9em;border-left:3px solid #2980B9;",
        i18n$t("RNA Velocity — import strict des matrices spliced/unspliced alignees sur l'objet Seurat. "),
        i18n$t("Aucune inference n'est effectuee dans l'application.")),

    radioButtons(ns("velocity_import_mode"), i18n$t("Source d'import"),
                  choices = setNames(c("rds", "mtx", "loom"),
                              c(.tr_plain("RDS combine (spliced + unspliced)"),
                                .tr_plain("Matrix Market (MTX + barcodes + features)"),
                                .tr_plain("Loom / H5AD (velocyto, scVelo)"))),
                  selected = "rds"),

    conditionalPanel(
      condition = "input.velocity_import_mode == 'rds'", ns = ns,
      fileInput(ns("velocity_rds_file"), i18n$t("Fichier RDS velocity (liste spliced/unspliced)"),
                accept = c(".rds", ".RDS", ".rda", ".RData"), width = "100%"),
      div(class = "small text-muted mb-2",
           i18n$t("Le RDS doit être une liste nommée avec spliced et unspliced. Noms inattendus bloqués."))
    ),

    conditionalPanel(
      condition = "input.velocity_import_mode == 'loom'", ns = ns,
      # Roadmap SC FUNCTION_TEST M0 : lecteur loom/h5ad (hdf5r, pre-computed).
      # Aucune inference : les couches spliced/unspliced du fichier sont
      # lues telles quelles et passent la meme validation que RDS/MTX.
      fileInput(ns("velocity_loom_file"), i18n$t("Fichier velocity .loom / .h5ad"),
                accept = c(".loom", ".h5ad"), width = "100%"),
      div(class = "small text-muted mb-2",
           i18n$t("Lecture stricte des couches spliced/unspliced (+ ambiguous si present). Aucune inference ; les symboles de gènes dupliqués sont rendus uniques (avertissement affiché)."))
    ),

    conditionalPanel(
      condition = "input.velocity_import_mode == 'mtx'", ns = ns,
      h6(i18n$t("Matrice spliced (MTX)"), style = "font-weight:bold;"),
      fileInput(ns("velocity_mtx_spliced"), i18n$t("Fichier .mtx (spliced)"), accept = c(".mtx", ".gz"), width = "100%"),
      fileInput(ns("velocity_mtx_spliced_barcodes"), i18n$t("Barcodes spliced"), accept = c(".tsv", ".txt", ".gz"), width = "100%"),
      fileInput(ns("velocity_mtx_spliced_features"), i18n$t("Features/genes spliced"), accept = c(".tsv", ".txt", ".gz"), width = "100%"),
      h6(i18n$t("Matrice unspliced (MTX)"), style = "font-weight:bold;"),
      fileInput(ns("velocity_mtx_unspliced"), i18n$t("Fichier .mtx (unspliced)"), accept = c(".mtx", ".gz"), width = "100%"),
      fileInput(ns("velocity_mtx_unspliced_barcodes"), i18n$t("Barcodes unspliced"), accept = c(".tsv", ".txt", ".gz"), width = "100%"),
      fileInput(ns("velocity_mtx_unspliced_features"), i18n$t("Features/genes unspliced"), accept = c(".tsv", ".txt", ".gz"), width = "100%"),
      h6(i18n$t("Matrice ambiguous (optionnel)"), style = "font-weight:bold;"),
      fileInput(ns("velocity_mtx_ambiguous"), i18n$t("Fichier .mtx (ambiguous)"), accept = c(".mtx", ".gz"), width = "100%"),
      fileInput(ns("velocity_mtx_ambiguous_barcodes"), i18n$t("Barcodes ambiguous"), accept = c(".tsv", ".txt", ".gz"), width = "100%"),
      fileInput(ns("velocity_mtx_ambiguous_features"), i18n$t("Features/genes ambiguous"), accept = c(".tsv", ".txt", ".gz"), width = "100%"),
      numericInput(ns("velocity_feature_column"),
                    i18n$t("Colonne d'identifiants dans features.tsv"),
                    value = 1, min = 1, step = 1, width = "100%"),
      div(class = "small text-muted mb-2",
           HTML(paste(i18n$t("Les fichiers 10x contiennent souvent gene_id (colonne 1) ET symbole "),
                 i18n$t("(colonne 2). La colonne choisie doit correspondre aux rownames de "),
                 i18n$t("l'objet Seurat — la premiere colonne n'est jamais supposee correcte."))))
    ),

    hr(),
    selectInput(
      ns("velocity_orientation"),
      i18n$t("Orientation des matrices"),
      choices = setNames(c("auto_strict", "genes_x_cells", "cells_x_genes"),
        c(.tr_plain("Détection automatique stricte"),
          .tr_plain("Genes x cellules"),
          .tr_plain("Cellules x genes"))),
      selected = "auto_strict"
    ),
    checkboxInput(
      ns("velocity_strip_cell_suffix"),
      i18n$t("Retirer les suffixes de barcodes (ex. -1, -2)"),
      value = FALSE
    ),
    checkboxInput(
      ns("velocity_strip_gene_version"),
      i18n$t("Retirer les suffixes de version Ensembl (ex. .1, .2)"),
      value = FALSE
    ),
    checkboxInput(
      ns("velocity_allow_low_overlap"),
      i18n$t("Forcer l'alignement malgre un recouvrement < 80%"),
      value = FALSE
    ),
    div(
      class = "small text-muted",
      i18n$t("Les corrections d'identifiants sont desactivees par defaut. "),
      " ",
      i18n$t("Une collision bloque la validation.")
    ),
    hr(),
    actionButton(ns("velocity_validate"), i18n$t("Valider les matrices velocity"),
                 class = "btn-primary w-100", icon = icon("check")),
    div(class = "small text-muted mt-1", textOutput(ns("velocity_status"))),
    hr(),
    # O13 (2026-10-27) : inférence scVelo OPT-IN — bouton dédié, SEULE voie
    # d'inférence du domaine, JAMAIS l'autopipeline ; bannière dédiée sur le
    # résultat produit (uiOutput conditionnel ci-dessous).
    actionButton(ns("velocity_scvelo_run"), i18n$t("Calculer les vecteurs (scVelo — opt-in)"),
                 class = "btn-warning w-100", icon = icon("paper-plane")),
    uiOutput(ns("velocity_scvelo_banner")),
    # Journal leger (2026-11-17, item #10) : place dans la sidebar sous le
    # bouton opt-in (preference operateur) au lieu de la carte de resultats.
    # Borne cote serveur (50 lignes) + hauteur max scrollable pour preserver
    # l'ergonomie de la sidebar.
    div(
      class = "mt-2",
      tags$strong(class = "small", i18n$t("Console de Log")),
      div(
        class = "small",
        style = "max-height:140px;overflow-y:auto;",
        verbatimTextOutput(ns("velocity_log"), placeholder = TRUE)
      )
    ),
    hr(),
    selectInput(ns("velocity_gene"), i18n$t("Gene pour phase portrait"), choices = character(0)),
    selectInput(ns("velocity_reduction"), i18n$t("Reduction pour vecteurs"), choices = c("umap", "pca")),
    # 2C (O28) : les 9 downloadButton empilés (constat 2C, packet round-7 :
    # `grep -c downloadButton` -> 9) sont regroupés en 2 dropdownButton compacts
    # (spéc §2.8) — « 📊 Graphiques » (phase portrait + embedding, PNG/PDF) et
    # « 💾 Données & Provenance » (provenance, résumé validation, alignement,
    # vecteurs par cellule, RDS complet). Les ID des handlers dl_velocity_* sont
    # INCHANGÉS : aucun contrat serveur ne bouge, seul l'habillage UI change.
    shinyWidgets::dropdownButton(
      inputId = ns("velocity_export_plots_menu"),
      label = i18n$t("📊 Exporter les Graphiques"),
      circle = FALSE, status = "info", width = "100%",
      div(class = "d-grid gap-2",
          downloadButton(ns("dl_velocity_phase_png"), i18n$t("Exporter phase portrait (PNG)"), class = "btn-sm btn-info w-100"),
          downloadButton(ns("dl_velocity_phase_pdf"), i18n$t("Exporter phase portrait (PDF)"), class = "btn-sm btn-info w-100"),
          downloadButton(ns("dl_velocity_embedding_png"), i18n$t("Exporter embedding (PNG)"), class = "btn-sm btn-info w-100"),
          downloadButton(ns("dl_velocity_embedding_pdf"), i18n$t("Exporter embedding (PDF)"), class = "btn-sm btn-info w-100"))),
    shinyWidgets::dropdownButton(
      inputId = ns("velocity_export_data_menu"),
      label = i18n$t("💾 Exporter les Données & Provenance"),
      circle = FALSE, status = "info", width = "100%",
      div(class = "d-grid gap-2",
          downloadButton(ns("dl_velocity_provenance"), i18n$t("Exporter provenance (CSV)"), class = "btn-sm btn-info w-100"),
          downloadButton(ns("dl_velocity_summary"), i18n$t("Exporter resume validation (CSV)"), class = "btn-sm btn-info w-100"),
          downloadButton(ns("dl_velocity_mapping"), i18n$t("Exporter alignement cellules/genes (CSV)"), class = "btn-sm btn-info w-100"),
          downloadButton(ns("dl_velocity_cell_vectors"), i18n$t("Exporter vecteurs par cellule (CSV)"), class = "btn-sm btn-info w-100"),
          downloadButton(ns("dl_velocity_result_rds"), i18n$t("Exporter resultat valide (RDS)"), class = "btn-sm btn-info w-100")))
  )
}

mod_sc_velocity_output_ui <- function(id) {
  ns <- NS(id)
  # CCC_NAV_DRAWER vague B jalon J3 (2026-11-12) : généralisation du pattern
  # tiroir optimisé du CCC (5f0a906 + bc9c09a) au panneau Velocity — candidat
  # C1 du packet gate 2 (GO opérateur vague B « GO J3 + J4 »). ⚠️ CE FICHIER
  # est EN VOL lane scvelo (fix opt-in loom, session concurrente) : le commit
  # J3 stage SES SEULS HUNKS (le hunk étranger loom reste NON commité).
  # JUMELLE : la carte .velocity_drawer_groups est redéclarée EN LOCAL dans
  # mod_sc_velocity_server (gel test-velocity-contract-freeze.R:152 — aucune
  # constante top-level) ; égalité verrouillée par test-mod-sc-velocity-ui.R.
  .velocity_drawer_groups <- list(
    "Vues velocity" = list(
      list(value = "velocity_tab_phase", label = "Phase portrait"),
      list(value = "velocity_tab_vectors", label = "Vecteurs UMAP"),
      list(value = "velocity_tab_coverage", label = "Couverture vecteurs"),
      list(value = "velocity_tab_qc", label = "QC alignement")))
  card(
    full_screen = TRUE,
    # Pastille ouvrant le tiroir (motif 2L CCC : dans l'en-tête de carte) :
    card_header(
      div(class = "d-flex justify-content-between align-items-center gap-2",
          tags$span(i18n$t("RNA Velocity — Phase portrait & vecteurs pre-calcules")),
          tags$button(type = "button",
                      class = "btn btn-sm btn-outline-secondary flex-shrink-0",
                      `data-bs-toggle` = "offcanvas",
                      `data-bs-target` = paste0("#", ns("velocity_drawer")),
                      `aria-controls` = ns("velocity_drawer"),
                      icon("list-ul"), " ", i18n$t("Navigation & filtres")))
    ),
    navset_tab(
      id = ns("velocity_out_tabs"),
      nav_panel(i18n$t("Phase portrait"), value = "velocity_tab_phase",
                plotOutput(ns("velocity_phase_plot"), height = "500px")),
      nav_panel(i18n$t("Vecteurs UMAP"), value = "velocity_tab_vectors",
                plotOutput(ns("velocity_embedding_plot"), height = "550px"),
                div(class = "small text-muted", textOutput(ns("velocity_vector_status")))),
      nav_panel(i18n$t("Couverture vecteurs"), value = "velocity_tab_coverage",
                plotOutput(ns("velocity_coverage_plot"), height = "450px")),
      nav_panel(i18n$t("QC alignement"), value = "velocity_tab_qc",
                plotOutput(ns("velocity_qc_plot"), height = "450px"))
    ),
    # Tiroir offcanvas (DEUXIÈME enfant de la carte — motif comm_drawer /
    # milo_drawer / sccoda_drawer / cross_drawer) :
    tags$div(
      id = ns("velocity_drawer"), class = "offcanvas offcanvas-end", tabindex = "-1",
      `aria-labelledby` = ns("velocity_drawer_title"),
      tags$div(class = "offcanvas-header",
               h5(id = ns("velocity_drawer_title"), class = "offcanvas-title",
                  icon("list-ul"), " ", i18n$t("Navigation Velocity")),
               # .tr_plain (strip span) — JAMAIS de $t() brut en contexte
               # attribut (BUG A11Y « Fermer », motif bc9c09a) :
               tags$button(type = "button", class = "btn-close",
                           `data-bs-dismiss` = "offcanvas",
                           `aria-label` = .tr_plain("Fermer"))),
      tags$div(class = "offcanvas-body",
        h6(i18n$t("État du résultat")),
        div(class = "small text-muted", `aria-live` = "polite",
            uiOutput(ns("velocity_drawer_summary"))),
        hr(),
        uiOutput(ns("velocity_drawer_nav"))
      )
    )
  )
}

# ── 2D (O28) : gène par défaut du phase portrait — heuristique ordonnée ──────
# À HEAD le sélecteur était rempli avec `selected = gene_choices[1]` (ordre
# alphabétique des gènes alignés) — constat 2D, packet round-7 ; spéc §2.9.
# Heuristique :
#   1. premier gène du top 10 des marqueurs différentiels significatifs
#      (p-adj < 0.05, tri croissant) présent dans les gènes alignés ;
#   2. sinon premier HVG (VariableFeatures) présent dans les gènes alignés ;
#   3. sinon gène de plus forte variance dans `unspliced` (matrice dense
#      uniquement — une matrice sparse ne doit JAMAIS être densifiée ici ;
#      variance calculée par blocs de lignes pour borner la mémoire) ;
#   4. fallback : gene_choices[1] (comportement d'antan).
# Fonction pure, testable headless (test-mod-sc-velocity.R).
.velocity_phase_default_gene <- function(gene_choices,
                                         markers_df = NULL,
                                         var_features = NULL,
                                         unspliced = NULL) {
  if (!length(gene_choices)) return(character(0))
  # Priorité 1 — marqueurs différentiels significatifs (top 10 par p-adj).
  if (is.data.frame(markers_df) && nrow(markers_df) > 0L &&
      all(c("gene", "p_val_adj") %in% names(markers_df))) {
    sig <- markers_df[!is.na(markers_df$p_val_adj) & markers_df$p_val_adj < 0.05, ,
                      drop = FALSE]
    sig <- sig[order(sig$p_val_adj), , drop = FALSE]
    hit <- intersect(as.character(head(sig$gene, 10L)), gene_choices)
    if (length(hit)) return(hit[1])
  }
  # Priorité 2 — HVG.
  hit <- intersect(as.character(var_features), gene_choices)
  if (length(hit)) return(hit[1])
  # Priorité 3 — variance maximale dans unspliced (matrice de base uniquement).
  if (is.matrix(unspliced) && nrow(unspliced) > 0L && ncol(unspliced) > 1L) {
    v <- tryCatch({
      m <- unspliced[intersect(rownames(unspliced), gene_choices), , drop = FALSE]
      # NB : jamais de `return(...)` dans ce bloc tryCatch — un return
      # s'échapperait de la FONCTION entière (essayé, mesuré) ; on rend la
      # valeur et le `any(ok)` en aval route vers le fallback.
      if (nrow(m)) {
        chunk <- 2000L
        n     <- nrow(m)
        out   <- numeric(n)
        for (start in seq.int(1L, n, by = chunk)) {
          idx <- seq.int(start, min(start + chunk - 1L, n))
          blk <- m[idx, , drop = FALSE]
          mu  <- rowMeans(blk)
          out[idx] <- rowSums((blk - mu)^2) / max(1L, ncol(blk) - 1L)
        }
        names(out) <- rownames(m)
        out
      } else numeric(0)
    }, error = function(e) numeric(0))
    ok <- !is.na(v)
    if (any(ok)) {
      return(names(v)[ok][which.max(v[ok])])
    }
  }
  # Fallback — comportement d'antan (premier gène des choix).
  gene_choices[1]
}

mod_sc_velocity_server <- function(id, global_data, shared_rv = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # ── i18n proxy ──────────────────────────────────────────────────────────
    .tr <- function(key) {
      tr <- isolate(global_data$i18n)
      if (is.null(tr)) return(key)
      tryCatch(.strip_i18n_html(tr$t(key)), error = function(e) key)
    }

    # ── helpers contrat Stage 8 ────────────────────────────────────────────
    # Libelle francais d'un etat de validite (fallback = code d'etat).
    .status_label <- function(st) {
      lab <- velocity_status_labels()
      if (is.character(st) && length(st) == 1L && st %in% names(lab)) lab[[st]] else st
    }

    # Fichiers sources declares (NOMS ORIGINAUX uniquement — jamais les
    # chemins locaux complets) + options de lecture, pour input_summary.
    .collect_input_files <- function(mode) {
      if (identical(mode, "rds")) {
        if (is.null(input$velocity_rds_file)) return(list())
        return(list(rds = input$velocity_rds_file$name))
      }
      if (identical(mode, "loom")) {
        if (is.null(input$velocity_loom_file)) return(list())
        return(list(loom = input$velocity_loom_file$name))
      }
      files <- list()
      for (slot in c("velocity_mtx_spliced", "velocity_mtx_spliced_barcodes",
                     "velocity_mtx_spliced_features", "velocity_mtx_unspliced",
                     "velocity_mtx_unspliced_barcodes", "velocity_mtx_unspliced_features",
                     "velocity_mtx_ambiguous", "velocity_mtx_ambiguous_barcodes",
                     "velocity_mtx_ambiguous_features")) {
        fi <- input[[slot]]
        if (!is.null(fi)) files[[slot]] <- fi$name
      }
      files$feature_column <- as.character(input$velocity_feature_column %||% 1L)
      files
    }

    # Plot d'erreur standardise (message de validation exploitable).
    .error_plot <- function(e) {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = conditionMessage(e), size = 5, colour = "red") +
        ggplot2::theme_void()
    }

    .velocity_has_reduction <- function(obj, sel_red) {
      if (is.null(obj) || is.null(sel_red) || !nzchar(sel_red)) return(FALSE)
      reds <- tryCatch(Seurat::Reductions(obj), error = function(e) names(obj@reductions))
      isTRUE(sel_red %in% reds)
    }

    .velocity_missing_reduction_msg <- function(sel_red) {
      sprintf(
        "Réduction '%s' absente de l'objet Seurat — exécutez d'abord le Pipeline SC (PCA / UMAP) avant de lancer la vélocité.",
        as.character(sel_red %||% "umap")
      )
    }

    # Coordonnees de base pour la vue embedding : embedding RDS partiel
    # prioritaire, sinon reduction Seurat restreinte aux cellules velocity
    # validees. NULL si introuvables (le demandeur affiche alors un message
    # exploitable — jamais de coordonnees substituees).
    .base_coords_for <- function(red) {
      if (!is.null(velocity_state$partial_embedding)) {
        return(velocity_state$partial_embedding)
      }
      obj <- global_data$sc_obj
      if (is.null(obj) || is.null(velocity_state$result) ||
          !.velocity_has_reduction(obj, red)) {
        return(NULL)
      }
      emb <- Seurat::Embeddings(obj, reduction = red)
      common_cells <- intersect(velocity_state$result$cell_names, rownames(emb))
      if (!length(common_cells)) return(NULL)
      emb[common_cells, seq_len(min(2L, ncol(emb))), drop = FALSE]
    }



    velocity_state <- reactiveValues(
      inputs = NULL,
      validated = NULL,
      result = NULL,
      alignment = NULL,
      umap_vectors = NULL,
      partial_embedding = NULL,
      embedding_alignment = NULL,
      vector_validation = NULL,
      plot_cache = list(),
      export_cache = list(),
      object_fingerprint = NULL,
      has_vectors = FALSE,
      optin_status = "idle",
      optin_error = NULL,
      status_msg = NULL,
      source_label = NULL
    )

    observeEvent(global_data$language, {
      updateRadioButtons(session, "velocity_import_mode",
        label = .tr("Source d'import"),
        choices = setNames(
          c("rds", "mtx", "loom"),
          c(.tr("RDS combine (spliced + unspliced)"), .tr("Matrix Market (MTX + barcodes + features)"),
            .tr("Loom / H5AD (velocyto, scVelo)"))
        ),
        selected = isolate(input$velocity_import_mode) %||% "rds"
      )
      updateSelectInput(session, "velocity_orientation",
        label = .tr("Orientation des matrices"),
        choices = setNames(
          c("auto_strict", "genes_x_cells", "cells_x_genes"),
          c(.tr("Détection automatique stricte"), .tr("Genes x cellules"), .tr("Cellules x genes"))
        ),
        selected = isolate(input$velocity_orientation) %||% "auto_strict"
      )
      updateCheckboxInput(session, "velocity_strip_cell_suffix", label = .tr("Retirer les suffixes de barcodes (ex. -1, -2)"))
      updateCheckboxInput(session, "velocity_strip_gene_version", label = .tr("Retirer les suffixes de version Ensembl (ex. .1, .2)"))
      updateCheckboxInput(session, "velocity_allow_low_overlap", label = .tr("Forcer l'alignement malgre un recouvrement < 80%"))
      updateActionButton(session, "velocity_validate", label = .tr("Valider les matrices velocity"))
      updateSelectizeInput(session, "velocity_gene", label = .tr("Gene pour phase portrait"), server = TRUE)
      updateSelectInput(session, "velocity_reduction", label = .tr("Reduction pour vecteurs"))
    }, ignoreInit = TRUE)
    velocity_status_rv <- reactiveVal(.tr("En attente d'import velocity..."))

    output$velocity_status <- renderText({ velocity_status_rv() })

    # Journal borne (2026-11-17) : 50 lignes max, horodatage local. Un seul
    # observeEvent sur le statut capture TOUTES les transitions de validation
    # et d'opt-in sans toucher aux 8 sites d'ecriture existants.
    velocity_log_rv <- reactiveVal(character(0))
    .velocity_log <- function(...) {
      ligne <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                       paste0(..., collapse = ""))
      velocity_log_rv(utils::tail(c(velocity_log_rv(), ligne), 50L))
      invisible(NULL)
    }
    observeEvent(velocity_status_rv(), {
      .velocity_log(velocity_status_rv())
    }, ignoreInit = FALSE)
    output$velocity_log <- renderText({ .vel_sanitize_text(paste(velocity_log_rv(), collapse = "\n")) })

    # ── CCC_NAV_DRAWER vague B jalon J3 (2026-11-12) : tiroir « Navigation » ──
    # Pattern 5f0a906 + bc9c09a, généralisé comme Milo/scCODA/Vues croisées.
    # JUMELLE de la déclaration output_ui (égalité verrouillée
    # test-mod-sc-velocity-ui.R) — aucune constante top-level (gel
    # test-velocity-contract-freeze.R:152).
    .velocity_drawer_groups <- list(
      "Vues velocity" = list(
        list(value = "velocity_tab_phase", label = "Phase portrait"),
        list(value = "velocity_tab_vectors", label = "Vecteurs UMAP"),
        list(value = "velocity_tab_coverage", label = "Couverture vecteurs"),
        list(value = "velocity_tab_qc", label = "QC alignement")))

    # Résumé d'état du tiroir : le statut VIVANT du panneau (MÊME source
    # velocity_status_rv — zéro divergence possible) + disponibilité des
    # vecteurs opt-in consommée de l'état (jamais inventée).
    output$velocity_drawer_summary <- renderUI({
      global_data$language  # re-rendu au changement de langue (motif du module)
      li <- function(n, lab) tags$li(tags$strong(n), " ", lab)
      vec_ok <- !is.null(velocity_state$umap_vectors)
      tags$ul(class = "list-unstyled mb-0",
              tags$li(velocity_status_rv()),
              li(if (vec_ok) "\u2713" else "\u2014", .tr("Vecteurs UMAP")))
    })

    # Navigation du tiroir rendue DYNAMIQUEMENT (motif P1+P2 comm_drawer) :
    output$velocity_drawer_nav <- renderUI({
      global_data$language  # re-rendu au changement de langue (motif du module)
      sel <- input$velocity_out_tabs %||% ""
      lapply(names(.velocity_drawer_groups), function(grp) {
        items <- .velocity_drawer_groups[[grp]]
        tags$div(class = "mb-3",
          h6(.tr(grp),
             tags$span(class = "badge text-bg-secondary ms-1", length(items))),
          lapply(items, function(pn) {
            active <- identical(sel, pn$value)
            actionButton(ns(paste0("velocity_jump_", pn$value)),
                         label = if (active)
                           tagList(tags$span(class = "text-primary", "\u25B8 "),
                                   .tr(pn$label))
                         else .tr(pn$label),
                         class = paste("btn btn-link btn-sm text-start p-0 d-block",
                                       if (active) "fw-bold"))
          }))
      })
    })

    # Sauts directs : un observateur par panneau (motif 2L CCC verbatim).
    lapply(.velocity_drawer_groups, function(grp) {
      lapply(grp, function(pn) {
        observeEvent(input[[paste0("velocity_jump_", pn$value)]], {
          nav_select("velocity_out_tabs", selected = pn$value)
        })
      })
    })

    output$velocity_vector_status <- renderText({
      parts <- character(0)
      res <- velocity_state$result
      if (!is.null(res) && !is.null(res$status)) {
        parts <- c(parts, paste0("Statut : ", .status_label(res$status)))
      }
      ali <- velocity_state$embedding_alignment
      if (!is.null(ali)) {
        parts <- c(parts, sprintf(
          "Coordonnees UMAP : %d/%d cellules couvertes (%d manquantes).",
          ali$n_embedding_matched, ali$n_velocity_cells, ali$n_embedding_missing
        ))
      }
      if (!is.null(velocity_state$umap_vectors)) {
        parts <- c(parts, sprintf(
          "%d vecteurs valides (%s).",
          nrow(velocity_state$umap_vectors),
          velocity_state$vector_validation$field %||% "?"
        ))
      } else {
        parts <- c(parts,
          "Aucun vecteur velocity pre-calcule valide — fleches desactivees.")
      }
      paste(parts, collapse = " ")
    })

    # Keep reduction choices in sync with Seurat object
    observeEvent(global_data$sc_obj, {
      obj <- global_data$sc_obj
      if (is.null(obj)) {
        updateSelectInput(session, "velocity_reduction", choices = character(0), selected = character(0))
        return()
      }
      reds <- names(obj@reductions)
      if (!length(reds)) reds <- c("umap", "pca")
      updateSelectInput(session, "velocity_reduction",
                        choices = reds,
                        selected = if ("umap" %in% reds) "umap" else reds[1])
    }, ignoreInit = TRUE)

    # Reset when Seurat object changes (BUG 9) — no velocity state may
    # survive from a previous Seurat object. If the new Seurat object has the
    # exact same cell/gene fingerprint (e.g. metadata/annotation update), keep
    # the validated velocity bundle intact instead of wiping it.
    observeEvent(global_data$sc_obj, {
      new_fp <- if (is.null(global_data$sc_obj)) NULL else {
        tryCatch(velocity_object_fingerprint(global_data$sc_obj), error = function(e) NULL)
      }
      if (!is.null(new_fp) && identical(new_fp, velocity_state$object_fingerprint) &&
          !is.null(velocity_state$result)) {
        .velocity_log(.tr("Objet Seurat mis à jour (empreinte identique) — résultats velocity conservés."))
        return()
      }
      velocity_state$inputs <- NULL
      velocity_state$validated <- NULL
      velocity_state$result <- NULL
      velocity_state$alignment <- NULL
      velocity_state$umap_vectors <- NULL
      velocity_state$partial_embedding <- NULL
      velocity_state$embedding_alignment <- NULL
      velocity_state$vector_validation <- NULL
      velocity_state$plot_cache <- list()
      velocity_state$export_cache <- list()
      velocity_state$object_fingerprint <- NULL
      velocity_state$has_vectors <- FALSE
      velocity_state$optin_status <- "idle"
      velocity_state$optin_error <- NULL
      velocity_state$source_label <- NULL
      velocity_optin_banner_rv(NULL)
      # Lane 3b : l'état publié pour le snapshot est périmé avec l'objet.
      shared_rv$velocity_panel_state <- NULL

      reset_msg <- .tr("Objet Seurat modifie : reimportez et validez les donnees velocity.")
      velocity_state$status_msg <- reset_msg
      velocity_status_rv(reset_msg)

      updateSelectizeInput(
        session,
        "velocity_gene",
        choices = character(0),
        selected = character(0),
        server = TRUE
      )

      updateCheckboxInput(
        session,
        "velocity_strip_cell_suffix",
        value = FALSE
      )

      updateCheckboxInput(
        session,
        "velocity_strip_gene_version",
        value = FALSE
      )

      updateCheckboxInput(
        session,
        "velocity_allow_low_overlap",
        value = FALSE
      )

      updateSelectInput(
        session,
        "velocity_orientation",
        selected = "auto_strict"
      )
    }, ignoreInit = TRUE)

    # Lane 3b (2026-11-14) : réhydratation du panneau Vélocité depuis le
    # snapshot — app.R pose shared_rv$velocity_restore (one-shot) APRÈS les
    # purges du cycle de chargement (epoch mod_sc + sc_obj ci-dessus, cf.
    # session_results_pending). Les datapath des fileInput (inputs) ne
    # voyagent JAMAIS (temporaires) ; les caches sont dérivés, recréés.
    observeEvent(shared_rv$velocity_restore, {
      pay <- shared_rv$velocity_restore
      shared_rv$velocity_restore <- NULL  # one-shot
      if (!is.list(pay)) return()
      # Le résultat consommé (Communication/viz/rapport) est restauré même
      # sans état de panneau (snapshots lane 3 v1).
      shared_rv$velocity_result <- pay$velocity_result
      st <- pay$panel_state
      if (is.null(st) || !is.list(st)) return()
      velocity_state$validated           <- st$validated
      velocity_state$result              <- st$result
      velocity_state$alignment           <- st$alignment
      velocity_state$umap_vectors        <- st$umap_vectors
      velocity_state$partial_embedding   <- st$partial_embedding
      velocity_state$embedding_alignment <- st$embedding_alignment
      velocity_state$vector_validation   <- st$vector_validation
      # garde d'obsolescence : fingerprint recalculé SUR l'objet restauré.
      velocity_state$object_fingerprint <-
        if (!is.null(global_data$sc_obj)) velocity_object_fingerprint(global_data$sc_obj) else NULL
      velocity_state$inputs <- NULL
      velocity_state$plot_cache <- list()
      velocity_state$export_cache <- list()
      velocity_state$optin_error <- NULL
      velocity_state$has_vectors <- !is.null(st$umap_vectors)
      velocity_state$optin_status <- st$optin_status %||%
        if (identical(st$vector_validation$field, "scvelo_optin") &&
            isTRUE(st$vector_validation$ok)) "done" else "idle"
      velocity_state$source_label <- st$source_label %||% NULL
      # Retour opérateur « phase portrait isn't » : le sélecteur de gène doit
      # être peuplé SINON req(input$velocity_gene) laisse le phase portrait
      # vide après reload — même helper que la validation.
      .velocity_refresh_gene_selector()
      restore_msg <- st$status_msg %||% .tr("État velocity restauré depuis la session sauvegardée — ni ré-inféré ni recalculé.")
      velocity_state$status_msg <- restore_msg
      velocity_status_rv(restore_msg)
      .velocity_publish_state()
    }, ignoreNULL = TRUE)

    # Helper to check fingerprint before plot/export — etat derive
    # stale_against_current_seurat_object (contrat Stage 8).
    .check_fingerprint <- function() {
      req(global_data$sc_obj)
      req(velocity_state$object_fingerprint)
      current_fp <- velocity_object_fingerprint(global_data$sc_obj)
      if (!identical(velocity_state$object_fingerprint, current_fp)) {
        stale_msg <- .status_label("stale_against_current_seurat_object")
        showNotification(stale_msg, type = "error", duration = 8)
        shiny::validate(shiny::need(FALSE, stale_msg))
      }
      invisible(TRUE)
    }

    # Lane 3b (2026-11-14) : publication de l'état du PANNEAU Vélocité vers
    # l'état partagé — le snapshot de session l'embarque (le panneau retrouve
    # ses matrices validées au reload, pas seulement le résultat consommé).
    .velocity_publish_state <- function() {
      probe <- .velocity_drive_probe_fields(velocity_state)
      velocity_state$has_vectors <- probe$has_vectors
      velocity_state$optin_status <- probe$optin_status
      if (is.null(shared_rv)) return(invisible(NULL))
      shared_rv$velocity_panel_state <- list(
        result              = velocity_state$result,
        validated           = velocity_state$validated,
        alignment           = velocity_state$alignment,
        umap_vectors        = velocity_state$umap_vectors,
        partial_embedding   = velocity_state$partial_embedding,
        embedding_alignment = velocity_state$embedding_alignment,
        vector_validation   = velocity_state$vector_validation,
        object_fingerprint  = velocity_state$object_fingerprint,
        has_vectors         = probe$has_vectors,
        optin_status        = probe$optin_status,
        status_msg          = velocity_state$status_msg,
        source_label        = velocity_state$source_label
      )
      invisible(NULL)
    }

    # Lane 3b (retour opérateur « phase portrait isn't ») : alimentation du
    # sélecteur de gène du phase portrait — UN SEUL helper partagé par la
    # validation et la réhydratation (sans lui, req(input$velocity_gene)
    # laisse le phase portrait vide après un reload : le sélecteur n'était
    # peuplé que par la validation). 2D (O28) : heuristique marqueurs → HVG →
    # variance unspliced → fallback.
    .velocity_refresh_gene_selector <- function() {
      validated <- velocity_state$validated
      if (is.null(validated) || is.null(validated$spliced)) return(invisible(NULL))
      gene_choices <- rownames(validated$spliced)
      default_gene <- .velocity_phase_default_gene(
        gene_choices,
        markers_df   = shared_rv$markers_data,
        var_features = tryCatch(VariableFeatures(global_data$sc_obj),
                                error = function(e) character(0)),
        unspliced    = validated$unspliced
      )
      updateSelectizeInput(session, "velocity_gene", choices = gene_choices,
                        selected = default_gene %||% gene_choices[1], server = TRUE)
      invisible(NULL)
    }

    # ── O13 (2026-10-27) : inférence scVelo OPT-IN — voie unique, explicite ──
    # Bouton dédié (JAMAIS l'autopipeline) ; bannière dédiée sur le résultat ;
    # moteur absent/échec ⇒ notification d'erreur classée, aucun repli silencieux.
    velocity_optin_banner_rv <- reactiveVal(NULL)
    output$velocity_scvelo_banner <- renderUI({
      lab <- velocity_optin_banner_rv()
      if (is.null(lab)) return(NULL)
      tags$div(class = "alert alert-warning", style = "font-size:0.82rem; margin-top:6px;", lab)
    })

    # ── Lane async opt-in (2026-11-14, §M priorité 3) : l'inférence scVelo ───
    # part en daemon — MÊME pattern que LIANA P1-6b (mod_sc_communication.R,
    # validé en production) : ExtendedTask + pool unique (ensure_app_daemons,
    # spatial_async.R) + mirai (moteur sourcé idempotent, chemins absolus) +
    # timeout + garde de ré-entrée + repli synchrone TRACÉ + observateur
    # d'achèvement. Avant : 10-90 s bloquants dans le process Shiny.
    # Le moteur scVelo (sc_velocity_scvelo_optin_infer) tourne sous reticulate ;
    # les daemons héritent de l'env du process principal (RETICULATE_PYTHON est
    # posé par le launcher lane 2) et chargent sc_velocity.R idempotemment.
    # Triage suite de fond 2026-11-14 : `.sc_project_root()` (session_checkpoint.R)
    # n'est PAS sourcé par tous les fichiers qui montent ce module (testServer
    # optin-handoff mesuré : « could not find function » à l'INIT du module).
    # Défensive, PAS une re-dérivation : en app, .sc_project_root existe
    # (session_checkpoint.R sourcé avant les modules) et le fallback getwd()
    # est la racine du projet (l'app y est lancée).
    .ts_sc_daemon_root <- function() {
      if (exists(".sc_project_root", where = globalenv(), inherits = FALSE)) {
        return(get(".sc_project_root", envir = globalenv())())
      }
      getwd()
    }
    ext_task_cls_vel <- if (exists("ExtendedTask", where = asNamespace("shiny"), inherits = FALSE)) {
      get("ExtendedTask", envir = asNamespace("shiny"))
    } else if (exists("ExtendedTask", where = asNamespace("bslib"), inherits = FALSE)) {
      get("ExtendedTask", envir = asNamespace("bslib"))
    } else {
      get0("ExtendedTask", envir = globalenv())
    }
    scvelo_core_deps <- c(
      file.path(.ts_sc_daemon_root(), "tools", "pkg_skeleton", "R", "sc_velocity.R"),
      file.path(.ts_sc_daemon_root(), "tools", "pkg_skeleton", "R", "error_state.R")
    )
    scvelo_task <- if (!is.null(ext_task_cls_vel) && requireNamespace("mirai", quietly = TRUE)) {
      t <- ext_task_cls_vel$new(function(canon, emb, params, deps = scvelo_core_deps) {
        ensure_app_daemons()
        mirai::mirai(
          {
            for (f in deps) if (file.exists(f)) source(f, local = FALSE, encoding = "UTF-8")
            sc_velocity_scvelo_optin_infer(canon, embedding = emb, params = params)
          },
          canon = canon, emb = emb, params = params, deps = scvelo_core_deps,
          .timeout = TS_SC_AUTO_PIPELINE_TIMEOUT_S * 1000L
        )
      })
      tryCatch(bslib::bind_task_button(t, "velocity_scvelo_run"), error = function(e) NULL)
      t
    } else {
      NULL
    }
    # Contexte d'invocation (réduction/canonique AU MOMENT du clic) —
    # l'achèvement async applique le résultat avec les MÊMES entrées, même si
    # l'opérateur a changé le dropdown entre-temps.
    scvelo_ctx <- reactiveValues(canon = NULL, sel_red = NULL)

    .vel_sanitize_text <- function(x) {
      if (!is.character(x) || length(x) == 0L) return(x)
      out <- if (exists(".ts_nul_free_text", mode = "function")) {
        .ts_nul_free_text(x)
      } else {
        vapply(x, function(s) {
          if (is.null(s) || !length(s)) return("")
          res <- tryCatch({
            r <- charToRaw(s)
            if (any(r == as.raw(0L))) rawToChar(r[r != as.raw(0L)]) else s
          }, error = function(e) "")
          enc2utf8(res)
        }, character(1L), USE.NAMES = FALSE)
      }
      vapply(out, function(s) {
        if (is.na(s) || !nzchar(s)) return(if (is.na(s)) "" else s)
        if (!validUTF8(s)) {
          s <- iconv(s, from = "UTF-8", to = "UTF-8", sub = "byte")
          if (is.na(s)) s <- ""
        }
        enc2utf8(s)
      }, character(1L), USE.NAMES = FALSE)
    }

    .vel_sanitize_payload <- function(obj) {
      if (is.character(obj)) {
        res <- .vel_sanitize_text(obj)
        names(res) <- if (!is.null(names(obj))) .vel_sanitize_text(names(obj)) else NULL
        return(res)
      }
      if (is.list(obj)) {
        for (i in seq_along(obj)) {
          obj[[i]] <- .vel_sanitize_payload(obj[[i]])
        }
        if (!is.null(names(obj))) names(obj) <- .vel_sanitize_text(names(obj))
      }
      obj
    }

    .vel_scvelo_fail <- function(e) {
      # P0.4 lane 1 (2026-11-13) : l'échec opt-in laisse un état PERSISTANT
      # (velocity_state$optin_error + statut) — même surface pour le repli
      # synchrone et l'achèvement async.
      msg_clean <- .vel_sanitize_text(conditionMessage(e))
      velocity_state$optin_error <- list(
        class = class(e)[[1L]],
        message = msg_clean,
        when = format(Sys.time(), "%H:%M:%S")
      )
      velocity_state$optin_status <- "error"
      message(sprintf("[velocity] inférence scVelo opt-in refusée [%s] : %s",
                      velocity_state$optin_error$class,
                      velocity_state$optin_error$message))
      err_status <- paste0(.tr("Inférence scVelo (opt-in) refusée :"), " ", msg_clean)
      velocity_state$status_msg <- err_status
      velocity_status_rv(err_status)
      .velocity_publish_state()
      showNotification(paste(.tr("Inférence scVelo (opt-in) refusée :"), msg_clean),
                       type = "error", duration = 12)
    }

    .apply_scvelo_result <- function(inferred, sel_red, canon) {
      uv <- validate_precomputed_velocity_vectors(
        vectors = inferred$embedding_vectors,
        cell_names = velocity_state$validated$cell_names,
        selected_reduction = sel_red,
        stored_reduction = sel_red
      )
      velocity_state$umap_vectors <- uv
      velocity_state$vector_validation <- list(field = "scvelo_optin", ok = TRUE)
      # FIX 2026-11-11 (démo velocity, chemin loom+opt-in) : un loom n'a pas de
      # réduction stockée (embedding_reduction = NA_character_) et le garde de
      # finalize_velocity_result() refusait 'NA' != sel_red — alors que les
      # vecteurs viennent d'être calculés dans CETTE session POUR sel_red
      # (embedding passé au moteur = Embeddings(sc_obj, sel_red)). Taguer le
      # canonique en conséquence est la sémantique correcte des vecteurs
      # opt-in (jamais précalculés, toujours pour la réduction courante).
      velocity_state$validated$embedding_reduction <- sel_red
      canonical2 <- finalize_velocity_result(
        validated = velocity_state$validated,
        input_mode = canon$input_summary$input_mode %||% "rds",
        input_files = canon$input_summary$files,
        seurat_obj = global_data$sc_obj,
        assay_used = NULL,
        requested_reduction = sel_red,
        embedding_alignment = velocity_state$embedding_alignment,
        velocity_vectors = uv,
        vector_validation = velocity_state$vector_validation,
        extra_warnings = c("O13 opt-in : vecteurs inférés dans l'app (scVelo) — jamais l'autopipeline."),
        analysis_id = "sc-velocity"
      )
      velocity_state$result <- canonical2
      shared_rv$velocity_result <- canonical2
      provenance_append(shared_rv, inferred$provenance)
      velocity_optin_banner_rv(.tr(inferred$banner))
      done_status <- .tr("Vecteurs inférés (scVelo opt-in) — affichage mis à jour.")
      velocity_state$status_msg <- done_status
      velocity_status_rv(done_status)
      # P0.4 lane 1 (2026-11-13) : succès => on efface toute trace d'un
      # échec opt-in antérieur (l'état honnête reflète la DERNIÈRE tentative).
      velocity_state$optin_error <- NULL
      velocity_state$optin_status <- "done"
      .velocity_publish_state()
      showNotification(.tr("Vecteurs calculés via scVelo (opt-in)."), type = "message", duration = 6)
    }

    drive_counter_validate <- shiny::reactiveVal(0L)
    drive_counter_scvelo   <- shiny::reactiveVal(0L)
    drive_trigger_validate <- shiny::reactive({
      dc <- drive_counter_validate()
      bc <- input$velocity_validate
      if (identical(dc, 0L) && (is.null(bc) || identical(as.integer(bc), 0L))) return(NULL)
      list(dc, bc)
    })
    drive_trigger_scvelo   <- shiny::reactive({
      dc <- drive_counter_scvelo()
      bc <- input$velocity_scvelo_run
      if (identical(dc, 0L) && (is.null(bc) || identical(as.integer(bc), 0L))) return(NULL)
      list(dc, bc)
    })

    observeEvent(drive_trigger_scvelo(), {
      canon <- velocity_state$result
      if (is.null(canon)) {
        showNotification(
          .tr("Validez d'abord les matrices velocity (aucun résultat canonique à inférer)."),
          type = "warning", duration = 8)
        return()
      }
      # garde de ré-entrée (P1-6b) : un calcul déjà en vol ne relance pas
      if (!is.null(scvelo_task) && identical(scvelo_task$status(), "running")) return(invisible())
      sel_red <- input$velocity_reduction %||% "umap"
      if (!.velocity_has_reduction(global_data$sc_obj, sel_red)) {
        .vel_scvelo_fail(errorCondition(
          .velocity_missing_reduction_msg(sel_red),
          class = c("sc_velocity_error", "error", "condition")
        ))
        return()
      }
      emb_all <- Seurat::Embeddings(global_data$sc_obj, reduction = sel_red)
      emb <- emb_all[canon$cell_names, seq_len(min(2L, ncol(emb_all))), drop = FALSE]
      params <- list(seed = 1L, n_pcs = 30L, n_neighbors = 30L, mode = "stochastic")
      pool_ok <- isTRUE(tryCatch(ensure_app_daemons(), error = function(e) FALSE))
      if (!is.null(scvelo_task) && pool_ok) {
        # ASYNC : le calcul Python part en daemon — la session reste réactive.
        scvelo_ctx$canon <- canon
        scvelo_ctx$sel_red <- sel_red
        velocity_state$optin_error <- NULL
        velocity_state$optin_status <- "running"
        run_status <- .tr("Inférence scVelo (opt-in) en daemon — la session reste réactive, comptez quelques dizaines de secondes…")
        velocity_state$status_msg <- run_status
        velocity_status_rv(run_status)
        .velocity_publish_state()
        showNotification(run_status, type = "message", duration = 6)
        scvelo_task$invoke(canon, emb, params)
      } else {
        # Repli synchrone TRACÉ (jamais muet — convention run_job) ; la
        # notification persistante lane 1 reste retirée par on.exit.
        velocity_state$optin_error <- NULL
        velocity_state$optin_status <- "running"
        .velocity_publish_state()
        showNotification(.tr("Worker mirai indisponible — exécution synchrone : l'interface restera bloquée pendant le calcul."),
                         type = "warning", duration = 8)
        nid_optin <- showNotification(
          .tr("Inférence scVelo en cours (opt-in) — calcul Python en tâche de fond du moteur, comptez quelques dizaines de secondes…"),
          type = "message", duration = NULL, closeButton = FALSE)
        on.exit(removeNotification(nid_optin), add = TRUE)
        tryCatch({
          inferred <- sc_velocity_scvelo_optin_infer(canon, embedding = emb, params = params)
          .apply_scvelo_result(inferred, sel_red, canon)
        }, error = function(e) {
          .vel_scvelo_fail(e)
        })
      }
    })

    # Achèvement async : application du résultat / échec persistant — mêmes
    # surfaces que le repli synchrone (statut, notification classée, optin_error).
    if (!is.null(scvelo_task)) {
      observeEvent(scvelo_task$status(), {
        st <- scvelo_task$status()
        if (identical(st, "running")) {
          velocity_state$optin_status <- "running"
          .velocity_publish_state()
        } else if (identical(st, "success")) {
          tryCatch({
            res_clean <- .vel_sanitize_payload(scvelo_task$result())
            .apply_scvelo_result(res_clean, scvelo_ctx$sel_red, scvelo_ctx$canon)
          }, error = function(e) {
            .vel_scvelo_fail(e)
          })
        } else if (identical(st, "error")) {
          err <- tryCatch(scvelo_task$result(), error = function(e) e)
          err_raw <- if (inherits(err, "condition")) {
            conditionMessage(err)
          } else {
            paste(as.character(unlist(.vel_sanitize_payload(err), use.names = FALSE)), collapse = " ")
          }
          err_msg <- enc2utf8(.vel_sanitize_text(paste(as.character(err_raw), collapse = " ")))
          if (length(err_msg) == 0L || is.na(err_msg[[1L]])) err_msg <- ""
          err_cls <- if (inherits(err, "error")) class(err) else c("error", "condition")
          .vel_scvelo_fail(errorCondition(err_msg[[1L]], class = err_cls))
        }
      })
    }

    # ── §H.6 (2026-11-14) : pipeline d'ingestion PARTAGÉ — zéro voie parallèle ──
    # Le corps complet du validate (lecture → validation → alignement →
    # canonique → publication) est factorisé ici ; le bouton humain et
    # l'importeur drive (import_path, module `import_velocity`) l'appellent
    # tous les deux, chacun avec SON read_fn : le bouton lit les datapaths du
    # fileInput, l'importeur lit le chemin CONFINÉ résolu par
    # `ts_resolve_data_path()` (mêmes lecteurs `read_velocity_*`, même chaîne
    # `validate_velocity_matrices()` — on contourne le widget, on ne le
    # falsifie pas, spec S5). Lève en cas d'échec ; l'appelant dispose du
    # tryCatch qui fixe le statut et notifie.
    .velocity_ingest <- function(mode, read_fn, source_files = NULL) {
      seurat_cells <- colnames(global_data$sc_obj)
      seurat_genes <- rownames(global_data$sc_obj)

      read <- read_fn()
      spliced <- read$spliced
      unspliced <- read$unspliced
      ambiguous <- read$ambiguous %||% NULL
      rds_meta <- read$meta %||% list()
      vector_field <- read$vector_field
      vector_matrix <- read$vector_matrix

      # 3. Read orientation, 4. suffix/version, 5. low overlap confirmation (BUG 5/6)
      orientation <- input$velocity_orientation %||% "auto_strict"
      strip_cell_suffix <- isTRUE(input$velocity_strip_cell_suffix)
      strip_gene_version <- isTRUE(input$velocity_strip_gene_version)
      allow_low_overlap <- isTRUE(input$velocity_allow_low_overlap)

      # 6. Validate matrices, 7. Align by barcode and gene ID.
      # Normalization applies to BOTH sides (velocity + Seurat barcodes).
      validated <- validate_velocity_matrices(
        spliced = spliced,
        unspliced = unspliced,
        ambiguous = ambiguous,
        seurat_cells = seurat_cells,
        seurat_genes = seurat_genes,
        orientation = orientation,
        strip_cell_suffix = strip_cell_suffix,
        strip_gene_version = strip_gene_version,
        allow_low_overlap = allow_low_overlap
      )

      # Zero aligned cells makes every downstream view (phase portrait,
      # embedding, export) empty without any visible explanation.
      if (validated$n_cells_matched == 0L) {
        stop(errorCondition(
          paste0(
            "Aucune cellule alignee entre les matrices velocity et ",
            "l'objet Seurat courant. Verifiez que l'objet charge ",
            "correspond au meme run Cell Ranger et activez le retrait ",
            "des suffixes (-1) si necessaire."
          ),
          class = "sc_velocity_error"
        ))
      }

      # BUG 10: preserve RDS metadata on the validated result (never in
      # Seurat meta.data — this stays module-local).
      validated$velocity_source <- rds_meta$velocity_source %||% "unknown"
      validated$velocity_method <- rds_meta$velocity_method %||% "precomputed"
      validated$input_orientation <- rds_meta$orientation %||% NA_character_
      validated$embedding_reduction <- rds_meta$embedding_reduction %||%
        rds_meta$stored_reduction %||%
        rds_meta$reduction %||% NA_character_
      validated$clusters <- rds_meta$clusters %||% NULL
      validated$umap_embedding <- rds_meta$umap_embedding %||% NULL

      # 8. Store validated result (structure plate, compatibilite)
      velocity_state$validated <- validated
      velocity_state$alignment <- list(cells = validated$cell_names, genes = validated$gene_names)

      # Reset per-import derived state before repopulating
      velocity_state$partial_embedding <- NULL
      velocity_state$embedding_alignment <- NULL
      velocity_state$vector_validation <- NULL
      velocity_state$umap_vectors <- NULL
      velocity_state$has_vectors <- FALSE
      velocity_state$optin_status <- "idle"
      velocity_state$optin_error <- NULL
      velocity_optin_banner_rv(NULL)

      sel_red <- input$velocity_reduction %||% "umap"
      orchestration_warnings <- character(0)
      if (is.null(validated$umap_embedding) &&
          !.velocity_has_reduction(global_data$sc_obj, sel_red)) {
        msg_red <- .velocity_missing_reduction_msg(sel_red)
        orchestration_warnings <- c(orchestration_warnings, msg_red)
        showNotification(msg_red, type = "warning", duration = 8)
      }

      # BUG 3/BUG 4: partial coordinate display via align_velocity_embedding()
      if (!is.null(validated$umap_embedding)) {
        ali <- tryCatch(
          align_velocity_embedding(
            embedding = validated$umap_embedding,
            cell_names = validated$cell_names,
            selected_reduction = sel_red,
            stored_reduction = validated$embedding_reduction
          ),
          error = function(e) {
            orchestration_warnings <<- c(orchestration_warnings,
              paste0("Embedding UMAP ignore : ", conditionMessage(e)))
            showNotification(paste("Embedding UMAP ignore :", conditionMessage(e)),
                             type = "warning", duration = 8)
            NULL
          }
        )
        if (!is.null(ali)) {
          velocity_state$partial_embedding <- ali$embedding
          velocity_state$embedding_alignment <- ali
          validated$embedding_alignment <- ali
        }
      }

      # BUG 8: strict vector validation — no arrows on any failure
      if (!is.null(vector_matrix)) {
        vec_stored_red <- validated$embedding_reduction
        if (!identical(vector_field, "umap_velocity") &&
            !identical(vector_field, "embedding_velocity") &&
            !identical(vector_field, "umap_vectors") &&
            !identical(vector_field, "vectors")) {
          vector_matrix <- NULL
        }
        if (!is.null(vector_matrix)) {
          tryCatch({
            velocity_state$umap_vectors <- validate_precomputed_velocity_vectors(
              vectors = vector_matrix,
              cell_names = validated$cell_names,
              selected_reduction = sel_red,
              stored_reduction = vec_stored_red
            )
            velocity_state$vector_validation <- list(field = vector_field, ok = TRUE)
          }, error = function(e) {
            velocity_state$umap_vectors <- NULL
            velocity_state$vector_validation <- list(field = vector_field, ok = FALSE,
                                                     reason = conditionMessage(e))
            orchestration_warnings <<- c(orchestration_warnings,
              paste0("Vecteurs velocity rejetes : ", conditionMessage(e)))
            showNotification(paste("Vecteurs velocity refuses :", conditionMessage(e)),
                             type = "error", duration = 8)
          })
        }
      }

      # 9. Store current object fingerprint
      velocity_state$object_fingerprint <- velocity_object_fingerprint(global_data$sc_obj)

      # Stage 8 (3B-1) : resultat canonique — contrat documente dans
      # R/sc/sc_velocity.R. La provenance est PRODUITE dans l'objet (regle 7
      # AGENTS.md) puis appendee a l'etat partage ; consolidation reservee au
      # rapport (macro-step 4F).
      resolved_files <- source_files %||% .collect_input_files(mode)
      file_label <- if (length(resolved_files)) {
        first_val <- resolved_files[[mode]] %||% resolved_files[[1L]]
        if (is.character(first_val) && length(first_val) == 1L && nzchar(first_val)) first_val else mode
      } else {
        mode
      }
      src_label <- sprintf("%s (%s)", file_label, mode)

      canonical <- finalize_velocity_result(
        validated           = validated,
        input_mode          = mode,
        input_files         = resolved_files,
        seurat_obj          = global_data$sc_obj,
        assay_used          = NULL,
        requested_reduction = sel_red,
        embedding_alignment = velocity_state$embedding_alignment,
        velocity_vectors    = velocity_state$umap_vectors,
        vector_validation   = velocity_state$vector_validation,
        extra_warnings      = orchestration_warnings,
        analysis_id         = "sc-velocity"
      )
      velocity_state$result <- canonical

      # Stage 17 (4F) : exposition ADDITIVE du resultat canonique au rapport
      # consolide (lecture seule ; velocity_state reste la reference du
      # panel 8 — aucun changement de comportement du panneau).
      shared_rv$velocity_result <- canonical

      provenance_append(shared_rv, canonical$provenance)

      # 10. Populate gene selector — 2D (O28) — helper partagé avec la
      # réhydratation lane 3b (zéro voie parallèle).
      .velocity_refresh_gene_selector()

      # 11. French status summary : statut canonique explicite (+ BUG 3
      # message when only coordinates exist). Validite technique uniquement —
      # aucune validite biologique n'est impliquee.
      msg <- sprintf(
        paste0(
          "Validation OK : %d cellules, %d genes. ",
          "Orientation : %s. Alignement cellules : %s. ",
          "Alignement genes : %s. Source : %s."
        ),
        validated$n_cells_matched,
        validated$n_genes_matched,
        validated$orientation,
        validated$cell_match_mode,
        validated$gene_match_mode,
        src_label
      )
      if (is.null(vector_matrix) && !is.null(velocity_state$partial_embedding)) {
        msg <- paste(
          msg,
          "Coordonnees UMAP disponibles, mais aucun vecteur velocity pre-calcule n'est fourni.",
          sep = "\n"
        )
      }
      msg <- paste0("[", .status_label(canonical$status), "]\n", msg)
      velocity_state$source_label <- src_label
      velocity_state$status_msg <- msg
      velocity_status_rv(msg)
      # Lane 3b : publication panneau (validate succès) — après mise à jour de
      # source_label / status_msg / has_vectors / optin_status.
      .velocity_publish_state()
      showNotification(.tr("Validation velocity reussie."), type = "message", duration = 4)
      invisible(NULL)
    }

    # ── §H.6 (2026-11-14) : importeur drive `import_path` (module
    # `import_velocity`) — la voie fileInput loom/h5ad/rds devient drivable
    # sans navigateur. Le chemin arrive CONFINÉ et NORMALISÉ — la validation
    # (confinement + allowlist d'extensions `exts`) est faite au SEAM par
    # `ts_drive_validate_velocity_path()` (motif import_sc : l'importeur
    # consomme des données validées, il ne les re-dérive pas). L'ingestion
    # passe par `.velocity_ingest()` — le pipeline EXACT du bouton humain
    # ci-dessous ; les options d'alignement (orientation, suffixes,
    # recouvrement) sont lues sur les widgets par `isolate()` : elles se
    # règlent par un `set_inputs` ANTÉRIEUR, seul ordre correct (motif
    # import_bulk, G3). Sonde d'état bornée (motif import_sc `drive_state`).
    drive_velocity_state <- function() {
      res <- tryCatch(shiny::isolate(velocity_state$result), error = function(e) NULL)
      obj <- tryCatch(shiny::isolate(global_data$sc_obj), error = function(e) NULL)
      probe <- tryCatch(
        shiny::isolate(.velocity_drive_probe_fields(velocity_state)),
        error = function(e) .velocity_drive_probe_fields(NULL)
      )
      status_msg <- tryCatch(shiny::isolate(velocity_state$status_msg), error = function(e) NULL)
      source_label <- tryCatch(shiny::isolate(velocity_state$source_label), error = function(e) NULL)
      list(module = TS_DRIVE_SC_VELOCITY_IMPORT_MODULE,
           has_result = !is.null(res),
           status = if (is.null(res)) NULL else res$status,
           n_cells_matched = if (is.null(res)) NULL else res$n_cells_matched,
           has_seurat = !is.null(obj),
           has_vectors = probe$has_vectors,
           optin_status = probe$optin_status,
           status_msg = status_msg,
           source_label = source_label)
    }
    # Garde défensive (motif ExtendedTask ci-dessus) : la couche drive
    # (drive_allowlist.R) est sourcée par global.R avant les modules — mais
    # des environnements de test minimalistes montent ce module SANS elle
    # (test-sc-velocity-optin-handoff mesuré : « could not find function
    # ts_drive_publish_importer » à l'INIT). Sans drive, pas d'importeur :
    # le panneau humain est inchangé.
    .vel_drive_ok <- exists("ts_drive_publish_importer", where = globalenv(),
                            inherits = FALSE) &&
      exists("ts_drive_registry", where = globalenv(), inherits = FALSE) &&
      exists("TS_DRIVE_SC_VELOCITY_IMPORT_MODULE", where = globalenv(),
             inherits = FALSE)
    if (.vel_drive_ok) ts_drive_publish_importer(global_data,
                              TS_DRIVE_SC_VELOCITY_IMPORT_MODULE,
                              function(request) {
      warnings <- character(0)
      path <- request$path
      ext <- tolower(tools::file_ext(path))
      mode <- if (identical(ext, "rds")) "rds" else "loom"
      # Le nom AFFICHÉ — jamais le chemin local complet (règle
      # .collect_input_files : noms originaux uniquement) ; `name` du
      # scénario est honoré (déjà trimmé par le validator), sinon basename —
      # la MÊME règle que le widget.
      display_name <- request$name %||% basename(path)
      source_files <- setNames(list(display_name), mode)
      read_fn <- function() {
        velocity_input <- if (identical(ext, "h5ad")) {
          read_velocity_h5ad(path)
        } else if (identical(ext, "loom")) {
          read_velocity_loom(path)
        } else {
          read_velocity_rds(path)
        }
        vf <- NULL; vm <- NULL
        for (f in c("umap_velocity", "embedding_velocity", "umap_vectors", "vectors")) {
          if (!is.null(velocity_input[[f]])) { vf <- f; vm <- velocity_input[[f]]; break }
        }
        list(spliced = velocity_input$spliced,
             unspliced = velocity_input$unspliced,
             ambiguous = velocity_input$ambiguous,
             meta = velocity_input,
             vector_field = vf, vector_matrix = vm)
      }
      tryCatch({
        .velocity_ingest(mode, read_fn, source_files = source_files)
        probe <- .velocity_drive_probe_fields(velocity_state)
        list(ok = TRUE, status = "applied", errors = character(0),
             warnings = warnings,
             has_vectors = probe$has_vectors,
             optin_status = probe$optin_status,
             source_label = velocity_state$source_label,
             status_msg = velocity_state$status_msg)
      }, error = function(e) {
        err_msg <- conditionMessage(e)
        err_status <- paste(.tr("Erreur validation velocity :"), err_msg)
        velocity_state$status_msg <- err_status
        velocity_status_rv(err_status)
        showNotification(err_status, type = "error", duration = 8)
        probe <- .velocity_drive_probe_fields(velocity_state)
        list(ok = FALSE, status = "error",
             errors = err_msg, warnings = warnings,
             has_vectors = probe$has_vectors,
             optin_status = probe$optin_status)
      })
    }, state = drive_velocity_state)

    drive_ready_validate <- function() {
      obj <- tryCatch(shiny::isolate(global_data$sc_obj), error = function(e) NULL)
      if (is.null(obj)) {
        return("no Seurat object loaded (global_data$sc_obj is NULL)")
      }
      mode <- tryCatch(shiny::isolate(input$velocity_import_mode), error = function(e) NULL) %||% "rds"
      f <- if (identical(mode, "loom")) {
        tryCatch(shiny::isolate(input$velocity_loom_file), error = function(e) NULL)
      } else {
        tryCatch(shiny::isolate(input$velocity_rds_file), error = function(e) NULL)
      }
      if (is.null(f)) {
        return(sprintf("no velocity file selected for mode '%s'", mode))
      }
      TRUE
    }
    drive_ready_scvelo <- function() {
      canon <- tryCatch(shiny::isolate(velocity_state$result), error = function(e) NULL)
      if (is.null(canon)) {
        return("no validated velocity result yet (run velocity_validate or import_velocity first)")
      }
      if (!is.null(scvelo_task) &&
          identical(tryCatch(scvelo_task$status(), error = function(e) ""), "running")) {
        return("scVelo inference is already running")
      }
      obj <- tryCatch(shiny::isolate(global_data$sc_obj), error = function(e) NULL)
      sel_red <- tryCatch(shiny::isolate(input$velocity_reduction), error = function(e) NULL) %||% "umap"
      if (!.velocity_has_reduction(obj, sel_red)) {
        return(sprintf("reduction '%s' is not available on the Seurat object", sel_red))
      }
      TRUE
    }
    if (exists("ts_drive_bind_button", mode = "function") &&
        exists("ts_drive_publish_token", mode = "function") &&
        exists("TS_DRIVE_BUTTONS")) {
      ts_drive_bind_button("sc-velocity-velocity_validate")
      ts_drive_bind_button("sc-velocity-velocity_scvelo_run")
      ts_drive_publish_token(global_data, "sc-velocity-velocity_validate",
                             drive_counter_validate, ready = drive_ready_validate,
                             state = drive_velocity_state)
      ts_drive_publish_token(global_data, "sc-velocity-velocity_scvelo_run",
                             drive_counter_scvelo, ready = drive_ready_scvelo,
                             state = drive_velocity_state)
    }

    # Validation button (BUG 1-8, 10)
    observeEvent(drive_trigger_validate(), {
      req(global_data$sc_obj)

      # P0.4 lane 1 (2026-11-13, demande opérateur) : même traitement que
      # l'opt-in — la validation lit les couches sur de grosses matrices
      # (quelques secondes sans aucun retour). Notification persistante
      # retirée par on.exit (succès COMME échec).
      nid_validate <- showNotification(
        .tr("Validation des matrices velocity en cours — lecture spliced/unspliced et alignement…"),
        type = "message", duration = NULL, closeButton = FALSE)
      on.exit(removeNotification(nid_validate), add = TRUE)

      tryCatch({
        # CHRYSALIS 2E : garde d'entree explicite — le req() ci-dessus laisse
        # l'observateur silencieux si l'objet est absent/non-Seurat ; l'assert
        # produit desormais un message FR clair dans le canal d'erreur.
        assert_seurat(global_data$sc_obj, context = "velocity")

        # §H.6 (2026-11-14) : le pipeline complet est PARTAGÉ
        # (.velocity_ingest ci-dessus, aussi consommé par l'importeur drive) ;
        # ce dispatch ne décide que COMMENT LIRE — depuis les datapaths du
        # fileInput. L'upload anonymise l'extension : loom/h5ad se dispatche
        # sur le nom original.
        mode <- input$velocity_import_mode %||% "rds"
        if (identical(mode, "rds")) {
          read_fn <- function() {
            req(input$velocity_rds_file)
            # read_velocity_rds() enforces the whitelist AND validates metadata
            velocity_input <- read_velocity_rds(input$velocity_rds_file$datapath)
            vf <- NULL; vm <- NULL
            # BUG 3/BUG 8: ONLY these fields are dx/dy vectors in embedding
            # space. umap_embedding is coordinates and is NEVER routed here.
            for (f in c("umap_velocity", "embedding_velocity", "umap_vectors", "vectors")) {
              if (!is.null(velocity_input[[f]])) { vf <- f; vm <- velocity_input[[f]]; break }
            }
            list(spliced = velocity_input$spliced,
                 unspliced = velocity_input$unspliced,
                 ambiguous = velocity_input$ambiguous,
                 meta = velocity_input,
                 vector_field = vf, vector_matrix = vm)
          }
        } else if (identical(mode, "loom")) {
          # Roadmap SC FUNCTION_TEST M0 : lecture .loom (velocyto) / .h5ad
          # (scVelo, layout old-anndata). Meme chaine de validation que RDS :
          # validate_velocity_matrices() s'applique tel quel dans
          # .velocity_ingest.
          read_fn <- function() {
            req(input$velocity_loom_file)
            lname <- input$velocity_loom_file$name %||% ""
            velocity_input <- if (grepl("\\.h5ad$", lname, ignore.case = TRUE)) {
              read_velocity_h5ad(input$velocity_loom_file$datapath)
            } else {
              read_velocity_loom(input$velocity_loom_file$datapath)
            }
            vf <- NULL; vm <- NULL
            for (f in c("umap_velocity", "embedding_velocity", "umap_vectors", "vectors")) {
              if (!is.null(velocity_input[[f]])) { vf <- f; vm <- velocity_input[[f]]; break }
            }
            list(spliced = velocity_input$spliced,
                 unspliced = velocity_input$unspliced,
                 ambiguous = velocity_input$ambiguous,
                 meta = velocity_input,
                 vector_field = vf, vector_matrix = vm)
          }
        } else {
          # MTX mode: require both spliced and unspliced triplets
          read_fn <- function() {
            req(input$velocity_mtx_spliced, input$velocity_mtx_spliced_barcodes, input$velocity_mtx_spliced_features,
                input$velocity_mtx_unspliced, input$velocity_mtx_unspliced_barcodes, input$velocity_mtx_unspliced_features)
            fc <- input$velocity_feature_column %||% 1L
            sp <- read_velocity_mtx(
              matrix_path = input$velocity_mtx_spliced$datapath,
              barcode_path = input$velocity_mtx_spliced_barcodes$datapath,
              feature_path = input$velocity_mtx_spliced_features$datapath,
              feature_column = fc
            )
            uns <- read_velocity_mtx(
              matrix_path = input$velocity_mtx_unspliced$datapath,
              barcode_path = input$velocity_mtx_unspliced_barcodes$datapath,
              feature_path = input$velocity_mtx_unspliced_features$datapath,
              feature_column = fc
            )
            # optional ambiguous triplet (only when all three files provided)
            amb <- NULL
            if (!is.null(input$velocity_mtx_ambiguous) &&
                !is.null(input$velocity_mtx_ambiguous_barcodes) &&
                !is.null(input$velocity_mtx_ambiguous_features)) {
              amb <- read_velocity_mtx(
                matrix_path = input$velocity_mtx_ambiguous$datapath,
                barcode_path = input$velocity_mtx_ambiguous_barcodes$datapath,
                feature_path = input$velocity_mtx_ambiguous_features$datapath,
                feature_column = fc
              )
            }
            list(spliced = sp, unspliced = uns, ambiguous = amb,
                 meta = list(), vector_field = NULL, vector_matrix = NULL)
          }
        }
        .velocity_ingest(mode, read_fn)


      }, error = function(e) {
        err_status <- paste(.tr("Erreur validation velocity :"), conditionMessage(e))
        velocity_state$status_msg <- err_status
        velocity_status_rv(err_status)
        showNotification(
          err_status,
          type = "error",
          duration = 8
        )
      })
    })

    # Phase portrait (Stage 9 : consommateur pur du resultat canonique ;
    # apercu avec cap declare et affiche)
    output$velocity_phase_plot <- renderPlot({
      req(velocity_state$result)
      req(input$velocity_gene)
      .check_fingerprint()
      if (ncol(velocity_state$result$spliced) == 0L) {
        return(ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
            label = "Aucune cellule alignee : revalidez les donnees velocity.",
            size = 5, colour = "red") +
          ggplot2::theme_void())
      }
      tryCatch({
        plot_velocity_phase_portrait(
          velocity_state$result, input$velocity_gene,
          data_scope = "preview", seurat_obj = global_data$sc_obj
        )
      }, error = function(e) .error_plot(e))
    })

    # Embedding plot (Stage 9) : les fleches proviennent UNIQUEMENT des
    # vecteurs valides embarques dans le resultat canonique ; coordonnees
    # RDS partielles prioritaires, sinon reduction Seurat. Reduction absente
    # => message de validation exploitable, jamais de substitution.
    output$velocity_embedding_plot <- renderPlot({
      req(velocity_state$result)
      .check_fingerprint()
      red <- input$velocity_reduction %||% "umap"

      base_coords <- .base_coords_for(red)
      if (is.null(base_coords)) {
        msg <- if (!.velocity_has_reduction(global_data$sc_obj, red)) {
          .velocity_missing_reduction_msg(red)
        } else {
          sprintf(paste0("Reduction '%s' introuvable sur l'objet Seurat ou ",
                         "aucune cellule velocity dans l'embedding : calculez ",
                         "la reduction ou choisissez-en une autre."), red)
        }
        shiny::validate(shiny::need(FALSE, msg))
      }

      tryCatch(
        plot_velocity_embedding(
          velocity_state$result, base_coords,
          selected_reduction = red, data_scope = "preview",
          seurat_obj = global_data$sc_obj
        ),
        error = function(e) .error_plot(e)
      )
    })

    # Vue densite/couverture des vecteurs valides (Stage 9)
    output$velocity_coverage_plot <- renderPlot({
      req(velocity_state$result)
      .check_fingerprint()
      tryCatch(
        plot_velocity_coverage(velocity_state$result,
                               seurat_obj = global_data$sc_obj),
        error = function(e) .error_plot(e)
      )
    })

    # Vue QC d'alignement (Stage 9)
    output$velocity_qc_plot <- renderPlot({
      req(velocity_state$result)
      .check_fingerprint()
      tryCatch(
        plot_velocity_alignment_qc(velocity_state$result,
                                   seurat_obj = global_data$sc_obj),
        error = function(e) .error_plot(e)
      )
    })

    # Provenance export
    output$dl_velocity_provenance <- downloadHandler(
      filename = function() paste0("velocity_provenance_", format(Sys.Date(), "%Y-%m-%d"), ".csv"),
      content = function(file) {
        req(velocity_state$result)
        req(input$velocity_gene)
        .check_fingerprint()
        df <- build_velocity_provenance_export(
          velocity_result = velocity_state$result,
          gene = input$velocity_gene,
          selected_reduction = input$velocity_reduction %||% NA_character_
        )
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    # ── Exports contrat Stage 8 ─────────────────────────────────────────────
    # Resume de validation (1 ligne : analysis_id, statut, reglages, compteurs,
    # avertissements, versions).
    output$dl_velocity_summary <- downloadHandler(
      filename = function() paste0("velocity_validation_summary_", format(Sys.Date(), "%Y-%m-%d"), ".csv"),
      content = function(file) {
        req(velocity_state$result)
        .check_fingerprint()
        df <- build_velocity_validation_summary(velocity_state$result)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    # Table d'alignement cellules/genes produite a la validation.
    output$dl_velocity_mapping <- downloadHandler(
      filename = function() paste0("velocity_alignment_mapping_", format(Sys.Date(), "%Y-%m-%d"), ".csv"),
      content = function(file) {
        req(velocity_state$result)
        .check_fingerprint()
        df <- build_velocity_alignment_mapping(velocity_state$result)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    # Resultat canonique complet (RDS) — identite, statut, reglages et
    # avertissements embarques dans l'objet ; analysis_id dans le nom de fichier.
    output$dl_velocity_result_rds <- downloadHandler(
      filename = function() {
        aid <- velocity_state$result$analysis_id %||% "sc-velocity"
        paste0("velocity_result_", aid, "_", format(Sys.Date(), "%Y-%m-%d"), ".rds")
      },
      content = function(file) {
        req(velocity_state$result)
        .check_fingerprint()
        saveRDS(velocity_state$result, file)
      }
    )

    # ── Exports figures Stage 9 ─────────────────────────────────────────────
    # Figures PNG/PDF en portee COMPLETE ("full") — la portee est ecrite dans
    # le sous-titre de la figure ; nom de fichier trace par analysis_id.
    output$dl_velocity_phase_png <- downloadHandler(
      filename = function() velocity_export_filename(velocity_state$result, "velocity_phase_portrait", "png"),
      content = function(file) {
        req(velocity_state$result); req(input$velocity_gene)
        .check_fingerprint()
        p <- plot_velocity_phase_portrait(
          velocity_state$result, input$velocity_gene,
          data_scope = "full", seurat_obj = global_data$sc_obj
        )
        ts_export_plot(file, p, width = 8, height = 6, dpi = 300)
      }
    )

    output$dl_velocity_phase_pdf <- downloadHandler(
      filename = function() velocity_export_filename(velocity_state$result, "velocity_phase_portrait", "pdf"),
      content = function(file) {
        req(velocity_state$result); req(input$velocity_gene)
        .check_fingerprint()
        p <- plot_velocity_phase_portrait(
          velocity_state$result, input$velocity_gene,
          data_scope = "full", seurat_obj = global_data$sc_obj
        )
        ts_export_plot(file, p, width = 8, height = 6)
      }
    )

    # Figure embedding en portee COMPLETE — partagee PNG/PDF. Coordonnees
    # introuvables => notification + annulation explicite (jamais de
    # substitution silencieuse).
    .export_embedding_figure <- function(file) {
      red <- input$velocity_reduction %||% "umap"
      coords <- .base_coords_for(red)
      if (is.null(coords)) {
        msg <- if (!.velocity_has_reduction(global_data$sc_obj, red)) {
          .velocity_missing_reduction_msg(red)
        } else {
          sprintf(paste0("Reduction '%s' introuvable ou aucune cellule ",
                         "velocity dans l'embedding : export annule."), red)
        }
        showNotification(msg, type = "error", duration = 8)
        req(FALSE)
      }
      p <- plot_velocity_embedding(
        velocity_state$result, coords,
        selected_reduction = red, data_scope = "full",
        seurat_obj = global_data$sc_obj
      )
      ts_export_plot(file, p, width = 8, height = 6, dpi = 300)
    }

    output$dl_velocity_embedding_png <- downloadHandler(
      filename = function() velocity_export_filename(velocity_state$result, "velocity_embedding", "png"),
      content = function(file) {
        req(velocity_state$result)
        .check_fingerprint()
        .export_embedding_figure(file)
      }
    )

    output$dl_velocity_embedding_pdf <- downloadHandler(
      filename = function() velocity_export_filename(velocity_state$result, "velocity_embedding", "pdf"),
      content = function(file) {
        req(velocity_state$result)
        .check_fingerprint()
        .export_embedding_figure(file)
      }
    )

    # Vecteurs/alignement par cellule (CSV) — depuis le resultat canonique.
    output$dl_velocity_cell_vectors <- downloadHandler(
      filename = function() velocity_export_filename(velocity_state$result, "velocity_cell_vectors", "csv"),
      content = function(file) {
        req(velocity_state$result)
        .check_fingerprint()
        df <- build_velocity_cell_vectors_export(velocity_state$result)
        utils::write.csv(df, file, row.names = FALSE)
      }
    )

    # Expose state for tests (optional)
    return(velocity_state)
  })
}
