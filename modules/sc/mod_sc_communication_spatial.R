# =============================================================================
# mod_sc_communication_spatial.R — onglet « Contexte spatial » du panneau
# Communication (roadmap avancée Phase 1 / V1.x-A)
# =============================================================================
# Orchestration uniquement : la science residue dans
# R/sc/sc_communication_spatial.R (contrat COMMUNICATION_SPATIAL_CONTRACT).
# Le module DECLAR explicitement la source de coordonnees (objet spatial
# courant = unites natives, ou reduction 2D = distances de projection NON
# physiques), le rayon (jamais infere) et les permutations (opt-in, seed
# fixee) ; la provenance est PRODUITE au calcul et appended a l'etat partage.
# Aucun recalcul de score, aucune modification du resultat canonique.
# =============================================================================

mod_sc_communication_spatial_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "alert alert-warning", style = "font-size:0.85em;padding:6px;",
        i18n$t("La distance est une contrainte spatiale — une compatibilité spatiale n'établit aucune communication biologique. ")),
    fluidRow(
      column(6, radioButtons(ns("spc_source"), i18n$t("Source de coordonnées"),
                             choices = setNames(c("reduction", "spatial_object"),
                                                c(.tr_plain("Projection 2D (réduction — distances NON physiques)"),
                                                  .tr_plain("Coordonnées spatiales (bundle spatial courant)"))),
                             selected = "reduction", width = "100%")),
      column(6,
             conditionalPanel(
               condition = "input.spc_source == 'reduction'", ns = ns,
               selectInput(ns("spc_reduction"), i18n$t("Réduction 2D"),
                           choices = character(0), width = "100%")),
             conditionalPanel(
               condition = "input.spc_source == 'spatial_object'", ns = ns,
               div(class = "small text-muted mt-4",
                   textOutput(ns("spc_spatial_info")))))),
    fluidRow(
      column(4, numericInput(ns("spc_radius"), i18n$t("Rayon (NA = non calculé)"),
                             value = NA, min = 0, width = "100%")),
      column(4, numericInput(ns("spc_nperm"), i18n$t("Permutations (0 = aucune)"),
                             value = 0, min = 0, max = 1000, step = 1, width = "100%")),
      column(4, div(class = "small text-muted mt-4",
                    i18n$t("Seed fixée : 1 (déterministe).")))),
    actionButton(ns("spc_compute"), i18n$t("Calculer le contexte spatial"),
                 class = "btn-primary w-100", icon = icon("ruler")),
    div(class = "small text-muted mt-1", textOutput(ns("spc_status"))),
    hr(),
    # P0.4 harmonisation (2026-11-14, demande opérateur) : plots AU-DESSUS de
    # la table — le sélecteur de paire pilote les plots, il les précède.
    # spc_table garde son scrollY CQF-1 (défile DANS son conteneur — épinglé
    # par test-sc-communication-spatial.R), NON converti au contrat non-fill.
    fluidRow(
      column(6, selectizeInput(ns("spc_pair"), i18n$t("Paire sender → receiver"),
                               choices = character(0), width = "100%")),
      column(6, downloadButton(ns("dl_spc_csv"), i18n$t("CSV (données du graphe)"), class = "btn-sm btn-outline-secondary mt-4"))),
    fluidRow(
      column(6, plotOutput(ns("spc_edges"), height = "360px"),
             fluidRow(
               column(6, downloadButton(ns("dl_spc_edges_png"), i18n$t("PNG"), class = "btn-sm btn-outline-secondary")),
               column(6, downloadButton(ns("dl_spc_edges_pdf"), i18n$t("PDF"), class = "btn-sm btn-outline-secondary")))),
      column(6, plotOutput(ns("spc_distances"), height = "360px"),
             fluidRow(
               column(6, downloadButton(ns("dl_spc_dist_png"), i18n$t("PNG"), class = "btn-sm btn-outline-secondary")),
               column(6, downloadButton(ns("dl_spc_dist_pdf"), i18n$t("PDF"), class = "btn-sm btn-outline-secondary"))))),
    DT::dataTableOutput(ns("spc_table"), height = "auto", fill = FALSE)
  )
}

mod_sc_communication_spatial_server <- function(id, comm_state, global_data, shared_rv) {
  moduleServer(id, function(input, output, session) {
    .tr <- function(key) {
      tr <- isolate(global_data$i18n)
      if (is.null(tr)) return(key)
      tryCatch(.strip_i18n_html(tr$t(key)), error = function(e) key)
    }
    .error_plot <- function(e) {
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = conditionMessage(e), size = 5, colour = "red") +
        ggplot2::theme_void()
    }
    spc_ctx <- reactiveVal(NULL)
    # PAS 1 (idéation saves 2026-11-17) : aide de publication de l'état léger.
    # ⚠️ JUMELLE de la même aide dans mod_sc_communication_velocity.R et
    # _trajectory.R — read-modify-write explicite (affectation imbriquée
    # rv$a$b <- x sur champ NULL non garantie par reactiveValues).
    .ctx_view_set <- function(drawer, fields) {
      cur <- comm_state$ctx_view
      if (!is.list(cur)) cur <- list()
      cur[[drawer]] <- fields
      comm_state$ctx_view <- cur
    }

    # R-3bis (idiome repo, cf. mod_bulk.R fix 2026-11-06) : une sortie
    # invalidée pendant que son panneau est caché ne livre JAMAIS son payload.
    # Les renders consomment le miroir de feuille active du drawer (publié par
    # le parent) — re-exécution à l'ouverture du sous-onglet.
    .ctx_tab_guard <- function() {
      req(state_get(shared_rv, "comm_tabs_current") %||% "" == "ccc_tab_spatial")
    }
    # Lane 3c (2026-11-14, retour live « white card ») : résultat source du
    # contexte courant — permet au reset de n'invalider que les contextes
    # PÉRIMÉS (voir l'observateur reset ci-dessous).
    spc_src <- reactiveVal(NULL)

    # Garde de peremption (resultat canonique vs objet courant) — garde legere,
    # le portage complet reste dans le module parent.
    .guard_stale <- function() {
      req(global_data$sc_obj, comm_state$result, comm_state$object_fingerprint)
      if (!identical(velocity_object_fingerprint(global_data$sc_obj),
                     comm_state$object_fingerprint)) {
        shiny::validate(shiny::need(FALSE, .tr("Résultat de communication périmé : réimportez.")))
      }
    }

    # Nouveau resultat importe -> contexte precedant perime (jamais affiche).
    # Lane 3c (retour live « white card », mesuré au log DIAG 2026-11-14) :
    # au reload, le parent fait clear (result <- NULL) PUIS restaure le
    # résultat — l'auto-calcul s'exécute entre les deux (il voit la valeur
    # COURANTE, déjà restaurée), et ce reset — dont le 1er run est sauté par
    # ignoreInit — arrivait APRÈS et effaçait un contexte calculé SUR le
    # résultat courant ; le garde anti-recalcul de l'auto ne re-déclenchant
    # pas sur spc_ctx, le panneau restait vide jusqu'au clic. Invalider
    # UNIQUEMENT si le résultat courant diffère de celui du contexte.
    observeEvent(comm_state$result, {
      if (identical(comm_state$result, spc_src())) return()
      spc_ctx(NULL)
      spc_src(NULL)
    }, ignoreInit = TRUE)

    # Redutions 2D disponibles (>= 2 colonnes), defaut umap — PAS 1 : la même
    # aide sert à l'observateur d'objet (peuplement) et à l'observateur de
    # restauration (garde de validité de la réduction sauvegardée).
    .spc_valid_reductions <- function() {
      obj <- global_data$sc_obj
      reds <- if (!is.null(obj)) names(obj@reductions) else character(0)
      reds[vapply(reds, function(rd) {
        em <- tryCatch(SeuratObject::Embeddings(obj, rd), error = function(e) NULL)
        !is.null(em) && ncol(em) >= 2L
      }, logical(1))]
    }
    observeEvent(global_data$sc_obj, {
      reds2 <- .spc_valid_reductions()
      updateSelectInput(session, "spc_reduction", choices = reds2,
                        selected = if ("umap" %in% reds2) "umap" else reds2[1])
    }, ignoreInit = FALSE)

    # PAS 1 (idéation saves 2026-11-17) : restauration des params de vue du
    # tiroir spatial — radio/numériques en direct ; la réduction n'est
    # ré-appliquée que si elle reste VALIDE sur l'objet courant (même garde
    # que le peuplement). Un champ absent du snapshot = inchangé.
    observeEvent(comm_state$ctx_view$spatial, {
      v <- comm_state$ctx_view$spatial
      if (!is.list(v)) return()
      if (!is.null(v$source) && !identical(v$source, input$spc_source)) {
        shiny::updateRadioButtons(session, "spc_source", selected = v$source)
      }
      if (!is.null(v$radius) && !is.na(v$radius) &&
          !identical(v$radius, input$spc_radius)) {
        shiny::updateNumericInput(session, "spc_radius", value = v$radius)
      }
      if (!is.null(v$nperm) && !identical(v$nperm, input$spc_nperm)) {
        shiny::updateNumericInput(session, "spc_nperm", value = v$nperm)
      }
      red <- v$reduction
      if (!is.null(red) && !identical(red, input$spc_reduction) &&
          red %in% .spc_valid_reductions()) {
        updateSelectInput(session, "spc_reduction", selected = red)
      }
    }, ignoreNULL = TRUE, ignoreInit = FALSE)

    # PAS 1 : publication de l'état léger — les params courants du tiroir.
    # Garde : un contexte établi (au calcul, ou après l'auto-calc Lane 3b qui
    # re-compose le contexte restauré). Sans contexte, RIEN n'est publié :
    # au restore, les observateurs d'entrée déclenchés par la ré-application
    # ne peuvent pas écraser l'état restauré.
    .spc_publish_view <- function() {
      if (is.null(spc_ctx())) return()
      .ctx_view_set("spatial", list(
        source    = input$spc_source,
        reduction = input$spc_reduction,
        radius    = input$spc_radius %||% NA_real_,
        nperm     = as.integer(input$spc_nperm %||% 0L),
        pair      = input$spc_pair %||% character(0)
      ))
    }
    observeEvent(input$spc_source,    .spc_publish_view(), ignoreInit = TRUE)
    observeEvent(input$spc_reduction, .spc_publish_view(), ignoreInit = TRUE)
    observeEvent(input$spc_radius,    .spc_publish_view(), ignoreInit = TRUE)
    observeEvent(input$spc_nperm,     .spc_publish_view(), ignoreInit = TRUE)
    observeEvent(input$spc_pair,      .spc_publish_view(), ignoreInit = TRUE)

    output$spc_spatial_info <- renderText({
      b <- global_data$spatial_obj
      n <- if (!is.null(b) && !is.null(b$coords)) nrow(b$coords) else 0L
      if (n == 0L) {
        .tr("Aucune coordonnée disponible dans le bundle spatial courant.")
      } else {
        sprintf(.tr("Bundle spatial courant : %d spot(s)/cellule(s) (coordonnées d'import, colonnes id/x/y)."), n)
      }
    })

    # Identites nommees depuis la colonne D'IMPORT (enregistree dans le
    # resultat canonique) sur l'objet courant.
    .identities <- function() {
      obj <- global_data$sc_obj
      col <- comm_state$result$identity_column %||% NA_character_
      if (is.null(obj) || is.na(col) || !col %in% colnames(obj@meta.data)) return(NULL)
      meta <- obj@meta.data
      if (!is.null(comm_state$context_cells)) {
        meta <- meta[rownames(meta) %in% as.character(comm_state$context_cells), , drop = FALSE]
      }
      setNames(as.character(meta[[col]]), rownames(meta))
    }

    .coordinates <- function() {
      if (identical(input$spc_source, "spatial_object")) {
        b <- global_data$spatial_obj
        coords <- b$coords %||% NULL
        if (is.null(coords) || !all(c("id", "x", "y") %in% names(coords))) {
          stop("Aucune coordonnée id/x/y disponible dans le bundle spatial courant.", call. = FALSE)
        }
        coords <- coords[, c("id", "x", "y"), drop = FALSE]
        if (!is.null(comm_state$context_cells)) {
          coords <- coords[as.character(coords$id) %in%
                            as.character(comm_state$context_cells), , drop = FALSE]
        }
        coords
      } else {
        obj <- global_data$sc_obj
        red <- input$spc_reduction
        if (is.null(obj) || is.null(red) || !nzchar(red)) {
          stop("Choisissez d'abord une réduction 2D.", call. = FALSE)
        }
        em <- SeuratObject::Embeddings(obj, red)[, 1:2, drop = FALSE]
        if (!is.null(comm_state$context_cells)) {
          em <- em[rownames(em) %in% as.character(comm_state$context_cells), , drop = FALSE]
        }
        colnames(em) <- c("x", "y")
        as.data.frame(em)
      }
    }

    .spc_compute_ctx <- function() {
      .guard_stale()
      tryCatch({
        idents <- .identities()
        if (is.null(idents)) {
          stop(sprintf(
            "Colonne d'identités '%s' introuvable sur l'objet courant — réimportez la communication.",
            comm_state$result$identity_column %||% NA_character_
          ), call. = FALSE)
        }
        ctx <- build_communication_spatial_context(
          comm_state$result,
          coordinates       = .coordinates(),
          cell_identities   = idents,
          identity_column   = comm_state$result$identity_column,
          coordinate_source = input$spc_source,
          coordinate_units  = if (identical(input$spc_source, "spatial_object")) "unités natives" else NA_character_,
          radius            = input$spc_radius %||% NA_real_,
          n_permutations    = as.integer(input$spc_nperm %||% 0L),
          perm_seed         = 1L
        )
        spc_ctx(ctx)
        spc_src(comm_state$result)
        # PAS 1 (idéation saves 2026-11-17) : params qui ONT PRODUIT ce
        # contexte publiés pour le snapshot (la paire suivra au round-trip
        # client via l'observateur input$spc_pair).
        .ctx_view_set("spatial", list(
          source    = input$spc_source,
          reduction = input$spc_reduction,
          radius    = input$spc_radius %||% NA_real_,
          nperm     = as.integer(input$spc_nperm %||% 0L)
        ))
        # Provenance PRODUITE au calcul (regle 7) — appended, jamais reconstruite.
        if (!is.null(shared_rv)) provenance_append(shared_rv, ctx$provenance)
        pairs <- paste(ctx$pair_table$sender_node, ctx$pair_table$receiver_node, sep = "\u2192")
        # PAS 1 : la paire sauvegardée est honorée si elle existe toujours
        # (restauration = même vue qu'avant la sauvegarde) ; sinon le défaut
        # historique (première paire, espace des libellés — inchangé).
        pair_vals <- paste(ctx$pair_table$sender_node, ctx$pair_table$receiver_node, sep = "\r")
        saved_pair <- if (is.list(comm_state$ctx_view$spatial)) comm_state$ctx_view$spatial$pair else NULL
        sel <- if (length(pair_vals) && length(saved_pair) && saved_pair %in% pair_vals) {
          saved_pair
        } else if (length(pairs)) pairs[1] else character(0)
        updateSelectizeInput(session, "spc_pair", choices = setNames(
          paste(ctx$pair_table$sender_node, ctx$pair_table$receiver_node, sep = "\r"), pairs
        ), selected = sel)
        showNotification(
          sprintf(.tr("Contexte spatial : %d paire(s) cartographiée(s) sur %d."),
                  ctx$qc$n_pairs_mapped, ctx$qc$n_pairs_total),
          type = "message", duration = 5
        )
      }, error = function(e) {
        spc_ctx(NULL)
        spc_src(NULL)
        showNotification(paste(.tr("Erreur contexte spatial :"), conditionMessage(e)),
                         type = "error", duration = 10)
      })
    }

    observeEvent(input$spc_compute, { .spc_compute_ctx() })

    # Lane 3b (retour opérateur, 2026-11-14) : AUTO-CALCUL du contexte spatial
    # — nécessite une réduction 2D sélectionnée (l'observer d'objet peuple les
    # choix) ; tant qu'elle manque, on attend (l'observer se re-déclenche quand
    # les choix arrivent). Le bouton reste (recalcul après changement de
    # source/rayon/permutations).
    observeEvent(c(comm_state$result, global_data$sc_obj, input$spc_reduction), {
      if (is.null(comm_state$result)) return()
      if (!is.null(spc_ctx())) return()
      red <- input$spc_reduction
      if (is.null(red) || !nzchar(red)) return()
      if (is.null(.identities())) return()
      .spc_compute_ctx()
    }, ignoreInit = FALSE)

    output$spc_status <- renderText({
      .ctx_tab_guard()
      ctx <- spc_ctx()
      req(ctx)
      sprintf(.tr("Statut : %s — %d cellule(s) cartographiée(s) ; %d paire(s) non cartographiée(s) (listées dans les avertissements)."),
              ctx$status, ctx$qc$n_cells_matched, ctx$qc$n_pairs_unmapped)
    })

    output$spc_table <- DT::renderDataTable({
      .ctx_tab_guard()
      ctx <- spc_ctx()
      req(ctx)
      # CQF-1 (fiche SC_COMM_QUICKFIXES_PROPOSAL_2026-11-02 §3) : scrollY
      # natif DT via extra_options — la table défile DANS son conteneur
      # 320px (UI, ligne 45) au lieu de déborder sur le fluidRow des plots.
      # scrollX reste hors défaut (2H/O28).
      ts_datatable(ctx$pair_table, page_length = 15L,
                   filename_base = "sc_comm_spatial_pair_table", filter = "none",
                   extra_options = list(scrollY = "260px", scrollCollapse = TRUE))
    })

    .selected_pair <- function() {
      ctx <- spc_ctx()
      k <- input$spc_pair
      if (is.null(ctx) || is.null(k) || !nzchar(k)) return(NULL)
      sp <- strsplit(k, "\r", fixed = TRUE)[[1L]]
      list(sender = sp[1L], receiver = sp[2L])
    }

    output$spc_edges <- renderPlot({
      .ctx_tab_guard()
      ctx <- spc_ctx()
      req(ctx)
      pr <- .selected_pair()
      req(pr)
      tryCatch(
        plot_communication_spatial_edges(
          ctx, pr$sender, pr$receiver,
          coordinates = .coordinates(), cell_identities = .identities()
        ),
        error = function(e) .error_plot(e)
      )
    })

    output$spc_distances <- renderPlot({
      .ctx_tab_guard()
      ctx <- spc_ctx()
      req(ctx)
      tryCatch(plot_communication_spatial_distance_summary(ctx), error = function(e) .error_plot(e))
    })

    output$dl_spc_csv <- downloadHandler(
      filename = function() {
        ctx <- spc_ctx()
        if (is.null(ctx)) return("communication_spatial_context.csv")
        paste0("communication_spatial_context_", ctx$analysis_id, "_", ctx$timestamp_utc, ".csv")
      },
      content = function(file) {
        ctx <- spc_ctx()
        req(ctx)
        .guard_stale()
        utils::write.csv(build_communication_spatial_export(ctx), file, row.names = FALSE)
      }
    )

    .fig_download <- function(output_id, kind, plot_fn) {
      output[[output_id]] <- downloadHandler(
        filename = function() {
          ctx <- spc_ctx()
          if (is.null(ctx)) return(paste0(kind, ".png"))
          paste0(kind, "_", ctx$analysis_id, "_", ctx$timestamp_utc, ".", tolower(tools::file_ext(output_id)))
        },
        content = function(file) {
          ctx <- spc_ctx()
          req(ctx)
          .guard_stale()
          ts_export_plot(file, plot_fn(ctx), width = 8, height = 6, dpi = 300)
        }
      )
    }
    .fig_download("dl_spc_edges_png", "communication_spatial_edges", function(ctx) {
      pr <- .selected_pair(); req(pr)
      plot_communication_spatial_edges(ctx, pr$sender, pr$receiver,
                                       coordinates = .coordinates(),
                                       cell_identities = .identities())
    })
    .fig_download("dl_spc_edges_pdf", "communication_spatial_edges", function(ctx) {
      pr <- .selected_pair(); req(pr)
      plot_communication_spatial_edges(ctx, pr$sender, pr$receiver,
                                       coordinates = .coordinates(),
                                       cell_identities = .identities())
    })
    .fig_download("dl_spc_dist_png", "communication_spatial_distances", function(ctx) {
      plot_communication_spatial_distance_summary(ctx)
    })
    .fig_download("dl_spc_dist_pdf", "communication_spatial_distances", function(ctx) {
      plot_communication_spatial_distance_summary(ctx)
    })

    invisible(NULL)
  })
}
