# =============================================================================
# mod_sc_communication_velocity.R — onglet « Contexte vélocité » du panneau
# Communication (roadmap avancée Phase 2 / V1.x-C)
# =============================================================================
# Orchestration uniquement : la science residue dans
# R/sc/sc_communication_velocity.R (contrat COMMUNICATION_VELOCITY_CONTRACT).
# Le module consomme le résultat velocity canonique (shared_rv$velocity_result,
# contrat Stage 10) et le seuil de recouvrement déclaré
# (TS_VELOCITY_OVERLAP_MIN de config/thresholds.R — jamais codé en dur).
# Sans vecteur validé : message explicite (état unavailable_no_vectors) —
# aucune substitution, aucune ré-inférence.
# =============================================================================

mod_sc_communication_velocity_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "alert alert-warning", style = "font-size:0.85em;padding:6px;",
        i18n$t("« Associé à l'état de vitesse » n'établit aucune causalité : la magnitude des vecteurs précalculés est une quantité descriptive par population.")),
    div(class = "small text-muted mb-2", textOutput(ns("vel_info"))),
    # 3F (O28, spéc §2.18 point 1) : garde PRÉVENTIVE — bandeau info + bouton
    # désactivé avant tout clic quand les vecteurs sont absents/inutilisables
    # (le constat 3F : la seule garde existante — stop :84 — était post-clic).
    uiOutput(ns("vel_guard")),
    fluidRow(
      column(6, numericInput(ns("vel_min_overlap"), i18n$t("Recouvrement minimal (fraction)"),
                             value = TS_VELOCITY_OVERLAP_MIN, min = 0, max = 1, step = 0.05,
                             width = "100%")),
      column(6, div(class = "small text-muted mt-4",
                    i18n$t("Défaut consommé de config/thresholds.R (TS_VELOCITY_OVERLAP_MIN).")))),
    actionButton(ns("vel_compute"), i18n$t("Calculer le contexte vélocité"),
                 class = "btn-primary w-100", icon = icon("tachometer-alt")),
    div(class = "small text-muted mt-1", textOutput(ns("vel_status"))),
    div(class = "small text-muted", textOutput(ns("vel_dir_summary"))),
    hr(),
    # P0.4 harmonisation (2026-11-14, demande opérateur) : plot AU-DESSUS de la
    # table + k1 complétion — hauteur fixe < contenu paginé = débord du wrapper
    # DataTables sur le voisin (mesuré live : conteneur ~400 px vs 807 px).
    plotOutput(ns("vel_plot"), height = "420px"),
    fluidRow(
      column(3, downloadButton(ns("dl_vel_png"), i18n$t("PNG"), class = "btn-sm btn-outline-secondary")),
      column(3, downloadButton(ns("dl_vel_pdf"), i18n$t("PDF"), class = "btn-sm btn-outline-secondary")),
      column(3, downloadButton(ns("dl_vel_csv"), i18n$t("CSV"), class = "btn-sm btn-outline-secondary"))),
    DT::dataTableOutput(ns("vel_table"), height = "auto", fill = FALSE)
  )
}

mod_sc_communication_velocity_server <- function(id, comm_state, global_data, shared_rv) {
  moduleServer(id, function(input, output, session) {
    .tr <- function(key) {
      tr <- isolate(global_data$i18n)
      if (is.null(tr)) return(key)
      tryCatch(.strip_i18n_html(tr$t(key)), error = function(e) key)
    }
    vel_ctx <- reactiveVal(NULL)
    observeEvent(vel_ctx(), {
      cur <- comm_state$ctx_results %||% list()
      cur$velocity <- vel_ctx()
      comm_state$ctx_results <- cur
    }, ignoreNULL = FALSE)
    # Lane 3c (2026-11-14, retour live « white card ») : résultat source du
    # contexte courant — permet au reset de n'invalider que les contextes
    # PÉRIMÉS (voir l'observateur reset ci-dessous).
    vel_src <- reactiveVal(NULL)

    # PAS 1 (idéation saves 2026-11-17) : publication de l'état léger de ce
    # tiroir dans comm_state$ctx_view (miroir snapshot par le parent).
    # ⚠️ JUMELLE de la même aide dans mod_sc_communication_spatial.R et
    # _trajectory.R — affectation imbriquée NON utilisée (rv$a$b <- x sur un
    # champ NULL n'est pas garantie par reactiveValues) : read-modify-write
    # explicite, une seule clé top-level = invalidation propre des observateurs.
    .ctx_view_set <- function(drawer, fields) {
      cur <- comm_state$ctx_view
      if (!is.list(cur)) cur <- list()
      cur[[drawer]] <- fields
      comm_state$ctx_view <- cur
    }

    # PAS 1 : restauration — ré-application des params de vue sauvegardés
    # (snapshot PAS 1). Numérique simple : mise à jour directe si elle diffère.
    observeEvent(comm_state$ctx_view$velocity, {
      v <- comm_state$ctx_view$velocity
      if (!is.list(v)) return()
      mo <- v$min_overlap
      if (!is.null(mo) && !is.na(mo) && !identical(mo, input$vel_min_overlap)) {
        shiny::updateNumericInput(session, "vel_min_overlap", value = mo)
      }
    }, ignoreNULL = TRUE, ignoreInit = FALSE)

    # R-3bis (idiome repo, cf. mod_bulk.R fix 2026-11-06 ; mesuré aussi ici
    # sur l'achèvement async CellChat) : une sortie invalidée pendant que son
    # panneau est caché ne livre JAMAIS son payload ensuite. Les renders de ce
    # panneau consomment le miroir de feuille active publié par le parent
    # (shared_rv$comm_tabs_current) — re-exécution à l'ouverture du sous-onglet,
    # livraison garantie. req() silencieux tant que la feuille n'est pas la nôtre.
    .ctx_tab_guard <- function() {
      req(state_get(shared_rv, "comm_tabs_current") %||% "" == "ccc_tab_velocite")
    }

    .guard_stale <- function() {
      req(global_data$sc_obj, comm_state$result, comm_state$object_fingerprint)
      if (!identical(velocity_object_fingerprint(global_data$sc_obj),
                     comm_state$object_fingerprint)) {
        shiny::validate(shiny::need(FALSE, .tr("Résultat de communication périmé : réimportez.")))
      }
    }

    observeEvent(comm_state$result, {
      # Lane 3c (retour live « white card », mesuré au log DIAG 2026-11-14) :
      # au reload, le parent fait clear (result <- NULL) PUIS restaure le
      # résultat — l'auto-calcul s'exécute entre les deux (il voit la valeur
      # COURANTE, déjà restaurée), et ce reset — dont le 1er run est sauté par
      # ignoreInit — arrivait APRÈS et effaçait un contexte calculé SUR le
      # résultat courant ; le garde anti-recalcul de l'auto ne re-déclenchant
      # pas sur vel_ctx, le panneau restait vide jusqu'au clic. Invalider
      # UNIQUEMENT si le résultat courant diffère de celui du contexte.
      if (identical(comm_state$result, vel_src())) return()
      vel_ctx(NULL)
      vel_src(NULL)
    }, ignoreInit = TRUE)

    output$vel_info <- renderText({
      vr <- if (!is.null(shared_rv)) state_get(shared_rv, "velocity_result") else NULL
      if (is.null(vr)) {
        .tr("Aucun résultat velocity importé : importez d'abord des vecteurs précalculés dans le panneau Vélocité (8.6).")
      } else {
        sprintf(.tr("Résultat velocity consommé : statut « %s », %d cellule(s) alignée(s). Les vecteurs ne sont jamais ré-inférés."),
                vr$status, length(vr$cell_names %||% character(0)))
      }
    })

    # ── 3F (O28, spéc §2.18 point 1) : garde préventive ──────────────────────
    # États utilisables du contrat velocity (sc_velocity.R:141-151) :
    # valid / valid_partial_embedding portent des vecteurs consommables ;
    # valid_no_vectors (objet valide, AUCUN vecteur) et les invalid_* n'en
    # portent pas. Préventif = bandeau + bouton désactivé AVANT le clic ; le
    # stop post-clic (:84) reste en seconde ligne de défense.
    .vel_vr <- reactive({
      if (!is.null(shared_rv)) state_get(shared_rv, "velocity_result") else NULL
    })
    .vel_ready <- function(vr) {
      st <- if (!is.null(vr)) vr$status %||% NA_character_ else NA_character_
      !is.null(vr) && st %in% c("valid", "valid_partial_embedding")
    }
    observeEvent(.vel_vr(), {
      vr <- .vel_vr()
      shiny::updateActionButton(session, "vel_compute", disabled = !.vel_ready(vr))
    })
    output$vel_guard <- renderUI({
      vr <- .vel_vr()
      if (.vel_ready(vr)) return(NULL)
      msg <- if (is.null(vr)) {
        .tr("Le calcul de contexte nécessite des vecteurs de vélocité précalculés complets — aucun résultat velocity importé.")
      } else if (identical(vr$status, "valid_no_vectors")) {
        .tr("Le calcul de contexte nécessite des vecteurs de vélocité précalculés complets (absents dans l'import actuel).")
      } else {
        sprintf(.tr("Statut velocity « %s » : vecteurs inutilisables — réimportez depuis le panneau Vélocité (8.6)."),
                as.character(vr$status %||% "inconnu"))
      }
      div(class = "alert alert-info", style = "font-size:0.85em;padding:6px;",
          icon("circle-info"), " ", msg)
    })

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

    # W2.2 §4.3 : directionnalité (alignement T_ij + corrélation vitesse-récepteur).
    # Descriptif uniquement (aucune causalité). T_ij n'est JAMAIS fabriquée ici :
    # uniquement lue si présente dans le résultat velocity ; toute absence => NULL.
    # Résultat publié dans session$userData$vel_directionality (aucun nouveau texte UI).
    vel_directionality <- reactiveVal(NULL)
    .vel_directional <- function(vr, idents) {
      if (!exists("compute_communication_velocity_directionality", mode = "function")) return(NULL)
      tryCatch({
        tm <- vr$transition_matrix
        if (is.null(tm) || !(is.matrix(tm) || inherits(tm, "Matrix"))) tm <- NULL
        rexpr <- tryCatch({
          tbl <- comm_state$result$canonical_table
          recs <- unique(as.character(tbl$receptor))
          obj <- global_data$sc_obj
          mat <- SeuratObject::GetAssayData(obj, layer = "data")
          recs <- intersect(recs[!is.na(recs)], rownames(mat))
          if (length(recs) > 0L) mat[recs, , drop = FALSE] else NULL
        }, error = function(e) NULL)
        has_vec <- !is.null(vr$velocity_vectors)
        if (is.null(tm) && (is.null(rexpr) || !has_vec)) return(NULL)
        compute_communication_velocity_directionality(
          comm_state$result, velocity_result = vr, cell_identities = idents,
          transition_matrix = tm, receptor_expr = rexpr
        )
      }, error = function(e) NULL)
    }

    .vel_compute_ctx <- function() {
      .guard_stale()
      tryCatch({
        vr <- if (!is.null(shared_rv)) state_get(shared_rv, "velocity_result") else NULL
        if (is.null(vr)) {
          stop("Aucun résultat velocity importé : utilisez d'abord le panneau Vélocité (8.6).", call. = FALSE)
        }
        idents <- .identities()
        if (is.null(idents)) {
          stop(sprintf("Colonne d'identités '%s' introuvable — réimportez la communication.",
                       comm_state$result$identity_column %||% NA), call. = FALSE)
        }
        # 3F (O28, spéc §2.18) : jauge de progression — le constat 3F mesurait
        # 0 withProgress dans ce module (calcul potentiellement long, gel
        # perçu de la session).
        ctx <- withProgress(message = .tr("Calcul du contexte vélocité..."), value = 0.4, {
          build_communication_velocity_context(
            comm_state$result,
            velocity_result  = vr,
            cell_identities  = idents,
            identity_column  = comm_state$result$identity_column,
            min_overlap      = input$vel_min_overlap %||% NA_real_
          )
        })
        vel_ctx(ctx)
        vel_src(comm_state$result)
        vel_directionality(.vel_directional(vr, idents))
        session$userData$vel_directionality <- vel_directionality()
        # PAS 1 (idéation saves 2026-11-17) : les params qui ONT PRODUIT ce
        # contexte sont publiés pour le snapshot — au restore, l'auto-calcul
        # (Lane 3b) recompose le MÊME contexte avec ces params.
        .ctx_view_set("velocity", list(min_overlap = input$vel_min_overlap %||% NA_real_))
        if (!is.null(shared_rv)) provenance_append(shared_rv, ctx$provenance)
        if (identical(ctx$status, "valid")) {
          showNotification(
            sprintf(.tr("Contexte vélocité : %d paire(s) cartographiée(s) — recouvrement %.1f%%."),
                    ctx$qc$n_pairs_mapped, 100 * ctx$qc$overlap_fraction),
            type = "message", duration = 5
          )
        } else {
          showNotification(
            sprintf(.tr("Contexte vélocité indisponible (statut : %s) — voir les avertissements."),
                    ctx$status),
            type = "warning", duration = 8
          )
        }
      }, error = function(e) {
        vel_ctx(NULL)
        vel_src(NULL)
        showNotification(paste(.tr("Erreur contexte vélocité :"), conditionMessage(e)),
                         type = "error", duration = 10)
      })
    }

    observeEvent(input$vel_compute, { .vel_compute_ctx() })

    # Lane 3b (retour opérateur, 2026-11-14) : AUTO-CALCUL du contexte dès que
    # les entrées sont réunies et qu'aucun contexte n'est établi — plus de
    # clic requis après un reload. Le bouton reste (recalcul explicite après
    # changement de paramètre). Garde 3F IDENTIQUE (.vel_ready) ; identités
    # absentes = pas d'auto (échec garanti sinon).
    observeEvent(c(comm_state$result, .vel_vr(), global_data$sc_obj), {
      if (is.null(comm_state$result)) return()
      if (!is.null(vel_ctx())) return()
      if (!.vel_ready(.vel_vr())) return()
      if (is.null(.identities())) return()
      .vel_compute_ctx()
    }, ignoreInit = FALSE)

    output$vel_status <- renderText({
      .ctx_tab_guard()
      ctx <- vel_ctx()
      req(ctx)
      if (identical(ctx$status, "valid")) {
        sprintf(.tr("Statut : %s — %d paire(s) cartographiée(s) ; recouvrement %.1f%%."),
                ctx$status, ctx$qc$n_pairs_mapped, 100 * ctx$qc$overlap_fraction)
      } else {
        paste(.tr("Statut :"), ctx$status, "—", utils::head(ctx$warnings, 1))
      }
    })

    # W2.2 §4.3 : synthèse directionnelle concise (descriptive, sans causalité,
    # aucune clé i18n : texte numérique/codes de statut uniquement). NULL => "".
    output$vel_dir_summary <- renderText({
      d <- vel_directionality()
      if (is.null(d) || !is.list(d)) return("")
      pt <- d$pair_table
      n_map <- if (is.data.frame(pt)) nrow(pt) else 0L
      if (n_map == 0L) return("")
      if (identical(d$transition_status, "valid")) {
        n_al <- sum(pt$transition_aligned %in% TRUE)
        tij <- sprintf("T_ij: %d/%d aligned", n_al, n_map)
      } else {
        tij <- sprintf("T_ij: n/a (%s)", d$transition_status %||% "NA")
      }
      n_cv <- sum(pt$cor_status == "valid", na.rm = TRUE)
      sprintf("%s ; rho: %d/%d valid (%s)", tij, n_cv, n_map, d$cor_status %||% "NA")
    })
    output$vel_table <- DT::renderDataTable({
      .ctx_tab_guard()
      ctx <- vel_ctx()
      req(ctx)
      dt <- ts_datatable(ctx$pair_table, page_length = 15L,
                         filename_base = "sc_comm_velocity_pair_table", filter = "none")
      num_cols <- intersect(
        c("mean_imported_score", "sender_velocity_mean", "sender_velocity_median",
          "receiver_velocity_mean", "receiver_velocity_median"),
        names(ctx$pair_table)
      )
      if (length(num_cols) > 0L) {
        dt <- DT::formatRound(dt, columns = num_cols, digits = 4)
      }
      dt
    })

    output$vel_plot <- renderPlot({
      .ctx_tab_guard()
      ctx <- vel_ctx()
      req(ctx)
      tryCatch(plot_communication_velocity_context(ctx), error = function(e) {
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = conditionMessage(e), size = 5, colour = "red") +
          ggplot2::theme_void()
      })
    })

    output$dl_vel_csv <- downloadHandler(
      filename = function() {
        ctx <- vel_ctx()
        if (is.null(ctx)) return("communication_velocity_context.csv")
        paste0("communication_velocity_context_", ctx$analysis_id, "_", ctx$timestamp_utc, ".csv")
      },
      content = function(file) {
        ctx <- vel_ctx()
        req(ctx)
        .guard_stale()
        utils::write.csv(build_communication_velocity_export(ctx), file, row.names = FALSE)
      }
    )

    .fig_download <- function(output_id, kind, ext) {
      output[[output_id]] <- downloadHandler(
        filename = function() {
          ctx <- vel_ctx()
          if (is.null(ctx)) return(paste0(kind, ".", ext))
          paste0(kind, "_", ctx$analysis_id, "_", ctx$timestamp_utc, ".", ext)
        },
        content = function(file) {
          ctx <- vel_ctx()
          req(ctx)
          .guard_stale()
          ts_export_plot(file, plot_communication_velocity_context(ctx),
                          width = 8, height = 6, dpi = 300)
        }
      )
    }
    .fig_download("dl_vel_png", "communication_velocity_context", "png")
    .fig_download("dl_vel_pdf", "communication_velocity_context", "pdf")

    invisible(NULL)
  })
}
