# =============================================================================
# mod_sc_da_cross.R — Abondance differentielle 4E-3 : vues croisées Milo x
# scCODA (Stage 16)
# =============================================================================
# Le module ORCHESTRE uniquement : il consomme les DEUX resultats canoniques
# exposes par les panels 8d (shared_rv$da_milo_result) et 8e
# (shared_rv$da_sccoda_result) et les vues pures de
# R/sc/sc_abundance_cross_views.R. AUCUN calcul : la comparaison est
# DESCRIPTIVE, sensible a la methode, sans p-value de consensus. Les options
# d'affichage sont enregistrees dans la provenance au moment de l'export.
# =============================================================================

mod_sc_da_cross_ui <- function(id) {
  ns <- NS(id)
  tagList(
    div(class = "alert alert-light",
        style = "font-size:0.9em;border-left:3px solid #2C3E50;",
        i18n$t("Vues croisées DESCRIPTIVES : Milo teste des VOISINAGES (régions d'un graphe), scCODA teste la COMPOSITION par échantillon — deux grandeurs différentes, jamais fusionnées. "),
        i18n$t("Aucune « p-value de consensus » n'est calculée.")),
    uiOutput(ns("cross_comparability_box")),
    actionButton(ns("cross_refresh"), i18n$t("Rafraîchir la comparaison"),
                 class = "btn-primary w-100", icon = icon("scale-balanced")),
    hr(),
    downloadButton(ns("dl_cross_concordance"), i18n$t("Exporter concordance par identité (CSV)"), class = "btn-sm btn-info w-100 mt-1")
  )
}

mod_sc_da_cross_output_ui <- function(id) {
  ns <- NS(id)
  # CCC_NAV_DRAWER vague B jalon J4 (2026-11-12) : généralisation du pattern
  # tiroir optimisé du CCC (5f0a906 + bc9c09a) au panneau Vues croisées (8f)
  # — candidat C4 du packet gate 2 (GO opérateur vague B après vérité-terrain
  # live ; J4 AVANT J3 : mod_sc_velocity.R EN VOL lane scvelo). JUMELLE : la
  # carte .cross_drawer_groups est redéclarée EN LOCAL dans
  # mod_sc_da_cross_server (gel test-da-cross-views-contract-freeze.R —
  # aucune constante top-level) ; égalité verrouillée par
  # test-mod-sc-da-cross-ui.R.
  .cross_drawer_groups <- list(
    "Vues croisées" = list(
      list(value = "cross_tab_concordance", label = "Concordance par identité"),
      list(value = "cross_tab_map", label = "Mappage voisinages → identités"),
      list(value = "cross_tab_compo", label = "Composition échantillons"),
      list(value = "cross_tab_effects", label = "Effets scCODA (rappel)")),
    "Limites" = list(
      list(value = "cross_tab_limites", label = "Limites")))
  card(
    full_screen = TRUE,
    # Pastille ouvrant le tiroir (motif 2L CCC : dans l'en-tête de carte) :
    card_header(
      div(class = "d-flex justify-content-between align-items-center gap-2",
          tags$span(i18n$t("Abondance différentielle — vues croisées (Milo × scCODA)")),
          tags$button(type = "button",
                      class = "btn btn-sm btn-outline-secondary flex-shrink-0",
                      `data-bs-toggle` = "offcanvas",
                      `data-bs-target` = paste0("#", ns("cross_drawer")),
                      `aria-controls` = ns("cross_drawer"),
                      icon("list-ul"), " ", i18n$t("Navigation & filtres")))
    ),
    navset_tab(
      id = ns("cross_out_tabs"),
      nav_panel(i18n$t("Concordance par identité"), value = "cross_tab_concordance",
                plotOutput(ns("cross_concordance"), height = "480px"),
                div(class = "m-1",
                    downloadButton(ns("dl_cross_concordance_png"), "PNG", class = "btn-sm btn-outline-secondary"),
                    downloadButton(ns("dl_cross_concordance_pdf"), "PDF", class = "btn-sm btn-outline-secondary")),
                DT::dataTableOutput(ns("cross_table"), height = "300px")),
      nav_panel(i18n$t("Mappage voisinages → identités"), value = "cross_tab_map",
                plotOutput(ns("cross_nhood_map"), height = "480px"),
                div(class = "m-1",
                    downloadButton(ns("dl_cross_nhood_png"), "PNG", class = "btn-sm btn-outline-secondary"),
                    downloadButton(ns("dl_cross_nhood_pdf"), "PDF", class = "btn-sm btn-outline-secondary"))),
      nav_panel(i18n$t("Composition échantillons"), value = "cross_tab_compo",
                plotOutput(ns("cross_compo"), height = "480px"),
                div(class = "m-1",
                    downloadButton(ns("dl_cross_compo_png"), "PNG", class = "btn-sm btn-outline-secondary"),
                    downloadButton(ns("dl_cross_compo_pdf"), "PDF", class = "btn-sm btn-outline-secondary"))),
      nav_panel(i18n$t("Effets scCODA (rappel)"), value = "cross_tab_effects",
                plotOutput(ns("cross_sccoda_effects"), height = "480px")),
      nav_panel(i18n$t("Limites"), value = "cross_tab_limites",
                div(class = "alert alert-danger m-3", style = "font-size:0.9em;",
                    i18n$t("Deux GRANDEURS DIFFÉRENTES : la médiane des logFC de voisinages (Milo, niveau voisinage) et l'effet compositionnel (scCODA, niveau échantillon) ne mesurent pas la même chose — la concordance est descriptive et exploratoire.")),
                div(class = "alert alert-warning m-3", style = "font-size:0.9em;",
                    i18n$t("Les cellules ne sont pas des réplicats biologiques : chaque vue porte le contexte de composition par échantillon. Un désaccord peut provenir des échelles, des effectifs ou d'une biologie captée différemment — rien ne tranche.")),
                div(class = "alert alert-light m-3", style = "font-size:0.9em;",
                    uiOutput(ns("cross_meanings"))))
    ),
    # Tiroir offcanvas (DEUXIÈME enfant de la carte — motif comm_drawer /
    # milo_drawer / sccoda_drawer) :
    tags$div(
      id = ns("cross_drawer"), class = "offcanvas offcanvas-end", tabindex = "-1",
      `aria-labelledby` = ns("cross_drawer_title"),
      tags$div(class = "offcanvas-header",
               h5(id = ns("cross_drawer_title"), class = "offcanvas-title",
                  icon("list-ul"), " ", i18n$t("Navigation Vues croisées")),
               # .tr_plain (strip span) — JAMAIS de $t() brut en contexte
               # attribut (BUG A11Y « Fermer », motif bc9c09a) :
               tags$button(type = "button", class = "btn-close",
                           `data-bs-dismiss` = "offcanvas",
                           `aria-label` = .tr_plain("Fermer"))),
      tags$div(class = "offcanvas-body",
        h6(i18n$t("État du résultat")),
        div(class = "small text-muted", `aria-live` = "polite",
            uiOutput(ns("cross_drawer_summary"))),
        hr(),
        uiOutput(ns("cross_drawer_nav"))
      )
    )
  )
}

mod_sc_da_cross_server <- function(id, global_data, shared_rv = NULL) {
  moduleServer(id, function(input, output, session) {
    # ── i18n proxy ──────────────────────────────────────────────────────────
    .tr <- function(key) {
      tr <- isolate(global_data$i18n)
      if (is.null(tr)) return(key)
      tryCatch(.strip_i18n_html(tr$t(key)), error = function(e) key)
    }

    cross_state <- reactiveValues(summary = NULL, provenance = NULL)
    .has_both <- function() {
      !is.null(shared_rv) && !is.null(shared_rv$da_milo_result) &&
        !is.null(shared_rv$da_sccoda_result)
    }

    # ── CCC_NAV_DRAWER vague B jalon J4 (2026-11-12) : tiroir « Navigation » ──
    # Pattern 5f0a906 + bc9c09a, généralisé comme Milo (6c8eb84) et scCODA
    # (12f181c). JUMELLE de la déclaration output_ui (égalité verrouillée
    # test-mod-sc-da-cross-ui.R) — aucune constante top-level (gel
    # test-da-cross-views-contract-freeze.R).
    .cross_drawer_groups <- list(
      "Vues croisées" = list(
        list(value = "cross_tab_concordance", label = "Concordance par identité"),
        list(value = "cross_tab_map", label = "Mappage voisinages → identités"),
        list(value = "cross_tab_compo", label = "Composition échantillons"),
        list(value = "cross_tab_effects", label = "Effets scCODA (rappel)")),
      "Limites" = list(
        list(value = "cross_tab_limites", label = "Limites")))

    # Résumé d'état du tiroir : disponibilité des DEUX résultats canoniques
    # (consommés tels quels via shared_rv — jamais recalculés) + verdict de
    # COMPARABILITÉ du contrat (.cross_summary → fully_comparable). Même clé
    # d'attente que la boîte de comparabilité du panneau (zéro divergence).
    output$cross_drawer_summary <- renderUI({
      global_data$language  # re-rendu au changement de langue (motif du module)
      li <- function(n, lab) tags$li(tags$strong(n), " ", lab)
      milo_ok <- !is.null(shared_rv) && !is.null(shared_rv$da_milo_result)
      sccoda_ok <- !is.null(shared_rv) && !is.null(shared_rv$da_sccoda_result)
      if (!milo_ok && !sccoda_ok) {
        return(tags$ul(class = "list-unstyled mb-0",
                       tags$li(.tr("Lancez Milo (8d) ET scCODA (8e) pour activer la lecture croisée."))))
      }
      cmp <- if (.has_both() && !is.null(.cross_summary())) {
        .cross_summary()$comparability$fully_comparable
      } else NA
      cmp_lab <- if (isTRUE(cmp)) "\u2713" else
        if (isFALSE(cmp)) .tr("réserves") else "\u2014"
      tags$ul(class = "list-unstyled mb-0",
              style = "column-count:2;column-gap:1.5rem;",
              li(if (milo_ok) "\u2713" else "\u2014", .tr("Milo (8d)")),
              li(if (sccoda_ok) "\u2713" else "\u2014", .tr("scCODA (8e)")),
              li(cmp_lab, .tr("Comparabilité")))
    })

    # Navigation du tiroir rendue DYNAMIQUEMENT (motif P1+P2 comm_drawer) :
    output$cross_drawer_nav <- renderUI({
      global_data$language  # re-rendu au changement de langue (motif du module)
      sel <- input$cross_out_tabs %||% ""
      lapply(names(.cross_drawer_groups), function(grp) {
        items <- .cross_drawer_groups[[grp]]
        tags$div(class = "mb-3",
          h6(.tr(grp),
             tags$span(class = "badge text-bg-secondary ms-1", length(items))),
          lapply(items, function(pn) {
            active <- identical(sel, pn$value)
            actionButton(ns(paste0("cross_jump_", pn$value)),
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
    lapply(.cross_drawer_groups, function(grp) {
      lapply(grp, function(pn) {
        observeEvent(input[[paste0("cross_jump_", pn$value)]], {
          nav_select("cross_out_tabs", selected = pn$value)
        })
      })
    })

    # ── Comparabilite (pure lecture des deux resultats canoniques) ────────
    .cross_summary <- eventReactive(c(input$cross_refresh, .has_both()), {
      if (!.has_both()) return(NULL)
      tryCatch(
        build_da_cross_method_summary(shared_rv$da_milo_result,
                                      shared_rv$da_sccoda_result),
        error = function(e) {
          showNotification(paste(.tr("Erreur vues croisées :"), conditionMessage(e)),
                           type = "error", duration = 10)
          NULL
        }
      )
    })

    output$cross_comparability_box <- renderUI({
      sm <- .cross_summary()
      if (is.null(sm)) {
        return(div(class = "alert alert-light py-1", style = "font-size:0.85em;",
                   .tr("Lancez Milo (8d) ET scCODA (8e) pour activer la lecture croisée.")))
      }
      if (isTRUE(sm$comparability$fully_comparable)) {
        div(class = "alert alert-success py-1", style = "font-size:0.85em;",
            .tr("Résultats comparables : même objet, même contraste, même colonne d'identité."))
      } else {
        div(class = "alert alert-warning py-1", style = "font-size:0.85em;",
            tags$b(.tr("Comparabilité partielle — lisez les réserves :")),
            tags$ul(lapply(sm$comparability$caveats, function(cv) tags$li(cv))))
      }
    })

    # ── Vues (consommatrices pures des deux resultats) ─────────────────────
    .concordance_plot <- function() {
      plot_da_cross_concordance(shared_rv$da_milo_result,
                                shared_rv$da_sccoda_result)
    }
    .nhood_plot <- function() {
      plot_da_cross_nhood_mapping(shared_rv$da_milo_result)
    }
    .compo_plot <- function() {
      plot_da_cross_sample_composition(shared_rv$da_milo_result,
                                       shared_rv$da_sccoda_result)
    }
    .sccoda_effects_plot <- function() {
      plot_sccoda_effects(shared_rv$da_sccoda_result)
    }

    output$cross_concordance <- renderPlot({
      req(.has_both()); .cross_summary(); .concordance_plot()
    })
    output$cross_nhood_map <- renderPlot({
      req(.has_both()); .cross_summary(); .nhood_plot()
    })
    output$cross_compo <- renderPlot({
      req(.has_both()); .cross_summary(); .compo_plot()
    })
    output$cross_sccoda_effects <- renderPlot({
      req(.has_both()); .cross_summary(); .sccoda_effects_plot()
    })

    output$cross_table <- DT::renderDataTable({
      sm <- .cross_summary()
      req(sm)
      ts_datatable(sm$concordance, page_length = 15L,
                   filename_base = "sc_da_concordance", filter = "none") |>
        DT::formatSignif(columns = c("milo_frac_significant",
                                     "milo_median_logfc", "sccoda_effect",
                                     "sccoda_inclusion_probability"),
                         digits = 4)
    })

    output$cross_meanings <- renderUI({
      sm <- .cross_summary()
      req(sm)
      tags$ul(
        tags$li(sm$method_meanings$milo),
        tags$li(sm$method_meanings$sccoda),
        tags$li(sm$method_meanings$comparison)
      )
    })

    # ── Exports (chaque export PRODUIT une entree de provenance figee) ─────
    .export_figure <- function(plot_fn, file) {
      ts_export_plot(file, plot_fn(), width = 8, height = 6.5, dpi = 300)
    }
    .append_provenance <- function() {
      sm <- .cross_summary()
      if (is.null(sm)) return(NULL)
      entry <- build_da_cross_provenance(
        shared_rv$da_milo_result, shared_rv$da_sccoda_result,
        options = list(view = "da_cross_module")
      )
      provenance_append(shared_rv, entry)
      entry
    }

    output$dl_cross_concordance <- downloadHandler(
      filename = function() da_cross_export_filename("da_cross_concordance", "csv"),
      content = function(file) {
        req(.has_both())
        .append_provenance()
        utils::write.csv(build_da_cross_concordance_export(
          shared_rv$da_milo_result, shared_rv$da_sccoda_result
        ), file, row.names = FALSE)
      }
    )
    output$dl_cross_concordance_png <- downloadHandler(
      filename = function() da_cross_export_filename("da_cross_concordance", "png"),
      content = function(file) { req(.has_both()); .append_provenance(); .export_figure(.concordance_plot, file) }
    )
    output$dl_cross_concordance_pdf <- downloadHandler(
      filename = function() da_cross_export_filename("da_cross_concordance", "pdf"),
      content = function(file) { req(.has_both()); .append_provenance(); .export_figure(.concordance_plot, file) }
    )
    output$dl_cross_nhood_png <- downloadHandler(
      filename = function() da_cross_export_filename("da_cross_nhood_mapping", "png"),
      content = function(file) { req(.has_both()); .append_provenance(); .export_figure(.nhood_plot, file) }
    )
    output$dl_cross_nhood_pdf <- downloadHandler(
      filename = function() da_cross_export_filename("da_cross_nhood_mapping", "pdf"),
      content = function(file) { req(.has_both()); .append_provenance(); .export_figure(.nhood_plot, file) }
    )
    output$dl_cross_compo_png <- downloadHandler(
      filename = function() da_cross_export_filename("da_cross_sample_composition", "png"),
      content = function(file) { req(.has_both()); .append_provenance(); .export_figure(.compo_plot, file) }
    )
    output$dl_cross_compo_pdf <- downloadHandler(
      filename = function() da_cross_export_filename("da_cross_sample_composition", "pdf"),
      content = function(file) { req(.has_both()); .append_provenance(); .export_figure(.compo_plot, file) }
    )

    return(cross_state)
  })
}
