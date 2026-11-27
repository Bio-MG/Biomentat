# =============================================================================
# mod_bulk_de_summary.R — Bulk Child 2: Resume Up/Down tab
# (Step-3.6 refactor — extracted from the monolithic mod_bulk_de.R)
# =============================================================================
# Barchart + table of significant Up/Down genes across ALL computed
# contrasts. Reuses the same helper as the HTML/PDF report
# (summarize_contrasts_updown(), helpers_bulk.R) — single source of truth —
# and recomputes LIVE on input$lfc_thresh/padj_thresh, with no DE recalculation.
# =============================================================================

.de_summary_server <- function(input, output, session, ns, global_data, shared_rv) {

  .tr <- function(key) {
    tr <- global_data$i18n
    if (is.null(tr)) return(key)
    tryCatch(.strip_i18n_html(tr$t(key)), error = function(e) key)
  }


  updown_summary <- reactive({
    req(length(shared_rv$contrasts) > 0)
    summarize_contrasts_updown(shared_rv$contrasts, lfc_thresh = input$lfc_thresh,
                               padj_thresh = input$padj_thresh,
                               active_contrast = shared_rv$active_contrast)
  })

  updown_plot <- reactive({
    global_data$language                     # i18n trigger
    plot_updown_barchart(updown_summary(), tr = .tr_fn(global_data),
                         palette = shared_rv$bulk_palette %||% "default",
                         manual_colors = shared_rv$volcano_role_colors)
  })

  output$plot_updown <- renderPlot({
    validate(need(length(shared_rv$contrasts) > 0,
                  .tr("Aucun contraste calculé — lancez l'Étape 2 (DE) d'abord.")))
    updown_plot()
  })

  output$dl_updown_png <- downloadHandler(
    filename = function() paste0("updown_summary_", Sys.Date(), ".png"),
    content  = function(file) ts_export_plot(file, updown_plot(), width = 8, height = 5.5, dpi = 300)
  )

  output$table_updown <- renderDT({
    global_data$language  # i18n
    # R-3bis (fix 2026-11-06) : lecture du miroir d'onglet (mod_bulk.R) —
    # ré-exécute ce render au moment où l'onglet devient VISIBLE, ce qui
    # débloque la DT « empoisonnée » par une invalidation pendant cachée
    # (mécanisme mesuré en app minimale : Shiny 1.14 + DT 0.34 + bslib).
    force(shared_rv$main_tabs_current)
    req(updown_summary())
    df <- updown_summary()
    df$Actif <- ifelse(df$actif, .tr("\u2192 actif"), "")
    df$actif <- NULL
    colnames(df) <- c(.tr("Contraste"), .tr("G\u00e8nes test\u00e9s"), .tr("Significatifs"),
                      "Up", "Down", .tr("Actif"))
    ts_datatable(df, page_length = 15L, filename_base = "bulk_de_summary",
                 filter = "none")
  })

  output$dl_updown_csv <- downloadHandler(
    filename = function() paste0("updown_summary_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(updown_summary(), file, row.names = FALSE)
  )

  # R-3 (fix 2026-11-06, clan BULK-7 « jamais rendu ») : mesuré live (probe
  # DOM), la DT Up/Down n'a JAMAIS reçu son payload (aucun script[data-for])
  # quand son invalidation a eu lieu dans un onglet inactif — la reprise de
  # suspension ne la re-flushe pas, alors que le plot du même reactive rend.
  # Idiome Shiny standard : désactiver la suspension de CETTE sortie — elle
  # calcule et envoie son payload même cachée ; le binding DT dessine (et
  # redessine à l'activation via les hooks de redimensionnement bslib).
  outputOptions(output, "table_updown", suspendWhenHidden = FALSE)
}


# ── Server: DTU (Differential Transcript Usage via bulk_dtu.R) ────────────────
.de_dtu_server <- function(input, output, session, ns, global_data, shared_rv) {
  .tr <- function(key) {
    tr <- global_data$i18n
    if (is.null(tr)) return(key)
    tryCatch(.strip_i18n_html(tr$t(key)), error = function(e) key)
  }

  dtu_res_rv <- reactiveVal(NULL)
  tx2gene_preview_rv <- reactiveVal(NULL)

  observeEvent(input$btn_build_tx2gene, {
    global_data$language
    counts <- shared_rv$filtered_counts %||% global_data$bulk_obj$counts
    if (is.null(counts)) {
      showNotification(.tr("Chargez des comptages de transcrits avant de g\u00e9n\u00e9rer tx2gene."),
                       type = "error", duration = 8)
      return()
    }
    src <- input$dtu_tx2gene_source %||% "auto"
    gtf_p <- if (!is.null(input$dtu_gtf_file)) input$dtu_gtf_file$datapath else NULL
    res <- tryCatch(
      build_tx2gene(counts = counts, source = src, gtf_path = gtf_p),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      showNotification(sprintf(.tr("G\u00e9n\u00e9ration tx2gene impossible : %s"), conditionMessage(res)),
                       type = "error", duration = 10)
      return()
    }
    global_data$bulk_obj$tx2gene <- res$table
    tx2gene_preview_rv(res)
    showNotification(
      sprintf(.tr("\u2713 Table tx2gene g\u00e9n\u00e9r\u00e9e (%d / %d transcrits mapp\u00e9s)."),
              res$n_mapped, res$n_input),
      type = "message", duration = 5
    )
  })

  output$dtu_tx2gene_preview <- renderText({
    global_data$language
    p <- tx2gene_preview_rv()
    if (is.null(p)) return("")
    sprintf(.tr("tx2gene (%s) : %d / %d transcrits mapp\u00e9s (%d non mapp\u00e9s)"),
            p$source, p$n_mapped, p$n_input, p$n_unmapped)
  })

  observeEvent(input$run_dtu, {
    global_data$language
    counts <- shared_rv$filtered_counts %||% global_data$bulk_obj$counts
    meta   <- global_data$bulk_obj$metadata
    cond   <- input$condition_col %||% shared_rv$active_condition_col
    if (is.null(counts) || is.null(meta) || is.null(cond) || !cond %in% colnames(meta)) {
      showNotification(.tr("Chargez des comptages de transcrits et s\u00e9lectionnez une colonne de condition avant de lancer DTU."),
                       type = "error", duration = 8)
      return()
    }
    tx2gene <- global_data$bulk_obj$tx2gene
    f_info <- input$dtu_tx2gene_file
    if (!is.null(f_info) && !is.null(f_info$datapath) && file.exists(f_info$datapath)) {
      sep <- if (grepl("\\.(tsv|txt)$", f_info$name %||% "", ignore.case = TRUE)) "\t" else ","
      tx2gene <- tryCatch(
        utils::read.delim(f_info$datapath, sep = sep, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) NULL
      )
    }
    if (is.null(tx2gene) || !is.data.frame(tx2gene)) {
      showNotification(.tr("Fournissez une table tx2gene (CSV/TSV avec colonnes transcript_id et gene_id)."),
                       type = "error", duration = 8)
      return()
    }
    groups <- factor(as.character(meta[[cond]]))
    res <- tryCatch(
      run_bulk_dtu(
        counts = counts,
        tx2gene = tx2gene,
        design_or_groups = groups,
        test_method = input$dtu_test_method %||% "Simes"
      ),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      showNotification(sprintf(.tr("Analyse DTU impossible : %s"), conditionMessage(res)),
                       type = "error", duration = 10)
      return()
    }
    dtu_res_rv(res)
    showNotification(
      sprintf(.tr("\u2713 DTU termin\u00e9 (%d g\u00e8nes, %d transcrits test\u00e9s)."),
              res$n_genes_tested, res$n_transcripts_tested),
      type = "message", duration = 5
    )
  })

  output$dtu_status <- renderText({
    global_data$language
    r <- dtu_res_rv()
    if (is.null(r)) {
      return(.tr("En attente \u2014 importez une table tx2gene et cliquez sur \u00ab Lancer l'analyse DTU \u00bb."))
    }
    n_sig <- if ("FDR" %in% colnames(r$gene_table)) sum(r$gene_table$FDR < 0.05, na.rm = TRUE) else 0L
    sprintf(.tr("DTU (%s) : %d g\u00e8nes test\u00e9s (%d FDR < 0.05) | %d transcrits test\u00e9s"),
            r$test_method, r$n_genes_tested, n_sig, r$n_transcripts_tested)
  })

  output$dtu_gene_table <- renderDT({
    r <- req(dtu_res_rv())
    ts_datatable(r$gene_table, page_length = 15L, filename_base = "bulk_dtu_genes")
  })

  output$dtu_tx_table <- renderDT({
    r <- req(dtu_res_rv())
    ts_datatable(r$transcript_table, page_length = 15L, filename_base = "bulk_dtu_transcripts")
  })
}
