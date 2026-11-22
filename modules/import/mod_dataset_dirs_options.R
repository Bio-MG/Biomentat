# =============================================================================
# modules/import/mod_dataset_dirs_options.R — Saved dataset discovery folders
# =============================================================================
# Shared UI & server helpers for Bulk, Single-Cell, and Spatial import modules
# so users can save/remove custom dataset discovery folders in persistent
# options (`ts_save_dataset_dir` / `ts_remove_dataset_dir` / `ts_datasets_dirs`).
# =============================================================================

dataset_dirs_options_ui <- function(ns, i18n) {
  div(
    class = "border rounded p-2 mb-2 bg-light",
    div(
      class = "d-flex align-items-center gap-2 mb-1",
      div(
        style = "flex: 1 1 auto;",
        textInput(
          ns("custom_dataset_dir"),
          label = NULL,
          placeholder = .tr_plain("Chemin d'un dossier de datasets...")
        )
      ),
      # span i18n$t (pas .tr_plain) : le shim JS bascule le libellé du
      # bouton en live au changement de langue.
      actionButton(
        ns("btn_save_dataset_dir"),
        i18n$t("Enregistrer ce dossier"),
        class = "btn-outline-secondary btn-sm",
        icon = icon("bookmark")
      )
    ),
    uiOutput(ns("saved_dataset_dirs_ui"))
  )
}

dataset_dirs_options_server <- function(input, output, session,
                                        modality = c("sc", "bulk", "spatial"),
                                        on_changed = NULL,
                                        language = NULL) {
  modality <- match.arg(modality)
  ns <- session$ns

  rv_saved_dirs <- reactiveVal(
    if (exists("ts_saved_dataset_dirs", mode = "function")) {
      ts_saved_dataset_dirs(modality = modality)
    } else {
      character(0)
    }
  )

  observeEvent(input$btn_save_dataset_dir, {
    p <- trimws(as.character(input$custom_dataset_dir %||% ""))
    if (!nzchar(p)) {
      showNotification(.tr_plain("Veuillez saisir un chemin de dossier."), type = "warning")
      return()
    }
    if (exists("ts_save_dataset_dir", mode = "function")) {
      updated <- tryCatch(
        ts_save_dataset_dir(p, modality = modality),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error")
          NULL
        }
      )
      if (!is.null(updated)) {
        rv_saved_dirs(updated)
        updateTextInput(session, "custom_dataset_dir", value = "")
        showNotification(.tr_plain("Dossier enregistré dans les options."), type = "message")
        if (is.function(on_changed)) on_changed()
      }
    }
  })

  observeEvent(input$btn_remove_dataset_dir, {
    sel <- trimws(as.character(input$saved_dataset_dir_select %||% ""))
    if (!nzchar(sel)) return()
    if (exists("ts_remove_dataset_dir", mode = "function")) {
      updated <- ts_remove_dataset_dir(sel, modality = modality)
      rv_saved_dirs(updated)
      showNotification(.tr_plain("Dossier retiré des options."), type = "message")
      if (is.function(on_changed)) on_changed()
    }
  })

  output$saved_dataset_dirs_ui <- renderUI({
    # i18n : re-render au changement de langue (les textes ci-dessous passent
    # par .tr_plain, en-session = langue courante).
    if (!is.null(language)) language()
    dirs <- rv_saved_dirs()
    if (length(dirs) == 0L) {
      return(tags$small(class = "text-muted", .tr_plain("Aucun dossier personnalisé enregistré.")))
    }
    div(
      class = "d-flex align-items-center gap-2 mt-1",
      div(
        style = "flex: 1 1 auto;",
        selectInput(
          ns("saved_dataset_dir_select"),
          label = NULL,
          choices = dirs,
          selected = dirs[[1L]]
        )
      ),
      actionButton(
        ns("btn_remove_dataset_dir"),
        .tr_plain("Retirer"),
        class = "btn-outline-danger btn-sm",
        icon = icon("trash")
      )
    )
  })

  rv_saved_dirs

  # i18n : le placeholder de textInput est un attribut HTML — invisible au
  # shim JS. Re-poussé côté serveur au changement de langue (motif app.R
  # 1458) ; le libellé du bouton save, lui, bascule via le span i18n$t.
  if (!is.null(language)) {
    observeEvent(language(), {
      updateTextInput(session, "custom_dataset_dir",
                      placeholder = .tr_plain("Chemin d'un dossier de datasets..."))
    }, ignoreInit = TRUE)
  }
}
