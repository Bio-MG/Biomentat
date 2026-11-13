# =============================================================================
# test-mod-bulk-i18n-quickwin.R — Quick-win I18N wiring (§2.1 & §2.2)
# =============================================================================

suppressPackageStartupMessages(library(shiny))

test_that("mod_bulk_network_ui uses reactive .tr() (data-key spans) for banner & threshold", {
  skip_if_not_installed("shiny.i18n")
  tr <- shiny.i18n::Translator$new(
    translation_json_path = file.path(ts_project_root(), "i18n/translation.json")
  )
  tr$set_translation_language("fr")
  shiny.i18n::usei18n(tr)

  e <- new.env(parent = globalenv())
  e$i18n <- tr
  e$I18N_AVAILABLE <- TRUE
  e$.strip_i18n_html <- function(x) {
    if (inherits(x, c("shiny.tag", "shiny.tag.list"))) {
      return(paste(as.character(x), collapse = ""))
    }
    as.character(x)[1L]
  }
  e$.tr_plain <- function(key) {
    res <- tryCatch(e$i18n$t(key), error = function(err) key)
    gsub("<[^>]+>", "", e$.strip_i18n_html(res))
  }
  e$.tr <- function(key) {
    tryCatch(e$i18n$t(key), error = function(err) key)
  }

  sys.source(file.path(ts_project_root(), "modules/bulk/mod_bulk_network.R"), envir = e)

  html <- as.character(e$mod_bulk_network_ui("net"))
  keys <- c(
    "Interactome dérivé des voies Reactome — ce n'est PAS un PPI.",
    "Une arête signifie « ces deux protéines participent à la même voie », pas « interaction physique démontrée ». Les nœuds gris (relais) sont des CO-MEMBRES DE VOIE prédits, pas des partenaires de liaison.",
    "Seuls les gènes dont la prime DÉPASSE strictement ce seuil deviennent des primes. 1,3 = p = 0,05 sur l'échelle -log10(padj)."
  )
  for (k in keys) {
    k_esc <- htmltools::htmlEscape(k, attribute = TRUE)
    expect_true(
      grepl(paste0('data-key="', k_esc, '"'), html, fixed = TRUE),
      info = paste("Expected shiny.i18n data-key span for:", k)
    )
  }

  tr$set_translation_language("en")
  expect_match(as.character(e$.tr(keys[1])), "Reactome pathways", fixed = TRUE)
  expect_match(as.character(e$.tr(keys[2])), "PREDICTED PATHWAY CO-MEMBERS", fixed = TRUE)
  expect_match(as.character(e$.tr(keys[3])), "STRICTLY EXCEEDS", fixed = TRUE)
})

test_that("mod_bulk_merge_server updates merge_datasets & merge_condition placeholders on language switch", {
  skip_if_not_installed("shiny.i18n")
  tr <- shiny.i18n::Translator$new(
    translation_json_path = file.path(ts_project_root(), "i18n/translation.json")
  )
  tr$set_translation_language("fr")
  shiny.i18n::usei18n(tr)

  captured <- list()
  e <- new.env(parent = globalenv())
  e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
  e$i18n <- tr
  e$.strip_i18n_html <- function(x) {
    txt <- paste(as.character(x), collapse = "")
    gsub("<[^>]+>", "", txt)
  }
  e$updateSelectizeInput <- function(session, inputId, label = NULL, choices = NULL,
                                     selected = NULL, options = list(), server = FALSE) {
    captured[[length(captured) + 1L]] <<- list(
      id = inputId, label = label, options = options
    )
    invisible(NULL)
  }
  e$updateCheckboxInput <- function(...) invisible(NULL)
  e$updateTextInput <- function(...) invisible(NULL)
  e$updateActionButton <- function(...) invisible(NULL)
  e$renderDT <- function(...) function(...) NULL
  e$DTOutput <- function(...) NULL
  e$navset_card_underline <- function(...) NULL
  e$nav_panel <- function(...) NULL

  sys.source(file.path(ts_project_root(), "modules/bulk/mod_bulk_merge.R"), envir = e)

  gd <- shiny::reactiveValues(
    language = "fr",
    i18n = tr,
    bulk_datasets = list()
  )

  shiny::testServer(e$mod_bulk_merge_server, args = list(global_data = gd), {
    captured <<- list()
    tr$set_translation_language("en")
    gd$language <- "en"
    session$flushReact()
  })

  ds_calls <- Filter(function(x) identical(x$id, "merge_datasets"), captured)
  cond_calls <- Filter(function(x) identical(x$id, "merge_condition"), captured)

  expect_true(length(ds_calls) >= 1L)
  expect_true(length(cond_calls) >= 1L)

  last_ds <- ds_calls[[length(ds_calls)]]
  last_cond <- cond_calls[[length(cond_calls)]]

  expect_identical(last_ds$options$placeholder, "Select datasets...")
  expect_identical(last_cond$options$placeholder, "None — blind correction if ComBat")
})
