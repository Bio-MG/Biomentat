# =============================================================================
# test-mod-bulk-network-i18n.R — Reactive i18n choices for mod_bulk_network
# =============================================================================

suppressPackageStartupMessages(library(shiny))

.make_bulk_network_i18n_env <- function(tr) {
  captured <- list()
  e <- new.env(parent = globalenv())
  e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
  e$i18n <- tr
  e$.strip_i18n_html <- function(x) {
    txt <- paste(as.character(x), collapse = "")
    gsub("<[^>]+>", "", txt)
  }
  e$.t_fmt <- function(template, ...) {
    vals <- list(...)
    for (nm in names(vals)) {
      template <- gsub(paste0("{", nm, "}"), format(vals[[nm]]), template, fixed = TRUE)
    }
    template
  }
  record_call <- function(fn, inputId, label = NULL, choices = NULL, selected = NULL) {
    captured[[length(captured) + 1L]] <<- list(
      fn = fn, id = inputId, label = label, choices = choices, selected = selected
    )
    invisible(NULL)
  }
  e$updateSelectInput <- function(session, inputId, label = NULL, choices = NULL,
                                  selected = NULL) {
    record_call("updateSelectInput", inputId, label, choices, selected)
  }
  e$updateRadioButtons <- function(session, inputId, label = NULL, choices = NULL,
                                   selected = NULL, inline = FALSE,
                                   choiceNames = NULL, choiceValues = NULL) {
    record_call("updateRadioButtons", inputId, label, choices, selected)
  }
  e$updateNumericInput <- function(...) invisible(NULL)
  e$updateActionButton <- function(...) invisible(NULL)
  e$renderDT <- function(...) function(...) NULL
  e$DTOutput <- function(...) NULL
  e$bulk_network_pcsf_params <- function(params) list(max_join_hops = 5L)

  sys.source(file.path(ts_project_root(), "modules/bulk/mod_bulk_network.R"), envir = e)
  list(env = e, get_captured = function() captured, clear = function() captured <<- list())
}

test_that("mod_bulk_network_server updates network_species, network_source, network_prize on init and language switch while preserving selected", {
  skip_if_not_installed("shiny.i18n")
  tr <- shiny.i18n::Translator$new(
    translation_json_path = file.path(ts_project_root(), "i18n/translation.json")
  )
  tr$set_translation_language("en")
  shiny.i18n::usei18n(tr)

  harness <- .make_bulk_network_i18n_env(tr)
  e <- harness$env

  gd <- shiny::reactiveValues(language = "en", i18n = tr)
  srv <- shiny::reactiveValues()

  shiny::testServer(
    e$mod_bulk_network_server,
    args = list(global_data = gd, shared_rv = srv),
    {
      session$flushReact()
      init_calls <- harness$get_captured()
      expect_true(length(init_calls) >= 3L)

      sp_init <- Filter(function(x) identical(x$id, "network_species"), init_calls)
      src_init <- Filter(function(x) identical(x$id, "network_source"), init_calls)
      prz_init <- Filter(function(x) identical(x$id, "network_prize"), init_calls)

      expect_true(length(sp_init) >= 1L)
      expect_true(length(src_init) >= 1L)
      expect_true(length(prz_init) >= 1L)

      last_sp_en <- sp_init[[length(sp_init)]]
      last_src_en <- src_init[[length(src_init)]]
      last_prz_en <- prz_init[[length(prz_init)]]

      expect_identical(unname(last_sp_en$choices), "hsapiens")
      expect_false(identical(names(last_sp_en$choices), "Humain (hsapiens)"))
      expect_identical(last_sp_en$selected, "hsapiens")

      expect_identical(unname(last_src_en$choices), c("all_sig", "up", "down"))
      expect_false(identical(names(last_src_en$choices),
                             c("Tous les significatifs", "Surexprimés (up)", "Sous-exprimés (down)")))
      expect_identical(last_src_en$selected, "all_sig")

      expect_identical(unname(last_prz_en$choices), c("padj", "lfc"))
      expect_false(identical(names(last_prz_en$choices),
                             c("-log10(padj) — recommandé", "|log2FC| — direction ignorée")))
      expect_identical(last_prz_en$selected, "padj")

      # Change user selections, then switch language to French and verify preservation
      session$setInputs(
        network_species = "hsapiens",
        network_source = "up",
        network_prize = "lfc"
      )
      harness$clear()
      tr$set_translation_language("fr")
      gd$language <- "fr"
      session$flushReact()

      fr_calls <- harness$get_captured()
      sp_fr <- Filter(function(x) identical(x$id, "network_species"), fr_calls)
      src_fr <- Filter(function(x) identical(x$id, "network_source"), fr_calls)
      prz_fr <- Filter(function(x) identical(x$id, "network_prize"), fr_calls)

      expect_identical(names(sp_fr[[length(sp_fr)]]$choices), "Humain (hsapiens)")
      expect_identical(sp_fr[[length(sp_fr)]]$selected, "hsapiens")

      expect_identical(names(src_fr[[length(src_fr)]]$choices),
                       c("Tous les significatifs", "Surexprimés (up)", "Sous-exprimés (down)"))
      expect_identical(src_fr[[length(src_fr)]]$selected, "up")

      expect_identical(names(prz_fr[[length(prz_fr)]]$choices),
                       c("-log10(padj) — recommandé", "|log2FC| — direction ignorée"))
      expect_identical(prz_fr[[length(prz_fr)]]$selected, "lfc")
    }
  )
})
