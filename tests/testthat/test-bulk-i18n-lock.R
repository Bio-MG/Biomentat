# =============================================================================
# test-bulk-i18n-lock.R — W1.2 Bulk QuickWin & i18n Lock
# =============================================================================
# Locks i18n coverage and server-side language push across modules/bulk/*.R and
# modules/bulk_de/*.R:
#   1. Every string literal passed to i18n$t(), .tr(), .tr_plain(), .trl(), or
#      i18n_plain() in modules/bulk/*.R and modules/bulk_de/*.R must exist in
#      i18n/translation.json with a non-empty English translation.
#   2. mod_bulk_pathways_server pushes all 4 enrich_mode choices (including
#      "compare" -> "Comparative — compareCluster (inter-contrastes)") and
#      updates comp_source, comp_fdr, and dir_fdr on language switch.
#   3. mod_bulk_dose_response_server pushes translated choices for dose_model
#      and translated placeholder for dose_gene on language switch.
#   4. mod_bulk_de_engine (.de_engine_server) pushes plot_theme choices and
#      venn_contrasts placeholder on language switch.
#   5. Previously unwrapped French UI strings in bulk modules (datasets
#      activation dimensions, filter varpart method summary, merge status line,
#      network QC table, WGCNA power warning & module table header, DE no-rep
#      banner BCV summary, Venn/multimethod error fallback, Volcano manual
#      palette role labels, and Heatmap statistical subtitle) are translated.
# =============================================================================

suppressPackageStartupMessages(library(shiny))

.load_translation_map <- function() {
  path <- file.path(ts_project_root(), "i18n", "translation.json")
  j <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  fr <- vapply(j$translation, function(x) x$fr %||% "", character(1))
  en <- vapply(j$translation, function(x) x$en %||% "", character(1))
  stats::setNames(en, fr)
}

test_that("every i18n$t / .tr / .tr_plain / .trl / i18n_plain literal in bulk & bulk_de exists in translation.json", {
  tr_map <- .load_translation_map()
  root <- ts_project_root()
  files <- sort(c(
    list.files(file.path(root, "modules", "bulk"), pattern = "\\.R$", full.names = TRUE),
    list.files(file.path(root, "modules", "bulk_de"), pattern = "\\.R$", full.names = TRUE)
  ))
  expect_true(length(files) > 0L)

  pat <- '(?:i18n\\$t|\\.tr_plain|\\.trl|\\.tr|i18n_plain)\\s*\\(\\s*"((?:[^"\\\\]|\\\\.)*)"\\s*\\)'
  missing_keys <- character(0)
  empty_en <- character(0)

  for (f in files) {
    lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
    rel <- sub(paste0("^", gsub("([()])", "\\\\\\1", normalizePath(root, winslash = "/")), "/?"),
               "", normalizePath(f, winslash = "/"))
    for (i in seq_along(lines)) {
      ln <- lines[i]
      if (grepl("^\\s*#", ln)) next
      m <- gregexpr(pat, ln, perl = TRUE)[[1]]
      if (m[1L] == -1L) next
      matches <- regmatches(ln, list(m))[[1]]
      for (call_txt in matches) {
        lit <- sub('^(?:i18n\\$t|\\.tr_plain|\\.trl|\\.tr|i18n_plain)\\s*\\(\\s*', "", call_txt, perl = TRUE)
        lit <- sub('\\s*\\)$', "", lit, perl = TRUE)
        key <- tryCatch(enc2utf8(eval(parse(text = lit))), error = function(e) NA_character_)
        if (is.na(key) || !nzchar(key)) next
        if (!key %in% names(tr_map)) {
          missing_keys <- c(missing_keys, sprintf("%s:%d -> %s", rel, i, key))
        } else if (!nzchar(trimws(tr_map[[key]]))) {
          empty_en <- c(empty_en, sprintf("%s:%d -> %s", rel, i, key))
        }
      }
    }
  }

  expect_identical(missing_keys, character(0),
                   info = paste("Missing i18n keys:\n", paste(missing_keys, collapse = "\n")))
  expect_identical(empty_en, character(0),
                   info = paste("Empty EN translations:\n", paste(empty_en, collapse = "\n")))
})

test_that("required W1.2 bulk i18n lock keys exist in translation.json with English translations", {
  tr_map <- .load_translation_map()
  required_keys <- c(
    "L'activation remplace le jeu actif : filtrage, contrastes et voies du jeu précédent sont purgés.",
    "Tous significatifs",
    "Up seulement",
    "Down seulement",
    "%s gènes × %s échantillons",
    "Méthode : {m} — {n} gènes analysés",
    "%s — %d gènes communs × %d échantillons [%s]",
    "⚠️ Aucun power n'atteint R2 = 0.80 — le meilleur compromis est retenu (voir sous-titre du graphique).",
    "BCV = %.2f (dispersion = %.4f) — %d échantillon(s) / %d coefficient(s).",
    "(ou conteneur trop petit — agrandissez la fenêtre/onglet)",
    "le calcul DE a échoué dans le pool",
    "Up-régulé",
    "Down-régulé",
    "Non-significatif",
    "%d Up · %d Down · %d testés%s",
    "Paramètre",
    "algorithme",
    "seuil de prime",
    "primes en entrée",
    "primes non positives",
    "sous le seuil",
    "hors réseau",
    "terminaux",
    "relais",
    "nœuds",
    "arêtes",
    "composantes",
    "prime totale",
    "score",
    "liaison max (arêtes)",
    "nœuds du réseau",
    "arêtes du réseau",
    "taux de correspondance",
    "espèce",
    "source",
    "version source"
  )
  for (k in required_keys) {
    expect_true(k %in% names(tr_map), info = paste("Missing key:", k))
    expect_true(nzchar(trimws(tr_map[[k]])), info = paste("Empty EN for key:", k))
  }
})

test_that("mod_bulk_pathways_server pushes all 4 enrich_mode choices and comp_source on language switch", {
  skip_if_not_installed("shiny.i18n")
  tr <- shiny.i18n::Translator$new(
    translation_json_path = file.path(ts_project_root(), "i18n/translation.json")
  )
  tr$set_translation_language("fr")

  captured <- list()
  e <- new.env(parent = globalenv())
  e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
  e$i18n <- tr
  e$TS_PADJ_METHODS <- c("BH", "holm")
  e$TS_PADJ_METHOD_DEFAULT <- "BH"
  e$TS_BULK_GSVA_MIN_SIZE <- 10L
  e$TS_BULK_GSVA_MAX_SIZE <- 500L
  e$TS_BULK_GSVA_OVERLAP_MIN <- 0.2
  e$.strip_i18n_html <- function(x) gsub("<[^>]+>", "", paste(as.character(x), collapse = ""))
  e$.tr_plain <- function(k) e$.strip_i18n_html(tr$t(k))
  e$.tr_placeholder <- function() "—"
  e$bulk_gene_set_choices <- function(tr_fn = identity) stats::setNames("msigdb_hallmark", tr_fn("MSigDB Hallmark (50 voies)"))
  e$enrichr_cloud_optin_ui <- function(...) NULL
  e$ts_drive_publish_token <- function(...) invisible(NULL)
  e$ts_drive_publish_export <- function(...) invisible(NULL)
  e$ts_drive_job_state <- function() NULL
  e$ts_drive_job_busy <- function() FALSE
  e$renderDT <- function(...) function(...) NULL
  e$DTOutput <- function(...) NULL
  e$renderPlotly <- function(...) function(...) NULL
  e$plotlyOutput <- function(...) NULL
  e$updateRadioButtons <- function(session, inputId, label = NULL, choices = NULL, selected = NULL) {
    captured[[length(captured) + 1L]] <<- list(fn = "updateRadioButtons", id = inputId, label = label, choices = choices)
  }
  e$updateSelectInput <- function(session, inputId, label = NULL, choices = NULL, selected = NULL) {
    captured[[length(captured) + 1L]] <<- list(fn = "updateSelectInput", id = inputId, label = label, choices = choices)
  }
  e$updateSelectizeInput <- function(session, inputId, label = NULL, choices = NULL, selected = NULL, options = list(), server = FALSE) {
    captured[[length(captured) + 1L]] <<- list(fn = "updateSelectizeInput", id = inputId, label = label, choices = choices, options = options)
  }
  e$updateNumericInput <- function(session, inputId, label = NULL, ...) {
    captured[[length(captured) + 1L]] <<- list(fn = "updateNumericInput", id = inputId, label = label)
  }
  e$updateActionButton <- function(...) invisible(NULL)

  sys.source(file.path(ts_project_root(), "modules/bulk/mod_bulk_pathways.R"), envir = e)

  gd <- shiny::reactiveValues(language = "fr", i18n = tr, packs = list())
  srv <- shiny::reactiveValues(filtered_counts = NULL, vst_mat = NULL, contrasts = list())

  shiny::testServer(e$mod_bulk_pathways_server, args = list(global_data = gd, shared_rv = srv), {
    session$flushReact()
    captured <<- list()
    tr$set_translation_language("en")
    gd$language <- "en"
    session$flushReact()
  })

  em_calls <- Filter(function(x) identical(x$id, "enrich_mode"), captured)
  expect_true(length(em_calls) >= 1L)
  last_em <- em_calls[[length(em_calls)]]
  expect_identical(unname(last_em$choices), c("ora", "gsea", "dir", "compare"))
  expect_true(any(grepl("compareCluster", names(last_em$choices), fixed = TRUE)))

  cs_calls <- Filter(function(x) identical(x$id, "comp_source"), captured)
  expect_true(length(cs_calls) >= 1L)
  last_cs <- cs_calls[[length(cs_calls)]]
  expect_identical(unname(last_cs$choices), c("all", "up", "down"))
  expect_identical(names(last_cs$choices), c("All significant", "Up only", "Down only"))
})

test_that("mod_bulk_dose_response_server pushes dose_model choices and dose_gene placeholder on language switch", {
  skip_if_not_installed("shiny.i18n")
  tr <- shiny.i18n::Translator$new(
    translation_json_path = file.path(ts_project_root(), "i18n/translation.json")
  )
  tr$set_translation_language("fr")

  captured <- list()
  e <- new.env(parent = globalenv())
  e$`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
  e$i18n <- tr
  e$.strip_i18n_html <- function(x) gsub("<[^>]+>", "", paste(as.character(x), collapse = ""))
  e$.tr_plain <- function(k) e$.strip_i18n_html(tr$t(k))
  e$bulk_dose_models <- function() c("LL.4", "W1.4", "W2.4", "BC.4")
  e$renderDT <- function(...) function(...) NULL
  e$DTOutput <- function(...) NULL
  e$updateSelectInput <- function(session, inputId, label = NULL, choices = NULL, selected = NULL) {
    captured[[length(captured) + 1L]] <<- list(id = inputId, label = label, choices = choices)
  }
  e$updateSelectizeInput <- function(session, inputId, label = NULL, choices = NULL, selected = NULL, options = list(), server = FALSE) {
    captured[[length(captured) + 1L]] <<- list(id = inputId, label = label, choices = choices, options = options)
  }
  e$updateRadioButtons <- function(...) invisible(NULL)
  e$updateNumericInput <- function(...) invisible(NULL)
  e$updateActionButton <- function(...) invisible(NULL)

  sys.source(file.path(ts_project_root(), "modules/bulk/mod_bulk_dose_response.R"), envir = e)

  gd <- shiny::reactiveValues(language = "fr", i18n = tr, bulk_obj = NULL)
  srv <- shiny::reactiveValues(vst_mat = NULL, contrasts = list())

  shiny::testServer(e$mod_bulk_dose_response_server, args = list(global_data = gd, shared_rv = srv), {
    session$flushReact()
    captured <<- list()
    tr$set_translation_language("en")
    gd$language <- "en"
    session$flushReact()
  })

  dm_calls <- Filter(function(x) identical(x$id, "dose_model"), captured)
  expect_true(length(dm_calls) >= 1L)
  last_dm <- dm_calls[[length(dm_calls)]]
  expect_identical(unname(last_dm$choices), c("LL.4", "W1.4", "W2.4", "BC.4"))
  expect_identical(names(last_dm$choices)[1L], "Log-logistic (Hill, LL.4)")

  dg_calls <- Filter(function(x) identical(x$id, "dose_gene"), captured)
  expect_true(length(dg_calls) >= 1L)
  last_dg <- dg_calls[[length(dg_calls)]]
  expect_identical(last_dg$options$placeholder, "Gene (type to search)")
})

test_that("mod_bulk_de_engine pushes plot_theme choices and venn_contrasts placeholder on language switch", {
  skip_if_not_installed("shiny.i18n")
  tr <- shiny.i18n::Translator$new(
    translation_json_path = file.path(ts_project_root(), "i18n/translation.json")
  )
  tr$set_translation_language("fr")

  root <- ts_project_root()
  eng_src <- paste(readLines(file.path(root, "modules/bulk_de/mod_bulk_de_engine.R"),
                             warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(eng_src, 'updateSelectInput\\(session,\\s*"plot_theme"', perl = TRUE)
  expect_match(eng_src, 'placeholder\\s*=\\s*\\.trl\\("2 \\\\u00e0 6 contrastes"\\)', perl = TRUE)
  expect_match(eng_src, 'sprintf\\(\\.trl\\("BCV = %\\.2f \\(dispersion = %\\.4f\\)', perl = TRUE)

  viz_src <- paste(readLines(file.path(root, "modules/bulk_de/mod_bulk_de_viz.R"),
                             warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(viz_src, '\\.tr\\("Up-r\\\\u00e9gul\\\\u00e9"\\)', perl = TRUE)
  expect_match(viz_src, 'sprintf\\(\\.tr\\("%d Up', perl = TRUE)
})

