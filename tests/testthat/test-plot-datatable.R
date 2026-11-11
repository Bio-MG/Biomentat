# =============================================================================
# test-plot-datatable.R — PLOT-S3 : wrapper DT harmonise ts_datatable()
# =============================================================================
# Objectif central (MISE A JOUR jalon DT-EXPORT 2026-09-14, contrat §6
# option B) :
#   - buttons : defaut = TS_DT_BUTTONS_DEFAULT (TRUE) — les tables de
#     resultats sont boutonnees ; les tables d'APERCU passent
#     buttons = FALSE explicitement (ni extensions ni dom alors poses)
#   - page_length est TOUJOURS OBLIGATOIRE dans la signature : la
#     normalisation a 15 se fait site par site, jamais par un defaut cache
#   - plus AUCUN appel direct a DT::datatable() dans modules/ et R/
#     (hors R/plotting/datatable.R lui-meme)
# =============================================================================
suppressWarnings(suppressPackageStartupMessages(library(DT)))

# D2 M4 (2026-10-25) : R/plotting vit desormais dans l'arbre du paquet
# (tools/pkg_skeleton/R/) — charge via load_project_pkg() (cf. helper-source.R).
load_project_pkg()

.tsd_df <- function() {
  data.frame(gene = c("A", "B", "C"), padj = c(0.0001, 0.02, 0.4),
             stringsAsFactors = FALSE)
}

# --- 1. Contrat du resolveur -------------------------------------------------

test_that("surface publique figee", {
  expect_setequal(
    ts_datatable_public_api(),
    c("ts_datatable", "ts_datatable_buttons",
      "ts_datatable_page_lengths", "ts_datatable_public_api",
      "ts_dt_future", "ts_dt_future_install", "ts_dt_future_restore")
  )
  expect_false(anyDuplicated(ts_datatable_public_api()) > 0)
})

test_that("page_length est OBLIGATOIRE — aucune valeur par defaut n'est neutre", {
  df <- .tsd_df()
  e <- tryCatch(ts_datatable(df), error = function(e) e)
  expect_s3_class(e, "plot_datatable_error")
  expect_identical(e$state, "invalid_page_length")
  expect_match(conditionMessage(e), "obligatoire")
  expect_match(conditionMessage(e), "10, 15, 8, 20 et 6")
  # il ne doit PAS y avoir de valeur par defaut dans la signature
  expect_true(is.name(formals(ts_datatable)$page_length) ||
              is.symbol(formals(ts_datatable)$page_length))
})

test_that("page_length invalide => invalid_page_length", {
  df <- .tsd_df()
  for (bad in list(NULL, NA, 0, -5, "15", c(10, 15))) {
    e <- tryCatch(ts_datatable(df, page_length = bad), error = function(e) e)
    expect_s3_class(e, "plot_datatable_error")
    expect_identical(e$state, "invalid_page_length")
  }
})

test_that("les tailles de page reelles de l'app sont toutes connues", {
  expect_setequal(ts_datatable_page_lengths(), c(6L, 8L, 10L, 15L, 20L))
})

test_that("donnees invalides => invalid_data ; matrix/liste sont coercés", {
  e <- tryCatch(ts_datatable("pas une table", page_length = 10),
                error = function(e) e)
  expect_s3_class(e, "plot_datatable_error")
  expect_identical(e$state, "invalid_data")

  m <- matrix(1:6, nrow = 2, dimnames = list(NULL, c("a", "b", "c")))
  expect_s3_class(ts_datatable(m, page_length = 10), "datatables")
  expect_s3_class(ts_datatable(list(a = 1:3), page_length = 10), "datatables")
})

test_that("buttons non logique => invalid_buttons", {
  df <- .tsd_df()
  for (bad in list("oui", NA, c(TRUE, FALSE))) {
    e <- tryCatch(ts_datatable(df, page_length = 10, buttons = bad),
                  error = function(e) e)
    expect_s3_class(e, "plot_datatable_error")
    expect_identical(e$state, "invalid_buttons")
  }
})

# --- 2. GARANTIE du defaut DT-EXPORT (contrat §6, option B) -------------------

test_that("appel nu = table de resultats boutonnee (defaut TS_DT_BUTTONS_DEFAULT)", {
  df <- .tsd_df()
  got <- ts_datatable(df, page_length = 15)
  expect_identical(got$x$options$dom, "Bfrtip")
  expect_identical(got$x$options$pageLength, 15)
  # 2H (O28, 2026-10-30) : scrollX n'est PLUS le defaut (arbitrage constat 2H —
  # voir bloc 4 ci-dessous) ; le filtre par colonne, lui, reste "top".
  expect_identical(got$x$options$scrollX, FALSE)
  expect_true("Buttons" %in% unlist(got$x$extensions))
  # aucun jeu de boutons nomme : DT applique ses boutons par defaut
  expect_null(got$x$options$buttons)
  expect_identical(got$x$filter, "top")
  expect_false(isTRUE(got$x$rownames))
})

test_that("le defaut buttons suit la constante TS_DT_BUTTONS_DEFAULT", {
  df <- .tsd_df()
  old <- if (exists("TS_DT_BUTTONS_DEFAULT", envir = globalenv())) {
    get("TS_DT_BUTTONS_DEFAULT", envir = globalenv())
  } else {
    NULL
  }
  assign("TS_DT_BUTTONS_DEFAULT", FALSE, envir = globalenv())
  on.exit({
    if (is.null(old)) {
      if (exists("TS_DT_BUTTONS_DEFAULT", envir = globalenv()))
        rm("TS_DT_BUTTONS_DEFAULT", envir = globalenv())
    } else {
      assign("TS_DT_BUTTONS_DEFAULT", old, envir = globalenv())
    }
  }, add = TRUE)
  got <- ts_datatable(df, page_length = 15)
  expect_null(got$x$options$dom)
  expect_true(is.null(got$x$extensions) || length(got$x$extensions) == 0L)
})

test_that("buttons = FALSE ne pose NI extensions NI dom (tables d'apercu)", {
  df <- .tsd_df()
  got <- ts_datatable(df, page_length = 8, buttons = FALSE)
  expect_null(got$x$options$dom)
  expect_true(is.null(got$x$extensions) || length(got$x$extensions) == 0L)
})

test_that("buttons = TRUE reproduit exactement le site pathways (dom + Buttons)", {
  df <- .tsd_df()
  ref <- DT::datatable(df, filter = "top", rownames = FALSE,
                       options = list(pageLength = 10, scrollX = TRUE,
                                      dom = "Bfrtip"),
                       extensions = "Buttons")
  # Le site pathways historique déclarait scrollX = TRUE EXPLICITEMENT (voir
  # les options du réf ci-dessus) : la reproduction passe donc scroll_x = TRUE.
  # Depuis le 2H (O28), le defaut du wrapper est FALSE — la règle « chaque site
  # garde ses valeurs écrites explicitement » n'est pas touchée.
  got <- ts_datatable(df, page_length = 10, buttons = TRUE, scroll_x = TRUE)
  expect_identical(got$x$options$dom, "Bfrtip")
  expect_identical(got$x$options$pageLength, 10)
  expect_identical(got$x$options$scrollX, TRUE)
  expect_true("Buttons" %in% unlist(got$x$extensions))
  # aucun jeu de boutons nomme : DT applique ses boutons par defaut
  expect_null(got$x$options$buttons)
  expect_identical(got$x$options$dom, ref$x$options$dom)
})

test_that("filename_base ajoute les 4 boutons nommes (forme recommandee)", {
  df <- .tsd_df()
  got <- ts_datatable(df, page_length = 10, buttons = TRUE,
                      filename_base = "pathways")
  expect_identical(got$x$options$dom, "Bfrtip")
  btns <- got$x$options$buttons
  expect_type(btns, "list")
  flat <- unlist(btns)
  expect_true("copy"  %in% flat)
  expect_true("print" %in% flat)
  expect_true("csv"   %in% flat)
  expect_true("excel" %in% flat)
  expect_true("pathways" %in% flat)
})

test_that("dom / extensions explicites ont la priorite sur buttons", {
  df <- .tsd_df()
  got <- ts_datatable(df, page_length = 10, buttons = TRUE, dom = "tip")
  expect_identical(got$x$options$dom, "tip")
})

test_that("extra_options : fusionnee avec priorite maximale, sinon invalid_extra_options", {
  df <- .tsd_df()
  got <- ts_datatable(df, page_length = 15, filename_base = "x",
                      extra_options = list(
                        language = list(search = "Filtrer :"),
                        lengthMenu = c(10, 15, 25)))
  expect_identical(got$x$options$language$search, "Filtrer :")
  expect_identical(got$x$options$lengthMenu, c(10, 15, 25))
  # les options posees par le wrapper restent presentes
  expect_identical(got$x$options$dom, "Bfrtip")
  # extra_options a la priorite sur le wrapper (ex. surcharger dom)
  got2 <- ts_datatable(df, page_length = 15, extra_options = list(dom = "t"))
  expect_identical(got2$x$options$dom, "t")
  # erreur classee si liste non nommee
  e <- tryCatch(ts_datatable(df, page_length = 15, extra_options = list(1, 2)),
                error = function(e) e)
  expect_s3_class(e, "plot_datatable_error")
  expect_identical(e$state, "invalid_extra_options")
})

test_that("les parametres de forme sont respectes et les ... transmis", {
  df <- .tsd_df()
  got <- ts_datatable(df, page_length = 20, filter = "none",
                      rownames = TRUE, scroll_x = FALSE,
                      selection = list(mode = "none"))
  expect_identical(got$x$options$pageLength, 20)
  expect_identical(got$x$options$scrollX, FALSE)
  expect_identical(got$x$filter, "none")
  # pas d'erreur => les `...` (selection, etc.) ont bien ete acceptes
  expect_s3_class(got, "datatables")
})

test_that("ts_datatable_buttons() est deterministe", {
  b <- ts_datatable_buttons("x")
  expect_identical(b, ts_datatable_buttons("x"))
  expect_identical(b[[1]], "copy")
  expect_identical(b[[3]]$filename, "x")
  expect_identical(b[[4]]$extend, "excel")
})

# --- 3. Anti-regression -------------------------------------------------------

.tsd_migrated <- c(
  "tools/pkg_skeleton/R/bulk_helpers.R",
  "tools/pkg_skeleton/R/pathway_helpers.R",
  "tools/pkg_skeleton/R/sc_helpers.R",
  "modules/bulk_de/mod_bulk_de_multimethod.R",
  "modules/bulk_de/mod_bulk_de_venn.R",
  "modules/spatial/mod_spatial_cluster.R"
)

test_that("PLOT-S3 : les 6 sites canoniques passent par ts_datatable()", {
  for (rel in .tsd_migrated) {
    path <- file.path(ts_project_root(), rel)
    code <- readLines(path, warn = FALSE, encoding = "UTF-8")
    code <- grep("^\\s*#", code, value = TRUE, invert = TRUE)
    src <- paste(code, collapse = "\n")
    # On cible precisement la forme canonique migree, pas n'importe quelle
    # table du fichier (certains fichiers contiennent d'autres tables, de forme
    # differente, volontairement laissees hors du perimetre PLOT-S3).
    pat <- paste0(
      "(?:DT::)?datatable\\s*\\(\\s*[A-Za-z0-9_.$]+\\s*,\\s*",
      "filter\\s*=\\s*.top.\\s*,\\s*rownames\\s*=\\s*FALSE\\s*,\\s*",
      "options\\s*=\\s*list\\(\\s*pageLength\\s*=\\s*[0-9]+\\s*,\\s*",
      "scrollX\\s*=\\s*TRUE")
    hits <- grepl(pat, src, perl = TRUE)
    expect_false(hits,
      paste0(rel, " : une table de forme canonique appelle encore datatable() ",
             "directement au lieu de ts_datatable()"))
  }
})

test_that("DT-EXPORT : les sites canoniques sont normalises (page 15 + boutons)", {
  # page_length = 15L + filename_base sur les tables de resultats
  expect_match(paste(readLines(file.path(ts_project_root(),
                                         "tools/pkg_skeleton/R/bulk_helpers.R")), collapse = "\n"),
               "ts_datatable\\(df_display, page_length = 15L, filename_base = \"bulk_degenes\"\\)")
  expect_match(paste(readLines(file.path(ts_project_root(),
                                         "tools/pkg_skeleton/R/pathway_helpers.R")), collapse = "\n"),
               "page_length = 15L")
  expect_match(paste(readLines(file.path(ts_project_root(),
                                         "modules/spatial/mod_spatial_cluster.R")), collapse = "\n"),
               "page_length = 15L")
})

test_that("DT-EXPORT : plus AUCUN appel datatable() direct hors du wrapper", {
  roots <- c(file.path(ts_project_root(), "modules"),
             file.path(ts_project_root(), "R"))
  offenders <- character(0)
  for (root in roots) {
    files <- list.files(root, pattern = "\\.R$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
    for (f in files) {
      if (normalizePath(f, winslash = "/") ==
          normalizePath(file.path(ts_project_root(), "R/plotting/datatable.R"),
                        winslash = "/")) {
        next
      }
      code <- readLines(f, warn = FALSE, encoding = "UTF-8")
      code <- grep("^\\s*#", code, value = TRUE, invert = TRUE)
      src <- paste(code, collapse = "\n")
      # (?<![A-Za-z0-9_]) : exclut ts_datatable( ; casse lowercase = appels reels
      if (grepl("(?<![A-Za-z0-9_])(?:DT::)?datatable\\s*\\(", src, perl = TRUE)) {
        offenders <- c(offenders, f)
      }
    }
  }
  expect_length(offenders, 0L)
  info <- if (length(offenders)) paste(offenders, collapse = ", ") else ""
  expect_identical(info, "")
})

# --- 4. 2H (O28) : arbitrage scrollX vs filter="top" --------------------------
# Constat 2H (packet round-7) : `scrollX = TRUE` + `filter = "top"` produisait
# le bug DataTables du double header / double filtre (le clonage du conteneur
# `thead` par le defilement horizontal duplique les boîtes de sélection du
# header flottant ; preuve datatable.R:105 à HEAD 88f9a43, aucun headerCallback ;
# spéc §2.19). Arbitrage retenu : le filtre par colonne — fonctionnalité
# porteuse d'ergonomie sur les tables de résultats — reste le defaut, le
# scrollX devient un opt-in explicite pour les tables réellement larges. Les
# deux correctifs JS/CSS de la spéc étant invérifiables hors navigateur,
# l'arbitrage retire la CAUSE du glitch (déterministe, vérifiable en R).

test_that("2H : appel nu = filter 'top' actif, scrollX OFF (arbitrage O28)", {
  df <- .tsd_df()
  got <- ts_datatable(df, page_length = 15)
  # L'ergonomie porteuse reste le defaut :
  expect_identical(got$x$filter, "top")
  # Le declencheur du double header quitte le defaut :
  expect_identical(got$x$options$scrollX, FALSE)
  # Le defaut vit dans la SIGNATURE (jamais un effet de bord d'une constante) :
  expect_identical(formals(ts_datatable)$scroll_x, FALSE)
})

test_that("2H : les tables larges optent explicitement pour le scrollX", {
  df <- .tsd_df()
  got <- ts_datatable(df, page_length = 10, scroll_x = TRUE)
  expect_identical(got$x$options$scrollX, TRUE)
  # Le reste du contrat est intact dans la même construction :
  expect_identical(got$x$filter, "top")
  expect_identical(got$x$options$pageLength, 10)
})
# --- PC (2026-11-09) : canal future opt-in ------------------------------------

test_that("PC : ts_dt_future() est FAUX par defaut (opt-in strict, zero changement)", {
  withr::with_options(list(TS_DT_FUTURE = NULL), {
    old_env <- Sys.getenv("TS_DT_FUTURE", unset = NA_character_)
    do.call(Sys.unsetenv, list("TS_DT_FUTURE"))
    on.exit(if (!is.na(old_env)) do.call(Sys.setenv, list(TS_DT_FUTURE = old_env)), add = TRUE)
    expect_false(ts_dt_future())
  })
})

test_that("PC : ts_dt_future() accepte option logique, chaine et variable d'env", {
  withr::with_options(list(TS_DT_FUTURE = TRUE), {
    old_env <- Sys.getenv("TS_DT_FUTURE", unset = NA_character_)
    do.call(Sys.unsetenv, list("TS_DT_FUTURE"))
    on.exit(if (!is.na(old_env)) do.call(Sys.setenv, list(TS_DT_FUTURE = old_env)), add = TRUE)
    expect_true(ts_dt_future())
  })
  withr::with_options(list(TS_DT_FUTURE = "1"), expect_true(ts_dt_future()))
  withr::with_options(list(TS_DT_FUTURE = "false"), expect_false(ts_dt_future()))
  withr::with_options(list(TS_DT_FUTURE = NULL), {
    withr::with_envvar(c(TS_DT_FUTURE = "1"), expect_true(ts_dt_future()))
    # une option explicite FALSE gagne sur l'env :
    withr::with_options(list(TS_DT_FUTURE = FALSE), {
      withr::with_envvar(c(TS_DT_FUTURE = "1"), expect_false(ts_dt_future()))
    })
  })
})

test_that("PC : ts_dt_future_install() pose le defaut du formal future, idempotent, et le test restaure", {
  skip_if_not_installed("future.mirai")
  on.exit(ts_dt_future_restore(), add = TRUE)

  expect_true(isTRUE(ts_dt_future_install()))
  got <- formals(DT::renderDT)$future
  expect_identical(got, quote(getOption("TS_DT_FUTURE", FALSE)))
  # idempotent :
  expect_true(isTRUE(ts_dt_future_install()))
  expect_identical(formals(DT::renderDT)$future, quote(getOption("TS_DT_FUTURE", FALSE)))

  # Le defaut pose EVALUE au defaut stock : option absente => FALSE (zero
  # changement), option TRUE => TRUE. Le formal reste surchargeable par un
  # appel explicite (site qui passe future = ... garde la priorite).
  withr::with_options(list(TS_DT_FUTURE = NULL), {
    expect_false(eval(formals(DT::renderDT)$future))
  })
  withr::with_options(list(TS_DT_FUTURE = TRUE), {
    expect_true(eval(formals(DT::renderDT)$future))
  })
})

test_that("PC : ts_dt_future_restore() revient au defaut stock FALSE", {
  skip_if_not_installed("future.mirai")
  ts_dt_future_install()
  expect_true(isTRUE(ts_dt_future_restore()))
  expect_identical(formals(DT::renderDT)$future, FALSE)
})

# --- CQF-3 (2026-11-11) : langue DT FR par defaut (additif contrat) -----------
# Constat mesure (triage UIUX Spatial 2026-11-08, RC1 ; triage SC 2026-11-09,
# RC6) : ts_datatable() ne pose JAMAIS `language` => 71 des 73 tables de l'app
# heritent des libelles anglais DataTables (Copy/Print/CSV/Excel, `Search:`,
# lengthMenu anglais) dans une app a defaut FR (I18N_DEFAULT_LANG = "fr").
# Fix additif : defaut FR literal dans le moteur (option TS_DT_LANGUAGE = voie
# i18n app-side, meme arbitrage O6 que TS_PACKS_REGISTRY), extra_options garde
# la priorite maximale (les 2 sites existants qui passent language = ... ne
# changent pas). dom / filter restent GELÉS (contrat §2H / §6).

test_that("CQF-3 : appel nu = language FR par defaut (search/lengthMenu/info/boutons)", {
  df <- .tsd_df()
  withr::with_options(list(TS_DT_LANGUAGE = NULL), {
    got <- ts_datatable(df, page_length = 15, filename_base = "x")
    lang <- got$x$options$language
    expect_identical(lang$search, "Recherche :")
    expect_identical(lang$lengthMenu, "Afficher _MENU_ lignes")
    expect_identical(lang$info, "Lignes _START_ à _END_ sur _TOTAL_")
    expect_identical(lang$infoEmpty, "Aucune ligne")
    expect_identical(lang$infoFiltered, "(filtrées à partir de _MAX_ lignes)")
    expect_identical(lang$emptyTable, "Aucune donnée disponible")
    expect_identical(lang$zeroRecords, "Aucun résultat trouvé")
    expect_identical(lang$processing, "Traitement...")
    expect_identical(lang$paginate$first, "Premier")
    expect_identical(lang$paginate$last, "Dernier")
    expect_identical(lang$paginate[["next"]], "Suivant")
    expect_identical(lang$paginate$previous, "Précédent")
    expect_identical(lang$buttons$copy, "Copier")
    expect_identical(lang$buttons$print, "Imprimer")
    expect_identical(lang$buttons$pageLength[[2]], "Tout afficher")
    # csv/excel : libelles identiques FR/EN — pas de cle i18n
    expect_identical(lang$buttons$csv, "CSV")
    expect_identical(lang$buttons$excel, "Excel")
  })
})

test_that("CQF-3 : l'option TS_DT_LANGUAGE remplace le defaut wholesale (voie i18n EN)", {
  df <- .tsd_df()
  withr::with_options(list(TS_DT_LANGUAGE = list(
    search = "Search:", info = "Showing _START_ to _END_ of _TOTAL_ entries",
    buttons = list(copy = "Copy", print = "Print", csv = "CSV", excel = "Excel")
  )), {
    got <- ts_datatable(df, page_length = 15)
    expect_identical(got$x$options$language$search, "Search:")
    expect_identical(got$x$options$language$buttons$copy, "Copy")
    # l'option remplace le defaut TEL QUEL (pas de fusion partielle silencieuse)
    expect_null(got$x$options$language$emptyTable)
  })
})

test_that("CQF-3 : TS_DT_LANGUAGE invalide => avertissement + repli FR honnete", {
  df <- .tsd_df()
  withr::with_options(list(TS_DT_LANGUAGE = "pas une liste"), {
    expect_warning(
      got <- ts_datatable(df, page_length = 15),
      "TS_DT_LANGUAGE"
    )
    expect_identical(got$x$options$language$search, "Recherche :")
  })
  withr::with_options(list(TS_DT_LANGUAGE = list(1, 2, 3)), {
    expect_warning(
      got2 <- ts_datatable(df, page_length = 15),
      "TS_DT_LANGUAGE"
    )
    expect_identical(got2$x$options$language$search, "Recherche :")
  })
})

test_that("CQF-3 : extra_options garde la priorite maximale sur le defaut language (fusion recursive)", {
  df <- .tsd_df()
  withr::with_options(list(TS_DT_LANGUAGE = NULL), {
    got <- ts_datatable(df, page_length = 15,
                        extra_options = list(language = list(search = "Filtrer :")))
    expect_identical(got$x$options$language$search, "Filtrer :")
    # les cles non surchargees restent FR par defaut (fusion recursive)
    expect_identical(got$x$options$language$emptyTable, "Aucune donnée disponible")
  })
})

test_that("CQF-3 : dom / filter du contrat restent GELÉS sous le defaut language", {
  df <- .tsd_df()
  withr::with_options(list(TS_DT_LANGUAGE = NULL), {
    got <- ts_datatable(df, page_length = 15)
    expect_identical(got$x$options$dom, "Bfrtip")
    expect_identical(got$x$filter, "top")
  })
})
