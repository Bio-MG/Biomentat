# =============================================================================
# test-lane2-nul-defense.R — Lane 2 (2026-11-14, QC P0.4 §H.0)
# =============================================================================
# CRASH MESURÉ 19:06:57 : `Error in execCallbacks : embedded nul in string:
# '[19:06:57] Mapping IDs : symboles dC)jC... b\0\024 ignorC).'` — UNE chaîne
# i18n mojibakée + NUL embarqué qui entre dans une mise à jour réactive (log
# sidebar, mod_bulk.R:624) tue TOUT le process app (toutes sessions).
#
# Cause racine double :
#   (a) tools/launch_dev_drive.R désarmait la recette documentée
#       (Sys.unsetenv("LC_ALL")) puis tentait Sys.setlocale("LC_CTYPE",
#       "fr_FR.UTF-8") — nom UNIX, échec silencieux sur Windows ;
#   (b) AUCUNE défense d'encodage : le NUL traverse .tr() -> log() ->
#       renderText -> Shiny sans jamais être démaillé.
#
# Fix (ce lot) :
#   (a) le launcher se SOIGNE lui-même (Sys.setlocale("LC_ALL",
#       "French_France.65001") réussit en session Windows — mesuré) et REFUSE
#       de démarrer si la locale n'est pas effective (override conscient :
#       options(ts.drive.allow_bad_locale = TRUE)) ;
#   (b) .ts_nul_free_text() (global.R) démaille les NULs par les octets
#       (charToRaw -> filtre 0x00 -> rawToChar) + enc2utf8, branché AU POINT
#       DE CHOIX .strip_i18n_html (tous les proxys .tr()/.tr_plain() de tous
#       les modules) ET au renderText du log auto-pipeline (le vecteur exact
#       du crash).
#
# NOTE honnête sur la falsification : le R pur ne peut PAS construire de
# chaîne à NUL (littéraux refusés au parse ; readChar tronque au 1er NUL —
# mesuré) — les NULs arrivent des conversions C. Le contrat est donc testé
# par : passthrough byte-exact des textes normaux + contrat des points de
# branchement + mécanisme octet-vérifié (le filtre 0x00 est le seul chemin
# qui peut retirer un NUL d'une chaîne R).
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")

test_that("lane2 : .ts_nul_free_text est un passthrough byte-exact des textes normaux", {
  s <- "Mapping IDs : symboles déjà détectés (ou type inconnu) — étape ignorée. ✓"
  out <- .ts_nul_free_text(s)
  expect_identical(out, s)
  expect_identical(charToRaw(out), charToRaw(s))
  # vecteur + chaînes vides/NA-like
  v <- c("", "Filtrage & VST...", "✓ 309 cellules × 35 000 gènes")
  outv <- .ts_nul_free_text(v)
  expect_identical(outv, v)
  # non-caractère : traversé tel quel (les proxys peuvent renvoyer des tags)
  expect_identical(.ts_nul_free_text(42L), 42L)
  expect_null(.ts_nul_free_text(NULL))
})

test_that("lane2 : le mécanisme octet (filtre 0x00) est le chemin de démaillage", {
  # Le seul moyen fiable de retirer un NUL d'une chaîne R est le passage par
  # les octets. Contrat : charToRaw ne perd aucune octet, le filtre 0x00
  # retire exactement les octets nuls, rawToChar des octets non-nuls est sûr.
  r <- charToRaw("abc")
  expect_false(any(r == as.raw(0L)))
  mixed <- as.raw(c(0x61, 0x00, 0x62, 0x00, 0x63))
  expect_identical(rawToChar(mixed[mixed != as.raw(0L)]), "abc")
})

test_that("lane2 : la défense est branchée au point de choix i18n et au log sidebar", {
  # le helper vit dans io_helpers.R (sourcé par tout harnais) :
  srci <- paste(readLines(file.path(ts_project_root(), "tools", "pkg_skeleton",
                                    "R", "io_helpers.R"),
                          warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(srci, ".ts_nul_free_text <- function", fixed = TRUE,
               info = "helper lane2 absent de io_helpers.R")
  # .strip_i18n_html (point de choix de TOUS les proxys .tr/.tr_plain) passe
  # par le démaillage :
  src <- paste(readLines(file.path(ts_project_root(), "global.R"),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, "return(.ts_nul_free_text(flat))", fixed = TRUE,
               info = ".strip_i18n_html ne démaille pas la branche tag")
  expect_match(src, "if (is.character(res)) return(.ts_nul_free_text(res))",
               fixed = TRUE,
               info = ".strip_i18n_html ne démaille pas la branche scalaire")
  # le vecteur exact du crash : le log auto-pipeline est démaillé au rendu
  srcb <- paste(readLines(file.path(ts_project_root(), "modules", "bulk",
                                    "mod_bulk.R"),
                          warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(srcb, "renderText({ .ts_nul_free_text(auto_log_rv()) })",
               fixed = TRUE,
               info = "le log auto-pipeline (vecteur du crash) n'est pas défendu")
})

test_that("lane2 : le launcher impose la recette locale et ne la désarme plus", {
  src <- paste(readLines(file.path(ts_project_root(), "tools",
                                   "launch_dev_drive.R"),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_identical(lengths(regmatches(src, gregexpr('Sys.unsetenv("LC_ALL")',
                                                    src, fixed = TRUE)))[1], 0L,
                   info = "le launcher DÉSARME encore la recette LC_ALL")
  expect_match(src, 'Sys.setlocale("LC_ALL", "French_France.65001")', fixed = TRUE,
               info = "auto-soin de la recette absent")
  expect_match(src, "ts.drive.allow_bad_locale", fixed = TRUE,
               info = "override conscient absent")
  expect_match(src, "embedded nul", fixed = TRUE,
               info = "le refus de boot n'explique pas le risque (crash §H.0)")
})
