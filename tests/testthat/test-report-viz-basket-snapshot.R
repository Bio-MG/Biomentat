# =============================================================================
# test-report-viz-basket-snapshot.R — P4.1b gap n°3 (2026-11-17)
# =============================================================================
# IDEATION_SAVE_PLOTS_2026-11-17.md, gap n°3 : le panier de viz du rapport
# SC (report_viz_list — configs légères écrites par mod_sc_viz « 📌 Ajouter
# au Rapport », PAS des objets plots) vivait déjà dans sc_shared mais
# n'était NI capturé au snapshot NI purgé au remplacement de monde.
#
# Câblage posé :
#   state.R : .ts_bump_all_epochs purge le panier — WORLD-SCOPED, les DEUX
#           chemins destructifs (load_session_file / confirm_reset) passent
#           par ce bump (la purge standard sc_purge_shared_results ne le
#           couvre pas : il n'est pas un résultat dérivé) ;
#   app.R : save (results$report_viz_list) + pending payload + ré-application
#           INCONDITIONNELLE (snapshot antérieur sans panier = panier VIDE,
#           jamais un panier d'un ancien monde mélangé au monde restauré ;
#           l'ordre purge → ré-application est garanti par l'ordre de
#           création des observateurs : pending créé après mod_sc_server).
#
# Le panier n'exige AUCUN hook module : le UI du panier (mod_sc.R) et le
# rendu rapport (param saved_viz_list) lisent shared_rv$report_viz_list
# réactivement — ré-appliqué = ré-affiché.
# =============================================================================

Sys.setenv(OMNIPATHR_LOGFILE = "none")

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/state.R")

suppressPackageStartupMessages({
  library(shiny)
})

.read_app <- function(rel) {
  paste(readLines(file.path(ts_project_root(), rel),
                  warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# CONTRAT ÉTAT — le bump d'epoch purgé le panier (world-scoped)
# ══════════════════════════════════════════════════════════════════════════════

test_that("p4.1b : le bump d'epoch purge le panier de viz (les deux chemins destructifs)", {
  gd <- reactiveValues(
    sc_obj_epoch = 0L, bulk_obj_epoch = 0L,
    sc_shared = reactiveValues(report_viz_list = list(viz1 = list(type = "dimplot"),
                                                      viz2 = list(type = "vlnplot"))))
  shiny::isolate(.ts_bump_all_epochs(gd))
  shiny::isolate({
    expect_identical(gd$sc_obj_epoch, 1L)
    expect_identical(gd$bulk_obj_epoch, 1L)
    # Le panier appartient au monde remplacé — vide après bump.
    expect_identical(gd$sc_shared$report_viz_list, list())
  })
})

test_that("p4.1b : le bump SANS sc_shared reste un no-op sûr (boot pré-init)", {
  gd <- reactiveValues(sc_obj_epoch = 2L, bulk_obj_epoch = 5L)
  shiny::isolate(expect_silent(.ts_bump_all_epochs(gd)))
  shiny::isolate({
    expect_identical(gd$sc_obj_epoch, 3L)
    expect_identical(gd$bulk_obj_epoch, 6L)
  })
})

# ══════════════════════════════════════════════════════════════════════════════
# CONTRATS SOURCE app.R — save + payload pending + ré-application
# ══════════════════════════════════════════════════════════════════════════════

test_that("p4.1b : le snapshot embarque le panier (save_session_btn)", {
  src <- .read_app("app.R")
  expect_match(src, "report_viz_list             = shared$report_viz_list,",
               fixed = TRUE,
               info = "le panier de viz n'est pas capturé dans le bloc results")
  # Le poids estimé est honnête : le panier compte dans `extra`.
  expect_match(src, "shared$report_viz_list, shared$provenance) else list()",
               fixed = TRUE,
               info = "le panier de viz ne compte pas dans le poids estimé")
})

test_that("p4.1b : la restauration ré-applique le panier INCONDITIONNELLEMENT", {
  src <- .read_app("app.R")
  # Le payload pending porte le champ (même absent du snapshot — NULL).
  expect_match(src, "report_viz_list = res$report_viz_list,",
               fixed = TRUE,
               info = "le payload pending ne porte pas le panier de viz")
  # Ré-application sans garde de présence : snapshot antérieur = panier VIDE
  # (le bump l'a déjà purgé — jamais de panier d'un ancien monde).
  expect_match(src, "shared$report_viz_list <- pay$report_viz_list %||% list()",
               fixed = TRUE,
               info = "la ré-application n'est pas inconditionnelle (vide si absent)")
})

test_that("p4.1b : le panier reste des CONFIGS légères (jamais d'objet plot)", {
  src <- .read_app("modules/sc/mod_sc_viz.R")
  # Le seul write du panier : current[[title]] <- cfg — une config d'UI.
  expect_match(src, "current[[title]] <- cfg", fixed = TRUE,
               info = "le panier n'écrit plus des configs (objet plot détecté ?)")
})
