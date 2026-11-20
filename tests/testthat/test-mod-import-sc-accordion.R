# =============================================================================
# test-mod-import-sc-accordion.R -- W1.2 #3b : accordion « Datasets découverts »
# doit être pilotable/expansible (id namespacé, multiple, panneau ré-ouvert au scan)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(DT)
  library(bsicons)
  library(shinyFiles)
  library(Seurat)
})
source_project_file("tools/pkg_skeleton/R/rdata_io.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("tools/pkg_skeleton/R/session_checkpoint.R")
source_project_file("tools/pkg_skeleton/R/biomentat_handoff.R")
source_project_file("tools/pkg_skeleton/R/dataset_directory.R")
source_project_file("tools/pkg_skeleton/R/sc_multi.R")
source_project_file("modules/import/mod_rdata_picker.R")
source_project_file("modules/import/mod_import_sc.R")

.sca_html <- function() {
  saved <- list(
    i18n = get0("i18n", envir = globalenv(), inherits = FALSE),
    .tr_plain = get0(".tr_plain", envir = globalenv(), inherits = FALSE),
    i18n_plain = get0("i18n_plain", envir = globalenv(), inherits = FALSE)
  )
  on.exit({
    for (nm in names(saved)) {
      if (is.null(saved[[nm]])) {
        rm(list = intersect(nm, ls(globalenv())), envir = globalenv())
      } else assign(nm, saved[[nm]], envir = globalenv())
    }
  }, add = TRUE)
  assign("i18n", list(t = function(x) x), envir = globalenv())
  assign(".tr_plain", function(x) x, envir = globalenv())
  assign("i18n_plain", function(x) x, envir = globalenv())
  paste(as.character(mod_import_sc_ui("import_sc")), collapse = "\n")
}

test_that("accordion SC porte un id namespacé (pilotable par accordion_panel_open)", {
  html <- .sca_html()
  expect_match(html, 'id="import_sc-acc_sources"', fixed = TRUE)
})

test_that("panel_datasets est un item accordion repliable avec cible unique", {
  html <- .sca_html()
  m <- regmatches(html, regexpr(
    '(?s)data-value="panel_datasets">.*?data-bs-target="#([^"]+)"', html, perl = TRUE))
  expect_length(m, 1L)
  target <- sub('(?s).*data-bs-target="#([^"]+)"$', "\\1", m, perl = TRUE)
  expect_equal(lengths(regmatches(html, gregexpr(paste0('id="', target, '"'), html, fixed = TRUE))), 1L)
  ids <- regmatches(html, gregexpr('id="[^"]+"', html))[[1]]
  expect_false(anyDuplicated(ids) > 0L)
})

test_that("le scan rouvre le panneau « Datasets découverts » via accordion_panel_open", {
  src <- paste(readLines(file.path(.ts_find_project_root(), "modules/import/mod_import_sc.R"), warn = FALSE, encoding = "UTF-8"),
               collapse = "\n")
  expect_match(src, 'accordion_panel_open\\(\\s*"acc_sources"\\s*,\\s*"panel_datasets"')
})
