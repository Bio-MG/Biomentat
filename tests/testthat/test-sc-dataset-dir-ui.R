# =============================================================================
# test-sc-dataset-dir-ui.R -- W1.3 : SC UX (btn_load_dir, .sc_utf8, discovery)
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

.scu_stub_ui <- function() {
  saved <- list(
    i18n = get0("i18n", envir = globalenv(), inherits = FALSE),
    .tr_plain = get0(".tr_plain", envir = globalenv(), inherits = FALSE),
    i18n_plain = get0("i18n_plain", envir = globalenv(), inherits = FALSE)
  )
  assign("i18n", list(t = function(x) x), envir = globalenv())
  assign(".tr_plain", function(x) x, envir = globalenv())
  assign("i18n_plain", function(x) x, envir = globalenv())
  saved
}

.scu_restore_ui <- function(saved) {
  for (nm in names(saved)) {
    if (is.null(saved[[nm]])) {
      rm(list = intersect(nm, ls(globalenv())), envir = globalenv())
    } else {
      assign(nm, saved[[nm]], envir = globalenv())
    }
  }
}

.scu_env <- function() {
  gd <- new.env(parent = emptyenv())
  gd$i18n <- NULL
  gd$drive_registry <- new.env(parent = emptyenv())
  gd$sc_datasets <- list()
  gd
}

.scu_make_seurat_rds <- function(path) {
  cnt <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L),
    j = c(1L, 1L, 2L, 2L),
    x = c(5, 8, 3, 9),
    dims = c(2L, 2L),
    dimnames = list(c("G1", "G2"), c("CELL1", "CELL2"))
  )
  obj <- SeuratObject::CreateSeuratObject(counts = cnt, project = "S")
  saveRDS(obj, path)
  path
}

test_that("mod_import_sc_ui : btn_load_dir initialement disabled et panel_datasets present", {
  saved <- .scu_stub_ui()
  on.exit(.scu_restore_ui(saved), add = TRUE)

  html <- paste(as.character(mod_import_sc_ui("import_sc")), collapse = "\n")
  expect_match(html, 'id="import_sc-btn_load_dir"', fixed = TRUE)
  expect_match(html, 'shinyjs-disabled', fixed = TRUE)
  expect_match(html, 'data-value="panel_datasets"', fixed = TRUE)
  expect_match(html, 'id="import_sc-btn_scan_datasets"', fixed = TRUE)
})

test_that(".sc_utf8 et add_log() convertissent les octets Latin-1 sans caracteres corrompus", {
  raw_latin1 <- iconv("erreur d'évaluation : accès refusé", from = "UTF-8", to = "latin1")
  Encoding(raw_latin1) <- "unknown"
  fixed <- .sc_utf8(raw_latin1)
  expect_true(all(validUTF8(fixed)))
  expect_identical(fixed, "erreur d'évaluation : accès refusé")

  gd <- .scu_env()
  shiny::testServer(mod_import_sc_server, args = list(global_data = gd), expr = {
    add_log(raw_latin1)
    lg <- shiny::isolate(logs())
    expect_true(all(validUTF8(lg)))
    expect_match(lg, "évaluation", fixed = TRUE)
  })
})

test_that("mod_import_sc_server : aucun scan au cold-boot, scan a la demande et double garde AC3", {
  root <- tempfile("sc_ds_ui_")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  .scu_make_seurat_rds(file.path(root, "valid_sc.rds"))
  dir.create(file.path(root, "empty_dir"))
  bad_dir <- file.path(root, "bad_10x")
  dir.create(bad_dir)
  writeLines("%%MatrixMarket", file.path(bad_dir, "matrix.mtx"))

  old_ds <- Sys.getenv("TRANSCRIPTO_DATASETS_DIR", unset = NA_character_)
  on.exit({
    if (is.na(old_ds)) Sys.unsetenv("TRANSCRIPTO_DATASETS_DIR") else Sys.setenv(TRANSCRIPTO_DATASETS_DIR = old_ds)
  }, add = TRUE)
  Sys.setenv(TRANSCRIPTO_DATASETS_DIR = root)

  gd <- .scu_env()
  shiny::testServer(mod_import_sc_server, args = list(global_data = gd), expr = {
    session$flushReact()
    expect_null(shiny::isolate(ds_scan_rv()))
    expect_null(shiny::isolate(global_data$sc_obj))
    expect_false(grepl("Scan du répertoire", shiny::isolate(logs()), fixed = TRUE))

    session$setInputs(btn_scan_datasets = 1L)
    session$flushReact()
    sc <- shiny::isolate(ds_scan_rv())
    expect_length(sc, 3L)
    expect_identical(sc[["valid_sc.rds"]]$status, "valid")
    expect_identical(sc[["empty_dir"]]$status, "empty")
    expect_identical(sc[["bad_10x"]]$status, "invalid")

    rows <- .sc_datasets_table_rows(sc, ns = function(x) x, tr = function(k) k)
    expect_null(shiny::tagGetAttribute(rows[["valid_sc.rds"]]$button, "disabled"))
    expect_false(is.null(shiny::tagGetAttribute(rows[["empty_dir"]]$button, "disabled")))
    expect_false(is.null(shiny::tagGetAttribute(rows[["bad_10x"]]$button, "disabled")))

    # Force click on invalid entry -> refused server-side
    bad_btn <- rows[["bad_10x"]]$btn_id
    do.call(session$setInputs, stats::setNames(list(1L), bad_btn))
    session$flushReact()
    expect_null(shiny::isolate(global_data$sc_obj))
    expect_match(shiny::isolate(logs()), "bad_10x", fixed = TRUE)

    # Click valid entry -> loads sc_obj and bumps epoch
    ok_btn <- rows[["valid_sc.rds"]]$btn_id
    do.call(session$setInputs, stats::setNames(list(1L), ok_btn))
    session$flushReact()
    expect_false(is.null(shiny::isolate(global_data$sc_obj)))
    expect_equal(ncol(shiny::isolate(global_data$sc_obj)), 2L)
    expect_identical(shiny::isolate(global_data$sc_obj_epoch), 1L)

    # btn_load_dir with empty sample_list() does not crash or modify sc_obj
    prev_epoch <- shiny::isolate(global_data$sc_obj_epoch)
    session$setInputs(btn_load_dir = 1L)
    session$flushReact()
    expect_identical(shiny::isolate(global_data$sc_obj_epoch), prev_epoch)
  })
})

test_that("mod_import_sc_server : uncompressed 10x_v3 et gzipped 10x_v2 se chargent sans modifier le dossier source", {
  root <- tempfile("sc_10x_stage_")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  mtx <- c("%%MatrixMarket matrix coordinate integer general", "%", "2 2 4",
           "1 1 5", "2 1 8", "1 2 3", "2 2 9")

  # 1. Uncompressed 10x_v3
  v3_dir <- file.path(root, "01_plain_v3")
  dir.create(v3_dir)
  writeLines(mtx, file.path(v3_dir, "matrix.mtx"))
  writeLines(c("CELL1-1", "CELL2-1"), file.path(v3_dir, "barcodes.tsv"))
  writeLines(c("ENSG1\tGene1\tGene Expression", "ENSG2\tGene2\tGene Expression"),
             file.path(v3_dir, "features.tsv"))
  v3_before <- sort(list.files(v3_dir))

  # 2. Gzipped 10x_v2
  v2_dir <- file.path(root, "02_gz_v2")
  dir.create(v2_dir)
  write_gz <- function(lines, path) {
    con <- gzfile(path, "wt")
    on.exit(close(con), add = TRUE)
    writeLines(lines, con)
  }
  write_gz(mtx, file.path(v2_dir, "matrix.mtx.gz"))
  write_gz(c("CELLA-1", "CELLB-1"), file.path(v2_dir, "barcodes.tsv.gz"))
  write_gz(c("ENSG1\tGeneA", "ENSG2\tGeneB"), file.path(v2_dir, "genes.tsv.gz"))
  v2_before <- sort(list.files(v2_dir))

  old_ds <- Sys.getenv("TRANSCRIPTO_DATASETS_DIR", unset = NA_character_)
  on.exit({
    if (is.na(old_ds)) Sys.unsetenv("TRANSCRIPTO_DATASETS_DIR") else Sys.setenv(TRANSCRIPTO_DATASETS_DIR = old_ds)
  }, add = TRUE)
  Sys.setenv(TRANSCRIPTO_DATASETS_DIR = root)

  gd <- .scu_env()
  shiny::testServer(mod_import_sc_server, args = list(global_data = gd), expr = {
    session$setInputs(btn_scan_datasets = 1L)
    session$flushReact()
    sc <- shiny::isolate(ds_scan_rv())
    expect_identical(sc[["01_plain_v3"]]$status, "valid")
    expect_identical(sc[["01_plain_v3"]]$format, "10x_v3")
    expect_identical(sc[["02_gz_v2"]]$status, "valid")
    expect_identical(sc[["02_gz_v2"]]$format, "10x_v2")

    rows <- .sc_datasets_table_rows(sc, ns = function(x) x, tr = function(k) k)

    # Load uncompressed 10x_v3
    do.call(session$setInputs, stats::setNames(list(1L), rows[["01_plain_v3"]]$btn_id))
    session$flushReact()
    obj3 <- shiny::isolate(global_data$sc_obj)
    expect_false(is.null(obj3))
    expect_equal(ncol(obj3), 2L)
    expect_equal(nrow(obj3), 2L)
    expect_identical(shiny::isolate(global_data$sc_obj_epoch), 1L)
    expect_identical(sort(list.files(v3_dir)), v3_before)

    # Load gzipped 10x_v2
    do.call(session$setInputs, stats::setNames(list(1L), rows[["02_gz_v2"]]$btn_id))
    session$flushReact()
    obj2 <- shiny::isolate(global_data$sc_obj)
    expect_false(is.null(obj2))
    expect_equal(ncol(obj2), 2L)
    expect_equal(nrow(obj2), 2L)
    expect_identical(shiny::isolate(global_data$sc_obj_epoch), 2L)
    expect_identical(sort(list.files(v2_dir)), v2_before)
  })
})

test_that(".charger_sc_dataset preserve orig.ident sur un Seurat .rds pre-fusionne", {
  root <- tempfile("sc_premerged_")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  cnt <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L, 1L, 2L),
    j = c(1L, 1L, 2L, 2L, 3L, 3L),
    x = c(5, 8, 3, 9, 4, 7),
    dims = c(2L, 3L),
    dimnames = list(c("G1", "G2"), c("C1", "C2", "C3"))
  )
  multi_obj <- SeuratObject::CreateSeuratObject(counts = cnt, project = "Cohort")
  multi_obj$orig.ident <- c("Ctrl_1", "Ctrl_1", "Treat_1")
  saveRDS(multi_obj, file.path(root, "pbmc_cohort.rds"))

  default_obj <- SeuratObject::CreateSeuratObject(counts = cnt, project = "SeuratProject")
  saveRDS(default_obj, file.path(root, "plain_default.rds"))

  old_ds <- Sys.getenv("TRANSCRIPTO_DATASETS_DIR", unset = NA_character_)
  on.exit({
    if (is.na(old_ds)) Sys.unsetenv("TRANSCRIPTO_DATASETS_DIR") else Sys.setenv(TRANSCRIPTO_DATASETS_DIR = old_ds)
  }, add = TRUE)
  Sys.setenv(TRANSCRIPTO_DATASETS_DIR = root)

  gd <- .scu_env()
  shiny::testServer(mod_import_sc_server, args = list(global_data = gd), expr = {
    session$setInputs(btn_scan_datasets = 1L)
    session$flushReact()
    sc <- shiny::isolate(ds_scan_rv())
    rows <- .sc_datasets_table_rows(sc, ns = function(x) x, tr = function(k) k)

    do.call(session$setInputs, stats::setNames(list(1L), rows[["pbmc_cohort.rds"]]$btn_id))
    session$flushReact()
    loaded_multi <- shiny::isolate(global_data$sc_obj)
    expect_identical(as.character(loaded_multi$orig.ident), c("Ctrl_1", "Ctrl_1", "Treat_1"))

    do.call(session$setInputs, stats::setNames(list(1L), rows[["plain_default.rds"]]$btn_id))
    session$flushReact()
    loaded_def <- shiny::isolate(global_data$sc_obj)
    expect_true(all(as.character(loaded_def$orig.ident) == "plain_default.rds"))
  })
})
