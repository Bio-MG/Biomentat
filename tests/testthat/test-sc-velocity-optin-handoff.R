# =============================================================================
# test-sc-velocity-optin-handoff.R — P0.4 lane 1 (2026-11-13)
# =============================================================================
# Leçon mesurée 2026-11-13 (P0.4 live, rapport QC/p04_live_check_2026-11-13
# §H.8) : l'opérateur a lancé l'inférence opt-in scVelo sur 4444 cellules ;
# scVelo a CALCULÉ (log serveur : 'velocity' + 'velocity_graph' + 'computing
# velocity embedding'), puis l'état est RESTÉ « valid_no_vectors » — la vue
# « Vecteurs UMAP » est restée sans flèches ET le contexte vélocité de
# Communication a refusé (« les vecteurs ne sont jamais ré-inférés »).
#
# CAUSES mesurées au source :
#   (1) sc_velocity_scvelo_optin_infer() convertit TOUT en dense :
#       as.matrix(result$spliced[...]) + np$array(t(spliced)) (2 copies
#       pleines côté R + 2 côté python) ET re-convertit la couche gène
#       COMPLÈTE adata$layers["velocity"] vers R (~1,2 Gio à 4444×35k) alors
#       qu'AUCUN consommateur ne lit gene_vectors (grep global mesuré) —
#       cause probable de l'échec mémoire de la session live, avalé par le
#       tryCatch du module (aucune trace console).
#   (2) l'échec opt-in ne laisse AUCUNE trace persistante : après la
#       notification de 12 s, l'opérateur voit « aucun vecteur fourni » sans
#       la raison ; les sessions headless n'ont RIEN au console.
#
# Ce lot (RED avant fix) :
#   (A) contrats source anti-coercition dense + gene_vectors OFF par défaut ;
#   (B) trace persistante de l'échec opt-in (velocity_state$optin_error +
#       statut portant la raison) — stub moteur absent ;
#   (C) bout-en-bout testServer avec le moteur scVelo RÉEL (skip si absent) :
#       validate → opt-in → shared_rv$velocity_result doit porter des
#       vecteurs (status « valid ») — la publication est le contrat rendu.
# =============================================================================

.velopt_genes <- paste0("gene", seq_len(40))
.velopt_cells <- paste0("cell", seq_len(40))

.velopt_counts <- function() {
  set.seed(11)
  sp <- matrix(rpois(length(.velopt_genes) * length(.velopt_cells), 5),
               nrow = length(.velopt_genes), ncol = length(.velopt_cells),
               dimnames = list(.velopt_genes, .velopt_cells))
  un <- matrix(rpois(length(.velopt_genes) * length(.velopt_cells), 2),
               nrow = length(.velopt_genes), ncol = length(.velopt_cells),
               dimnames = list(.velopt_genes, .velopt_cells))
  list(spliced = Matrix::Matrix(sp, sparse = TRUE),
       unspliced = Matrix::Matrix(un, sparse = TRUE))
}

.velopt_seurat_with_umap <- function() {
  d <- .velopt_counts()
  seu <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = as.matrix(d$spliced)))
  emb <- matrix(stats::runif(length(.velopt_cells) * 2),
                nrow = length(.velopt_cells), ncol = 2,
                dimnames = list(.velopt_cells, c("UMAP_1", "UMAP_2")))
  seu[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = emb, key = "UMAP_", assay = "RNA")
  seu
}

.velopt_python_available <- function() {
  py <- file.path(ts_project_root(), "python_env_sccoda", "Scripts", "python.exe")
  if (!file.exists(py)) return(FALSE)
  if (!identical(Sys.getenv("RETICULATE_PYTHON"), "")) {
    return(TRUE) # déjà configuré (launcher / .Renviron)
  }
  Sys.setenv(RETICULATE_PYTHON = normalizePath(py, winslash = "/"))
  TRUE
}

# ---------------------------------------------------------------------------
# (A) Contrats source — anti-coercition dense + gene_vectors OFF par défaut
# ---------------------------------------------------------------------------
test_that("P0.4 lane1 : l'inférence opt-in ne dense-ifie plus les matrices (source)", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "tools/pkg_skeleton/R/sc_velocity.R"),
                         warn = FALSE), collapse = "\n")
  # plus de conversion dense R-side des matrices validées avant AnnData
  expect_false(grepl("as\\.matrix\\(result\\$spliced", src))
  expect_false(grepl("as\\.matrix\\(result\\$unspliced", src))
  # plus d'injection dense np$array côté python — entrée scipy sparse
  expect_false(grepl("np\\$array\\(t\\(spliced\\)\\)", src))
  expect_false(grepl("np\\$array\\(t\\(unspl\\)\\)", src))
  expect_true(grepl("csc_matrix|csr_matrix|scipy\\.sparse", src))
})

test_that("P0.4 lane1 : la couche gène complète n'est plus re-convertie vers R (source)", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "tools/pkg_skeleton/R/sc_velocity.R"),
                         warn = FALSE), collapse = "\n")
  # le py_to_r de la couche gène ENTIÈRE (1,2 Gio à l'échelle réelle, zéro
  # consommateur mesuré) est RETIRÉ — plus de branche conditionnelle, plus de
  # paramètre mort : l'embedding seul (cellules × 2) est converti
  expect_false(grepl('adata\\$layers\\["velocity"\\]', src))
  expect_false(grepl("return_gene_vectors", src))
  expect_false(grepl("gene_vectors", src))
})

# ---------------------------------------------------------------------------
# (B) Trace persistante de l'échec opt-in (moteur stubbé absent)
# ---------------------------------------------------------------------------
test_that("P0.4 lane1 : échec opt-in => optin_error persistant + statut avec raison", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("Seurat")
  source_project_file("tools/pkg_skeleton/R/io_helpers.R")
  source_project_file("tools/pkg_skeleton/R/validation.R")
  source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
  source_project_file("modules/sc/mod_sc_velocity.R")

  # stub : moteur absent — assignment dans GLOBALENV (l'environnement lexical
  # du module sourcé par helper-source.R) ; test_that évaluerait dans une
  # closure locale invisible du module. Restauration à la sortie.
  stub_orig <- if (exists("sc_velocity_scvelo_available", envir = globalenv()))
    get("sc_velocity_scvelo_available", envir = globalenv()) else NULL
  assign("sc_velocity_scvelo_available", function() FALSE, envir = globalenv())
  on.exit({
    if (is.null(stub_orig)) rm("sc_velocity_scvelo_available", envir = globalenv())
    else assign("sc_velocity_scvelo_available", stub_orig, envir = globalenv())
  }, add = TRUE)

  seu <- .velopt_seurat_with_umap()
  rds_path <- tempfile(fileext = ".rds")
  on.exit(unlink(rds_path), add = TRUE)
  d <- .velopt_counts()
  saveRDS(list(spliced = as.matrix(d$spliced),
               unspliced = as.matrix(d$unspliced)), rds_path)

  gd <- shiny::reactiveValues(sc_obj = seu, language = "fr", i18n = NULL)
  shared <- shiny::reactiveValues()
  shiny::testServer(mod_sc_velocity_server, args = list(global_data = gd,
                                                        shared_rv = shared), {
    session$setInputs(
      velocity_import_mode = "rds",
      velocity_rds_file = data.frame(name = "vel.rds", datapath = rds_path,
                                     stringsAsFactors = FALSE),
      velocity_orientation = "auto_strict",
      velocity_strip_cell_suffix = FALSE,
      velocity_strip_gene_version = FALSE,
      velocity_allow_low_overlap = FALSE,
      velocity_reduction = "umap"
    )
    session$setInputs(velocity_validate = 1L)
    session$flushReact()
    expect_identical(velocity_state$result$status, "valid_no_vectors")

    session$setInputs(velocity_scvelo_run = 1L)
    session$flushReact()

    # RED (avant fix) : l'échec ne laissait AUCUNE trace persistante
    expect_false(is.null(velocity_state$optin_error))
    expect_true(is.character(velocity_state$optin_error$message) &&
                  nzchar(velocity_state$optin_error$message))
    expect_match(velocity_status_rv() %||% "", "refus", fixed = TRUE)
    # l'état consommé par Communication ne prétend toujours pas aux vecteurs
    expect_identical(shared$velocity_result$status, "valid_no_vectors")
  })
})

# ---------------------------------------------------------------------------
# (C) Bout-en-bout avec le moteur scVelo RÉEL — publication du résultat
# ---------------------------------------------------------------------------
test_that("P0.4 lane1 : opt-in réel publie les vecteurs (status valid + shared_rv)", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("Seurat")
  skip_if_not(.velopt_python_available(), "venv python_env_sccoda absent")
  source_project_file("tools/pkg_skeleton/R/io_helpers.R")
  source_project_file("tools/pkg_skeleton/R/validation.R")
  source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
  source_project_file("modules/sc/mod_sc_velocity.R")

  seu <- .velopt_seurat_with_umap()
  rds_path <- tempfile(fileext = ".rds")
  on.exit(unlink(rds_path), add = TRUE)
  d <- .velopt_counts()
  saveRDS(list(spliced = as.matrix(d$spliced),
               unspliced = as.matrix(d$unspliced)), rds_path)

  gd <- shiny::reactiveValues(sc_obj = seu, language = "fr", i18n = NULL)
  shared <- shiny::reactiveValues()
  shiny::testServer(mod_sc_velocity_server, args = list(global_data = gd,
                                                        shared_rv = shared), {
    session$setInputs(
      velocity_import_mode = "rds",
      velocity_rds_file = data.frame(name = "vel.rds", datapath = rds_path,
                                     stringsAsFactors = FALSE),
      velocity_orientation = "auto_strict",
      velocity_strip_cell_suffix = FALSE,
      velocity_strip_gene_version = FALSE,
      velocity_allow_low_overlap = FALSE,
      velocity_reduction = "umap"
    )
    session$setInputs(velocity_validate = 1L)
    session$flushReact()
    expect_identical(velocity_state$result$status, "valid_no_vectors")

    session$setInputs(velocity_scvelo_run = 1L)
    session$flushReact()

    res <- shared$velocity_result
    expect_false(is.null(res))
    expect_identical(res$status, "valid")
    expect_false(is.null(res$velocity_vectors))
    expect_identical(rownames(res$velocity_vectors),
                     as.character(velocity_state$validated$cell_names))
    # l'état moteur est resté propre
    expect_null(velocity_state$optin_error)
  })
})
