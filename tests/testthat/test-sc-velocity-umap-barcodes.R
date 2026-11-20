# =============================================================================
# test-sc-velocity-umap-barcodes.R — W1.1 items #3 & #4
#   - Item #4 : normalisation prefixes/suffixes barcodes + indice mismatch
#   - Item #3 : garde FR exploitable si reduction 'umap' absente de Seurat
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_velocity.R")

test_that("normalize_velocity_cell_barcodes strips colon and underscore prefixes when strip_prefix = TRUE", {
  ids <- c("WB:AAACCCGGGTTT", "S1_AAAGGGCCCTTT", "S1_cell1")
  # Default strip_prefix = FALSE preserves raw IDs
  expect_identical(normalize_velocity_cell_barcodes(ids), ids)

  res <- normalize_velocity_cell_barcodes(ids, strip_prefix = TRUE)
  expect_identical(res, c("AAACCCGGGTTT", "AAAGGGCCCTTT", "S1_cell1"))
})

test_that("normalize_velocity_cell_barcodes strips velocyto trailing x with strip_suffix = TRUE", {
  ids <- c("AAACCCGGGTTT-1x", "AAAGGGCCCTTTx", "AAATTTTTTTTT-2", "AAACCCAAAGGG")
  res <- normalize_velocity_cell_barcodes(ids, strip_suffix = TRUE)
  expect_identical(res, c("AAACCCGGGTTT", "AAAGGGCCCTTT", "AAATTTTTTTTT", "AAACCCAAAGGG"))
})

test_that("normalize_velocity_cell_barcodes strips both prefix and suffix together", {
  ids <- c("WB:AAACCCGGGTTT-1x", "S2_AAAGGGCCCTTT-1")
  res <- normalize_velocity_cell_barcodes(ids, strip_suffix = TRUE, strip_prefix = TRUE)
  expect_identical(res, c("AAACCCGGGTTT", "AAAGGGCCCTTT"))
})

test_that("normalize_velocity_cell_barcodes detects collisions when strip_prefix = TRUE", {
  ids <- c("S1:AAACCCGGGTTT", "S2:AAACCCGGGTTT")
  err <- expect_error(
    normalize_velocity_cell_barcodes(ids, strip_prefix = TRUE),
    "collision"
  )
  expect_true(inherits(err, "velocity_validation_error"))
  expect_identical(err$state, "invalid_input")
})

test_that("validate_velocity_matrices aligns prefixed/suffixed barcodes when strip_cell_prefix = TRUE", {
  skip_if_not_installed("Matrix")
  genes <- paste0("G", 1:4)
  vel_cells <- c("WB:AAACCCGGGTTT-1x", "WB:AAAGGGCCCTTT-1x", "WB:AAATTTTTTTTT-1x")
  seu_cells <- c("S1_AAACCCGGGTTT-1", "S1_AAAGGGCCCTTT-1", "S1_AAATTTTTTTTT-1")

  sp <- Matrix::Matrix(1, nrow = 4, ncol = 3, sparse = TRUE,
                       dimnames = list(genes, vel_cells))
  un <- Matrix::Matrix(2, nrow = 4, ncol = 3, sparse = TRUE,
                       dimnames = list(genes, vel_cells))

  res <- validate_velocity_matrices(
    spliced = sp,
    unspliced = un,
    seurat_cells = seu_cells,
    seurat_genes = genes,
    orientation = "genes_x_cells",
    strip_cell_suffix = TRUE,
    strip_cell_prefix = TRUE
  )
  expect_identical(res$n_cells_matched, 3L)
  expect_identical(res$cell_names, seu_cells)
  expect_true(isTRUE(res$strip_cell_prefix))
})

test_that("validate_velocity_matrices gives actionable prefix/suffix hint on 0 raw matches", {
  skip_if_not_installed("Matrix")
  genes <- paste0("G", 1:4)
  vel_cells <- c("WB:AAACCCGGGTTT-1x", "WB:AAAGGGCCCTTT-1x", "WB:AAATTTTTTTTT-1x")
  seu_cells <- c("S1_AAACCCGGGTTT-1", "S1_AAAGGGCCCTTT-1", "S1_AAATTTTTTTTT-1")

  sp <- Matrix::Matrix(1, nrow = 4, ncol = 3, sparse = TRUE,
                       dimnames = list(genes, vel_cells))
  un <- Matrix::Matrix(2, nrow = 4, ncol = 3, sparse = TRUE,
                       dimnames = list(genes, vel_cells))

  err <- expect_error(
    validate_velocity_matrices(
      spliced = sp,
      unspliced = un,
      seurat_cells = seu_cells,
      seurat_genes = genes,
      orientation = "genes_x_cells",
      strip_cell_suffix = FALSE,
      strip_cell_prefix = FALSE
    ),
    "3 barcode\\(s\\) correspondent"
  )
  expect_identical(err$state, "invalid_cell_alignment")
  expect_match(conditionMessage(err), "prefixes/suffixes", fixed = TRUE)
})

test_that("mod_sc_velocity_server emits actionable FR UMAP guard when umap reduction is absent", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("Seurat")
  source_project_file("tools/pkg_skeleton/R/io_helpers.R")
  source_project_file("tools/pkg_skeleton/R/validation.R")
  source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
  source_project_file("modules/sc/mod_sc_velocity.R")

  genes <- paste0("G", 1:4)
  cells <- paste0("Cell", 1:6)
  set.seed(1)
  sp <- matrix(rpois(24, 5), nrow = 4, ncol = 6, dimnames = list(genes, cells))
  un <- matrix(rpois(24, 2), nrow = 4, ncol = 6, dimnames = list(genes, cells))
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = sp))
  expect_false("umap" %in% Seurat::Reductions(seu))

  rds_path <- tempfile(fileext = ".rds")
  on.exit(unlink(rds_path), add = TRUE)
  saveRDS(list(spliced = sp, unspliced = un), rds_path)

  gd <- shiny::reactiveValues(sc_obj = seu, language = "fr", i18n = NULL)
  srv <- shiny::reactiveValues()
  expected_msg <- "Réduction 'umap' absente de l'objet Seurat — exécutez d'abord le Pipeline SC (PCA / UMAP) avant de lancer la vélocité."

  shiny::testServer(mod_sc_velocity_server, args = list(global_data = gd, shared_rv = srv), {
    session$setInputs(
      velocity_import_mode = "rds",
      velocity_rds_file = data.frame(name = "vel.rds", datapath = rds_path, stringsAsFactors = FALSE),
      velocity_orientation = "auto_strict",
      velocity_strip_cell_suffix = FALSE,
      velocity_strip_gene_version = FALSE,
      velocity_allow_low_overlap = FALSE,
      velocity_reduction = "umap"
    )
    session$setInputs(velocity_validate = 1L)
    expect_false(is.null(velocity_state$result))
    expect_true(any(grepl(expected_msg, velocity_state$result$warnings, fixed = TRUE)))

    session$setInputs(velocity_scvelo_run = 1L)
    expect_match(velocity_status_rv() %||% "", expected_msg, fixed = TRUE)
    expect_identical(velocity_state$optin_error$message, expected_msg)
  })
})
