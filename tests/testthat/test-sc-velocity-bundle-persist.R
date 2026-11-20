# =============================================================================
# test-sc-velocity-bundle-persist.R — W1.1 #14 (§2c)
# =============================================================================
# Velocity bundle persistence & fingerprint rehydration across sc_obj updates:
#   - When velocity result/vectors are imported or computed, cache a bundle
#     (fingerprint, result, umap_vectors, panel_state) outside velocity_state.
#   - When global_data$sc_obj changes to an object with a DIFFERENT fingerprint,
#     clear velocity_state and shared_rv$velocity_panel_state, while keeping the
#     cached velocity_bundle for checkpoint restore.
#   - When global_data$sc_obj is updated or restored to an object with a MATCHING
#     fingerprint, rehydrate velocity_state (result, umap_vectors,
#     object_fingerprint) and shared_rv$velocity_panel_state instead of wiping.
# =============================================================================

skip_if_not_installed("Matrix")
skip_if_not_installed("SeuratObject")

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/validation.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/provenance.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
source_project_file("modules/sc/mod_sc_velocity.R")

.vbp_make_seurat <- function(genes, cells) {
  m <- Matrix::Matrix(
    seq_len(length(genes) * length(cells)),
    nrow = length(genes),
    ncol = length(cells),
    sparse = TRUE,
    dimnames = list(genes, cells)
  )
  SeuratObject::CreateSeuratObject(counts = m)
}

.vbp_make_panel <- function(genes, cells, fp) {
  vec <- matrix(
    seq_len(length(cells) * 2L),
    nrow = length(cells),
    ncol = 2L,
    dimnames = list(cells, c("dx", "dy"))
  )
  mat <- matrix(1, nrow = length(genes), ncol = length(cells),
                dimnames = list(genes, cells))
  canon <- list(
    status = "valid",
    cell_names = cells,
    gene_names = genes,
    velocity_vectors = vec,
    params = list(reduction = "umap"),
    object_identity = list(fingerprint = fp),
    provenance = list(list(analysis_id = "sc-velocity", method = "velocity"))
  )
  list(
    canonical = canon,
    vectors = vec,
    panel_state = list(
      result              = canon,
      validated           = list(cell_names = cells, spliced = mat, unspliced = mat),
      alignment           = list(cells = cells, genes = genes),
      umap_vectors        = vec,
      partial_embedding   = NULL,
      embedding_alignment = NULL,
      vector_validation   = list(field = "vectors", ok = TRUE),
      object_fingerprint  = fp,
      status_msg          = "Validation OK",
      source_label        = "velocity.rds (rds)"
    )
  )
}

test_that("bundle is cached outside velocity_state and survives differing sc_obj reset", {
  genes <- paste0("gene", seq_len(5))
  cells_a <- paste0("cell", seq_len(4))
  cells_b <- paste0("cell", seq_len(3))
  obj_a <- .vbp_make_seurat(genes, cells_a)
  obj_b <- .vbp_make_seurat(genes, cells_b)
  fp_a <- velocity_object_fingerprint(obj_a)
  fix <- .vbp_make_panel(genes, cells_a, fp_a)

  gd <- shiny::reactiveValues(sc_obj = obj_a, language = "fr", i18n = NULL)
  shared <- shiny::reactiveValues()

  shiny::testServer(mod_sc_velocity_server,
                    args = list(global_data = gd, shared_rv = shared), {
    shared$velocity_restore <- list(
      velocity_result = fix$canonical,
      panel_state     = fix$panel_state
    )
    session$flushReact()

    # Bundle cached on shared_rv with matching fingerprint
    expect_false(is.null(shared$velocity_bundle))
    expect_identical(shared$velocity_bundle$fingerprint, fp_a)
    expect_identical(shared$velocity_bundle$result$status, "valid")
    expect_identical(shared$velocity_bundle$umap_vectors, fix$vectors)
    expect_true(is.list(shared$velocity_bundle$panel_state))

    # Switch to a different Seurat object -> velocity_state cleared, bundle kept
    gd$sc_obj <- obj_b
    session$flushReact()

    expect_null(velocity_state$result)
    expect_null(velocity_state$umap_vectors)
    expect_null(velocity_state$object_fingerprint)
    expect_null(shared$velocity_panel_state)

    expect_false(is.null(shared$velocity_bundle))
    expect_identical(shared$velocity_bundle$fingerprint, fp_a)
  })
})

test_that("checkpoint restore to matching fingerprint rehydrates velocity_state and panel_state", {
  genes <- paste0("gene", seq_len(5))
  cells_a <- paste0("cell", seq_len(4))
  cells_b <- paste0("cell", seq_len(3))
  obj_a <- .vbp_make_seurat(genes, cells_a)
  obj_a_meta <- obj_a
  obj_a_meta$cluster <- "C1"
  obj_b <- .vbp_make_seurat(genes, cells_b)
  fp_a <- velocity_object_fingerprint(obj_a)
  fix <- .vbp_make_panel(genes, cells_a, fp_a)

  gd <- shiny::reactiveValues(sc_obj = obj_a, language = "fr", i18n = NULL)
  shared <- shiny::reactiveValues()

  shiny::testServer(mod_sc_velocity_server,
                    args = list(global_data = gd, shared_rv = shared), {
    shared$velocity_restore <- list(
      velocity_result = fix$canonical,
      panel_state     = fix$panel_state
    )
    session$flushReact()

    # Switch away to obj_b (wipes velocity_state)
    gd$sc_obj <- obj_b
    session$flushReact()
    expect_null(velocity_state$result)
    expect_null(shared$velocity_panel_state)

    # Restore checkpoint with same fingerprint (obj_a_meta) -> rehydrates!
    gd$sc_obj <- obj_a_meta
    session$flushReact()

    expect_identical(velocity_state$object_fingerprint, fp_a)
    expect_identical(velocity_state$result$status, "valid")
    expect_identical(velocity_state$umap_vectors, fix$vectors)
    expect_false(is.null(shared$velocity_panel_state))
    expect_identical(shared$velocity_panel_state$object_fingerprint, fp_a)
    expect_identical(shared$velocity_result$status, "valid")
    expect_match(output$velocity_vector_status, "vecteurs valides")
  })
})

test_that("metadata update with same fingerprint preserves velocity_state and bundle", {
  genes <- paste0("gene", seq_len(5))
  cells_a <- paste0("cell", seq_len(4))
  obj_a <- .vbp_make_seurat(genes, cells_a)
  obj_a_meta <- obj_a
  obj_a_meta$celltype <- "T_cell"
  fp_a <- velocity_object_fingerprint(obj_a)
  fix <- .vbp_make_panel(genes, cells_a, fp_a)

  gd <- shiny::reactiveValues(sc_obj = obj_a, language = "fr", i18n = NULL)
  shared <- shiny::reactiveValues()

  shiny::testServer(mod_sc_velocity_server,
                    args = list(global_data = gd, shared_rv = shared), {
    shared$velocity_restore <- list(
      velocity_result = fix$canonical,
      panel_state     = fix$panel_state
    )
    session$flushReact()

    gd$sc_obj <- obj_a_meta
    session$flushReact()

    expect_identical(velocity_state$object_fingerprint, fp_a)
    expect_identical(velocity_state$result$status, "valid")
    expect_identical(velocity_state$umap_vectors, fix$vectors)
    expect_false(is.null(shared$velocity_panel_state))
    expect_identical(shared$velocity_bundle$fingerprint, fp_a)
  })
})

test_that("module-level bundle rehydrates even when shared_rv is NULL", {
  genes <- paste0("gene", seq_len(5))
  cells_a <- paste0("cell", seq_len(4))
  cells_b <- paste0("cell", seq_len(3))
  obj_a <- .vbp_make_seurat(genes, cells_a)
  obj_b <- .vbp_make_seurat(genes, cells_b)
  fp_a <- velocity_object_fingerprint(obj_a)
  fix <- .vbp_make_panel(genes, cells_a, fp_a)

  gd <- shiny::reactiveValues(sc_obj = obj_a, language = "fr", i18n = NULL)

  shiny::testServer(mod_sc_velocity_server,
                    args = list(global_data = gd, shared_rv = NULL), {
    velocity_state$validated <- fix$panel_state$validated
    velocity_state$result <- fix$canonical
    velocity_state$umap_vectors <- fix$vectors
    velocity_state$vector_validation <- fix$panel_state$vector_validation
    velocity_state$object_fingerprint <- fp_a
    velocity_state$status_msg <- "Validation OK"
    velocity_state$source_label <- "velocity.rds (rds)"
    .velocity_publish_state()
    session$flushReact()

    gd$sc_obj <- obj_b
    session$flushReact()
    expect_null(velocity_state$result)

    gd$sc_obj <- obj_a
    session$flushReact()
    expect_identical(velocity_state$object_fingerprint, fp_a)
    expect_identical(velocity_state$result$status, "valid")
    expect_identical(velocity_state$umap_vectors, fix$vectors)
  })
})
