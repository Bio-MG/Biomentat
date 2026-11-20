# =============================================================================
# test-ccc-lr-expression.R — W2.2 §11.8 Gem #3 : vérification descriptive de
# l'expression single-cell L-R (build_communication_lr_expression +
# plot_communication_lr_expression)
# =============================================================================

.lr_make_mock_obj <- function() {
  # 6 cellules : 3 CD4 T (sender), 3 B (receiver)
  # Gènes : LIG1, LIG2, REC1, REC2, ZERO_G
  mat <- matrix(
    c(
      # C1   C2   C3   C4   C5   C6
      2.0, 4.0, 0.0, 0.0, 0.0, 0.0, # LIG1 : CD4 T 2/3 (66.67%), mean 2.0, nz 3.0 ; B 0/3 (0%), mean 0, nz NA
      1.0, 1.0, 1.0, 3.0, 0.0, 0.0, # LIG2 : CD4 T 3/3 (100%), mean 1.0, nz 1.0 ; B 1/3 (33.33%), mean 1.0, nz 3.0
      0.0, 0.0, 0.0, 2.0, 6.0, 0.0, # REC1 : CD4 T 0/3 (0%) ; B 2/3 (66.67%), mean 8/3, nz 4.0
      0.0, 0.0, 0.0, 4.0, 4.0, 4.0, # REC2 : CD4 T 0/3 (0%) ; B 3/3 (100%), mean 4.0, nz 4.0
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0  # ZEROG : 0 partout
    ),
    nrow = 5L,
    ncol = 6L,
    byrow = TRUE,
    dimnames = list(
      c("LIG1", "LIG2", "REC1", "REC2", "ZEROG"),
      paste0("Cell_", 1:6)
    )
  )
  meta <- data.frame(
    cell_type = c("CD4 T", "CD4 T", "CD4 T", "B", "B", "B"),
    row.names = paste0("Cell_", 1:6),
    stringsAsFactors = FALSE
  )
  list(expr = mat, meta.data = meta)
}

test_that("build_communication_lr_expression and plot_communication_lr_expression are in communication_views_public_api()", {
  api <- communication_views_public_api()
  expect_true("build_communication_lr_expression" %in% api)
  expect_true("plot_communication_lr_expression" %in% api)
})

test_that("build_communication_lr_expression validates inputs with communication_import_error", {
  obj <- .lr_make_mock_obj()

  expect_error(
    build_communication_lr_expression(NULL, "CD4 T", "B", "LIG1", "REC1", "cell_type"),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_lr_expression("not_an_obj", "CD4 T", "B", "LIG1", "REC1", "cell_type"),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_lr_expression(obj, "CD4 T", "B", "LIG1", "REC1", "missing_col"),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_lr_expression(obj, "", "B", "LIG1", "REC1", "cell_type"),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_lr_expression(obj, "CD4 T", NA_character_, "LIG1", "REC1", "cell_type"),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_lr_expression(obj, "CD4 T", "B", "___", "REC1", "cell_type"),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_lr_expression(obj, "CD4 T", "B", "LIG1", character(0), "cell_type"),
    class = "communication_import_error"
  )
})

test_that("build_communication_lr_expression computes exact single-gene metrics and epistemic note", {
  obj <- .lr_make_mock_obj()
  res <- build_communication_lr_expression(
    seurat_obj = obj,
    sender = "CD4 T",
    receiver = "B",
    ligand = "LIG1",
    receptor = "REC1",
    identity_column = "cell_type"
  )

  expect_true(is.list(res))
  expect_true(is.data.frame(res$table))
  expect_identical(
    colnames(res$table),
    c("gene", "role", "population", "n_cells", "pct_expressing",
      "mean_expr", "mean_expr_nonzero", "gene_found")
  )
  expect_identical(nrow(res$table), 4L)
  expect_identical(
    res$epistemic_note,
    "Vérification descriptive d'expression L-R — la co-expression ne prouve pas une interaction physique ni causale"
  )

  # LIG1 chez Sender (CD4 T): c(2, 4, 0)
  row_lig_s <- res$table[res$table$gene == "LIG1" & res$table$population == "Sender", ]
  expect_identical(row_lig_s$role, "Ligand")
  expect_identical(row_lig_s$n_cells, 3L)
  expect_true(row_lig_s$gene_found)
  expect_equal(row_lig_s$pct_expressing, 100 * 2 / 3, tolerance = 1e-6)
  expect_equal(row_lig_s$mean_expr, 2.0, tolerance = 1e-6)
  expect_equal(row_lig_s$mean_expr_nonzero, 3.0, tolerance = 1e-6)

  # LIG1 chez Receiver (B): c(0, 0, 0) -> pct = 0, mean = 0, mean_nonzero = NA_real_
  row_lig_r <- res$table[res$table$gene == "LIG1" & res$table$population == "Receiver", ]
  expect_equal(row_lig_r$pct_expressing, 0.0)
  expect_equal(row_lig_r$mean_expr, 0.0)
  expect_true(is.na(row_lig_r$mean_expr_nonzero))

  # REC1 chez Receiver (B): c(2, 6, 0)
  row_rec_r <- res$table[res$table$gene == "REC1" & res$table$population == "Receiver", ]
  expect_identical(row_rec_r$role, "Receptor")
  expect_equal(row_rec_r$pct_expressing, 100 * 2 / 3, tolerance = 1e-6)
  expect_equal(row_rec_r$mean_expr, 8 / 3, tolerance = 1e-6)
  expect_equal(row_rec_r$mean_expr_nonzero, 4.0, tolerance = 1e-6)

  expect_identical(res$qc$n_sender_cells, 3L)
  expect_identical(res$qc$n_receiver_cells, 3L)
  expect_identical(res$qc$n_ligand_subunits, 1L)
  expect_identical(res$qc$n_receptor_subunits, 1L)
  expect_identical(res$qc$n_genes_found, 2L)
  expect_identical(res$qc$n_genes_missing, 0L)
  expect_true(res$qc$all_receptor_subunits_expressed)
})

test_that("build_communication_lr_expression parses multi-subunit complexes deterministically", {
  obj <- .lr_make_mock_obj()
  res <- build_communication_lr_expression(
    seurat_obj = obj,
    sender = "CD4 T",
    receiver = "B",
    ligand = " LIG1 + LIG2 , LIG1 ",
    receptor = "REC1_REC2",
    identity_column = "cell_type"
  )

  expect_identical(res$params$ligand_subunits, c("LIG1", "LIG2"))
  expect_identical(res$params$receptor_subunits, c("REC1", "REC2"))
  expect_identical(res$qc$n_ligand_subunits, 2L)
  expect_identical(res$qc$n_receptor_subunits, 2L)
  expect_identical(nrow(res$table), 8L)
  expect_true(res$qc$all_receptor_subunits_expressed)

  # Si un récepteur est non exprimé (ZEROG), all_receptor_subunits_expressed == FALSE
  res_zero <- build_communication_lr_expression(
    seurat_obj = obj,
    sender = "CD4 T",
    receiver = "B",
    ligand = "LIG1",
    receptor = "REC1:ZEROG",
    identity_column = "cell_type"
  )
  expect_false(res_zero$qc$all_receptor_subunits_expressed)
})

test_that("build_communication_lr_expression handles autocrine (sender == receiver)", {
  obj <- .lr_make_mock_obj()
  res <- build_communication_lr_expression(
    seurat_obj = obj,
    sender = "CD4 T",
    receiver = "CD4 T",
    ligand = "LIG1",
    receptor = "LIG2",
    identity_column = "cell_type"
  )

  expect_true(res$params$is_autocrine)
  expect_identical(res$qc$n_sender_cells, 3L)
  expect_identical(res$qc$n_receiver_cells, 3L)
  expect_identical(nrow(res$table), 4L)
  s_lig1 <- res$table[res$table$gene == "LIG1" & res$table$population == "Sender", ]
  r_lig1 <- res$table[res$table$gene == "LIG1" & res$table$population == "Receiver", ]
  expect_equal(s_lig1$mean_expr, r_lig1$mean_expr)
})

test_that("build_communication_lr_expression never fabricates values for missing genes", {
  obj <- .lr_make_mock_obj()
  res <- build_communication_lr_expression(
    seurat_obj = obj,
    sender = "CD4 T",
    receiver = "B",
    ligand = "LIG1_MISSINGL",
    receptor = "REC1_MISSINGR",
    identity_column = "cell_type"
  )

  expect_identical(res$qc$n_genes_found, 2L)
  expect_identical(res$qc$n_genes_missing, 2L)
  expect_false(res$qc$all_receptor_subunits_expressed)

  miss <- res$table[res$table$gene %in% c("MISSINGL", "MISSINGR"), ]
  expect_identical(nrow(miss), 4L)
  expect_true(all(!miss$gene_found))
  expect_true(all(is.na(miss$pct_expressing)))
  expect_true(all(is.na(miss$mean_expr)))
  expect_true(all(is.na(miss$mean_expr_nonzero)))
})

test_that("build_communication_lr_expression extracts expression from a real Seurat object", {
  skip_if_not_installed("SeuratObject")
  obj_mock <- .lr_make_mock_obj()
  sp_mat <- Matrix::Matrix(obj_mock$expr, sparse = TRUE)
  so <- SeuratObject::CreateSeuratObject(
    counts = sp_mat,
    meta.data = obj_mock$meta.data
  )
  # Peuple la couche data explicitement (compatible Seurat v4/v5)
  so <- tryCatch(
    SeuratObject::SetAssayData(so, layer = "data", new.data = sp_mat),
    error = function(e) SeuratObject::SetAssayData(so, slot = "data", new.data = sp_mat)
  )

  res <- build_communication_lr_expression(
    seurat_obj = so,
    sender = "CD4 T",
    receiver = "B",
    ligand = "LIG1",
    receptor = "REC1_REC2",
    identity_column = "cell_type",
    assay = "RNA",
    layer = "data"
  )

  expect_identical(res$qc$n_sender_cells, 3L)
  expect_identical(res$qc$n_receiver_cells, 3L)
  expect_identical(res$qc$n_genes_found, 3L)
  expect_true(res$qc$all_receptor_subunits_expressed)
  expect_identical(res$params$assay, "RNA")
})

test_that("plot_communication_lr_expression renders ggplot DotPlot and empty-state fallback", {
  expect_error(
    plot_communication_lr_expression(NULL),
    class = "communication_import_error"
  )

  obj <- .lr_make_mock_obj()
  res <- build_communication_lr_expression(
    seurat_obj = obj,
    sender = "CD4 T",
    receiver = "B",
    ligand = "LIG1_LIG2",
    receptor = "REC1_REC2",
    identity_column = "cell_type"
  )
  p <- plot_communication_lr_expression(res)
  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 1L)
  expect_match(p$labels$caption %||% "", "co-expression ne prouve pas", fixed = TRUE)

  # Tous gènes absents -> fallback .communication_empty_view_plot()
  res_missing <- build_communication_lr_expression(
    seurat_obj = obj,
    sender = "CD4 T",
    receiver = "B",
    ligand = "NO_SUCH_L",
    receptor = "NO_SUCH_R",
    identity_column = "cell_type"
  )
  p_empty <- plot_communication_lr_expression(res_missing)
  expect_s3_class(p_empty, "ggplot")
  expect_match(p_empty$labels$subtitle %||% "", "Aucun", fixed = TRUE)
})
