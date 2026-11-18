# =============================================================================
# test-sc-gse96583-3x2-chain.R — W3.3 real 3x2 GSE96583 lifecycle + DA chain
# + post-4584db4/125814d live mirai daemon UTF-8 IPC check
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_multi.R")
source_project_file("tools/pkg_skeleton/R/sc_abundance_design.R")
source_project_file("tools/pkg_skeleton/R/sc_abundance_milo.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
source_project_file("tools/pkg_skeleton/R/spatial_async.R")

.gse96583_by_sample_root <- function() {
  file.path(
    dirname(ts_project_root()),
    "transcriptoshiny_test_datasets",
    "single_cell",
    "GSE96583_human_PBMC_SLE_8donor_2cond",
    "converted",
    "by_sample"
  )
}

test_that("W3.3 : GSE96583 3x2 lifecycle (6 triplets -> sc_multi -> merge -> DA design -> Milo)", {
  root <- .gse96583_by_sample_root()
  skip_if_not(dir.exists(root), message = "GSE96583 by_sample directory absent")
  skip_if_not_installed("Seurat")
  suppressPackageStartupMessages(library(Seurat))

  samples <- c("1015_ctrl", "1015_stim", "1256_ctrl", "1256_stim", "1488_ctrl", "1488_stim")
  sample_dirs <- file.path(root, samples)
  skip_if_not(all(dir.exists(sample_dirs)), message = "All 6 3x2 GSE96583 folders required")

  datasets <- list()
  objs <- list()
  withr::with_seed(20261118L, {
    for (i in seq_along(samples)) {
      s <- samples[i]
      d <- sample_dirs[i]
      mat <- suppressWarnings(Seurat::Read10X(data.dir = d))
      # Deterministic 60-cell stratified subsample per replicate for fast CI/BOOST run
      keep_cells <- colnames(mat)[seq_len(min(60L, ncol(mat)))]
      sub_mat <- mat[, keep_cells, drop = FALSE]
      obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = sub_mat, project = s))
      obj$orig.ident <- s

      meta_csv <- utils::read.csv(file.path(d, "metadata.csv"), stringsAsFactors = FALSE)
      m <- match(colnames(obj), meta_csv$barcode)
      obj$donor     <- as.character(meta_csv$donor[m])
      obj$condition <- as.character(meta_csv$condition[m])
      obj$celltype  <- as.character(meta_csv$celltype[m])
      obj$multiplets <- as.character(meta_csv$multiplets[m])

      datasets <- sc_multi_register(datasets, label = s, obj = obj,
                                    relation = "standalone", producer = "import")
      objs[[s]] <- obj
    }
  })

  sm <- sc_multi_summary(datasets)
  expect_identical(nrow(sm), 6L)
  expect_setequal(sm$label, samples)

  merged <- suppressWarnings(merge(objs[[1L]], y = objs[-1L], add.cell.ids = samples))
  merged <- SeuratObject::JoinLayers(merged)
  expect_identical(length(unique(merged$orig.ident)), 6L)
  expect_true(all(c("control", "stim") %in% unique(merged$condition)))
  expect_true("doublet" %in% unique(merged$multiplets))

  # Normalize + PCA for Milo
  merged <- suppressWarnings({
    o <- Seurat::NormalizeData(merged, verbose = FALSE)
    o <- Seurat::FindVariableFeatures(o, nfeatures = 300L, verbose = FALSE)
    o <- Seurat::ScaleData(o, features = Seurat::VariableFeatures(o), verbose = FALSE)
    Seurat::RunPCA(o, features = Seurat::VariableFeatures(o), npcs = 15L, verbose = FALSE)
  })

  val <- validate_da_design(
    merged@meta.data,
    sample_id = "orig.ident",
    condition = "condition",
    replicate_id = "donor",
    identity = "celltype"
  )
  da_design <- finalize_da_design_result(val, seurat_obj = merged)
  expect_true(da_design$status %in% c("valid", "valid_with_warnings"))
  expect_identical(nrow(da_design$sample_summary), 6L)
  expect_identical(nrow(da_design$condition_summary), 2L)

  if (isTRUE(tryCatch(milo_available(), error = function(e) FALSE))) {
    milo_res <- suppressWarnings(run_milo_da(
      merged,
      da_design_result = da_design,
      reduction = "pca",
      target_condition = "stim",
      reference_condition = "control",
      use_batch = FALSE,
      identity_col = "celltype",
      k = 10L,
      prop = 0.2,
      d = 10L
    ))
    expect_no_error(assert_milo_result(milo_res, seurat_obj = merged))
    expect_true(nrow(milo_res$da_table) >= 1L)
  }
})

test_that("Post-4584db4 / 125814d : daemon mirai preload UTF-8 retourne un label valide sans NUL embarqué", {
  skip_if_not_installed("mirai")
  on.exit(stop_spatial_daemons(), add = TRUE)
  init_spatial_daemons(
    n_daemons = 1L,
    source_files = c("tools/pkg_skeleton/R/sc_velocity.R")
  )
  m <- mirai::mirai(list(
    banner = sc_velocity_optin_banner_label(),
    utf8   = isTRUE(tryCatch(l10n_info()$`UTF-8`, error = function(e) FALSE))
  ))
  res <- mirai::call_mirai(m)$data
  expect_false(inherits(res, "miraiError"))
  expect_true(is.character(res$banner) && nzchar(res$banner))
  expect_true(validUTF8(res$banner))
  expect_false(any(charToRaw(res$banner) == as.raw(0L)))
  expect_true(isTRUE(res$utf8))
})
