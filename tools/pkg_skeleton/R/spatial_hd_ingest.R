# =============================================================================
# tools/pkg_skeleton/R/spatial_hd_ingest.R — Spatial HD / MERFISH Ingest (BK-10)
# =============================================================================
# Pure domain logic (zero Shiny symbols). Supports:
#   1. Platform detection for Visium HD (binned_outputs / tissue_positions) and
#      MERSCOPE / MERFISH (cell_by_gene.csv + cell_metadata.csv).
#   2. Sparse ingestion of MERSCOPE/MERFISH cell-by-gene + coordinate metadata
#      with optional Blank control probe filtering.
#   3. Deterministic sparse bin pooling for Visium HD (e.g., 2um -> 8um) via
#      sparse indicator matrix multiplication without dense materialization.
#
# Dependencies: base R, utils, stats, Matrix.
# =============================================================================

#' Internal typed error constructor for spatial HD ingestion
#'
#' @param msg Character scalar error message.
#' @return An error condition of class `"spatial_hd_ingest_error"`.
#' @keywords internal
.spatial_hd_stop <- function(msg) {
  errorCondition(msg, class = "spatial_hd_ingest_error")
}

#' Detect high-resolution spatial transcriptomics directory format
#'
#' @param dir_path Character path to an input spatial dataset directory.
#' @return Named list with `platform` (`"merscope_merfish"`, `"visium_hd"`, or
#'   `"unknown"`), `bin_sizes` (character vector for Visium HD), and `files`.
#' @export
detect_spatial_hd_format <- function(dir_path) {
  if (is.null(dir_path) || !is.character(dir_path) || length(dir_path) != 1L ||
      is.na(dir_path) || !nzchar(trimws(dir_path)) || !dir.exists(dir_path)) {
    stop(.spatial_hd_stop("'dir_path' must point to an existing directory."))
  }
  dir_norm <- normalizePath(dir_path, winslash = "/", mustWork = FALSE)

  cbg <- file.path(dir_norm, "cell_by_gene.csv")
  cmeta <- file.path(dir_norm, "cell_metadata.csv")
  if (file.exists(cbg) && file.exists(cmeta)) {
    return(list(
      platform = "merscope_merfish",
      bin_sizes = character(0L),
      files = list(counts = cbg, metadata = cmeta)
    ))
  }

  binned_dir <- file.path(dir_norm, "binned_outputs")
  tp_pq <- file.path(dir_norm, "tissue_positions.parquet")
  tp_csv <- file.path(dir_norm, "tissue_positions.csv")
  if (dir.exists(binned_dir) || file.exists(tp_pq) || file.exists(tp_csv)) {
    bin_sizes <- if (dir.exists(binned_dir)) {
      sort(basename(list.dirs(binned_dir, full.names = TRUE, recursive = FALSE)))
    } else {
      character(0L)
    }
    return(list(
      platform = "visium_hd",
      bin_sizes = bin_sizes,
      files = list(
        binned_dir = if (dir.exists(binned_dir)) binned_dir else NA_character_,
        tissue_positions = if (file.exists(tp_pq)) tp_pq else if (file.exists(tp_csv)) tp_csv else NA_character_
      )
    ))
  }

  list(platform = "unknown", bin_sizes = character(0L), files = list())
}

#' Read MERSCOPE / MERFISH cell_by_gene.csv and cell_metadata.csv
#'
#' @param counts_path Path to `cell_by_gene.csv` (cells in rows, genes in columns).
#' @param metadata_path Path to `cell_metadata.csv` with spatial coordinates.
#' @param drop_blanks Logical; drop `^Blank[-_]` control probes when `TRUE`.
#' @param min_counts Minimum total counts per cell to retain (default `1L`).
#' @return Named list with sparse `counts` (`dgCMatrix`, genes x cells),
#'   `coords` (`data.frame(id, x, y)`), `metadata`, `n_cells`, `n_genes`,
#'   `n_blanks_removed`, `platform = "merscope_merfish"`, and `provenance`.
#' @export
read_merscope_merfish <- function(counts_path,
                                  metadata_path,
                                  drop_blanks = TRUE,
                                  min_counts = 1L) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(.spatial_hd_stop("Package 'Matrix' is required for read_merscope_merfish()."))
  }
  if (is.null(counts_path) || !file.exists(counts_path)) {
    stop(.spatial_hd_stop(sprintf("MERSCOPE/MERFISH counts file not found: %s", as.character(counts_path))))
  }
  if (is.null(metadata_path) || !file.exists(metadata_path)) {
    stop(.spatial_hd_stop(sprintf("MERSCOPE/MERFISH metadata file not found: %s", as.character(metadata_path))))
  }

  raw_cnt <- utils::read.csv(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
  raw_meta <- utils::read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(raw_cnt) < 1L || ncol(raw_cnt) < 2L) {
    stop(.spatial_hd_stop("counts_path must contain at least 1 cell row and 1 gene column."))
  }
  if (nrow(raw_meta) < 1L || ncol(raw_meta) < 2L) {
    stop(.spatial_hd_stop("metadata_path must contain at least 1 cell row and coordinate columns."))
  }

  extract_ids <- function(df) {
    id_candidates <- intersect(c("EntityID", "cell", "cell_id", "id", ""), colnames(df))
    if (length(id_candidates) > 0L) {
      col_nm <- id_candidates[1L]
      ids <- as.character(df[[col_nm]])
      df <- df[, setdiff(colnames(df), col_nm), drop = FALSE]
    } else {
      ids <- as.character(df[[1L]])
      df <- df[, -1L, drop = FALSE]
    }
    list(ids = ids, df = df)
  }

  cnt_p <- extract_ids(raw_cnt)
  meta_p <- extract_ids(raw_meta)

  cnt_ids <- cnt_p$ids
  meta_ids <- meta_p$ids
  common_ids <- intersect(cnt_ids, meta_ids)
  if (length(common_ids) < 1L) {
    stop(.spatial_hd_stop("No overlapping cell IDs found between counts_path and metadata_path."))
  }

  cnt_df <- cnt_p$df[match(common_ids, cnt_ids), , drop = FALSE]
  meta_df <- meta_p$df[match(common_ids, meta_ids), , drop = FALSE]
  rownames(meta_df) <- common_ids

  x_col <- intersect(c("center_x", "x", "X"), colnames(meta_df))
  y_col <- intersect(c("center_y", "y", "Y"), colnames(meta_df))
  if (length(x_col) == 0L || length(y_col) == 0L) {
    stop(.spatial_hd_stop("metadata_path must contain coordinate columns ('center_x'/'center_y' or 'x'/'y')."))
  }

  n_blanks_removed <- 0L
  if (isTRUE(drop_blanks)) {
    is_blank <- grepl("^Blank[-_]", colnames(cnt_df), ignore.case = TRUE)
    n_blanks_removed <- as.integer(sum(is_blank))
    if (n_blanks_removed > 0L) {
      cnt_df <- cnt_df[, !is_blank, drop = FALSE]
    }
  }
  if (ncol(cnt_df) < 1L) {
    stop(.spatial_hd_stop("No gene columns remain after blank probe filtering."))
  }

  num_mat <- as.matrix(cnt_df)
  storage.mode(num_mat) <- "double"
  mat_sparse <- methods::as(Matrix::Matrix(t(num_mat), sparse = TRUE), "dgCMatrix")
  rownames(mat_sparse) <- colnames(cnt_df)
  colnames(mat_sparse) <- common_ids

  keep_cells <- Matrix::colSums(mat_sparse) >= as.numeric(min_counts)
  mat_sparse <- mat_sparse[, keep_cells, drop = FALSE]
  meta_df <- meta_df[keep_cells, , drop = FALSE]
  kept_ids <- colnames(mat_sparse)

  coords <- data.frame(
    id = kept_ids,
    x = as.numeric(meta_df[[x_col[1L]]]),
    y = as.numeric(meta_df[[y_col[1L]]]),
    stringsAsFactors = FALSE
  )

  list(
    counts = mat_sparse,
    coords = coords,
    metadata = meta_df,
    n_cells = as.integer(ncol(mat_sparse)),
    n_genes = as.integer(nrow(mat_sparse)),
    n_blanks_removed = n_blanks_removed,
    platform = "merscope_merfish",
    provenance = list(
      counts_path = counts_path,
      metadata_path = metadata_path,
      drop_blanks = isTRUE(drop_blanks),
      min_counts = as.numeric(min_counts)
    )
  )
}

#' Pool Visium HD high-resolution bins into coarser spatial super-bins
#'
#' @param counts_mat Sparse or dense gene-by-bin count matrix (`colnames` = bin IDs).
#' @param coords Data.frame with `id`, `x`, `y` matching `colnames(counts_mat)`.
#' @param bin_factor Integer pooling factor >= 2 along each axis (default `4L`,
#'   e.g., 2 um -> 8 um).
#' @return Named list with pooled sparse `counts` (`dgCMatrix`), `coords`
#'   (`data.frame(id, x, y, n_subbins)`), `bin_factor`, `n_input_bins`,
#'   `n_output_bins`, and `platform = "visium_hd"`.
#' @export
aggregate_visium_hd_bins <- function(counts_mat, coords, bin_factor = 4L) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(.spatial_hd_stop("Package 'Matrix' is required for aggregate_visium_hd_bins()."))
  }
  if (is.null(counts_mat) || length(dim(counts_mat)) != 2L || is.null(colnames(counts_mat))) {
    stop(.spatial_hd_stop("'counts_mat' must be a 2D matrix with column names."))
  }
  if (!is.data.frame(coords) || !all(c("id", "x", "y") %in% colnames(coords))) {
    stop(.spatial_hd_stop("'coords' must be a data.frame with columns 'id', 'x', 'y'."))
  }
  bf <- suppressWarnings(as.integer(bin_factor))
  if (length(bf) != 1L || is.na(bf) || bf < 2L) {
    stop(.spatial_hd_stop("'bin_factor' must be an integer >= 2."))
  }

  idx <- match(colnames(counts_mat), coords$id)
  if (anyNA(idx)) {
    stop(.spatial_hd_stop("All columns of 'counts_mat' must be present in 'coords$id'."))
  }
  x_v <- as.numeric(coords$x[idx])
  y_v <- as.numeric(coords$y[idx])
  if (any(!is.finite(x_v)) || any(!is.finite(y_v))) {
    stop(.spatial_hd_stop("'coords' x and y must be finite numeric values."))
  }

  step_x <- stats::median(diff(sort(unique(x_v))))
  step_y <- stats::median(diff(sort(unique(y_v))))
  if (!is.finite(step_x) || step_x <= 0) step_x <- 1
  if (!is.finite(step_y) || step_y <= 0) step_y <- 1

  gx <- as.integer(floor((x_v - min(x_v)) / (step_x * bf) + 1e-9))
  gy <- as.integer(floor((y_v - min(y_v)) / (step_y * bf) + 1e-9))
  super_key <- paste0("hd_bin_", gx, "_", gy)
  super_fac <- factor(super_key, levels = unique(super_key))
  super_levels <- levels(super_fac)

  G <- Matrix::sparseMatrix(
    i = seq_along(super_fac),
    j = as.integer(super_fac),
    x = 1,
    dims = c(length(super_fac), length(super_levels))
  )
  sp_mat <- if (inherits(counts_mat, "dgCMatrix")) counts_mat else methods::as(Matrix::Matrix(counts_mat, sparse = TRUE), "dgCMatrix")
  pooled <- methods::as(sp_mat %*% G, "dgCMatrix")
  colnames(pooled) <- super_levels
  rownames(pooled) <- rownames(counts_mat)

  n_sub <- as.integer(tabulate(as.integer(super_fac), nbins = length(super_levels)))
  cx <- as.numeric(tapply(x_v, super_fac, mean))
  cy <- as.numeric(tapply(y_v, super_fac, mean))
  pooled_coords <- data.frame(
    id = super_levels,
    x = cx,
    y = cy,
    n_subbins = n_sub,
    stringsAsFactors = FALSE
  )

  list(
    counts = pooled,
    coords = pooled_coords,
    bin_factor = bf,
    n_input_bins = as.integer(ncol(counts_mat)),
    n_output_bins = as.integer(ncol(pooled)),
    platform = "visium_hd"
  )
}

#' Load a MERSCOPE / MERFISH directory into a Seurat Spatial object
#'
#' @param dir_path Directory containing `cell_by_gene.csv` and `cell_metadata.csv`.
#' @param sample_name Character sample identifier.
#' @param drop_blanks Logical; drop `^Blank[-_]` control probes when `TRUE`.
#' @param min_counts Minimum total counts per cell to retain.
#' @return A `Seurat` object with `coord_x`/`coord_y` metadata and
#'   `attr(obj, "ts_spatial_mode") == "merscope"`.
#' @export
load_spatial_merscope <- function(dir_path,
                                  sample_name = "MERSCOPE_Sample",
                                  drop_blanks = TRUE,
                                  min_counts = 1L) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop(.spatial_hd_stop("Package 'Seurat' is required for load_spatial_merscope()."))
  }
  det <- detect_spatial_hd_format(dir_path)
  if (!identical(det$platform, "merscope_merfish")) {
    stop(.spatial_hd_stop(sprintf(
      "Directory does not contain MERSCOPE/MERFISH files (cell_by_gene.csv + cell_metadata.csv): %s",
      as.character(dir_path)
    )))
  }
  res <- read_merscope_merfish(
    counts_path = det$files$counts,
    metadata_path = det$files$metadata,
    drop_blanks = isTRUE(drop_blanks),
    min_counts = as.numeric(min_counts)
  )
  obj <- Seurat::CreateSeuratObject(counts = res$counts, assay = "Spatial", project = sample_name)
  obj$orig.ident <- sample_name
  obj$coord_x <- res$coords$x
  obj$coord_y <- res$coords$y
  attr(obj, "ts_spatial_mode") <- "merscope"
  attr(obj, "ts_merscope_blanks_removed") <- res$n_blanks_removed
  obj
}

#' Pool a Visium HD Seurat object into coarser spatial super-bins
#'
#' @param seurat_obj A `Seurat` object loaded from Visium HD.
#' @param bin_factor Integer pooling factor (`<= 1L` returns `seurat_obj` unchanged).
#' @return A `Seurat` object with pooled sparse counts and `coord_x`/`coord_y`/`n_subbins`.
#' @export
apply_visium_hd_bin_pooling <- function(seurat_obj, bin_factor = 1L) {
  bf <- suppressWarnings(as.integer(bin_factor))
  if (length(bf) != 1L || is.na(bf) || bf <= 1L) {
    return(seurat_obj)
  }
  if (!requireNamespace("Seurat", quietly = TRUE) || !requireNamespace("SeuratObject", quietly = TRUE)) {
    stop(.spatial_hd_stop("Packages 'Seurat' and 'SeuratObject' are required for apply_visium_hd_bin_pooling()."))
  }
  assay <- Seurat::DefaultAssay(seurat_obj)
  counts_mat <- SeuratObject::LayerData(seurat_obj, assay = assay, layer = "counts")
  coords <- get_spatial_coords(seurat_obj)
  if (is.null(coords) || nrow(coords) == 0L) {
    stop(.spatial_hd_stop("Cannot pool Visium HD bins: spatial coordinates are missing."))
  }
  pooled <- aggregate_visium_hd_bins(counts_mat, coords, bin_factor = bf)
  proj <- tryCatch(seurat_obj@project.name, error = function(e) "VisiumHD_Pooled")
  obj <- Seurat::CreateSeuratObject(counts = pooled$counts, assay = assay, project = proj)
  obj$orig.ident <- proj
  obj$coord_x <- pooled$coords$x
  obj$coord_y <- pooled$coords$y
  obj$n_subbins <- pooled$coords$n_subbins
  attr(obj, "ts_spatial_mode") <- "visium_hd_pooled"
  attr(obj, "ts_hd_pool_factor") <- pooled$bin_factor
  obj
}
