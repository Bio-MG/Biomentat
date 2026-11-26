# =============================================================================
# tools/build_demo_datasets.R — Build compact real offline demo .rda bundles
# =============================================================================
# Generates:
#   - tools/pkg_skeleton/inst/extdata/demo/demo_bulk.rda    (GSE164073 2000x6)
#   - tools/pkg_skeleton/inst/extdata/demo/demo_sc.rda      (GSE96583 1500x400)
#   - tools/pkg_skeleton/inst/extdata/demo/demo_spatial.rda (stxBrain 1500x400)
# Usage: Rscript tools/build_demo_datasets.R
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SeuratObject)
})

out_dir <- file.path("tools", "pkg_skeleton", "inst", "extdata", "demo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

.extract_counts <- function(obj, assay = "RNA") {
  mat <- tryCatch(
    SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(mat) || nrow(mat) == 0L) {
    mat <- SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  }
  methods::as(mat, "dgCMatrix")
}

.select_top_genes <- function(mat, n_target = 1500L, mt_pattern = "^MT-") {
  rs <- Matrix::rowSums(mat)
  detected <- names(rs)[rs > 0]
  mt_genes <- detected[grepl(mt_pattern, detected)]
  non_mt <- setdiff(detected, mt_genes)
  rem <- max(0L, n_target - length(mt_genes))
  if (length(non_mt) > rem) {
    # Score by log-mean + variance proxy on sparse counts
    m_sub <- mat[non_mt, , drop = FALSE]
    mu <- Matrix::rowMeans(m_sub)
    v <- Matrix::rowMeans(m_sub^2) - mu^2
    score <- log1p(mu) + 0.5 * log1p(pmax(v, 0))
    ord <- order(score, decreasing = TRUE)
    non_mt <- non_mt[ord[seq_len(rem)]]
  }
  keep <- c(mt_genes, non_mt)
  mat[keep, , drop = FALSE]
}

# ── 1. Bulk Demo (GSE164073_cornea_3v3) ─────────────────────────────────────
build_bulk_demo <- function() {
  counts_path <- file.path("data_sc", "datasets", "GSE164073_cornea_3v3", "counts.tsv")
  ss_path <- file.path("data_sc", "datasets", "GSE164073_cornea_3v3", "samplesheet.tsv")
  if (!file.exists(counts_path) || !file.exists(ss_path)) {
    stop("Missing GSE164073_cornea_3v3 source files.")
  }
  df <- utils::read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
  ss <- utils::read.delim(ss_path, check.names = FALSE, stringsAsFactors = FALSE)

  gene_col <- colnames(df)[1L]
  genes <- trimws(as.character(df[[gene_col]]))
  valid <- !is.na(genes) & nzchar(genes) & !duplicated(genes)
  df <- df[valid, , drop = FALSE]
  genes <- genes[valid]

  sample_ids <- as.character(ss$sample_id)
  mat <- as.matrix(df[, sample_ids, drop = FALSE])
  storage.mode(mat) <- "integer"
  rownames(mat) <- genes

  keep_expr <- rowSums(mat, na.rm = TRUE) >= 20L
  mat <- mat[keep_expr, , drop = FALSE]

  mock_idx <- which(ss$condition == "mock")
  cov2_idx <- which(ss$condition == "CoV2")
  m_mock <- rowMeans(mat[, mock_idx, drop = FALSE])
  m_cov2 <- rowMeans(mat[, cov2_idx, drop = FALSE])
  score <- abs(log2((m_cov2 + 1) / (m_mock + 1))) * log10(rowMeans(mat) + 1)
  n_keep <- min(2000L, nrow(mat))
  ord <- order(score, decreasing = TRUE)
  mat <- mat[ord[seq_len(n_keep)], , drop = FALSE]

  meta <- data.frame(
    sample_id = sample_ids,
    condition = factor(ss$condition, levels = c("mock", "CoV2")),
    tissue = as.character(ss$tissue),
    replicate = factor(ss$replicate),
    row.names = sample_ids,
    stringsAsFactors = FALSE
  )

  demo_data <- list(
    counts = mat,
    metadata = meta,
    provenance = list(
      source = "GSE164073_cornea_3v3",
      modality = "bulk",
      organism = "Homo sapiens",
      real_data = TRUE
    )
  )
  out_file <- file.path(out_dir, "demo_bulk.rda")
  save(demo_data, file = out_file, compress = "xz")
  message(sprintf("Saved %s (%d x %d, %.1f KB)",
                  out_file, nrow(mat), ncol(mat), file.info(out_file)$size / 1024))
}

# ── 2. Single-Cell Demo (GSE96583 PBMC pbmc2cond_test.rds) ───────────────────
build_sc_demo <- function() {
  rds_path <- file.path("data_sc", "datasets", "pbmc2cond_test.rds")
  if (!file.exists(rds_path)) stop("Missing pbmc2cond_test.rds source file.")
  obj <- readRDS(rds_path)
  md <- obj@meta.data

  set.seed(42L)
  ctrl_cells <- rownames(md)[md$condition == "control"]
  stim_cells <- rownames(md)[md$condition == "stim"]
  pick <- c(
    sample(ctrl_cells, min(200L, length(ctrl_cells))),
    sample(stim_cells, min(200L, length(stim_cells)))
  )

  mat <- .extract_counts(obj, assay = "RNA")[, pick, drop = FALSE]
  mat <- .select_top_genes(mat, n_target = 1500L, mt_pattern = "^MT-")
  mat <- methods::as(mat, "dgCMatrix")

  md_sub <- md[pick, , drop = FALSE]
  meta <- data.frame(
    sample_id = as.character(md_sub$sample_id),
    condition = as.character(md_sub$condition),
    donor = as.character(md_sub$donor),
    celltype = as.character(md_sub$celltype),
    cluster = as.character(md_sub$cluster),
    row.names = pick,
    stringsAsFactors = FALSE
  )

  demo_data <- list(
    counts = mat,
    metadata = meta,
    provenance = list(
      source = "GSE96583_pbmc_2cond",
      modality = "sc",
      organism = "Homo sapiens",
      real_data = TRUE
    )
  )
  out_file <- file.path(out_dir, "demo_sc.rda")
  save(demo_data, file = out_file, compress = "xz")
  message(sprintf("Saved %s (%d x %d, %.1f KB)",
                  out_file, nrow(mat), ncol(mat), file.info(out_file)$size / 1024))
}

# ── 3. Spatial Demo (stxBrain.SeuratData anterior1 Visium) ──────────────────
build_spatial_demo <- function() {
  if (!requireNamespace("stxBrain.SeuratData", quietly = TRUE)) {
    stop("stxBrain.SeuratData package is not installed.")
  }
  br <- SeuratData::LoadData("stxBrain", type = "anterior1")
  tc <- Seurat::GetTissueCoordinates(br)
  cn <- tolower(colnames(tc))
  col_x <- if ("imagecol" %in% cn) colnames(tc)[cn == "imagecol"][1L] else if ("x" %in% cn) colnames(tc)[cn == "x"][1L] else colnames(tc)[2L]
  col_y <- if ("imagerow" %in% cn) colnames(tc)[cn == "imagerow"][1L] else if ("y" %in% cn) colnames(tc)[cn == "y"][1L] else colnames(tc)[1L]
  spot_ids <- if ("cell" %in% cn) as.character(tc[[colnames(tc)[cn == "cell"][1L]]]) else rownames(tc)

  coords_all <- data.frame(
    id = spot_ids,
    x = as.numeric(tc[[col_x]]),
    y = as.numeric(tc[[col_y]]),
    fov = "anterior1",
    stringsAsFactors = FALSE
  )
  coords_all <- coords_all[match(colnames(br), coords_all$id), , drop = FALSE]

  # Contiguous 400-spot crop around spatial centroid
  cx <- stats::median(coords_all$x, na.rm = TRUE)
  cy <- stats::median(coords_all$y, na.rm = TRUE)
  dist2 <- (coords_all$x - cx)^2 + (coords_all$y - cy)^2
  keep_idx <- order(dist2)[seq_len(min(400L, nrow(coords_all)))]
  pick <- coords_all$id[keep_idx]
  coords <- coords_all[keep_idx, , drop = FALSE]
  rownames(coords) <- NULL

  mat <- .extract_counts(br, assay = "Spatial")[, pick, drop = FALSE]
  mat <- .select_top_genes(mat, n_target = 1500L, mt_pattern = "^mt-")
  mat <- methods::as(mat, "dgCMatrix")

  meta <- data.frame(
    sample_id = rep("anterior1", length(pick)),
    region = ifelse(coords$x <= stats::median(coords$x), "Zone_A", "Zone_B"),
    row.names = pick,
    stringsAsFactors = FALSE
  )

  demo_data <- list(
    counts = mat,
    coords = coords,
    metadata = meta,
    histology = NULL,
    provenance = list(
      source = "stxBrain_anterior1_visium",
      modality = "spatial",
      organism = "Mus musculus",
      real_data = TRUE
    )
  )
  out_file <- file.path(out_dir, "demo_spatial.rda")
  save(demo_data, file = out_file, compress = "xz")
  message(sprintf("Saved %s (%d x %d, %.1f KB)",
                  out_file, nrow(mat), ncol(mat), file.info(out_file)$size / 1024))
}

build_bulk_demo()
build_sc_demo()
build_spatial_demo()
