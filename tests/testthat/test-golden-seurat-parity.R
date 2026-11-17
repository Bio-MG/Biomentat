# =============================================================================
# test-golden-seurat-parity.R — W3.1 GOLDEN / parité Seurat & transforms affines
# =============================================================================
# 1. Parité numérique entre run_sc_manual_pipeline_core() (sc_pipeline.R) et
#    la recette Seurat de référence appelée directement sur la même fixture
#    déterministe (NormalizeData -> FindVariableFeatures -> ScaleData ->
#    RunPCA -> FindNeighbors -> FindClusters + FindMarkers).
# 2. Parité géométrique D_4 (8 symétries du groupe diédral) sur
#    apply_coord_orientation() (spatial_io.R) : involutions, rotations
#    90/180/270 deg, préservation des distances euclidiennes et de l'emprise.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_bpcells.R")
source_project_file("tools/pkg_skeleton/R/sc_pipeline.R")
source_project_file("tools/pkg_skeleton/R/spatial_io.R")

.golden_sc_fixture <- function(seed = 20261117L) {
  withr::with_seed(seed, {
    ng <- 150L
    nc <- 60L
    counts <- matrix(stats::rpois(ng * nc, lambda = 4), nrow = ng, ncol = nc)
    # Two clear biological states (cells 1..30 vs 31..60)
    counts[1:25, 1:30]  <- counts[1:25, 1:30] + stats::rpois(25L * 30L, lambda = 18)
    counts[26:50, 31:60] <- counts[26:50, 31:60] + stats::rpois(25L * 30L, lambda = 18)
    rownames(counts) <- paste0("Gene", seq_len(ng))
    colnames(counts) <- paste0("Cell", seq_len(nc))
    suppressWarnings(Seurat::CreateSeuratObject(counts = counts, project = "golden_sc"))
  })
}

test_that("GOLDEN Seurat : run_sc_manual_pipeline_core() est en parité avec la recette Seurat directe", {
  skip_if_not_installed("Seurat")
  suppressPackageStartupMessages(library(Seurat))
  obj_in <- .golden_sc_fixture()
  seed_use <- 989L
  pca_dim <- 10L
  clust_res <- 0.5

  res <- suppressWarnings(run_sc_manual_pipeline_core(
    obj_in,
    params = list(
      qc_min_gene = 10,
      qc_max_gene = 5000,
      qc_mt = 50,
      doublets_enabled = FALSE,
      run_ambient = FALSE,
      norm_method = "log",
      pca_dim = pca_dim,
      clust_res = clust_res,
      cluster_algo = 1L,
      reduction_method = "umap",
      seed = seed_use
    )
  ))
  expect_true(isTRUE(res$ok),
              info = paste("step:", res$step, "err:",
                           tryCatch(conditionMessage(res$error), error = function(e) "")))
  eng <- res$obj

  # Reference Seurat recipe directly on the same cells
  ref <- suppressWarnings({
    o <- subset(obj_in, cells = colnames(eng))
    o <- Seurat::NormalizeData(o, normalization.method = "LogNormalize",
                               scale.factor = 10000, verbose = FALSE)
    o <- Seurat::FindVariableFeatures(o, selection.method = "vst",
                                      nfeatures = min(2000L, nrow(o)), verbose = FALSE)
    o <- Seurat::ScaleData(o, features = Seurat::VariableFeatures(o), verbose = FALSE)
    npcs_ref <- min(pca_dim, ncol(o) - 1L, length(Seurat::VariableFeatures(o)) - 1L)
    o <- withr::with_seed(seed_use,
      Seurat::RunPCA(o, features = Seurat::VariableFeatures(o),
                     npcs = npcs_ref, verbose = FALSE))
    o <- Seurat::FindNeighbors(o, dims = seq_len(npcs_ref), verbose = FALSE)
    o <- withr::with_seed(seed_use,
      Seurat::FindClusters(o, resolution = clust_res, algorithm = 1L, verbose = FALSE))
    o
  })

  # 1. Normalized data parity
  eng_data <- as.matrix(SeuratObject::LayerData(eng, assay = "RNA", layer = "data"))
  ref_data <- as.matrix(SeuratObject::LayerData(ref, assay = "RNA", layer = "data"))
  expect_equal(eng_data, ref_data, tolerance = 1e-12)

  # 2. Variable features parity
  expect_setequal(Seurat::VariableFeatures(eng), Seurat::VariableFeatures(ref))

  # 3. PCA embeddings parity (up to column sign ambiguity on singular vectors)
  eng_pca <- Seurat::Embeddings(eng, "pca")[, seq_len(pca_dim), drop = FALSE]
  ref_pca <- Seurat::Embeddings(ref, "pca")[, seq_len(pca_dim), drop = FALSE]
  for (j in seq_len(ncol(eng_pca))) {
    s <- if (sum(eng_pca[, j] * ref_pca[, j]) < 0) -1 else 1
    expect_equal(s * eng_pca[, j], ref_pca[, j], tolerance = 1e-8)
  }

  # 4. Cluster assignments & FindMarkers parity
  expect_identical(as.character(eng$seurat_clusters), as.character(ref$seurat_clusters))
  Seurat::Idents(eng) <- eng$seurat_clusters
  Seurat::Idents(ref) <- ref$seurat_clusters
  lv <- levels(Seurat::Idents(eng))
  if (length(lv) >= 2L) {
    mk_eng <- Seurat::FindMarkers(eng, ident.1 = lv[1L], ident.2 = lv[2L],
                                  only.pos = FALSE, min.pct = 0.1,
                                  logfc.threshold = 0, verbose = FALSE)
    mk_ref <- Seurat::FindMarkers(ref, ident.1 = lv[1L], ident.2 = lv[2L],
                                  only.pos = FALSE, min.pct = 0.1,
                                  logfc.threshold = 0, verbose = FALSE)
    expect_identical(rownames(mk_eng), rownames(mk_ref))
    expect_equal(mk_eng$avg_log2FC, mk_ref$avg_log2FC, tolerance = 1e-10)
    expect_equal(mk_eng$p_val, mk_ref$p_val, tolerance = 1e-10)
  }
})

test_that("GOLDEN Spatial D4 : groupe diédral complet (rotations 90/180/270, isométrie & emprise)", {
  set.seed(42L)
  coords <- data.frame(
    id = paste0("spot_", 1:25),
    x = stats::runif(25, 110, 350),
    y = stats::runif(25, 40, 190),
    stringsAsFactors = FALSE
  )
  d0 <- as.matrix(stats::dist(coords[, c("x", "y")]))

  # All 8 D4 combinations preserve pairwise Euclidean distances (isometry)
  combos <- expand.grid(swap = c(FALSE, TRUE), fx = c(FALSE, TRUE), fy = c(FALSE, TRUE))
  for (i in seq_len(nrow(combos))) {
    tr <- apply_coord_orientation(
      coords,
      swap_xy = combos$swap[i],
      flip_x = combos$fx[i],
      flip_y = combos$fy[i],
      img_width = 500, img_height = 400
    )
    expect_equal(as.matrix(stats::dist(tr[, c("x", "y")])), d0, tolerance = 1e-10)
  }

  # 90-degree CW rotation (swap_xy = TRUE, flip_x = TRUE) applied 4 times = identity
  rot90 <- function(c) apply_coord_orientation(c, swap_xy = TRUE, flip_x = TRUE,
                                               img_width = 400, img_height = 400)
  c4 <- rot90(rot90(rot90(rot90(coords))))
  expect_equal(c4$x, coords$x, tolerance = 1e-12)
  expect_equal(c4$y, coords$y, tolerance = 1e-12)

  # 180-degree rotation (flip_x = TRUE, flip_y = TRUE) is an involution even without image extent
  c180 <- apply_coord_orientation(coords, flip_x = TRUE, flip_y = TRUE)
  expect_equal(range(c180$x), range(coords$x), tolerance = 1e-12)
  expect_equal(range(c180$y), range(coords$y), tolerance = 1e-12)
  c360 <- apply_coord_orientation(c180, flip_x = TRUE, flip_y = TRUE)
  expect_equal(c360$x, coords$x, tolerance = 1e-12)
  expect_equal(c360$y, coords$y, tolerance = 1e-12)
})
