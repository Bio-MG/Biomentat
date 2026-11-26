# =============================================================================
# tools/pkg_skeleton/R/demo_loader.R — P3.3: 1-click offline demo loader
# =============================================================================
# Pure-R offline demo dataset loader for Bulk, Single-Cell, and Spatial
# modalities. Loads compact real biological subsets bundled under
# inst/extdata/demo/ with deterministic synthetic fallback if absent.
# =============================================================================

.demo_stop <- function(msg) {
  stop(errorCondition(msg, class = "demo_loader_error"))
}

.demo_rda_path <- function(modality) {
  fname <- paste0("demo_", modality, ".rda")
  candidates <- c(
    file.path("tools", "pkg_skeleton", "inst", "extdata", "demo", fname),
    file.path("..", "..", "tools", "pkg_skeleton", "inst", "extdata", "demo", fname),
    file.path("..", "tools", "pkg_skeleton", "inst", "extdata", "demo", fname)
  )
  pkg_path <- tryCatch(
    system.file("extdata", "demo", fname, package = "biomentat"),
    error = function(e) ""
  )
  if (nzchar(pkg_path)) candidates <- c(candidates, pkg_path)
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0L) hit[1L] else NULL
}

#' List available offline 1-click demo datasets
#'
#' @return Data.frame describing available demo modalities.
biomentat_demo_manifest <- function() {
  data.frame(
    modality = c("bulk", "sc", "spatial"),
    title = c(
      "Bulk RNA-seq Demo — Human Cornea SARS-CoV-2 (GSE164073, 2000 genes x 6 samples)",
      "Single-Cell Demo — Human PBMC IFN-beta (GSE96583, 1500 genes x 400 cells)",
      "Spatial Visium Demo — Mouse Brain Anterior (stxBrain, 1500 genes x 400 spots)"
    ),
    offline_only = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

#' Load an offline demo dataset for a given modality
#'
#' @param modality Character, one of `"bulk"`, `"sc"`, or `"spatial"`.
#' @param seed Integer RNG seed for fallback reproducibility.
#' @return List containing `counts`, `metadata`, `provenance` (and `coords` for spatial).
load_biomentat_demo <- function(modality = "bulk", seed = 42L) {
  if (!is.character(modality) || length(modality) != 1L ||
      !(modality %in% c("bulk", "sc", "spatial"))) {
    .demo_stop("load_biomentat_demo() : modality doit etre 'bulk', 'sc' ou 'spatial'.")
  }

  rda_path <- .demo_rda_path(modality)
  if (!is.null(rda_path)) {
    e <- new.env(parent = emptyenv())
    loaded <- load(rda_path, envir = e)
    if ("demo_data" %in% loaded) {
      return(e$demo_data)
    }
  }

  set.seed(seed)
  if (identical(modality, "bulk")) {
    ng <- 200L; ns <- 8L
    counts <- matrix(stats::rpois(ng * ns, lambda = 40), nrow = ng, ncol = ns)
    counts[1:20, 5:8] <- counts[1:20, 5:8] + 120L
    rownames(counts) <- paste0("Gene-", seq_len(ng))
    colnames(counts) <- paste0("Sample_", seq_len(ns))
    meta <- data.frame(
      sample_id = colnames(counts),
      condition = factor(rep(c("CTL", "TRT"), each = 4L), levels = c("CTL", "TRT")),
      batch = factor(rep(c("B1", "B2"), times = 4L)),
      row.names = colnames(counts),
      stringsAsFactors = FALSE
    )
    return(list(
      counts = counts,
      metadata = meta,
      provenance = list(source = "demo_loader_fallback", modality = "bulk", seed = seed, real_data = FALSE)
    ))
  }

  if (identical(modality, "sc")) {
    ng <- 100L; nc <- 60L
    mat <- matrix(stats::rpois(ng * nc, lambda = 2), nrow = ng, ncol = nc)
    mat[1:15, 1:30] <- mat[1:15, 1:30] + 12L
    mat[16:30, 31:60] <- mat[16:30, 31:60] + 12L
    rownames(mat) <- paste0("Gene-", seq_len(ng))
    colnames(mat) <- paste0("Cell_", seq_len(nc))
    meta <- data.frame(
      sample_id = rep(c("S1", "S2", "S3", "S4"), each = 15L),
      condition = rep(c("CTL", "TRT"), each = 30L),
      cluster = rep(c("Cluster_0", "Cluster_1"), each = 30L),
      row.names = colnames(mat),
      stringsAsFactors = FALSE
    )
    return(list(
      counts = Matrix::Matrix(mat, sparse = TRUE),
      metadata = meta,
      provenance = list(source = "demo_loader_fallback", modality = "sc", seed = seed, real_data = FALSE)
    ))
  }

  ng <- 80L; np <- 50L
  mat <- matrix(stats::rpois(ng * np, lambda = 3), nrow = ng, ncol = np)
  mat[1:10, 1:25] <- mat[1:10, 1:25] + 15L
  rownames(mat) <- paste0("Gene-", seq_len(ng))
  colnames(mat) <- paste0("Spot_", seq_len(np))
  coords <- data.frame(
    id = colnames(mat),
    x = stats::runif(np, 0, 100),
    y = stats::runif(np, 0, 100),
    fov = "demo_fov",
    stringsAsFactors = FALSE
  )
  meta <- data.frame(
    sample_id = rep(c("Sec1", "Sec2"), each = 25L),
    condition = rep(c("CTL", "TRT"), each = 25L),
    niche = rep(c("Niche_A", "Niche_B"), each = 25L),
    row.names = colnames(mat),
    stringsAsFactors = FALSE
  )
  list(
    counts = Matrix::Matrix(mat, sparse = TRUE),
    coords = coords,
    metadata = meta,
    provenance = list(source = "demo_loader_fallback", modality = "spatial", seed = seed, real_data = FALSE)
  )
}
