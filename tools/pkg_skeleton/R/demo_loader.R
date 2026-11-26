# =============================================================================
# tools/pkg_skeleton/R/demo_loader.R — P3.3: 1-click offline demo loader
# =============================================================================
# Pure-R offline demo dataset loader for Bulk, Single-Cell, and Spatial
# modalities. Loads compact real biological subsets bundled under
# inst/extdata/demo/.
# =============================================================================

.demo_stop <- function(msg) {
  stop(errorCondition(msg, class = "demo_loader_error"))
}

.demo_rda_path <- function(modality) {
  fname <- paste0("demo_", modality, ".rda")
  candidates <- character(0)
  if (exists(".sc_project_root", mode = "function")) {
    root <- tryCatch(.sc_project_root(), error = function(e) "")
    if (is.character(root) && length(root) == 1L && nzchar(root)) {
      candidates <- c(
        candidates,
        file.path(root, "tools", "pkg_skeleton", "inst", "extdata", "demo", fname)
      )
    }
  }
  candidates <- c(
    candidates,
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
#' @param seed Integer RNG seed (kept for signature compatibility).
#' @return List containing `counts`, `metadata`, `provenance` (and `coords` for spatial).
load_biomentat_demo <- function(modality = "bulk", seed = 42L) {
  if (!is.character(modality) || length(modality) != 1L ||
      !(modality %in% c("bulk", "sc", "spatial"))) {
    .demo_stop("load_biomentat_demo() : modality doit etre 'bulk', 'sc' ou 'spatial'.")
  }

  rda_path <- .demo_rda_path(modality)
  if (is.null(rda_path)) {
    .demo_stop(sprintf("Bundle démo introuvable pour '%s' (%s).",
                       modality, paste0("demo_", modality, ".rda")))
  }
  e <- new.env(parent = emptyenv())
  loaded <- load(rda_path, envir = e)
  if (!("demo_data" %in% loaded)) {
    .demo_stop(sprintf("Objet 'demo_data' absent dans %s.", rda_path))
  }
  e$demo_data
}
