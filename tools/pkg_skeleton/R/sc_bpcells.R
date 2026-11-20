# =============================================================================
# helpers_sc_bpcells.R — Disk-backed (BPCells) pipeline support (Step-3.7A)
# =============================================================================
# "Intelligent" BPCells integration for the FULL SC pipeline (not just import):
#   - Large objects get their raw counts converted to an on-disk BPCells
#     matrix BEFORE Normalisation/PCA, so the pipeline never has to hold the
#     full genes x cells matrix in RAM.
#   - ScaleData() is redirected through smart_scale_data(), which restricts
#     scaling to VariableFeatures by default -- Seurat's own ScaleData()
#     default (features=NULL) scales EVERY gene, which is the single biggest
#     RAM trap in a routine pipeline, independent of whether counts are
#     disk-backed or not (a dense 30-40k-gene x N-cell matrix is the actual
#     memory blow-up, not the sparse raw counts).
#   - ensure_genes_scaled() extends scale.data coverage on demand for
#     DoHeatmap() calls on arbitrary marker genes that fall outside
#     VariableFeatures, without ever re-scaling the whole assay.
#
# Depends on: Seurat (v5, Assay5), BPCells (optional -- checked at call time,
# app degrades gracefully to the standard in-memory pipeline if absent).
# =============================================================================

#' Is the BPCells package installed?
.bpcells_available <- function() requireNamespace("BPCells", quietly = TRUE)

#' Current storage backend of a Seurat object's RNA counts layer.
#' @param obj Seurat object.
#' @return "disk" (BPCells/DelayedArray-backed) or "memory" (dense/dgCMatrix).
sc_backend_status <- function(obj) {
  if (is.null(obj) || !"RNA" %in% names(obj@assays)) return("memory")
  mat <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) tryCatch(GetAssayData(obj, assay = "RNA", slot = "counts"),
                                 error = function(e2) NULL)
  )
  if (is.null(mat)) return("memory")
  if (inherits(mat, "IterableMatrix") || inherits(mat, "DelayedMatrix")) "disk" else "memory"
}

#' Convert a Seurat object's raw counts to an on-disk BPCells matrix.
#'
#' Writes counts to `dir` (a fresh tempdir() subdirectory by default) and
#' rebuilds a lean Seurat object pointing at that on-disk matrix, carrying
#' over the COMPLETE cell metadata (`obj@meta.data`, design columns included).
#' Any existing normalization/PCA/UMAP/clusters is
#' lost -- same constraint as ID mapping (mod_sc_mapping.R): this MUST run
#' BEFORE Normalisation, right after QC.
#'
#' No-op (returns obj unchanged, `already_disk=TRUE`) if the object is
#' already disk-backed, so callers can call this unconditionally.
#'
#' @param obj Seurat object (raw counts in the RNA assay).
#' @param dir Target directory. NULL (default, recommended) uses a fresh
#'   tempfile() path -- callers should register it for cleanup via
#'   `session$onSessionEnded(function() unlink(dir, recursive=TRUE))`, since
#'   each conversion leaves one folder behind in tempdir() for the session's
#'   lifetime.
#' @return list(object, dir, n_cells, n_genes, already_disk). `dir` is NULL
#'   when already_disk is TRUE (nothing new was written, nothing to clean up).
convert_seurat_to_bpcells <- function(obj, dir = NULL) {
  if (!.bpcells_available())
    stop(errorCondition("Package 'BPCells' non installé. Installez-le via remotes::install_github('bnprks/BPCells/r').", class = "sc_bpcells_error"))

  if (sc_backend_status(obj) == "disk") {
    return(list(object = obj, dir = NULL, n_cells = ncol(obj), n_genes = nrow(obj), already_disk = TRUE))
  }

  if (is.null(dir)) dir <- tempfile(pattern = "bpcells_pipeline_")

  counts <- tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
  if (!inherits(counts, "dgCMatrix")) counts <- methods::as(counts, "CsparseMatrix")

  BPCells::write_matrix_dir(mat = counts, dir = dir, overwrite = TRUE)
  disk_mat <- BPCells::open_matrix_dir(dir = dir)

  # Preserve the COMPLETE cell metadata: design columns (batch, condition,
  # patient id, ...) are declared AFTER conversion (Stage-13 design panel),
  # so keeping only `orig.ident` would silently destroy the user's declared
  # design. CreateSeuratObject() recomputes nCount_RNA/nFeature_RNA from the
  # new counts; normalization/PCA/UMAP/cluster slots are still dropped
  # (conversion MUST run before normalisation).
  new_obj <- CreateSeuratObject(counts = disk_mat, meta.data = obj@meta.data,
                                project = Project(obj))

  list(object = new_obj, dir = dir, n_cells = ncol(new_obj), n_genes = nrow(new_obj), already_disk = FALSE)
}

#' ScaleData() restricted to a RAM-safe gene set (Step-3.7A default).
#'
#' Seurat's `ScaleData(features=NULL)` scales EVERY gene by default -- for a
#' 30-40k-gene assay this densifies the full genes x cells matrix into RAM
#' even though only the ~2000 variable features are ever used downstream
#' (PCA). This wrapper scales VariableFeatures(obj) plus any `extra_features`
#' the caller explicitly needs (e.g. marker genes about to be Heatmap'd),
#' and nothing else -- works identically for in-memory and disk-backed assays.
#'
#' @param obj Seurat object with NormalizeData()+FindVariableFeatures() already run.
#' @param extra_features Additional gene symbols to include.
#' @param vars.to.regress Optionnel (jalon SF, 2026-11-09) : colonnes metadata
#'   passées à `ScaleData(vars.to.regress = ...)` — ex. "CC.Difference"
#'   (S − G2M, défaut fiche §3-SF) ou c("S.Score","G2M.Score") (option
#'   complète). NULL (défaut) = comportement historique strict. L'appelant
#'   NE doit PAS l'utiliser sur un objet BPCells on-disk (régression ScaleData
#'   non garantie disk-backed — garde côté pipeline).
#' @return Seurat object with scale.data covering VariableFeatures ∪ extra_features.
smart_scale_data <- function(obj, extra_features = character(0),
                             vars.to.regress = NULL) {
  var_feat <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
  if (length(var_feat) == 0) {
    # FIX P1.1 Lead G (audit externe #4 §4, 2026-11-12) : le repli silencieux
    # historique `var_feat <- rownames(obj)` scalait TOUS les gènes — précisé-
    # ment quand la sélection HVG amont est absente ou a échoué — sous une
    # étiquette « RAM-safe ». Politique bornée explicite : sélection HVG vst
    # (plafond Seurat 2000) exécutée ICI et TRACÉE par un warning ; si elle ne
    # produit rien, erreur actionable — jamais de scaling du génome entier
    # silencieux.
    warning("smart_scale_data() : VariableFeatures vide (FindVariableFeatures non exécuté en amont ?) — sélection HVG bornée explicite (vst, nfeatures=2000) au lieu du scaling du génome entier.",
            call. = FALSE)
    obj <- tryCatch(
      Seurat::FindVariableFeatures(obj, selection.method = "vst",
                                   nfeatures = 2000L, verbose = FALSE),
      error = function(e) obj)
    var_feat <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
    if (length(var_feat) == 0) {
      stop(errorCondition(
        paste0("smart_scale_data() : aucune feature variable sélectionnable ",
               "(HVG vst vide) — refus de scaler tout le génome. Exécutez ",
               "FindVariableFeatures() en amont ou vérifiez l'objet (gènes ",
               "constants / objet dégénéré)."),
        class = "sc_bpcells_error"))
    }
  }
  feats <- unique(c(var_feat, intersect(extra_features, rownames(obj))))
  if (is.null(vars.to.regress)) {
    ScaleData(obj, features = feats, verbose = FALSE)
  } else {
    ScaleData(obj, features = feats, vars.to.regress = vars.to.regress,
              verbose = FALSE)
  }
}

#' Ensure a set of genes are covered by scale.data before a DoHeatmap() call,
#' extending coverage (VariableFeatures ∪ already-scaled ∪ requested genes)
#' rather than ever re-scaling the FULL assay just to plot a handful of
#' marker genes that happen to fall outside VariableFeatures.
#'
#' Note: Seurat's ScaleData() replaces (not merges) the scale.data slot for
#' the `features=` it's given, so this recomputes the small combined set each
#' time it's called rather than being a true incremental cache -- still cheap
#' (a few thousand genes at most, never the full 30-40k-gene assay), which is
#' the actual RAM-safety property being preserved here.
#'
#' @param obj Seurat object (NormalizeData() already run).
#' @param genes Gene symbols the caller is about to plot.
#' @return Seurat object, scale.data extended to cover `genes` if needed.
ensure_genes_scaled <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj))
  if (!length(genes)) return(obj)
  scaled_now <- tryCatch(rownames(GetAssayData(obj, layer = "scale.data")),
                         error = function(e) character(0))
  if (all(genes %in% scaled_now)) return(obj)
  ScaleData(obj, features = unique(c(scaled_now, genes)), verbose = FALSE)
}

#' Re-orient a disk-backed (BPCells) object's counts to row-major (gene-major)
#' storage before FindAllMarkers()/find_correlated_genes() (Step-3.8).
#'
#' Seurat/BPCells warns: "Column-major order detected; FindMarkers requires
#' row-major order. Consider first running BPCells::transpose_storage_order()
#' to avoid repeated transpositions." -- BPCells' default write orientation
#' (from convert_seurat_to_bpcells()) is cell-major (column-major), matching
#' how counts are read during import/QC/normalization; marker tests instead
#' scan gene-by-gene (row-major), so without this step every FindAllMarkers
#' call silently re-transposes internally, repeatedly, which is the actual
#' reason marker/correlation calculations can crawl on a disk-backed dataset
#' even after RAM-safety cell subsampling.
#'
#' Deliberately meant to be called AFTER subsample_seurat_for_analysis() by
#' its callers (mod_sc_markers.R, mod_sc_corr.R, mod_sc.R auto-pipeline) --
#' transposing the FULL multi-million-cell matrix would itself be a slow full
#' disk rewrite; transposing the already-capped subsample (typically a few
#' hundred thousand cells at most) is comparatively cheap and one-time.
#'
#' Fully defensive: no-op (returns `obj` unchanged) if the object isn't
#' disk-backed, if BPCells is unavailable, or if the transpose call itself
#' fails for any reason (API surface not guaranteed stable across BPCells
#' versions) -- callers never need their own tryCatch around this.
#'
#' @param obj Seurat object, ideally already cell-count-capped.
#' @param dir Target directory for the transposed copy. NULL uses tempfile().
#' @return list(object, dir, transposed). `dir` is NULL when transposed==FALSE
#'   (nothing written, nothing to clean up); when TRUE, caller should schedule
#'   `unlink(dir, recursive=TRUE)` (e.g. via session$onSessionEnded()).
optimize_bpcells_for_markers <- function(obj, dir = NULL) {
  no_op <- list(object = obj, dir = NULL, transposed = FALSE)
  if (!.bpcells_available() || sc_backend_status(obj) != "disk") return(no_op)

  mat <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (is.null(mat) || !inherits(mat, "IterableMatrix")) return(no_op)

  if (is.null(dir)) dir <- tempfile(pattern = "bpcells_rowmajor_")

  result <- tryCatch({
    row_major <- BPCells::transpose_storage_order(mat)
    BPCells::write_matrix_dir(mat = row_major, dir = dir, overwrite = TRUE)
    disk_mat  <- BPCells::open_matrix_dir(dir = dir)

    meta_keep <- obj@meta.data
    new_obj   <- CreateSeuratObject(counts = disk_mat, meta.data = meta_keep, project = Project(obj))
    list(object = new_obj, dir = dir, transposed = TRUE)
  }, error = function(e) NULL)

  if (is.null(result)) return(no_op)
  result
}

#' Detect and normalize legacy AnnData `.h5ad` sparse layout & orientation.
#'
#' Modern AnnData (>= 0.8) records `encoding-type` (`"csr_matrix"`,
#' `"csc_matrix"`, or `"array"`) on `/X` and layers, whereas legacy `.h5ad`
#' files often omit `encoding-type` and store sparse `(indptr, indices, data)`
#' groups in CSR (or CSC) orientation with `(n_obs, n_vars)` shape. R/Seurat
#' and BPCells expect `(n_vars, n_obs)` (genes x cells) orientation.
#'
#' @param encoding_type Character scalar from HDF5 `encoding-type` attribute,
#'   or `NULL`/`NA`/`""` when absent (legacy `.h5ad`).
#' @param dims Optional numeric/integer vector of length 2 (`c(nrow, ncol)`).
#' @param n_obs Optional non-negative integer count of observations (cells).
#' @param n_vars Optional non-negative integer count of variables (genes).
#' @return A list with elements `layout` (`"csr"`, `"csc"`, or `"dense"`),
#'   `needs_transpose` (logical), `legacy_fallback` (logical), and `note`
#'   (character scalar).
sc_normalize_h5ad_sparse_layout <- function(encoding_type = NULL,
                                            dims = NULL,
                                            n_obs = NULL,
                                            n_vars = NULL) {
  .stop_layout <- function(msg) {
    stop(errorCondition(msg, class = "sc_bpcells_error"))
  }

  .validate_count <- function(val, label) {
    if (is.null(val)) return(NULL)
    if (!is.numeric(val) || length(val) != 1L || is.na(val) || !is.finite(val) || val < 0) {
      .stop_layout(sprintf("sc_normalize_h5ad_sparse_layout() : '%s' invalide.", label))
    }
    as.integer(val)
  }

  missing_enc <- is.null(encoding_type) ||
    length(encoding_type) == 0L ||
    (length(encoding_type) == 1L && (is.na(encoding_type) || !nzchar(trimws(as.character(encoding_type)))))

  if (missing_enc) {
    legacy_fallback <- TRUE
    layout <- "csr"
  } else {
    if (!is.character(encoding_type) || length(encoding_type) != 1L) {
      .stop_layout("sc_normalize_h5ad_sparse_layout() : 'encoding_type' doit être une chaîne scalaire ou NULL.")
    }
    legacy_fallback <- FALSE
    enc_norm <- tolower(trimws(encoding_type))
    layout <- switch(
      enc_norm,
      csr_matrix = "csr",
      csr = "csr",
      csc_matrix = "csc",
      csc = "csc",
      array = "dense",
      dense = "dense",
      .stop_layout(sprintf("sc_normalize_h5ad_sparse_layout() : encoding_type non supporté : '%s'.", encoding_type))
    )
  }

  if (!is.null(dims)) {
    if (!is.numeric(dims) || length(dims) != 2L || anyNA(dims) || any(!is.finite(dims)) || any(dims < 0)) {
      .stop_layout("sc_normalize_h5ad_sparse_layout() : 'dims' doit être un vecteur numérique de longueur 2 (>= 0).")
    }
    dims <- as.integer(dims)
  }
  n_obs <- .validate_count(n_obs, "n_obs")
  n_vars <- .validate_count(n_vars, "n_vars")

  needs_transpose <- TRUE
  orientation_note <- "default AnnData (obs x var) orientation assumed; transpose to (genes x cells)"

  if (!is.null(dims) && !is.null(n_obs) && !is.null(n_vars)) {
    if (n_obs != n_vars && identical(dims, c(n_vars, n_obs))) {
      needs_transpose <- FALSE
      orientation_note <- "already in (genes x cells) orientation; no transpose needed"
    } else if (identical(dims, c(n_obs, n_vars))) {
      needs_transpose <- TRUE
      orientation_note <- if (n_obs == n_vars) {
        "square matrix (n_obs == n_vars); assuming AnnData (obs x var) orientation"
      } else {
        "AnnData (obs x var) orientation detected; transpose to (genes x cells)"
      }
    } else {
      .stop_layout(sprintf(
        "sc_normalize_h5ad_sparse_layout() : dims (%d, %d) incompatible avec n_obs=%d, n_vars=%d.",
        dims[1L], dims[2L], n_obs, n_vars
      ))
    }
  }

  prefix <- if (legacy_fallback) {
    sprintf("legacy .h5ad fallback (missing encoding-type -> '%s')", layout)
  } else {
    sprintf("layout '%s'", layout)
  }

  list(
    layout = layout,
    needs_transpose = needs_transpose,
    legacy_fallback = legacy_fallback,
    note = paste0(prefix, ": ", orientation_note, ".")
  )
}


