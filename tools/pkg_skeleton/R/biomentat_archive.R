# =============================================================================
# tools/pkg_skeleton/R/biomentat_archive.R — P3.3: Portable .biomentat archive
# =============================================================================
# Pure-R self-contained .biomentat bundle writer/reader carrying manifest.json
# (biomentat-exchange/v1 QC fields + SHA-256 integrity hash) and payload.rds.
# Uses base R gzip tar so it works cross-platform with zero external binaries.
# =============================================================================

.archive_stop <- function(msg) {
  stop(errorCondition(msg, class = "biomentat_archive_error"))
}

#' Create a portable .biomentat archive with QC manifest and SHA-256 checksum
#'
#' @param payload Non-null R object (typically a list/handoff bundle) to serialize.
#' @param out_path Destination file path ending in `.biomentat`.
#' @param qc_metrics Named list of QC summary fields (e.g. `n_genes`, `n_samples`).
#' @param metadata Named list of descriptive metadata (e.g. `modality`, `organism`).
#' @return Manifest list written into the archive.
create_biomentat_archive <- function(payload,
                                     out_path,
                                     qc_metrics = list(),
                                     metadata = list()) {
  if (is.null(payload)) {
    .archive_stop("create_biomentat_archive() : payload ne peut pas etre NULL.")
  }
  if (!is.character(out_path) || length(out_path) != 1L || !nzchar(out_path) ||
      !grepl("\\.biomentat$", out_path, ignore.case = TRUE)) {
    .archive_stop("create_biomentat_archive() : out_path doit se terminer par .biomentat.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE) ||
      !requireNamespace("digest", quietly = TRUE)) {
    .archive_stop("create_biomentat_archive() : les packages 'jsonlite' et 'digest' sont requis.")
  }

  stage_dir <- tempfile("biomentat_arch_")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)

  payload_file <- file.path(stage_dir, "payload.rds")
  saveRDS(payload, payload_file)
  sha <- digest::digest(payload_file, algo = "sha256", file = TRUE)

  manifest <- list(
    schema = "biomentat-archive/v1",
    exchange_schema = "biomentat-exchange/v1",
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    payload_sha256 = sha,
    qc_metrics = qc_metrics,
    metadata = metadata
  )
  manifest_file <- file.path(stage_dir, "manifest.json")
  jsonlite::write_json(manifest, manifest_file, auto_unbox = TRUE, pretty = TRUE)

  out_abs <- normalizePath(out_path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(out_abs), recursive = TRUE, showWarnings = FALSE)
  status <- withr::with_dir(stage_dir, {
    utils::tar(tarfile = out_abs, files = c("manifest.json", "payload.rds"),
               compression = "gzip", tar = "internal")
  })
  if (!identical(status, 0L) || !file.exists(out_abs)) {
    .archive_stop("create_biomentat_archive() : echec de la creation de l'archive tar.gz.")
  }
  manifest
}

#' Read and verify a portable .biomentat archive
#'
#' @param archive_path Path to a `.biomentat` archive.
#' @param verify_hash Logical; if TRUE, verify `payload.rds` SHA-256 against manifest.
#' @return List with `manifest` and `payload`.
read_biomentat_archive <- function(archive_path, verify_hash = TRUE) {
  if (!is.character(archive_path) || length(archive_path) != 1L || !file.exists(archive_path)) {
    .archive_stop("read_biomentat_archive() : fichier archive introuvable.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE) ||
      !requireNamespace("digest", quietly = TRUE)) {
    .archive_stop("read_biomentat_archive() : les packages 'jsonlite' et 'digest' sont requis.")
  }

  stage_dir <- tempfile("biomentat_unarch_")
  dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage_dir, recursive = TRUE), add = TRUE)

  utils::untar(archive_path, exdir = stage_dir, tar = "internal")
  manifest_file <- file.path(stage_dir, "manifest.json")
  payload_file <- file.path(stage_dir, "payload.rds")
  if (!file.exists(manifest_file) || !file.exists(payload_file)) {
    .archive_stop("read_biomentat_archive() : archive incomplete (manifest.json ou payload.rds manquant).")
  }

  manifest <- jsonlite::read_json(manifest_file, simplifyVector = TRUE)
  if (isTRUE(verify_hash)) {
    actual_sha <- digest::digest(payload_file, algo = "sha256", file = TRUE)
    if (!identical(actual_sha, manifest$payload_sha256)) {
      .archive_stop("read_biomentat_archive() : echec de verification SHA-256 de payload.rds.")
    }
  }

  list(manifest = manifest, payload = readRDS(payload_file))
}
