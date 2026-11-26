# =============================================================================
# tools/pkg_skeleton/R/dataset_directory.R — PoC « répertoire de datasets »
# (option O3, moteur pur)
# =============================================================================
# Fiche (fait foi) :
# docs/proposals/PROPOSAL_DATASET_DIRECTORY_POC_2026-11-02.md
#   §1 AC1 (états valid|warn|invalid|empty, raisons = vocabulaire v1 tel quel),
#      AC2 (strictement lecture seule, jamais au cold-boot),
#      AC5 (RED-first, freeze test) ;
#   §3.1 (racine 3 niveaux, motif ts_session_checkpoint_dir(),
#         session_checkpoint.R:120) ;
#   §3.2 (moteur pur : UNE profondeur de sous-répertoires, délégation
#         biomentat_handoff_validate_drop() — ZÉRO nouvelle règle) ;
#   §5 (échec d'un drop isolé — jamais d'erreur propagée ; le scan n'écrit
#       rien, l'état applicatif reste dans global_data/ASR — vue passive).
#
# Ce que ce fichier N'EST PAS (fiche §2/§7) :
#   * PAS le contrat v2 (SC/spatial .h5ad hors PoC) ;
#   * PAS un watcher filesystem (O4 différé) ;
#   * PAS une ingestion : aucun appel au scan au cold-boot (AC2) — le
#     branchement UI bulk est le jalon P1, l'outil drive passif
#     transcripto_drive_datasets le jalon P2 ; tant que P1 n'est pas livré,
#     biomentat_handoff.R reste « étranger non câblé » dans app.R (bloc 03) —
#     la dépendance runtime du scan est résolue À L'APPEL, jamais au
#     chargement.
#
# Purity rails (règle maison, motif session_checkpoint.R:9-22) :
#   * NO shiny, NO reactive, NO session — ce fichier doit sourcer() partout où
#     tourne un process R pur ;
#   * STRICTEMENT LECTURE SEULE (AC2) : le scan n'écrit RIEN. La seule
#     écriture tolérée dans la chaîne est validation.json, écrite par le
#     moteur v1 biomentat_handoff_write_validation() UNIQUEMENT quand un
#     out_dir explicite est fourni (biomentat_handoff.R:755-757) — le scan ne
#     passe JAMAIS out_dir ;
#   * ZÉRO nouvelle règle de validation (AC1) : chaque sous-répertoire est
#     délégué tel quel au moteur v1 biomentat_handoff_validate_drop()
#     (biomentat_handoff.R:634, alias biomentat_validate_drop) ; les raisons
#     retournées sont les états v1 (biomentat_handoff_states(),
#     biomentat_handoff.R:104-123) repris verbatim ;
#   * la racine n'est JAMAIS créée ici (motif checkpoint : création au moment
#     de l'écriture — il n'y a pas d'écriture) ; une racine absente donne une
#     table vide, jamais une erreur ;
#   * helpers transverses réutilisés du voisin session_checkpoint.R (zéro
#     duplication) : .sc_project_root (:96, remontée app.R) et
#     .sc_scalar_text (:81, garde scalaire NA/vide -> NULL).
# =============================================================================

# ── Résolution de la racine des datasets (fiche §3.1) ────────────────────────

#' Resolve the datasets root directory (fiche PoC §3.1).
#'
#' Resolution order — motif `ts_session_checkpoint_dir()`
#' (session_checkpoint.R:120), arbitrage racine tranché (opérateur « continue »
#' + best-judgment délégué, jalon P0) :
#' (1) explicit `dir` argument (tests / appel programmatique) ;
#' (2) env `TRANSCRIPTO_DATASETS_DIR` (knob dédié) ;
#' (3) PREMIÈRE entrée du knob opérateur existant `TRANSCRIPTO_DRIVE_DATA_DIR`
#'     + "/datasets" (même lecture multi-chemins que le checkpoint,
#'     session_checkpoint.R:128-135) ;
#' (4) `<project root>/<data_folder>/datasets`, where `<data_folder>` depends on
#'     `modality` (`"sc"` -> `"data_sc"`, `"bulk"` -> `"data_bulk"`,
#'     `"spatial"` -> `"data_spatial"`), the project root resolved by the SAME
#'     app.R ascent as the checkpoint (`.sc_project_root()`,
#'     session_checkpoint.R:96 — reused, zero duplication).
#' Le répertoire n'est PAS créé ici (appelant lecture-seule, AC2) — aucun
#' effet de bord disque, jamais tempdir() seul.
#'
#' @param dir Racine explicite ; prioritaire quand non NULL/vide.
#' @param modality Modalité cible (`"sc"`, `"bulk"`, `"spatial"` ; défaut `"sc"`).
#' @return Chemin absolu (slashes avant), non créé.
ts_datasets_dir <- function(dir = NULL, modality = c("sc", "bulk", "spatial")) {
  modality <- match.arg(modality)
  d_arg <- .sc_scalar_text(dir)
  if (!is.null(d_arg)) {
    d_arg <- trimws(d_arg)
    if (nzchar(d_arg)) {
      return(normalizePath(d_arg, winslash = "/", mustWork = FALSE))
    }
  }
  env <- trimws(Sys.getenv("TRANSCRIPTO_DATASETS_DIR", ""))
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  data_env <- Sys.getenv("TRANSCRIPTO_DRIVE_DATA_DIR", "")
  if (nzchar(trimws(data_env))) {
    first <- trimws(strsplit(data_env, .Platform$path.sep, fixed = TRUE)[[1]][1])
    if (nzchar(first)) {
      return(normalizePath(file.path(first, "datasets"),
                           winslash = "/", mustWork = FALSE))
    }
  }
  data_folder <- switch(modality,
    sc      = "data_sc",
    bulk    = "data_bulk",
    spatial = "data_spatial"
  )
  normalizePath(file.path(.sc_project_root(), data_folder, "datasets"),
                winslash = "/", mustWork = FALSE)
}

# ── Dossiers personnalisés enregistrés (options persistantes par modalité) ──

#' Resolve the JSON config path storing user-saved dataset discovery folders.
#'
#' @param config_path Optional explicit config file path.
#' @return Character scalar path (not created on read).
ts_dataset_dirs_config_path <- function(config_path = NULL) {
  explicit <- .sc_scalar_text(config_path)
  if (!is.null(explicit) && nzchar(trimws(explicit))) {
    return(normalizePath(trimws(explicit), winslash = "/", mustWork = FALSE))
  }
  env_cfg <- trimws(Sys.getenv("BIOMENTAT_DATASET_DIRS_CONFIG", ""))
  if (nzchar(env_cfg)) {
    return(normalizePath(env_cfg, winslash = "/", mustWork = FALSE))
  }
  normalizePath(
    file.path(tools::R_user_dir("biomentat", "config"), "dataset_dirs.json"),
    winslash = "/", mustWork = FALSE
  )
}

.ts_read_dataset_dirs_cfg <- function(config_path = NULL) {
  p <- ts_dataset_dirs_config_path(config_path)
  if (!file.exists(p)) return(list())
  raw <- tryCatch(jsonlite::fromJSON(p, simplifyVector = TRUE), error = function(e) list())
  if (!is.list(raw)) list() else raw
}

.ts_write_dataset_dirs_cfg <- function(cfg, config_path = NULL) {
  p <- ts_dataset_dirs_config_path(config_path)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(cfg, p, auto_unbox = TRUE, pretty = TRUE)
  invisible(p)
}

#' List user-saved dataset discovery folders for a modality.
#'
#' @param modality Target modality (`"sc"`, `"bulk"`, `"spatial"`).
#' @param config_path Optional explicit config file path.
#' @return Character vector of normalized folder paths.
ts_saved_dataset_dirs <- function(modality = c("sc", "bulk", "spatial"),
                                  config_path = NULL) {
  modality <- match.arg(modality)
  cfg <- .ts_read_dataset_dirs_cfg(config_path)
  vals <- cfg[[modality]]
  if (is.null(vals) || length(vals) == 0L) return(character(0))
  vals <- trimws(as.character(unlist(vals, use.names = FALSE)))
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0L) return(character(0))
  unique(normalizePath(vals, winslash = "/", mustWork = FALSE))
}

#' Save a dataset discovery folder into persistent options for a modality.
#'
#' @param path Folder path to save.
#' @param modality Target modality (`"sc"`, `"bulk"`, `"spatial"`).
#' @param config_path Optional explicit config file path.
#' @return Updated character vector of saved folder paths for `modality`.
ts_save_dataset_dir <- function(path,
                                modality = c("sc", "bulk", "spatial"),
                                config_path = NULL) {
  modality <- match.arg(modality)
  p <- .sc_scalar_text(path)
  if (is.null(p) || !nzchar(trimws(p))) {
    stop(errorCondition("ts_save_dataset_dir() : chemin de dossier requis.",
                        class = "dataset_dir_error"))
  }
  p_norm <- normalizePath(trimws(p), winslash = "/", mustWork = FALSE)
  cfg <- .ts_read_dataset_dirs_cfg(config_path)
  cur <- ts_saved_dataset_dirs(modality = modality, config_path = config_path)
  cfg[[modality]] <- as.list(unique(c(cur, p_norm)))
  .ts_write_dataset_dirs_cfg(cfg, config_path = config_path)
  ts_saved_dataset_dirs(modality = modality, config_path = config_path)
}

#' Remove a saved dataset discovery folder for a modality.
#'
#' @param path Folder path to remove.
#' @param modality Target modality (`"sc"`, `"bulk"`, `"spatial"`).
#' @param config_path Optional explicit config file path.
#' @return Updated character vector of saved folder paths for `modality`.
ts_remove_dataset_dir <- function(path,
                                  modality = c("sc", "bulk", "spatial"),
                                  config_path = NULL) {
  modality <- match.arg(modality)
  p <- .sc_scalar_text(path)
  if (is.null(p) || !nzchar(trimws(p))) {
    return(ts_saved_dataset_dirs(modality = modality, config_path = config_path))
  }
  p_norm <- normalizePath(trimws(p), winslash = "/", mustWork = FALSE)
  cfg <- .ts_read_dataset_dirs_cfg(config_path)
  cur <- ts_saved_dataset_dirs(modality = modality, config_path = config_path)
  kept <- cur[cur != p_norm]
  cfg[[modality]] <- as.list(kept)
  .ts_write_dataset_dirs_cfg(cfg, config_path = config_path)
  kept
}

#' Resolve all dataset discovery roots (default + saved folders) for a modality.
#'
#' @param dir Optional explicit root(s) overriding default + saved.
#' @param modality Target modality (`"sc"`, `"bulk"`, `"spatial"`).
#' @param config_path Optional explicit config file path.
#' @return Character vector of normalized root paths.
.ts_unified_datasets_candidates <- function(modality = c("sc", "bulk", "spatial"),
                                            project_root = .sc_project_root()) {
  modality <- match.arg(modality)
  data_folder <- switch(modality,
    sc      = "data_sc",
    bulk    = "data_bulk",
    spatial = "data_spatial"
  )
  cands <- character(0)
  if (!is.null(project_root) && length(project_root) == 1L && nzchar(trimws(project_root))) {
    cands <- c(
      cands,
      file.path(project_root, "datasets", data_folder),
      file.path(project_root, "datasets", data_folder, "datasets")
    )
  }
  env_root <- trimws(Sys.getenv("TRANSCRIPTO_DATASETS_DIR", ""))
  if (nzchar(env_root)) {
    cands <- c(
      cands,
      file.path(env_root, data_folder),
      file.path(env_root, data_folder, "datasets")
    )
  }
  cands <- unique(normalizePath(cands, winslash = "/", mustWork = FALSE))
  cands[dir.exists(cands)]
}

ts_datasets_dirs <- function(dir = NULL,
                             modality = c("sc", "bulk", "spatial"),
                             config_path = NULL) {
  modality <- match.arg(modality)
  if (!is.null(dir) && length(dir) > 0L) {
    vals <- trimws(as.character(dir))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    if (length(vals) > 0L) {
      return(unique(normalizePath(vals, winslash = "/", mustWork = FALSE)))
    }
  }
  unique(c(
    ts_datasets_dir(modality = modality),
    .ts_unified_datasets_candidates(modality = modality),
    ts_saved_dataset_dirs(modality = modality, config_path = config_path)
  ))
}

.ts_scan_multi_roots <- function(roots, scan_one_fn) {
  if (is.null(roots) || length(roots) == 0L) return(list())
  roots <- unique(trimws(as.character(roots)))
  roots <- roots[!is.na(roots) & nzchar(roots)]
  out <- list()
  for (r in roots) {
    part <- scan_one_fn(r)
    if (length(part) == 0L) next
    for (nm in names(part)) {
      key <- if (nm %in% names(out)) sprintf("%s [%s]", nm, basename(r)) else nm
      out[[key]] <- part[[nm]]
    }
  }
  out
}

# ── Scanner passif (fiche §3.2 — UNE profondeur, délégation v1) ──────────────

#' Frozen status vocabulary of the scan table (AC1).
#'
#' Motif `biomentat_handoff_states()` (biomentat_handoff.R:104) : le vocabulaire
#' est exposé pour être figé par le test T5 (motif test-report-contract-freeze.R).
#' "empty" = sous-répertoire sans rien à valider (aucun fichier, aucun
#' sous-répertoire) ; sa raison "empty_directory" est la SEULE raison hors
#' vocabulaire v1 (documentée fiche §1 AC1 / §3.1 — un dépot vide n'est pas
#' une validation échouée).
ts_dataset_dir_scan_states <- function() {
  c("valid", "warn", "invalid", "empty")
}

# Chemin absolu à slashes avant, sans exiger l'existence (mustWork = FALSE).
.ts_datasets_path <- function(...) {
  normalizePath(file.path(...), winslash = "/", mustWork = FALSE)
}

# Champ scalaire du contrat de sortie : NA_character_ quand indisponible
# (schéma figé T5 — jamais NULL, pour un consumeur UI/JSON homogène).
.ts_datasets_field <- function(x) {
  v <- .sc_scalar_text(x)
  if (is.null(v)) NA_character_ else v
}

# Lecture du manifeste v1 UNE fois, par délégation au lecteur EXISTANT
# biomentat_handoff_read_manifest() (biomentat_handoff.R:142) — zéro nouvelle
# règle de parsing ; tout échec -> NULL, jamais d'erreur propagée (§5). Le
# résultat sert le résumé ET la résolution des chemins (pas de double I/O).
.ts_datasets_read_manifest <- function(drop_dir) {
  tryCatch(biomentat_handoff_read_manifest(drop_dir), error = function(e) NULL)
}

.ts_bulk_detect_loose_pair <- function(drop_dir) {
  files <- sort(list.files(drop_dir, full.names = TRUE, recursive = FALSE))
  files <- files[!dir.exists(files)]
  tab_files <- files[grepl("\\.(tsv|csv|txt|tab)(\\.gz)?$", basename(files), ignore.case = TRUE)]
  if (length(tab_files) == 0L) return(NULL)
  bnames <- tolower(basename(tab_files))
  is_meta <- grepl("(meta|sample|coldata|design|pheno|clin)", bnames)
  meta_files <- tab_files[is_meta]
  non_meta <- tab_files[!is_meta]
  mat_files <- non_meta[grepl("(count|matrix|expr|raw|umi|read|gene)", tolower(basename(non_meta)))]
  mat_path <- if (length(mat_files) > 0L) {
    mat_files[1L]
  } else if (length(non_meta) > 0L && (length(meta_files) > 0L || length(non_meta) == 1L)) {
    non_meta[1L]
  } else {
    NULL
  }
  if (is.null(mat_path)) return(NULL)
  list(
    matrix = .ts_datasets_path(mat_path),
    samplesheet = if (length(meta_files) > 0L) .ts_datasets_path(meta_files[1L]) else NA_character_
  )
}

# Scan d'UN sous-répertoire (UNE profondeur — fiche §3.2/§5) : résumé du
# manifeste + chemins résolus, mtime, puis délégation de la validation au
# moteur v1. Échec d'un drop = entrée invalid isolée, les autres continuent.
.ts_datasets_scan_one <- function(drop_dir, validate = TRUE, allow_loose = FALSE) {
  dir_norm <- .ts_datasets_path(drop_dir)
  man <- .ts_datasets_read_manifest(drop_dir)

  entry <- list(
    status = NA_character_,
    raison = NA_character_,
    manifest = list(
      dataset_id = .ts_datasets_field(man$dataset_id),
      species = .ts_datasets_field(man$species),
      feature_id_type = .ts_datasets_field(man$feature_id_type),
      matrix_type = .ts_datasets_field(man$matrix_type)
    ),
    paths = list(
      dir = dir_norm,
      manifest = .ts_datasets_path(drop_dir, "manifest.json"),
      matrix = NA_character_,
      samplesheet = NA_character_
    ),
    mtime = file.mtime(dir_norm),
    checks = list()
  )
  # chemins de matrice/samplesheet résolus depuis le manifeste v1 (chemins
  # RELATIFS exigés par le contrat v1 — ici on ne fait qu'assembler, la
  # conformité des chemins reste la validation déléguée)
  if (!is.null(man)) {
    mp <- .sc_scalar_text(man$matrix_path)
    sp <- .sc_scalar_text(man$samplesheet_path)
    if (!is.null(mp)) entry$paths$matrix <- .ts_datasets_path(drop_dir, mp)
    if (!is.null(sp)) entry$paths$samplesheet <- .ts_datasets_path(drop_dir, sp)
  }

  # sous-répertoire sans rien à valider -> "empty" (AC1 ; seule raison hors v1)
  has_files <- length(list.files(drop_dir, all.files = FALSE, no.. = TRUE)) > 0L
  has_subdirs <- length(list.dirs(drop_dir, full.names = FALSE,
                                  recursive = FALSE)) > 0L
  if (!has_files && !has_subdirs) {
    entry$status <- "empty"
    entry$raison <- "empty_directory"
    return(entry)
  }

  # validate = FALSE : listing brut (table rapide pour UI/drive), SANS
  # validation — status/raison restent NA, résumé + chemins portés
  if (!isTRUE(validate)) {
    return(entry)
  }

  # Délégation v1 — ZÉRO nouvelle règle (AC1). out_dir jamais fourni : le
  # scan n'écrit RIEN (AC2 ; la seule écriture tolérée de la chaîne,
  # validation.json, n'est déclenchée que par un out_dir explicite,
  # biomentat_handoff.R:755-757).
  res <- tryCatch(
    biomentat_handoff_validate_drop(drop_dir = drop_dir, out_dir = NULL),
    error = function(e) NULL
  )
  if (is.null(res)) {
    # filet §5 : un drop qui ferait tomber le moteur lui-même est isolé —
    # raison = l'état catch-all du moteur v1 ("file_corrupted",
    # biomentat_handoff.R:727-728), jamais une erreur propagée
    entry$status <- "invalid"
    entry$raison <- "file_corrupted"
    return(entry)
  }

  entry$checks <- if (is.null(res$checks)) list() else res$checks
  if (isTRUE(res$valid)) {
    entry$status <- if (length(res$warnings) > 0L) "warn" else "valid"
  } else {
    entry$status <- "invalid"
    # raison = PREMIER code d'erreur v1, repris tel quel (AC1) — l'ordre des
    # erreurs est celui du moteur v1 (étapes manifeste -> matrice -> samplesheet)
    first_err <- if (length(res$errors) > 0L) res$errors[[1L]]$code else NULL
    entry$raison <- .ts_datasets_field(first_err)
    if (isTRUE(allow_loose) && identical(entry$raison, "manifest_missing")) {
      pair <- .ts_bulk_detect_loose_pair(drop_dir)
      if (!is.null(pair)) {
        entry$paths$matrix <- pair$matrix
        entry$paths$samplesheet <- pair$samplesheet
        if (!is.na(pair$samplesheet) && nzchar(pair$samplesheet)) {
          entry$status <- "valid"
          entry$raison <- NA_character_
          entry$manifest$matrix_type <- "counts+metadata"
        } else {
          entry$status <- "warn"
          entry$raison <- "metadata_missing"
          entry$manifest$matrix_type <- "counts_only"
        }
      }
    }
  }
  # dataset_id : la parole est au manifeste VALIDÉ par le moteur v1
  validated_id <- .sc_scalar_text(res$dataset_id)
  if (!is.null(validated_id)) entry$manifest$dataset_id <- validated_id
  # enrichissement du résumé depuis les checks v1 (déclaré vs détecté)
  if (!is.null(entry$checks$species_and_identifiers)) {
    spc <- .sc_scalar_text(entry$checks$species_and_identifiers$declared_species)
    fid <- .sc_scalar_text(entry$checks$species_and_identifiers$declared_id_type)
    if (!is.null(spc)) entry$manifest$species <- spc
    if (!is.null(fid)) entry$manifest$feature_id_type <- fid
  }
  if (!is.null(entry$checks$matrix_semantics)) {
    mtp <- .sc_scalar_text(entry$checks$matrix_semantics$matrix_type)
    if (!is.null(mtp)) entry$manifest$matrix_type <- mtp
  }
  entry
}

#' Scan the datasets root directory (fiche PoC §3.2 — passive, read-only).
#'
#' Lists ONE level of subdirectories of `datasets_dir` and delegates each one
#' to the existing v1 engine `biomentat_handoff_validate_drop()`
#' (biomentat_handoff.R:634) — the scan adds ZERO validation rule (AC1). The
#' result is a NAMED LIST keyed by the v1 `dataset_id` (fallback: the
#' subdirectory basename when no manifest is readable), each entry carrying:
#' `status` (valid | warn | invalid | empty — AC1, `ts_dataset_dir_scan_states()`),
#' `raison` (v1 vocabulary verbatim — first v1 error code; "empty_directory"
#' for an empty subdirectory; NA when valid or unvalidated), `manifest`
#' (dataset_id/species/feature_id_type/matrix_type summary), `paths`
#' (dir/manifest/matrix/samplesheet — resolved from the v1 manifest),
#' `mtime` (subdirectory mtime) and `checks` (the v1 checks, inherited
#' verbatim). With `validate = FALSE` the table is a bare listing:
#' status/raison stay NA, the manifest summary and paths are still carried.
#'
#' STRICT READ-ONLY (AC2) : nothing is written anywhere — the engine's only
#' write path (validation.json) requires an explicit `out_dir` which the scan
#' NEVER passes ; a failing subdirectory is isolated (fiche §5, never a
#' propagated error) ; a missing or non-directory root yields an empty table
#' (the scan never blocks and never creates anything).
#'
#' @param datasets_dir Root to scan (default: `ts_datasets_dir()`).
#' @param validate Delegate each subdirectory to the v1 engine (default TRUE).
#' @return Named list of per-dataset entries (see details); empty list for an
#'   empty or missing root.
ts_dataset_dir_scan <- function(datasets_dir = ts_datasets_dirs(), validate = TRUE, allow_loose = FALSE) {
  .ts_scan_multi_roots(datasets_dir, function(one_root) {
    root <- .sc_scalar_text(one_root)
    if (is.null(root) || !dir.exists(root)) {
      return(list())
    }
    old_coll <- Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", old_coll), add = TRUE)
    subdirs <- sort(list.dirs(root, full.names = TRUE, recursive = FALSE))

    out <- list()
    for (d in subdirs) {
      entry <- .ts_datasets_scan_one(d, validate = validate, allow_loose = allow_loose)
      key <- .sc_scalar_text(entry$manifest$dataset_id)
      if (is.null(key)) key <- basename(d)
      out[[key]] <- entry
    }
    out
  })
}

# ── Scanner passif Single-Cell (10X / .rds / .h5 / .h5ad / .loom) ───────────

# Détecte un répertoire 10X direct (filtered_feature_bc_matrix.h5 ou triplet
# matrix.mtx(.gz) + barcodes.tsv(.gz) + features.tsv(.gz)|genes.tsv(.gz)).
# Lecture seule stricte (ne crée jamais features.tsv.gz). Retourne NULL si le
# dossier ne contient ni h5 10x ni matrix.mtx(.gz).
.ts_sc_classify_10x_dir <- function(d) {
  h5 <- file.path(d, "filtered_feature_bc_matrix.h5")
  if (file.exists(h5)) {
    sz <- file.size(h5)
    if (is.na(sz) || sz <= 0) {
      return(list(valid = FALSE, format = NA_character_, raison = "empty_file"))
    }
    return(list(valid = TRUE, format = "10x_h5", raison = NA_character_))
  }
  has_mat <- any(file.exists(file.path(d, c("matrix.mtx", "matrix.mtx.gz"))))
  if (!has_mat) return(NULL)
  has_bc <- any(file.exists(file.path(d, c("barcodes.tsv", "barcodes.tsv.gz"))))
  has_feat <- any(file.exists(file.path(d, c("features.tsv", "features.tsv.gz"))))
  has_gene <- any(file.exists(file.path(d, c("genes.tsv", "genes.tsv.gz"))))
  if (has_bc && (has_feat || has_gene)) {
    fmt <- if (has_feat) "10x_v3" else "10x_v2"
    return(list(valid = TRUE, format = fmt, raison = NA_character_))
  }
  list(valid = FALSE, format = NA_character_, raison = "incomplete_10x_triplet")
}

# Détecte un drop Bulk biomentat-exchange/v1 afin de l'isoler proprement comme
# invalid (raison = "non_sc_format") lorsque data_sc/datasets en contient.
.ts_sc_is_bulk_drop <- function(d) {
  has_man <- file.exists(file.path(d, "manifest.json"))
  has_cnt <- file.exists(file.path(d, "counts.tsv"))
  has_sht <- file.exists(file.path(d, "samplesheet.tsv"))
  if (has_man) {
    man <- .ts_datasets_read_manifest(d)
    sv <- if (is.list(man)) .sc_scalar_text(man$schema_version) else NULL
    if (!is.null(sv) && grepl("biomentat-exchange", sv, fixed = TRUE)) {
      return(TRUE)
    }
  }
  isTRUE(has_cnt && has_sht)
}

.ts_sc_scan_one_dir <- function(d) {
  d_norm <- .ts_datasets_path(d)
  mt <- file.mtime(d_norm)
  has_files <- length(list.files(d, all.files = FALSE, no.. = TRUE)) > 0L
  kids <- sort(list.dirs(d, full.names = TRUE, recursive = FALSE))
  if (!has_files && length(kids) == 0L) {
    return(list(
      status = "empty", raison = "empty_directory", format = NA_character_,
      n_samples = 0L,
      paths = list(dir = d_norm, file = NA_character_, sample_dirs = character(0)),
      mtime = mt
    ))
  }

  tenx <- .ts_sc_classify_10x_dir(d)
  if (!is.null(tenx)) {
    if (isTRUE(tenx$valid)) {
      return(list(
        status = "valid", raison = NA_character_, format = tenx$format,
        n_samples = 1L,
        paths = list(dir = d_norm, file = NA_character_, sample_dirs = d_norm),
        mtime = mt
      ))
    }
    return(list(
      status = "invalid", raison = tenx$raison, format = NA_character_,
      n_samples = 0L,
      paths = list(dir = d_norm, file = NA_character_, sample_dirs = character(0)),
      mtime = mt
    ))
  }

  if (length(kids) > 0L) {
    valid_kids <- character(0)
    for (k in kids) {
      k_cls <- .ts_sc_classify_10x_dir(k)
      if (!is.null(k_cls) && isTRUE(k_cls$valid)) {
        valid_kids <- c(valid_kids, .ts_datasets_path(k))
      }
    }
    if (length(valid_kids) > 0L) {
      return(list(
        status = "valid", raison = NA_character_, format = "10x_multi",
        n_samples = length(valid_kids),
        paths = list(dir = d_norm, file = NA_character_, sample_dirs = valid_kids),
        mtime = mt
      ))
    }
  }

  if (.ts_sc_is_bulk_drop(d)) {
    return(list(
      status = "invalid", raison = "non_sc_format", format = "bulk_exchange",
      n_samples = 0L,
      paths = list(dir = d_norm, file = NA_character_, sample_dirs = character(0)),
      mtime = mt
    ))
  }

  rds <- sort(list.files(d, pattern = "\\.rds$", ignore.case = TRUE, full.names = TRUE))
  if (length(rds) == 1L && !dir.exists(rds[1L])) {
    sz <- file.size(rds[1L])
    if (is.na(sz) || sz <= 0) {
      return(list(
        status = "invalid", raison = "empty_file", format = NA_character_,
        n_samples = 0L,
        paths = list(dir = d_norm, file = .ts_datasets_path(rds[1L]), sample_dirs = character(0)),
        mtime = mt
      ))
    }
    return(list(
      status = "valid", raison = NA_character_, format = "rds_dir",
      n_samples = 1L,
      paths = list(dir = d_norm, file = .ts_datasets_path(rds[1L]), sample_dirs = character(0)),
      mtime = mt
    ))
  }

  list(
    status = "invalid", raison = "no_sc_data_found", format = NA_character_,
    n_samples = 0L,
    paths = list(dir = d_norm, file = NA_character_, sample_dirs = character(0)),
    mtime = mt
  )
}

#' Scan the Single-Cell datasets root directory (passive, strictly read-only).
#'
#' Lists ONE depth of subdirectories (10x_h5, 10x_v3, 10x_v2, 10x_multi, rds_dir,
#' empty_directory, incomplete_10x_triplet, Bulk exchange non_sc_format, or
#' no_sc_data_found) and top-level SC files (.rds, .h5, .h5ad, .loom).
#' Missing or non-directory root yields `list()`, never creates directories or
#' modifies files.
#'
#' @param datasets_dir Root to scan (default: `ts_datasets_dir(modality = "sc")`).
#' @return Named list of discovered SC dataset entries.
ts_sc_dataset_dir_scan <- function(datasets_dir = ts_datasets_dirs(modality = "sc")) {
  .ts_scan_multi_roots(datasets_dir, function(one_root) {
    root <- .sc_scalar_text(one_root)
    if (is.null(root) || !dir.exists(root)) {
      return(list())
    }
    old_coll <- Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", old_coll), add = TRUE)

    subdirs <- sort(list.dirs(root, full.names = TRUE, recursive = FALSE))
    top_files <- sort(list.files(root, full.names = TRUE, recursive = FALSE))
    top_files <- top_files[!dir.exists(top_files) &
                             grepl("\\.(rds|h5|h5ad|loom)$", basename(top_files),
                                   ignore.case = TRUE)]

    out <- list()
    for (d in subdirs) {
      out[[basename(d)]] <- .ts_sc_scan_one_dir(d)
    }
    for (f in top_files) {
      f_norm <- .ts_datasets_path(f)
      sz <- file.size(f)
      ext <- tolower(tools::file_ext(f))
      ok <- !is.na(sz) && sz > 0
      out[[basename(f)]] <- list(
        status = if (ok) "valid" else "invalid",
        raison = if (ok) NA_character_ else "empty_file",
        format = ext,
        n_samples = if (ok) 1L else 0L,
        paths = list(dir = NA_character_, file = f_norm, sample_dirs = character(0)),
        mtime = file.mtime(f_norm)
      )
    }
    out
  })
}

# ── Scanner Spatial (data_spatial/datasets) ───────────────────────────────────

.ts_spatial_scan_one_dir <- function(d) {
  d_norm <- .ts_datasets_path(d)
  mt <- file.mtime(d_norm)
  fnames <- list.files(d, full.names = FALSE)
  has_visium <- dir.exists(file.path(d, "spatial")) ||
    dir.exists(file.path(d, "outs", "spatial")) ||
    dir.exists(file.path(d, "binned_outputs")) ||
    any(grepl("^(filtered|raw)_feature_bc_matrix", fnames))
  if (has_visium) {
    vis_dir <- if (dir.exists(file.path(d, "outs", "spatial")) && !dir.exists(file.path(d, "spatial"))) {
      .ts_datasets_path(d, "outs")
    } else {
      d_norm
    }
    return(list(
      status = "valid",
      raison = NA_character_,
      format = "visium",
      paths = list(dir = vis_dir, file = NA_character_),
      mtime = mt
    ))
  }
  if (any(grepl("cell_feature_matrix", fnames)) &&
      (any(grepl("^cells\\.(parquet|csv(\\.gz)?)$", fnames)) ||
       dir.exists(file.path(d, "cell_feature_matrix")))) {
    return(list(
      status = "valid",
      raison = NA_character_,
      format = "xenium",
      paths = list(dir = d_norm, file = NA_character_),
      mtime = mt
    ))
  }
  if (any(grepl("exprMat", fnames, ignore.case = TRUE)) &&
      any(grepl("metadata", fnames, ignore.case = TRUE))) {
    return(list(
      status = "valid",
      raison = NA_character_,
      format = "cosmx",
      paths = list(dir = d_norm, file = NA_character_),
      mtime = mt
    ))
  }
  if (any(grepl("cell_by_gene", fnames, ignore.case = TRUE)) &&
      any(grepl("cell_metadata", fnames, ignore.case = TRUE))) {
    return(list(
      status = "valid",
      raison = NA_character_,
      format = "merscope",
      paths = list(dir = d_norm, file = NA_character_),
      mtime = mt
    ))
  }
  if (length(fnames) >= 2L &&
      any(grepl("(bead|coord|location|spatial)", fnames, ignore.case = TRUE)) &&
      any(grepl("\\.(csv|tsv|txt|h5|rds)(\\.gz)?$", fnames, ignore.case = TRUE))) {
    return(list(
      status = "valid",
      raison = NA_character_,
      format = "slideseq",
      paths = list(dir = d_norm, file = NA_character_),
      mtime = mt
    ))
  }
  files <- sort(list.files(d, pattern = "\\.(rds|h5ad|h5)$", ignore.case = TRUE, full.names = TRUE))
  files <- files[!dir.exists(files)]
  if (length(files) > 0L) {
    f1 <- .ts_datasets_path(files[1L])
    sz <- file.size(f1)
    ok <- !is.na(sz) && sz > 0
    return(list(
      status = if (ok) "valid" else "invalid",
      raison = if (ok) NA_character_ else "empty_file",
      format = tolower(tools::file_ext(f1)),
      paths = list(dir = d_norm, file = f1),
      mtime = mt
    ))
  }
  list(
    status = "invalid",
    raison = "no_spatial_files",
    format = "unknown",
    paths = list(dir = d_norm, file = NA_character_),
    mtime = mt
  )
}

#' Scan a Spatial datasets directory (1 depth, read-only).
#'
#' @param datasets_dir Root to scan (default: `ts_datasets_dir(modality = "spatial")`).
#' @return Named list of discovered Spatial dataset entries.
ts_spatial_dataset_dir_scan <- function(datasets_dir = ts_datasets_dirs(modality = "spatial")) {
  .ts_scan_multi_roots(datasets_dir, function(one_root) {
    root <- .sc_scalar_text(one_root)
    if (is.null(root) || !dir.exists(root)) {
      return(list())
    }
    old_coll <- Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", old_coll), add = TRUE)

    subdirs <- sort(list.dirs(root, full.names = TRUE, recursive = FALSE))
    top_files <- sort(list.files(root, full.names = TRUE, recursive = FALSE))
    top_files <- top_files[!dir.exists(top_files) &
                             grepl("\\.(rds|h5ad|h5)$", basename(top_files), ignore.case = TRUE)]

    out <- list()
    for (d in subdirs) {
      out[[basename(d)]] <- .ts_spatial_scan_one_dir(d)
    }
    for (f in top_files) {
      f_norm <- .ts_datasets_path(f)
      sz <- file.size(f)
      ext <- tolower(tools::file_ext(f))
      ok <- !is.na(sz) && sz > 0
      out[[basename(f)]] <- list(
        status = if (ok) "valid" else "invalid",
        raison = if (ok) NA_character_ else "empty_file",
        format = ext,
        paths = list(dir = NA_character_, file = f_norm),
        mtime = file.mtime(f_norm)
      )
    }
    out
  })
}
