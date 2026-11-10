# =============================================================================
# tools/pkg_skeleton/R/packs.R — PK-1 : pilier 1 « Système de packs » (socle)
# =============================================================================
# Fiche (fait foi) :
#   * docs/proposals/PROPOSAL_PACKS_ECOSYSTEM_2026-11-09.md
#       §3 PK-1 (lecteur manifest.json PUR — validation stricte, refus
#           honnêtes ; scan au démarrage ; registre session) ;
#       §4/D2 (emplacement : data_sc/packs/<pack>/manifest.json — zone
#           utilisateur, même philosophie que datasets/ ; config/ reste
#           code-app) ;
#       §4/D6 (discipline renv : les packs embarquent les DONNÉES, zéro
#           nouveau package — ce fichier n'appelle que jsonlite, déjà au
#           lock) ;
#       §5 (non-goals respectés ici : PAS d'outil MCP, PAS d'API cloud,
#           PAS de bump clusterProfiler, PAS de refonte des modules pathway) ;
#       §6 (RED→GREEN : fixture pack valide/invalide — refus honnêtes :
#           manifest manquant, type inconnu, fichier GMT absent) ;
#   * docs/ideation/PACKS_ECOSYSTEM_ENRICHR_PALETTES.md
#       §3 (schéma manifest « libraries » : type/id/name/version/source/
#           downloaded/attribution/files ; « Déposer un dossier = installer
#           une capacité. Aucun code applicatif à modifier » ; invalid
#           manifest → pack skipped avec message visible, JAMAIS un crash ;
#           lecture seule côté app) ;
#       §7.1 (schéma manifest « palette » : hex embarqués + attribution +
#           license) ;
#       §8 (data_sc/packs/<pack>/manifest.json) ;
#       §9 (champ attribution/license obligatoire pour les palettes — ne pas
#           embarquer sans vérifier).
#
# Périmètre PK-1 STRICT (ce que ce fichier N'EST PAS) :
#   * PAS l'ORA sur GMT (enricher() = PK-2, pilier 2) ;
#   * PAS la couche de convergence (PK-3) ;
#   * PAS le branchement palettes -> ts_theme/sélecteurs (PK-4, pilier 5) ;
#   * PAS le mode API Enrichr (PK-5 opt-in).
#   Ce fichier est le SOCLE lecteur : type « libraries » et « palette » sont
#   validés ici (le schéma des manifests est la matière de PK-1), leur
#   CONSOMMATION applicative vient dans les jalons suivants.
#
# Purity rails (règle maison, motif sc_multi.R / dataset_directory.R) :
#   * NO shiny, NO reactive, NO session — ce fichier doit sourcer() partout
#     où tourne un process R pur ;
#   * STRICTEMENT LECTURE SEULE (motif AC2 dataset_directory.R) : le lecteur
#     et le scan n'écrivent RIEN sur disque — data_sc/packs/ appartient à
#     l'utilisateur (app read-only side, fiche idéation §3) ;
#   * la racine n'est JAMAIS créée : une racine absente donne un résultat
#     vide, jamais une erreur (motif ts_dataset_dir_scan) ;
#   * helpers transverses réutilisés (zéro duplication) : .sc_project_root()
#     et .sc_scalar_text() du voisin session_checkpoint.R (:81/:96) — même
#     motif de réutilisation que dataset_directory.R:43-45 ;
#   * registre session PUR (motif sc_multi.R) : packs_register()/remove()
#     retournent la NOUVELLE liste nommée invisiblement et ne mutent JAMAIS
#     un global — l'appelant assigne global_data$packs lui-même ;
#   * un échec de lecture est un refus ISOLÉ porteur d'un état classé
#     packs_error + d'un message qui dit POURQUOI (fiche §6) — le scan
#     n'propage jamais d'erreur (fiche idéation §3/§9 : skipped with a
#     visible status message, never a crash).
#
# Consommateurs :
#   * app.R (bloc 03 pour le source(), server() pour le scan de session) ;
#   * PK-2 (enricher() sur les GMT des packs « libraries »), PK-4 (palettes
#     -> ts_theme), PK-3 (convergence) — via l'API pure ci-dessous.
# =============================================================================

#' Public API surface (gelée par tests/testthat/test-packs-manifest.R)
#' @return Character vector of exported function names.
packs_public_api <- function() {
  c(
    "packs_public_api", "packs_error_states", "packs_types",
    "packs_dir", "pack_read_manifest", "pack_scan_directory",
    "packs_register", "packs_remove", "packs_get", "packs_summary"
  )
}

#' Frozen error states (fiche PK-1 §6 — chaque refus porte un état + un pourquoi)
#' @return Character vector of `state` attribute values.
packs_error_states <- function() {
  c("invalid_input", "missing_manifest", "malformed_json", "missing_field",
    "invalid_field", "unknown_type", "missing_file", "duplicate_id",
    "unknown_id", "scan_failed")
}

#' Pack types implemented by THIS reader (fiche idéation §3)
#'
#' Le schéma autorise des types futurs (« template de rapport, … ») : un type
#' non implémenté est un REFUS honnête (unknown_type), jamais un chargement
#' silencieux — le pack apparaît dans les refus du scan avec le message.
#' @return Character vector of implemented type values.
packs_types <- function() {
  c("libraries", "palette")
}

#' Classed error constructor (motif .sc_multi_stop, sc_multi.R:65)
#' @param msg French error message (says WHY).
#' @param state One of packs_error_states().
#' @return An errorCondition with class `packs_error`.
.packs_stop <- function(msg, state) {
  errorCondition(msg, class = "packs_error", state = state)
}

# ── Résolution de la racine des packs (fiche D2 — motif ts_datasets_dir) ─────

#' Resolve the packs root directory (décision D2 : data_sc/packs/).
#'
#' Resolution order — motif `ts_datasets_dir()` (dataset_directory.R:68),
#' même arbitrage racine que le PoC datasets :
#' (1) explicit `dir` argument (tests / appel programmatique) ;
#' (2) env `TRANSCRIPTO_PACKS_DIR` (knob dédié) ;
#' (3) PREMIÈRE entrée du knob opérateur existant `TRANSCRIPTO_DRIVE_DATA_DIR`
#'     + "/packs" (même lecture multi-chemins que le checkpoint) ;
#' (4) `<project root>/data_sc/packs` — décision D2 de la fiche §4.
#' Le répertoire n'est PAS créé ici (appelant lecture-seule) — aucun effet de
#' bord disque, jamais tempdir() seul.
#'
#' @param dir Racine explicite ; prioritaire quand non NULL/vide.
#' @return Chemin absolu (slashes avant), non créé.
packs_dir <- function(dir = NULL) {
  if (!is.null(.sc_scalar_text(dir))) {
    return(normalizePath(dir, winslash = "/", mustWork = FALSE))
  }
  env <- Sys.getenv("TRANSCRIPTO_PACKS_DIR", "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = FALSE))
  }
  data_env <- Sys.getenv("TRANSCRIPTO_DRIVE_DATA_DIR", "")
  if (nzchar(data_env)) {
    first <- strsplit(data_env, .Platform$path.sep, fixed = TRUE)[[1]][1]
    if (nzchar(first)) {
      return(normalizePath(file.path(first, "packs"),
                           winslash = "/", mustWork = FALSE))
    }
  }
  normalizePath(file.path(.sc_project_root(), "data_sc", "packs"),
                winslash = "/", mustWork = FALSE)
}

# ── Helpers de coercion honnête (fromJSON simplifyVector = FALSE) ────────────

# Champ scalaire texte : character(1) non-NA, trimmé ; NULL sinon.
# (fromJSON(simplifyVector=FALSE) rend les scalaires JSON tels quels :
# "1" -> "1", 1 -> 1L — un champ déclaré texte qui arrive non-texte est un
# refus invalid_field, jamais une coercion silencieuse.)
.packs_scalar_str <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x)) return(NULL)
  trimws(x)
}

# Vecteur de chaînes depuis un tableau JSON : character vector trimmé ;
# NULL si la valeur n'est pas un tableau JSON de chaînes (objet JSON,
# scalaires non-texte…) — l'appelant traduit en invalid_field.
.packs_str_vec <- function(x) {
  if (is.null(x)) return(character(0))
  if (is.character(x)) {
    v <- as.character(x)
  } else if (is.list(x)) {
    # un objet JSON (liste NOMMÉE) n'est pas un tableau de chaînes
    if (!is.null(names(x)) && any(nzchar(names(x)))) return(NULL)
    ok <- vapply(x, function(e) is.character(e) && length(e) == 1L &&
                   !is.na(e), logical(1))
    if (!all(ok)) return(NULL)
    v <- as.character(unlist(x, use.names = FALSE))
  } else {
    return(NULL)
  }
  if (any(is.na(v))) return(NULL)
  v <- trimws(v)
  if (any(!nzchar(v))) return(NULL)
  v
}

# Champ optionnel -> NA_character_ quand absent (schéma homogène consommateurs).
.packs_opt_field <- function(x) {
  v <- .packs_scalar_str(x)
  if (is.null(v)) NA_character_ else v
}

# ── Lecteur PUR d'UN manifest.json (fiche §3/§6 — cœur de PK-1) ──────────────

# Ordre de validation (déterministe, documenté) :
#   argument -> existence du fichier -> parse JSON -> forme d'objet ->
#   type -> id -> name -> version -> attribution -> optionnels (source,
#   downloaded) -> license (requise pour palette, fiche §9) -> files
#   (requis non vide pour libraries ; existence vérifiée dès déclaré) ->
#   palettes (requises pour palette, hex validés).
pack_read_manifest <- function(manifest_path) {
  if (!is.character(manifest_path) || length(manifest_path) != 1L ||
      is.na(manifest_path) || !nzchar(trimws(manifest_path))) {
    stop(.packs_stop(
      "Chemin du manifest invalide : attendu le chemin d'UN fichier manifest.json de pack (chaîne unique non vide).",
      state = "invalid_input"))
  }
  if (dir.exists(manifest_path)) {
    stop(.packs_stop(
      sprintf("Chemin du manifest invalide : « %s » est un répertoire — attendu le chemin du FICHIER manifest.json du pack.",
              manifest_path),
      state = "invalid_input"))
  }
  if (!file.exists(manifest_path)) {
    stop(.packs_stop(
      sprintf("Manifest de pack introuvable : « %s » — un pack est un dossier contenant un manifest.json (fiche packs : data_sc/packs/<pack>/manifest.json).",
              manifest_path),
      state = "missing_manifest"))
  }

  # Parse : jsonlite est DÉJÀ au lock (D6 — zéro nouvelle dépendance). Les
  # structures restent brutes (simplifyVector = FALSE) pour que la validation
  # juge les types DÉCLARÉS, pas une coercion silencieuse.
  raw_txt <- paste(readLines(manifest_path, warn = FALSE, encoding = "UTF-8"),
                   collapse = "\n")
  obj <- tryCatch(
    jsonlite::fromJSON(raw_txt, simplifyVector = FALSE),
    error = function(e) e
  )
  if (inherits(obj, "error")) {
    stop(.packs_stop(
      sprintf("Manifest de pack illisible (JSON malformé) : « %s » — %s.",
              manifest_path, conditionMessage(obj)),
      state = "malformed_json"))
  }
  if (!is.list(obj) ||
      (length(obj) > 0L && (is.null(names(obj)) ||
        any(!nzchar(names(obj)))))) {
    stop(.packs_stop(
      sprintf("Manifest de pack malformé : « %s » — un objet JSON { … } est attendu (pas un tableau ni un scalaire).",
              manifest_path),
      state = "malformed_json"))
  }

  pack_dir <- normalizePath(dirname(manifest_path), winslash = "/",
                            mustWork = FALSE)

  # -- type (fiche idéation §3 : "libraries" | "palette" | <futur>) ----------
  if (is.null(obj$type)) {
    stop(.packs_stop(
      sprintf("Champ requis « type » absent du manifest « %s » — attendu libraries ou palette (fiche packs).",
              manifest_path),
      state = "missing_field"))
  }
  type_v <- .packs_scalar_str(obj$type)
  if (is.null(type_v)) {
    stop(.packs_stop(
      sprintf("Champ « type » malformé dans « %s » : chaîne unique attendue, reçu %s.",
              manifest_path, paste(deparse(obj$type), collapse = "")),
      state = "invalid_field"))
  }
  if (!type_v %in% packs_types()) {
    stop(.packs_stop(
      sprintf("Type de pack inconnu : « %s » (manifest « %s ») — types pris en charge : %s. Un type futur attend son implémentation dans le code applicatif.",
              type_v, manifest_path, paste(packs_types(), collapse = ", ")),
      state = "unknown_type"))
  }

  # -- id (clé du registre session ; garde sc_multi_check_label) -------------
  if (is.null(obj$id)) {
    stop(.packs_stop(
      sprintf("Champ requis « id » absent du manifest « %s » — l'id est la clé du registre des packs.",
              manifest_path),
      state = "missing_field"))
  }
  id_v <- .packs_scalar_str(obj$id)
  if (is.null(id_v) || !nzchar(id_v)) {
    stop(.packs_stop(
      sprintf("Champ « id » malformé dans « %s » : chaîne unique non vide attendue.",
              manifest_path),
      state = "invalid_field"))
  }
  if (nchar(id_v) > 80L) {
    stop(.packs_id_too_long(id_v, manifest_path))
  }
  if (grepl("[[:cntrl:]]", id_v)) {
    stop(.packs_stop(
      "Champ « id » invalide : caractères de contrôle (saut de ligne, tabulation…) interdits.",
      state = "invalid_field"))
  }

  # -- name / version / attribution (champs requis du schéma §3) -------------
  name_v <- .packs_require_text(obj$name, "name", manifest_path)
  ver_v <- .packs_require_text(obj$version, "version", manifest_path)
  attr_v <- .packs_require_text(obj$attribution, "attribution", manifest_path)

  # -- champs optionnels (présents => malformés refusés, pas ignorés) --------
  src_v <- .packs_opt_field(obj$source)
  if (!is.null(obj$source) && is.na(src_v)) {
    stop(.packs_stop(
      sprintf("Champ « source » malformé dans « %s » : chaîne unique non vide attendue (URL ou référence de provenance).",
              manifest_path),
      state = "invalid_field"))
  }
  dld_v <- NA_character_
  if (!is.null(obj$downloaded)) {
    dld_v <- .packs_scalar_str(obj$downloaded)
    if (is.null(dld_v) || !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", dld_v) ||
        is.na(as.Date(dld_v))) {
      stop(.packs_stop(
        sprintf("Champ « downloaded » malformé dans « %s » : date ISO AAAA-MM-JJ attendue (ex. 2026-10-31), reçu %s.",
                manifest_path, paste(deparse(obj$downloaded), collapse = "")),
        state = "invalid_field"))
    }
  }

  # -- license : requise pour « palette » (fiche §9 — ne pas embarquer des
  #    hex sans vérifier), optionnelle ailleurs ------------------------------
  lic_v <- NA_character_
  if (!is.null(obj$license)) {
    lic_v <- .packs_require_text(obj$license, "license", manifest_path)
  } else if (identical(type_v, "palette")) {
    stop(.packs_stop(
      sprintf("Champ requis « license » absent du manifest « %s » — un pack de palettes embarque des hex redistribués : attribution + license obligatoires (fiche packs §9).",
              manifest_path),
      state = "missing_field"))
  }

  # -- files (GMT) : requis non vide pour « libraries » (fiche §3) ;
  #    dès qu'ils sont déclarés, l'existence est vérifiée (refus honnête
  #    « fichier GMT absent », fiche PK-1 §6) ---------------------------------
  files_v <- character(0)
  if (!is.null(obj$files)) {
    files_v <- .packs_str_vec(obj$files)
    if (is.null(files_v)) {
      stop(.packs_stop(
        sprintf("Champ « files » malformé dans « %s » : tableau JSON de noms de fichiers (chaînes non vides) attendu.",
                manifest_path),
        state = "invalid_field"))
    }
    if (identical(type_v, "libraries") && length(files_v) == 0L) {
      stop(.packs_stop(
        sprintf("Champ « files » vide dans « %s » — un pack de type libraries déclare au moins un fichier GMT de librairie.",
                manifest_path),
        state = "invalid_field"))
    }
    for (f in files_v) {
      f_path <- file.path(pack_dir, f)
      if (!file.exists(f_path)) {
        stop(.packs_stop(
          sprintf("Fichier déclaré absent du pack : « %s » introuvable dans le dossier du pack %s (manifest « %s »).",
                  f, pack_dir, manifest_path),
          state = "missing_file"))
      }
    }
  } else if (identical(type_v, "libraries")) {
    stop(.packs_stop(
      sprintf("Champ requis « files » absent du manifest « %s » — un pack de type libraries déclare ses fichiers GMT (fiche packs).",
              manifest_path),
      state = "missing_field"))
  }

  # -- palettes : requises pour « palette » (fiche §7.1, hex embarqués) ------
  palettes_v <- NULL
  if (identical(type_v, "palette")) {
    if (is.null(obj$palettes)) {
      stop(.packs_stop(
        sprintf("Champ requis « palettes » absent du manifest « %s » — un pack de type palette embarque ses hex (fiche packs §7.1).",
                manifest_path),
        state = "missing_field"))
    }
    if (!is.list(obj$palettes) || length(obj$palettes) == 0L ||
        is.null(names(obj$palettes)) || any(!nzchar(names(obj$palettes)))) {
      stop(.packs_stop(
        sprintf("Champ « palettes » malformé dans « %s » : objet JSON non vide { nom_palette: [hex, …] } attendu.",
                manifest_path),
        state = "invalid_field"))
    }
    palettes_v <- lapply(names(obj$palettes), function(pnm) {
      pal <- .packs_str_vec(obj$palettes[[pnm]])
      if (is.null(pal) || length(pal) == 0L) {
        stop(.packs_stop(
          sprintf("Palette « %s » malformée dans « %s » : tableau JSON non vide de chaînes hex attendu.",
                  pnm, manifest_path),
          state = "invalid_field"))
      }
      bad <- !grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", pal)
      if (any(bad)) {
        stop(.packs_stop(
          sprintf("Palette « %s » : valeur(s) hex invalide(s) dans « %s » — attendu #RRGGBB (ou #RRGGBBAA), reçu « %s ».",
                  pnm, manifest_path,
                  paste(pal[bad], collapse = ", ")),
          state = "invalid_field"))
      }
      pal
    })
    names(palettes_v) <- names(obj$palettes)
  }

  entry <- list(
    id          = id_v,
    type        = type_v,
    name        = name_v,
    version     = ver_v,
    attribution = attr_v,
    source      = src_v,
    downloaded  = dld_v,
    license     = lic_v,
    files       = files_v,
    palettes    = palettes_v,
    dir         = pack_dir,
    manifest    = normalizePath(manifest_path, winslash = "/",
                                mustWork = FALSE)
  )
  class(entry) <- "biomentat_pack"
  entry
}

# Refus id trop long — helper séparé pour garder pack_read_manifest lisible.
.packs_id_too_long <- function(id_v, manifest_path) {
  .packs_stop(
    sprintf("Champ « id » trop long (%d caractères) dans « %s » — 80 caractères maximum : « %s… ».",
            nchar(id_v), manifest_path, substr(id_v, 1L, 40L)),
    state = "invalid_field")
}

# Champ texte requis : absent -> missing_field nommant le champ ; présent non
# textuel/vide -> invalid_field attendu-vs-reçu.
.packs_require_text <- function(x, field, manifest_path) {
  if (is.null(x)) {
    stop(.packs_stop(
      sprintf("Champ requis « %s » absent du manifest « %s ».",
              field, manifest_path),
      state = "missing_field"))
  }
  v <- .packs_scalar_str(x)
  if (is.null(v) || !nzchar(v)) {
    stop(.packs_stop(
      sprintf("Champ « %s » malformé dans « %s » : chaîne unique non vide attendue, reçu %s.",
              field, manifest_path, paste(deparse(x), collapse = "")),
      state = "invalid_field"))
  }
  v
}

# ── Scan passif d'UNE profondeur (fiche §8 — motif ts_dataset_dir_scan) ──────

#' Scan the packs root directory (fiche packs §3/§8 — passive, read-only).
#'
#' Lists ONE level of subdirectories of `packs_dir` and reads each
#' `<subdir>/manifest.json` via `pack_read_manifest()` — the scan adds ZERO
#' validation rule. Result structure:
#'   * `packs`     : named list of valid packs (key = pack `id`), each stamped
#'     `registered_at` (moment du scan = enregistrement de session) ;
#'   * `refusals`  : list of isolated refusals — chaque entrée porte `dir`,
#'     `manifest`, `state` (packs_error_states()) et `message` (le POURQUOI) ;
#'     un manifest invalide = pack SKIPPÉ avec message, JAMAIS un crash
#'     (fiche idéation §3/§9) ; un id déjà pris par un pack précédent du scan
#'     = refus duplicate_id (un registre ne porte pas deux packs de même clé) ;
#'   * `dir`       : la racine réellement scannée ;
#'   * `scanned_at`: horodatage du scan.
#'
#' STRICT READ-ONLY (motif AC2 dataset_directory.R) : rien n'est écrit nulle
#' part — une racine absente ou non-répertoire donne un résultat VIDE, jamais
#' une erreur, et la racine n'est jamais créée. L'ordre des sous-dossiers est
#' canonique (LC_COLLATE="C" pendant le listage, motif tools/sources_manifest.R).
#'
#' @param packs_root Root to scan (default: `packs_dir()`). NOTE : le paramètre
#'   NE s'appelle PAS packs_dir — il masquerait la fonction packs_dir() dans
#'   sa propre expression par défaut (référence récursive de promesse, mesuré
#'   au live-run PK-1) ; motif ts_dataset_dir_scan(datasets_dir=…).
#' @return list(packs =, refusals =, dir =, scanned_at =).
pack_scan_directory <- function(packs_root = packs_dir()) {
  root <- .sc_scalar_text(packs_root)
  if (is.null(root) || !dir.exists(root)) {
    # racine absente : AUCUN pack installé — résultat vide, jamais une erreur,
    # et la racine n'est PAS créée (lecture seule)
    return(list(
      packs = list(),
      refusals = list(),
      dir = normalizePath(if (is.null(root)) packs_root else root,
                          winslash = "/", mustWork = FALSE),
      scanned_at = Sys.time()
    ))
  }
  # ordre canonique indépendant de la locale (motif tools/sources_manifest.R)
  old_coll <- Sys.setlocale("LC_COLLATE", "C")
  on.exit(Sys.setlocale("LC_COLLATE", old_coll), add = TRUE)
  subdirs <- sort(list.dirs(root, full.names = TRUE, recursive = FALSE))

  packs <- list()
  refusals <- list()
  now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  for (d in subdirs) {
    d_norm <- normalizePath(d, winslash = "/", mustWork = FALSE)
    man_path <- file.path(d, "manifest.json")
    if (!file.exists(man_path)) {
      refusals[[length(refusals) + 1L]] <- list(
        dir = d_norm,
        manifest = normalizePath(man_path, winslash = "/", mustWork = FALSE),
        state = "missing_manifest",
        message = sprintf("Aucun manifest.json dans le dossier de pack « %s » — un pack est un dossier contenant un manifest.json (fiche packs : data_sc/packs/<pack>/manifest.json).",
                          d_norm)
      )
      next
    }
    pack <- tryCatch(pack_read_manifest(man_path), error = function(e) e)
    if (inherits(pack, "error")) {
      refusals[[length(refusals) + 1L]] <- list(
        dir = d_norm,
        manifest = normalizePath(man_path, winslash = "/", mustWork = FALSE),
        state = if (!is.null(pack$state) && pack$state %in% packs_error_states())
          pack$state else "scan_failed",
        message = conditionMessage(pack)
      )
      next
    }
    if (pack$id %in% names(packs)) {
      refusals[[length(refusals) + 1L]] <- list(
        dir = d_norm,
        manifest = normalizePath(man_path, winslash = "/", mustWork = FALSE),
        state = "duplicate_id",
        message = sprintf("Id de pack déjà pris : « %s » (dossier %s) est déjà porté par un pack précédent du scan — un registre ne peut pas porter deux packs de même clé ; renommez l'id dans le manifest.",
                          pack$id, d_norm)
      )
      next
    }
    pack$registered_at <- now
    packs[[pack$id]] <- pack
  }
  list(
    packs = packs,
    refusals = refusals,
    dir = normalizePath(root, winslash = "/", mustWork = FALSE),
    scanned_at = Sys.time()
  )
}

# ── Registre session PUR (motif sc_multi.R — JAMAIS de mutation de global) ───

# Conteneur registre : NULL ou liste nommée (motif .sc_multi_check_datasets).
.packs_check_registry <- function(registry) {
  if (is.null(registry)) return(invisible(TRUE))
  if (!is.list(registry) ||
      (length(registry) > 0L && (is.null(names(registry)) ||
        any(!nzchar(names(registry)))))) {
    stop(.packs_stop(
      "Conteneur de packs invalide : liste nommée attendue (clé = pack id).",
      state = "invalid_input"))
  }
  invisible(TRUE)
}

#' Register (or update) a pack (pure — returns the new named list)
#' @param registry Current container (NULL or named list, key = pack id).
#' @param pack A pack produced by pack_read_manifest() (class biomentat_pack).
#' @param overwrite Allow replacing an existing id (registered_at kept).
#' @return The new named list (invisible) — the caller assigns it to
#'   global_data$packs themselves (jamais de mutation de global ici).
packs_register <- function(registry, pack, overwrite = FALSE) {
  .packs_check_registry(registry)
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop(.packs_stop("overwrite doit être TRUE ou FALSE.",
                     state = "invalid_input"))
  }
  if (!inherits(pack, "biomentat_pack") || is.null(pack$id) ||
      !is.character(pack$id) || length(pack$id) != 1L || !nzchar(pack$id)) {
    stop(.packs_stop(
      "Pack invalide : produisez-le via pack_read_manifest() (objet biomentat_pack portant un id).",
      state = "invalid_input"))
  }
  exists_id <- pack$id %in% names(registry)
  if (exists_id && !isTRUE(overwrite)) {
    stop(.packs_stop(
      sprintf("Un pack d'id « %s » est déjà enregistré — renommez l'id dans son manifest ou autorisez explicitement le remplacement (overwrite = TRUE).",
              pack$id),
      state = "duplicate_id"))
  }
  now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  entry <- pack
  entry$registered_at <- if (exists_id && !is.null(registry[[pack$id]]$registered_at))
    registry[[pack$id]]$registered_at else now
  registry[[pack$id]] <- entry
  invisible(registry)
}

#' Remove a pack by id (pure — returns the new named list)
#' @param registry Current container.
#' @param id Pack id to remove.
#' @return The new named list (invisible).
packs_remove <- function(registry, id) {
  .packs_check_registry(registry)
  id_v <- .packs_scalar_str(id)
  if (is.null(id_v) || !nzchar(id_v)) {
    stop(.packs_stop("Id de pack invalide : chaîne unique non vide attendue.",
                     state = "invalid_input"))
  }
  if (!id_v %in% names(registry)) {
    stop(.packs_stop(
      sprintf("Pack « %s » introuvable dans le registre des packs.", id_v),
      state = "unknown_id"))
  }
  registry[[id_v]] <- NULL
  invisible(registry)
}

#' Get one pack by id (pure)
#' @param registry Current container.
#' @param id Pack id to fetch.
#' @return The pack entry (invisible).
packs_get <- function(registry, id) {
  .packs_check_registry(registry)
  id_v <- .packs_scalar_str(id)
  if (is.null(id_v) || !nzchar(id_v)) {
    stop(.packs_stop("Id de pack invalide : chaîne unique non vide attendue.",
                     state = "invalid_input"))
  }
  if (!id_v %in% names(registry)) {
    stop(.packs_stop(
      sprintf("Pack « %s » introuvable dans le registre des packs.", id_v),
      state = "unknown_id"))
  }
  invisible(registry[[id_v]])
}

#' Registry summary (frozen columns — consommateurs UI des jalons PK-2+)
#' @param registry Current container.
#' @return data.frame (0 rows when empty).
packs_summary <- function(registry) {
  .packs_check_registry(registry)
  cols <- c("id", "type", "name", "version", "n_files", "n_palettes",
            "attribution", "registered_at")
  if (is.null(registry) || length(registry) == 0L) {
    return(data.frame(
      id = character(0), type = character(0), name = character(0),
      version = character(0), n_files = integer(0), n_palettes = integer(0),
      attribution = character(0), registered_at = character(0),
      stringsAsFactors = FALSE
    ))
  }
  rows <- lapply(names(registry), function(nm) {
    e <- registry[[nm]]
    data.frame(
      id = nm,
      type = if (is.null(e$type)) NA_character_ else as.character(e$type),
      name = if (is.null(e$name)) NA_character_ else as.character(e$name),
      version = if (is.null(e$version)) NA_character_ else as.character(e$version),
      n_files = as.integer(length(e$files)),
      n_palettes = as.integer(length(e$palettes)),
      attribution = if (is.null(e$attribution)) NA_character_ else as.character(e$attribution),
      registered_at = if (is.null(e$registered_at)) NA_character_ else as.character(e$registered_at),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
