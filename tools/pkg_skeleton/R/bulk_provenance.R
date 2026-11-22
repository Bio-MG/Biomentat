# =============================================================================
# tools/pkg_skeleton/R/bulk_provenance.R — Manifeste de provenance Bulk RNA-seq (roadmap Bulk
# V2, Milestone 1) — PREREQUIS des analyses GSVA/WGCNA/survie.
# =============================================================================
# PRINCIPE (AGENTS.md regle 6) : la provenance est PRODUITE (a l'import, puis a
# chaque declaration utilisateur) et CONSOLIDEE pour l'affichage/l'export —
# jamais reconstruite apres coup. Ce fichier complete R/core/provenance.R
# (entrees d'analyse par etape, reutilisees telles quelles par les modules) :
# il couvre le manifeste de JEU DE DONNEES `bulk_obj$provenance` :
#   - genome de reference (GRCh38, GRCm39, ...),
#   - release d'annotation,
#   - type de normalisation (counts bruts / VST / TMM log2-CPM ...),
#   - session R (version, plateforme, packages cles).
#
# Fonctions pures, testables hors Shiny. Aucune dependance dure : les versions
# de packages sont relevees uniquement pour les packages installes.
# Erreurs : classed `bulk_provenance_error` (messages francais, call. = FALSE).
# =============================================================================

#' Surface publique figee du domaine provenance bulk (gel par test de freeze)
bulk_provenance_public_api <- function() {
  c("bulk_provenance_public_api",
    "bulk_provenance_known_genomes",
    "bulk_provenance_known_normalizations",
    "bulk_build_provenance",
    "bulk_ensure_provenance",
    "bulk_update_provenance",
    "bulk_provenance_dataframe",
    "bulk_provenance_session_packages")
}

#' Genomes de reference proposes dans l'UI (liste ouverte — "Autre" libre)
bulk_provenance_known_genomes <- function() {
  c("GRCh38", "GRCh37", "GRCm39", "GRCm38", "GRCz11", "Rnor_6.0", "Sus_scrofa", "Autre / inconnu")
}

#' Types de normalisation proposes dans l'UI
bulk_provenance_known_normalizations <- function() {
  c("Counts bruts", "VST (DESeq2)", "rlog (DESeq2)", "TMM log2-CPM (edgeR)",
    "log2-CPM", "TPM", "FPKM", "Autre / inconnu")
}

#' Versions des packages cles du pipeline bulk (installes uniquement)
#'
#' Releve nominal (jamais une re-derivation) : seule la version constatee est
#' enregistree. Les packages absents sont omis (la cle n'existe pas).
#'
#' Le chargement des namespaces est un simple CONSTAT de version : la fonction
#' ne doit laisser AUCUN etat global derriere elle (voir le correctif
#' `Matrix.warnDeprecatedCoerce` ci-dessous).
#'
#' @return Vecteur character nomme, au minimum vide.
bulk_provenance_session_packages <- function() {
  candidates <- c("DESeq2", "edgeR", "limma", "GSVA", "GSEABase", "WGCNA",
                  "variancePartition", "decoupleR", "msigdbr", "survival",
                  "survminer", "ComplexHeatmap", "mirai")
  # Effet de bord tiers neutralise : GSVA::.onLoad() fait
  # options(Matrix.warnDeprecatedCoerce = 2) SANS le restaurer. Charger GSVA
  # ici (uniquement pour LIRE sa version) escaladait donc toute depreciation
  # Matrix ulterieure en ERREUR FATALE (Matrix::Matrix.DeprecatedCoerce force
  # alors options(warn = 2L)). Casse reelle mesuree : Milo
  # (miloR::calcNhoodDistance -> as(<dgTMatrix>, "dgCMatrix")) echouait parce
  # que ce chargement laissait l'option derriere lui, globalement, pour tout le
  # reste de la session. On releve l'option avant le chargement et on la
  # restaure a l'identique (absente -> retiree, valeur NULL).
  mwd_before <- getOption("Matrix.warnDeprecatedCoerce")
  on.exit(options(Matrix.warnDeprecatedCoerce = mwd_before), add = TRUE)
  out <- character(0)
  for (p in candidates) {
    if (requireNamespace(p, quietly = TRUE)) {
      out[[p]] <- tryCatch(as.character(utils::packageVersion(p)),
                           error = function(e) NA_character_)
    }
  }
  out[!is.na(out)]
}

#' Construire le manifeste de provenance d'un jeu bulk (PRODUIT a l'import)
#'
#' @param reference_genome Character ou NA — ex. "GRCh38". NULL => NA (non declare).
#' @param annotation_release Character ou NA — ex. "Ensembl 110", "GENCODE v44".
#' @param normalization Character ou NA — ex. "Counts bruts", "VST (DESeq2)".
#' @param project Character ou NA — nom de projet (optionnel).
#' @param counts Matrice counts (genes x échantillons) ou NULL — enregistre
#'   uniquement les dimensions, jamais la matrice.
#' @param metadata data.frame métadonnées ou NULL — enregistre n et noms de colonnes.
#' @return Liste plate `type = "bulk_provenance"`.
bulk_build_provenance <- function(reference_genome = NULL, annotation_release = NULL,
                                  normalization = NULL, project = NULL,
                                  counts = NULL, metadata = NULL) {
  .as_decl <- function(x) {
    if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(as.character(x)))) {
      return(NA_character_)
    }
    trimws(as.character(x))
  }

  dims_counts <- tryCatch(dim(counts), error = function(e) NULL)
  dims_counts <- if (is.null(dims_counts) || length(dims_counts) != 2L) c(NA_integer_, NA_integer_) else as.integer(dims_counts)

  meta_cols <- tryCatch(if (is.data.frame(metadata)) colnames(metadata) else NULL,
                        error = function(e) NULL)

  list(
    type               = "bulk_provenance",
    reference_genome   = .as_decl(reference_genome),
    annotation_release = .as_decl(annotation_release),
    normalization      = .as_decl(normalization),
    project            = .as_decl(project),
    n_genes            = dims_counts[1],
    n_samples          = dims_counts[2],
    metadata_columns   = if (is.null(meta_cols)) character(0) else meta_cols,
    r_version          = paste0(R.version$major, ".", R.version$minor),
    r_version_string   = R.version.string,
    r_platform         = R.version$platform,
    r_running          = R.version$running %||% NA_character_,
    locale             = Sys.getlocale("LC_COLLATE"),
    key_packages       = bulk_provenance_session_packages(),
    created_utc        = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    updated_utc        = NA_character_,
    history            = list()
  )
}

#' Garantir la presence du manifeste sur un bulk_obj (creation paresseuse)
#'
#' Tolere les imports anterieurs a cette fonctionnalite : si `bulk_obj$provenance`
#' est absent ou corrompu, un manifeste est PRODUIT maintenant (horodatage du
#' moment de la creation, jamais antidate). Le bulk_obj est retourne INCHANGE
#' si un manifeste valide existe deja.
#'
#' @param bulk_obj Liste bulk (champs counts/metadata optionnels ici).
#' @param project Nom de projet utilise si un manifeste doit etre cree.
#' @return Le bulk_obj (invisible), avec `$provenance` garantit valide.
bulk_ensure_provenance <- function(bulk_obj, project = NULL) {
  if (is.null(bulk_obj)) {
    stop(errorCondition(
      "bulk_ensure_provenance() : bulk_obj est NULL — importez d'abord un jeu bulk.",
      class = "bulk_provenance_error", state = "invalid_input"))
  }
  if (!is.list(bulk_obj)) {
    stop(errorCondition(sprintf(
      "bulk_ensure_provenance() : un objet bulk (liste) est requis (recu : %s).",
      paste(class(bulk_obj), collapse = "/")),
      class = "bulk_provenance_error", state = "invalid_input"))
  }
  prov <- bulk_obj$provenance
  if (is.list(prov) && identical(prov$type, "bulk_provenance")) {
    return(invisible(bulk_obj))
  }
  bulk_obj$provenance <- bulk_build_provenance(
    normalization = "Counts bruts",
    project       = project,
    counts        = if (is.matrix(bulk_obj$counts)) bulk_obj$counts else NULL,
    metadata      = bulk_obj$metadata
  )
  invisible(bulk_obj)
}

#' Mettre a jour les declarations du manifeste (PRODUCTION d'une nouvelle etat)
#'
#' Les champs declares par l'utilisateur (genome, release, normalisation) sont
#' remplaces ; l'ancienne valeur est preservee dans `history` (jamais perdue) ;
#`  `updated_utc` horodate. Le manifeste lui-meme n'est JAMAIS reconstruit :
#' seuls les champs cites et l'horodatage changent.
#'
#' @param bulk_obj Liste bulk (bulk_ensure_provenance() est appele d'abord).
#' @param reference_genome,annotation_release,normalization Nouvelles declarations
#'   (NULL = champ inchange).
#' @return Le bulk_obj (invisible).
bulk_update_provenance <- function(bulk_obj, reference_genome = NULL,
                                   annotation_release = NULL, normalization = NULL) {
  bulk_obj <- bulk_ensure_provenance(bulk_obj)
  prov     <- bulk_obj$provenance

  .as_decl <- function(x) {
    if (is.null(x) || length(x) != 1L || is.na(x) || !nzchar(trimws(as.character(x)))) {
      return(NA_character_)
    }
    trimws(as.character(x))
  }

  updates <- list(
    reference_genome   = .as_decl(reference_genome),
    annotation_release = .as_decl(annotation_release),
    normalization      = .as_decl(normalization)
  )
  changed <- FALSE
  for (nm in names(updates)) {
    if (!is.na(updates[[nm]]) && !identical(prov[[nm]], updates[[nm]])) {
      prov$history <- c(prov$history, list(list(
        field = nm, previous = prov[[nm]], new = updates[[nm]],
        at_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
      )))
      prov[[nm]] <- updates[[nm]]
      changed <- TRUE
    }
  }
  if (changed) prov$updated_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  bulk_obj$provenance <- prov
  invisible(bulk_obj)
}

#' Aplatir le manifeste pour affichage / export (CONSOLIDATION)
#'
#' @param prov Manifeste `bulk_provenance` (ou NULL — renvoie un cadre vide).
#' @return data.frame 2 colonnes (Champ, Valeur), une ligne par champ scalaire ;
#'   les packages cles sont compactes en une ligne "nom=version; ...".
bulk_provenance_dataframe <- function(prov) {
  if (!is.list(prov) || !identical(prov$type, "bulk_provenance")) {
    return(data.frame(Champ = character(0), Valeur = character(0),
                      stringsAsFactors = FALSE))
  }
  pkg_str <- if (length(prov$key_packages) == 0L) "" else
    paste(vapply(names(prov$key_packages),
                 function(k) sprintf("%s=%s", k, prov$key_packages[[k]]),
                 character(1)), collapse = "; ")

  rows <- list(
    "Projet"                 = prov$project,
    "Génome de référence"    = prov$reference_genome,
    "Release d'annotation"   = prov$annotation_release,
    "Normalisation"          = prov$normalization,
    "Gènes × échantillons"   = if (!is.na(prov$n_genes) && !is.na(prov$n_samples))
                                  paste0(prov$n_genes, " × ", prov$n_samples) else NA_character_,
    "Version R"              = prov$r_version_string,
    "Plateforme"             = prov$r_platform,
    "Packages clés"          = pkg_str,
    "Créé (UTC)"             = prov$created_utc,
    "Mis à jour (UTC)"       = prov$updated_utc
  )
  data.frame(
    Champ  = names(rows),
    Valeur = vapply(rows, function(v) if (is.null(v) || is.na(v)) "" else as.character(v), character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Calculer l'empreinte SHA-256 complete d'un fichier de jeu de donnees (BK-08)
#'
#' Verifie que `path` designe un fichier regulier existant (non repertoire),
#' puis calcule son condensat cryptographique SHA-256 complet via `digest`.
#'
#' @param path Chemin unique vers un fichier existant.
#' @return Liste nommee `path`, `size_bytes`, `sha256`, `algo`.
#' @export
ts_dataset_sha256 <- function(path) {
  if (is.null(path) || !is.character(path) || length(path) != 1L ||
      is.na(path) || !nzchar(trimws(path))) {
    stop(errorCondition(
      "ts_dataset_sha256() : 'path' doit etre une chaine de caracteres non vide de longueur 1.",
      class = "biomentat_provenance_error"
    ))
  }
  if (!file.exists(path) || dir.exists(path) || isTRUE(file.info(path)$isdir)) {
    stop(errorCondition(
      sprintf("ts_dataset_sha256() : fichier introuvable ou repertoire invalide (%s).", path),
      class = "biomentat_provenance_error"
    ))
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop(errorCondition(
      "ts_dataset_sha256() : le package 'digest' est requis pour le calcul SHA-256.",
      class = "biomentat_provenance_error"
    ))
  }

  sha <- digest::digest(file = path, algo = "sha256")
  list(
    path       = normalizePath(path, winslash = "/", mustWork = FALSE),
    size_bytes = as.numeric(file.info(path)$size),
    sha256     = sha,
    algo       = "sha256"
  )
}

#' Construire des plis de validation croisee sans fuite de groupe (BK-08)
#'
#' Assigne les echantillons en `k` plis de sorte que tous les echantillons
#' partageant le meme `group_ids` (ex. donneur/patient) appartiennent exactement
#' au meme pli (fuite de groupe nulle). Si `strata` est fourni, les groupes sont
#' repartis en round-robin stratifie par strate.
#'
#' @param sample_ids Vecteur unique d'identifiants d'echantillons.
#' @param group_ids Vecteur d'identifiants de groupe (meme longueur que `sample_ids`).
#' @param strata Vecteur optionnel de strates (meme longueur que `sample_ids`).
#' @param k Nombre de plis (>= 2, <= nombre de groupes uniques).
#' @param seed Graine aleatoire optionnelle (NULL pour ne pas fixer la graine).
#' @return Liste nommee `fold_id`, `folds`, `k`, `n_samples`, `n_groups`,
#'   `leakage_verified`.
#' @export
ts_cv_leakage_free_folds <- function(sample_ids, group_ids, strata = NULL,
                                     k = 5L, seed = 42L) {
  if (is.null(sample_ids) || length(sample_ids) == 0L || anyNA(sample_ids) ||
      anyDuplicated(sample_ids) > 0L) {
    stop(errorCondition(
      "ts_cv_leakage_free_folds() : 'sample_ids' doit etre un vecteur non vide, sans NA ni doublons.",
      class = "biomentat_provenance_error"
    ))
  }
  n_samples <- length(sample_ids)
  if (is.null(group_ids) || length(group_ids) != n_samples || anyNA(group_ids)) {
    stop(errorCondition(
      "ts_cv_leakage_free_folds() : 'group_ids' doit avoir la meme longueur que 'sample_ids' et aucun NA.",
      class = "biomentat_provenance_error"
    ))
  }
  if (!is.null(strata) && (length(strata) != n_samples || anyNA(strata))) {
    stop(errorCondition(
      "ts_cv_leakage_free_folds() : 'strata' doit avoir la meme longueur que 'sample_ids' et aucun NA.",
      class = "biomentat_provenance_error"
    ))
  }
  if (is.null(k) || !is.numeric(k) || length(k) != 1L || is.na(k) ||
      k < 2 || k != as.integer(k)) {
    stop(errorCondition(
      "ts_cv_leakage_free_folds() : 'k' doit etre un entier >= 2.",
      class = "biomentat_provenance_error"
    ))
  }
  k <- as.integer(k)

  grp_chr <- as.character(group_ids)
  uniq_groups <- unique(grp_chr)
  n_groups <- length(uniq_groups)
  if (n_groups < k) {
    stop(errorCondition(
      sprintf("ts_cv_leakage_free_folds() : k (%d) depasse le nombre de groupes uniques (%d).",
              k, n_groups),
      class = "biomentat_provenance_error"
    ))
  }

  if (!is.null(strata)) {
    str_chr <- as.character(strata)
    grp_strata <- tapply(str_chr, grp_chr, function(x) length(unique(x)))
    if (any(grp_strata > 1L)) {
      stop(errorCondition(
        "ts_cv_leakage_free_folds() : un meme groupe (group_ids) appartient a plusieurs strates.",
        class = "biomentat_provenance_error"
      ))
    }
  }

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop(errorCondition(
        "ts_cv_leakage_free_folds() : 'seed' doit etre un scalaire numerique ou NULL.",
        class = "biomentat_provenance_error"
      ))
    }
    has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (has_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  group_fold <- integer(n_groups)
  names(group_fold) <- uniq_groups

  if (is.null(strata)) {
    shuffled <- sample(uniq_groups, size = n_groups, replace = FALSE)
    group_fold[shuffled] <- ((seq_len(n_groups) - 1L) %% k) + 1L
  } else {
    str_chr <- as.character(strata)
    grp_to_stratum <- vapply(uniq_groups, function(g) str_chr[match(g, grp_chr)], character(1))
    uniq_strata <- unique(str_chr)
    offset <- 0L
    for (st in uniq_strata) {
      st_groups <- uniq_groups[grp_to_stratum == st]
      shuffled_st <- if (length(st_groups) > 1L) sample(st_groups, size = length(st_groups), replace = FALSE) else st_groups
      idx <- seq_along(shuffled_st) - 1L + offset
      group_fold[shuffled_st] <- (idx %% k) + 1L
      offset <- offset + length(shuffled_st)
    }
  }

  fold_id <- unname(as.integer(group_fold[grp_chr]))
  folds <- lapply(seq_len(k), function(f) sample_ids[fold_id == f])
  names(folds) <- paste0("fold_", seq_len(k))

  fold_groups <- lapply(seq_len(k), function(f) unique(grp_chr[fold_id == f]))
  if (sum(lengths(fold_groups)) != n_groups) {
    stop(errorCondition(
      "ts_cv_leakage_free_folds() : fuite de groupe detectee entre les plis.",
      class = "biomentat_provenance_error"
    ))
  }

  list(
    fold_id          = fold_id,
    folds            = folds,
    k                = k,
    n_samples        = n_samples,
    n_groups         = n_groups,
    leakage_verified = TRUE
  )
}
