# =============================================================================
# provenance.R — Manifeste de provenance des analyses (CHRYSALIS 2C)
# =============================================================================
# PRINCIPE (AGENTS.md, regle 7) : la provenance est PRODUITE a chaque etape
# d'analyse — new_provenance_entry() au moment ou le calcul reussit,
# provenance_append() immediatement apres — et seulement CONSOLIDEE au moment
# du rapport (provenance_to_dataframe()). Elle n'est JAMAIS reconstruite apres
# coup. Aucune logique de rendu de rapport ici : uniquement la structure de
# donnees et les appenders.
#
# D2 M5 lot 1 (2026-10-25) : ce fichier demenage R/core/ -> tools/pkg_skeleton/R/
# (deplacement physique seul, app.R le source() a la meme position d'ordre ;
# cf. cartographie d2_reports_cartography.md §1).
#
# Dependances :
#   - state.R (state_get/state_set) pour l'appending dans l'etat
#     partage ; chargé paresseusement ci-dessous si absent (convention des
#     shims bulk_helpers.R / sc_state.R). Le socle vit dans le meme arbre
#     paquet depuis ce lot.
#   - 'digest' (optionnel, requireNamespace) pour l'empreinte du jeu de
#     donnees ; sans lui, dataset_hash reste NA.
#   - Aucune dependance directe Seurat/DESeq2 : les versions logicielles sont
#     collectees dynamiquement pour les packages installes uniquement.
# =============================================================================

if (!exists("state_get", envir = globalenv(), mode = "function") ||
    !exists("state_set", envir = globalenv(), mode = "function")) {
  # M5 lot 1 : state.R vit dans l'arbre du paquet — le shim suit le chemin.
  for (.p in c("tools/pkg_skeleton/R/state.R",
               file.path(getwd(), "tools", "pkg_skeleton", "R", "state.R"))) {
    if (file.exists(.p)) { sys.source(.p, envir = globalenv()); break }
  }
}

# Packages dont la version est enregistree quand ils sont installes (les
# cles scientifiques du roadmap ; la liste peut s'etendre sans breaking change
# puisque seuls les packages presents figurent dans le manifeste).
.provenance_pkg_candidates <- c(
  "Seurat", "SeuratObject", "DESeq2", "edgeR", "limma",
  "clusterProfiler", "slingshot", "scran", "BPCells", "mirai", "digest",
  "liana", "OmnipathR"
)

#' Collecter les versions logicielles (R + packages installes parmi les cles)
#'
#' @return Vecteur character nomme, toujours au moins l'entree "R".
.provenance_versions <- function() {
  v <- c(R = paste0(R.version$major, ".", R.version$minor))
  for (pkg in .provenance_pkg_candidates) {
    path <- tryCatch(find.package(pkg, quiet = TRUE), error = function(e) character(0))
    if (length(path) == 1L && nzchar(path)) {
      version <- tryCatch(
        as.character(utils::packageVersion(pkg)),
        error = function(e) NA_character_
      )
      if (!is.na(version)) v[pkg] <- version
    }
  }
  v
}

#' Identite d'un jeu de donnees : nom, empreinte, dimensions
#'
#' - NULL                      -> tout NA (etape sans dataset identifiable)
#' - chaine de caracteres      -> nom direct (ex. accession GEO), pas d'empreinte
#' - objet a `dim()` (matrice, data.frame, Seurat) -> dimensions + empreinte :
#'     * objet < ~10 Mo  : empreinte EXACTE (digest de l'objet serialise)
#'     * objet volumineux: empreinte SUPERFICIELLE (dims + classe + taille),
#'       hash_exact = FALSE — on ne serialise jamais des gigaoctets juste
#'       pour la provenance.
#'
#' @param dataset Objet ou nom de jeu de donnees.
#' @return list(name=, hash=, hash_exact=, dims=) — dims toujours longueur 2.
.provenance_dataset_identity <- function(dataset) {
  na_dims <- c(NA_integer_, NA_integer_)
  if (is.null(dataset)) {
    return(list(name = NA_character_, hash = NA_character_, hash_exact = NA,
                dims = na_dims))
  }
  if (is.character(dataset) && length(dataset) == 1L) {
    return(list(name = dataset, hash = NA_character_, hash_exact = NA,
                dims = na_dims))
  }

  dims <- tryCatch(dim(dataset), error = function(e) NULL)
  dims <- if (is.null(dims) || length(dims) == 0L || all(is.na(dims))) na_dims else dims
  if (length(dims) == 1L) dims <- c(dims, NA_integer_)
  dims <- as.integer(dims[seq_len(2L)])

  # Nom lisible pour un objet Seurat (nom du projet) ; NA sinon.
  nm <- tryCatch(
    if (inherits(dataset, "Seurat")) as.character(dataset@project.name) else NA_character_,
    error = function(e) NA_character_
  )

  hash <- NA_character_
  hash_exact <- NA
  if (requireNamespace("digest", quietly = TRUE)) {
    size_bytes <- tryCatch(as.numeric(utils::object.size(dataset)), error = function(e) 0)
    if (is.finite(size_bytes) && size_bytes > 0 && size_bytes < 1e7) {
      hash <- digest::digest(dataset, algo = "xxhash64")
      hash_exact <- TRUE
    } else {
      hash <- digest::digest(
        list(dims = dims, classe = class(dataset)[1], taille_Mo = round(size_bytes / 1e6, 3)),
        algo = "xxhash64"
      )
      hash_exact <- FALSE
    }
  }

  list(name = nm, hash = hash, hash_exact = hash_exact, dims = dims)
}

#' Construire une entree de provenance (PRODUITE a l'etape d'analyse)
#'
#' @param analysis_id Identifiant unique de l'etape (ex. "sc-trajectory-1").
#' @param method Methode appliquee (ex. "slingshot", "deseq2").
#' @param parameters Liste nommee des parametres REELLEMENT utilises — ne pas
#'   inventer de parametres, enregistrer uniquement ce qui existait deja.
#' @param dataset Objet (matrice/data.frame/Seurat) ou nom du jeu de donnees.
#' @param cells_used,cells_excluded Cellules utilisees/exclues (scalaire
#'   numerique ou vecteur d'identifiants ; NULL si non applicable).
#' @param seed Graine (NULL si non fixee).
#' @param warnings Vecteur character (vide si aucune).
#' @param timestamp POSIXct, horodatage de PRODUCTION.
#' @return Liste plate documentee.
new_provenance_entry <- function(analysis_id, method, parameters = list(),
                                 dataset = NULL, cells_used = NULL,
                                 cells_excluded = NULL, seed = NULL,
                                 warnings = character(0),
                                 timestamp = Sys.time()) {
  analysis_id <- as.character(analysis_id)
  if (length(analysis_id) != 1L || is.na(analysis_id) || !nzchar(analysis_id)) {
    stop("new_provenance_entry() : 'analysis_id' requis (chaine non vide).", call. = FALSE)
  }
  method <- as.character(method)
  if (length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("new_provenance_entry() : 'method' requise (chaine non vide).", call. = FALSE)
  }

  ds <- .provenance_dataset_identity(dataset)

  list(
    analysis_id     = analysis_id,
    method          = method,
    timestamp       = timestamp,
    dataset_name    = ds$name,
    dataset_hash    = ds$hash,
    hash_exact      = ds$hash_exact,
    dataset_dims    = ds$dims,
    parameters      = as.list(parameters),
    versions        = .provenance_versions(),
    cells_used      = cells_used,
    cells_excluded  = cells_excluded,
    seed            = seed,
    warnings        = as.character(warnings)
  )
}

#' Ajouter une entree au champ `provenance` de l'etat (PRODUITE, jamais reconstruite)
#'
#' @param state Etat partage (reactiveValues/environment, cf. R/core/state.R).
#'   NULL est tolere et ignore silencieusement : un module sans etat partage
#'   ne doit pas planter parce qu'il produit sa provenance.
#' @param entry Entree issue de new_provenance_entry().
#' @return invisible(index de l'entree), invisible(NULL) si state NULL.
provenance_append <- function(state, entry) {
  if (is.null(state)) return(invisible(NULL))
  if (!is.list(entry) || is.null(entry$analysis_id) || is.null(entry$method)) {
    stop("provenance_append() : 'entry' doit provenir de new_provenance_entry().", call. = FALSE)
  }
  cur <- state_get(state, "provenance")
  if (is.null(cur)) cur <- list()
  # NB : un data.frame EST un is.list() — refuse explicitement, sinon
  # l'appending [[<- exploserait avec un message obscur.
  if (!is.list(cur) || is.data.frame(cur)) {
    stop("provenance_append() : le champ 'provenance' de l'etat n'est pas une liste — etat corrompu.", call. = FALSE)
  }
  entry$index <- length(cur) + 1L
  cur[[entry$index]] <- entry
  state_set(state, "provenance", cur)
  invisible(entry$index)
}

#' Aplatir toutes les entrees de provenance pour consommation par le rapport
#'
#' CONSOLIDATION uniquement (regle 7) : cette fonction lit ce qui a ete
#' PRODUIT, ne deduit rien, ne complete rien.
#'
#' @param state Etat partage (ou NULL).
#' @return data.frame une ligne par entree, colonnes stables (vide mais
#'   correctement typee si aucune entree).
provenance_to_dataframe <- function(state) {
  entries <- if (is.null(state)) NULL else state_get(state, "provenance")

  .fmt_vec  <- function(v) paste(v, collapse = "; ")
  .params_str <- function(p) {
    if (length(p) == 0L) return("")
    .fmt_vec(vapply(names(p), function(k) {
      sprintf("%s=%s", k, paste(format(p[[k]]), collapse = ","))
    }, character(1)))
  }
  .versions_str <- function(v) {
    if (length(v) == 0L) return("")
    .fmt_vec(vapply(names(v), function(k) sprintf("%s=%s", k, v[[k]]), character(1)))
  }
  .count_str <- function(x) {
    if (is.null(x)) ""
    else if (length(x) == 1L) format(x)
    else as.character(length(x))
  }
  .ts_str <- function(ts) {
    tryCatch(format(ts, "%Y-%m-%d %H:%M:%S"), error = function(e) as.character(ts))
  }

  cols <- c("analysis_id", "method", "timestamp", "dataset_name", "dataset_hash",
            "hash_exact", "dataset_n_rows", "dataset_n_cols", "parameters",
            "versions", "cells_used", "cells_excluded", "seed", "warnings")
  empty_df <- as.data.frame(
    setNames(rep(list(character(0)), length(cols)), cols),
    stringsAsFactors = FALSE
  )
  if (is.null(entries) || length(entries) == 0L) return(empty_df)

  rows <- lapply(entries, function(e) {
    data.frame(
      stringsAsFactors = FALSE,
      analysis_id     = as.character(e$analysis_id %||% ""),
      method          = as.character(e$method %||% ""),
      timestamp       = .ts_str(e$timestamp),
      dataset_name    = as.character(e$dataset_name %||% ""),
      dataset_hash    = as.character(e$dataset_hash %||% ""),
      hash_exact      = as.character(e$hash_exact %||% ""),
      dataset_n_rows  = as.character(e$dataset_dims[[1]] %||% ""),
      dataset_n_cols  = as.character(e$dataset_dims[[2]] %||% ""),
      parameters      = .params_str(e$parameters),
      versions        = .versions_str(e$versions),
      cells_used      = .count_str(e$cells_used),
      cells_excluded  = .count_str(e$cells_excluded),
      seed            = .count_str(e$seed),
      warnings        = .fmt_vec(e$warnings %||% character(0))
    )
  })
  do.call(rbind, rows)
}

#' #15 biomentat-methods-drafter : rédige une section Matériel & Méthodes
#' déterministe à partir de la provenance enregistrée (aucune hallucination).
#'
#' @param state_or_entries État partagé (`state_get(..., "provenance")`) ou
#'   liste directe d'entrées `new_provenance_entry()`.
#' @param lang Langue de sortie (`"fr"` ou `"en"`).
#' @param dataset_label Nom lisible du jeu de données (optionnel).
#' @return Liste `list(markdown, text, n_steps, packages_cited)`.
#' @export
build_methods_draft <- function(state_or_entries,
                                lang = c("fr", "en"),
                                dataset_label = NULL) {
  lang <- match.arg(lang)
  entries <- if (is.null(state_or_entries)) {
    list()
  } else if (is.data.frame(state_or_entries)) {
    if ("analysis_id" %in% colnames(state_or_entries) && nrow(state_or_entries) > 0L) {
      lapply(seq_len(nrow(state_or_entries)), function(i) {
        row <- state_or_entries[i, , drop = FALSE]
        list(
          analysis_id = as.character(row$analysis_id[1L]),
          method = if ("method" %in% colnames(row)) as.character(row$method[1L]) else NA_character_,
          timestamp_utc = if ("timestamp_utc" %in% colnames(row)) as.character(row$timestamp_utc[1L]) else NA_character_,
          parameters = list(),
          dataset = list(),
          warnings = character(0)
        )
      })
    } else {
      list()
    }
  } else if (inherits(state_or_entries, "reactivevalues") || is.environment(state_or_entries)) {
    state_get(state_or_entries, "provenance") %||% list()
  } else if (is.list(state_or_entries) &&
             all(vapply(state_or_entries, function(x) is.list(x) && !is.null(x$analysis_id), logical(1)))) {
    state_or_entries
  } else {
    state_get(state_or_entries, "provenance") %||% list()
  }

  if (length(entries) == 0L) {
    msg <- if (identical(lang, "en")) {
      "No analysis steps recorded in provenance."
    } else {
      "Aucune \u00e9tape d'analyse enregistr\u00e9e dans la provenance."
    }
    return(list(markdown = msg, text = msg, n_steps = 0L, packages_cited = character(0)))
  }

  df <- provenance_to_dataframe(
    if (inherits(state_or_entries, "reactivevalues") || is.environment(state_or_entries)) {
      state_or_entries
    } else {
      list(provenance = entries)
    }
  )
  pkgs <- character(0)
  for (e in entries) {
    v <- e$versions
    if (length(v) > 0L && !is.null(names(v))) {
      for (nm in names(v)) pkgs[[nm]] <- as.character(v[[nm]])
    }
  }
  pkg_str <- if (length(pkgs) > 0L) {
    paste(sprintf("%s (%s)", names(pkgs), unname(pkgs)), collapse = ", ")
  } else {
    paste0("R (", R.version$major, ".", R.version$minor, ")")
  }

  hdr <- if (identical(lang, "en")) "## Materials and Methods" else "## Mat\u00e9riel et M\u00e9thodes"
  intro <- if (identical(lang, "en")) {
    sprintf(
      "Analyses were performed with Biomentat using %s%s.",
      pkg_str,
      if (!is.null(dataset_label) && nzchar(dataset_label)) sprintf(" on dataset '%s'", dataset_label) else ""
    )
  } else {
    sprintf(
      "Les analyses ont \u00e9t\u00e9 r\u00e9alis\u00e9es avec Biomentat \u00e0 l'aide de %s%s.",
      pkg_str,
      if (!is.null(dataset_label) && nzchar(dataset_label)) sprintf(" sur le jeu \u00ab %s \u00bb", dataset_label) else ""
    )
  }

  bullets <- vapply(seq_len(nrow(df)), function(i) {
    parts <- c(
      sprintf("**%s** (`%s`)", df$analysis_id[i], df$method[i]),
      if (nzchar(df$parameters[i])) paste0(if (identical(lang, "en")) "parameters: " else "param\u00e8tres : ", df$parameters[i]),
      if (nzchar(df$cells_used[i])) paste0(if (identical(lang, "en")) "cells used: " else "cellules utilis\u00e9es : ", df$cells_used[i]),
      if (nzchar(df$seed[i])) paste0("seed: ", df$seed[i])
    )
    paste0("- ", paste(parts[!is.na(parts) & nzchar(parts)], collapse = " ; "), ".")
  }, character(1L))

  md <- paste(c(hdr, "", intro, "", bullets), collapse = "\n")
  txt <- gsub("[`*]", "", md)
  list(
    markdown = md,
    text = txt,
    n_steps = nrow(df),
    packages_cited = pkgs
  )
}
