# =============================================================================
# tools/pkg_skeleton/R/sc_communication_velocity_directional.R — W2.2 §4.3 :
# directionnalité vélocité CCC (alignement T_ij et corrélation vitesse-récepteur)
# =============================================================================
# Consommateur PUR du résultat canonique communication (Stage 11) et du
# résultat velocity canonique (Stage 10) :
#   1. Alignement des probabilités de transition cellule->cellule T_ij pour
#      chaque paire (sender_node, receiver_node) par rapport au fond hors-diagonale ;
#   2. Corrélation de Spearman intra-population (cellules receiver) entre la
#      magnitude des vecteurs de vitesse et l'expression du/des récepteur(s).
#
# GARDE-FOU ÉPISTÉMIQUE STRICT :
#   Association descriptive uniquement (AUCUNE causalité). Aucune matrice T_ij
#   ni aucun vecteur de vitesse n'est ré-inféré ou fabriqué : toute donnée
#   absente retourne un statut explicite et des valeurs NA. La table canonique
#   du résultat communication parent n'est jamais modifiée.
# =============================================================================

.vdir_epistemic_note <- function() {
  paste0(
    "Association descriptive uniquement (AUCUNE causalité) : l'alignement des ",
    "probabilités de transition T_ij et la corrélation de Spearman ",
    "vitesse-récepteur ne démontrent pas que le ligand pilote la transition ",
    "cellulaire."
  )
}

.vdir_stop <- function(message, state = "invalid_input") {
  stop(errorCondition(
    message, state = state, class = "communication_context_error"
  ))
}

.vdir_validate_identities <- function(cell_identities) {
  if (!is.character(cell_identities) || is.null(names(cell_identities)) ||
      anyNA(names(cell_identities)) || !all(nzchar(names(cell_identities))) ||
      anyNA(cell_identities) || !all(nzchar(cell_identities))) {
    .vdir_stop(paste0(
      "compute_communication_velocity_directionality() : 'cell_identities' ",
      "doit être un vecteur character nommé (cellule -> identité), sans valeur vide."
    ))
  }
  if (anyDuplicated(names(cell_identities))) {
    .vdir_stop(paste0(
      "compute_communication_velocity_directionality() : noms de cellules ",
      "dupliqués dans 'cell_identities'."
    ))
  }
  invisible(TRUE)
}

.vdir_validate_transition_matrix <- function(transition_matrix) {
  is_mat <- is.matrix(transition_matrix) || inherits(transition_matrix, "Matrix")
  if (!is_mat) {
    .vdir_stop(paste0(
      "compute_communication_velocity_directionality() : 'transition_matrix' ",
      "doit être une matrice numérique carrée (cellules x cellules)."
    ))
  }
  if (nrow(transition_matrix) != ncol(transition_matrix) || nrow(transition_matrix) < 2L) {
    .vdir_stop(paste0(
      "compute_communication_velocity_directionality() : 'transition_matrix' ",
      "doit être une matrice carrée d'au moins 2x2 cellules."
    ))
  }
  rn <- rownames(transition_matrix)
  cn <- colnames(transition_matrix)
  if (is.null(rn) || is.null(cn) || anyNA(rn) || anyNA(cn) ||
      !all(nzchar(rn)) || !all(nzchar(cn)) ||
      anyDuplicated(rn) || anyDuplicated(cn)) {
    .vdir_stop(paste0(
      "compute_communication_velocity_directionality() : 'transition_matrix' ",
      "doit porter des rownames et colnames uniques non vides."
    ))
  }
  vals <- if (inherits(transition_matrix, "Matrix")) {
    as.numeric(as.matrix(transition_matrix))
  } else {
    transition_matrix
  }
  if (!is.numeric(vals)) {
    .vdir_stop(paste0(
      "compute_communication_velocity_directionality() : 'transition_matrix' ",
      "doit contenir des valeurs numériques."
    ))
  }
  if (any(vals < 0, na.rm = TRUE)) {
    .vdir_stop(paste0(
      "compute_communication_velocity_directionality() : valeurs négatives ",
      "interdites dans 'transition_matrix'."
    ))
  }
  invisible(TRUE)
}

.vdir_validate_receptor_expr <- function(receptor_expr) {
  if (is.numeric(receptor_expr) && is.null(dim(receptor_expr))) {
    nms <- names(receptor_expr)
    if (is.null(nms) || anyNA(nms) || !all(nzchar(nms)) || anyDuplicated(nms)) {
      .vdir_stop(paste0(
        "compute_communication_velocity_directionality() : un vecteur ",
        "'receptor_expr' doit être nommé par identifiants de cellules uniques."
      ))
    }
    return(list(mode = "vector", data = receptor_expr))
  }
  is_mat <- is.matrix(receptor_expr) || inherits(receptor_expr, "Matrix")
  if (is_mat) {
    rn <- rownames(receptor_expr)
    cn <- colnames(receptor_expr)
    if (is.null(rn) || is.null(cn) || anyNA(rn) || anyNA(cn) ||
        !all(nzchar(rn)) || !all(nzchar(cn)) || anyDuplicated(cn)) {
      .vdir_stop(paste0(
        "compute_communication_velocity_directionality() : une matrice ",
        "'receptor_expr' (gènes x cellules) doit avoir rownames et colnames valides."
      ))
    }
    if (is.matrix(receptor_expr) && !is.numeric(receptor_expr)) {
      .vdir_stop(paste0(
        "compute_communication_velocity_directionality() : la matrice ",
        "'receptor_expr' doit être numérique."
      ))
    }
    return(list(mode = "matrix", data = receptor_expr))
  }
  .vdir_stop(paste0(
    "compute_communication_velocity_directionality() : 'receptor_expr' doit ",
    "être un vecteur numérique nommé ou une matrice gènes x cellules."
  ))
}

.vdir_extract_receptor_vec <- function(rec_spec, receptors, r_cells) {
  if (identical(rec_spec$mode, "vector")) {
    common <- intersect(r_cells, names(rec_spec$data))
    return(list(vec = rec_spec$data[common], matched = "(vector)"))
  }
  m <- rec_spec$data
  common_cells <- intersect(r_cells, colnames(m))
  if (length(common_cells) == 0L) {
    return(list(vec = numeric(0), matched = character(0)))
  }
  rec_clean <- unique(receptors[!is.na(receptors) & nzchar(receptors)])
  rn <- rownames(m)
  matched <- character(0)
  for (rc in rec_clean) {
    if (rc %in% rn) {
      matched <- c(matched, rc)
    } else {
      parts <- unlist(strsplit(rc, "[_+:]", perl = TRUE))
      matched <- c(matched, intersect(parts, rn))
    }
  }
  matched <- unique(matched)
  if (length(matched) == 0L) {
    return(list(vec = numeric(0), matched = character(0)))
  }
  sub_m <- as.matrix(m[matched, common_cells, drop = FALSE])
  vec <- colMeans(sub_m, na.rm = TRUE)
  names(vec) <- common_cells
  list(vec = vec, matched = matched)
}

.vdir_spearman_within_pop <- function(mag, rec_vec, r_cells) {
  common <- intersect(r_cells, intersect(names(mag), names(rec_vec)))
  if (length(common) < 3L) {
    return(list(cor = NA_real_, n = as.integer(length(common)),
                status = "insufficient_cells"))
  }
  v <- as.numeric(mag[common])
  e <- as.numeric(rec_vec[common])
  ok <- is.finite(v) & is.finite(e)
  n_ok <- sum(ok)
  if (n_ok < 3L) {
    return(list(cor = NA_real_, n = as.integer(n_ok),
                status = "insufficient_cells"))
  }
  if (stats::sd(v[ok]) == 0 || stats::sd(e[ok]) == 0) {
    return(list(cor = NA_real_, n = as.integer(n_ok),
                status = "zero_variance"))
  }
  rho <- suppressWarnings(unname(stats::cor(v[ok], e[ok], method = "spearman")))
  list(cor = rho, n = as.integer(n_ok), status = "valid")
}

.vdir_resolve_nodes <- function(tbl) {
  sm <- tbl$sender_mapped %||% rep(NA_character_, nrow(tbl))
  rm <- tbl$receiver_mapped %||% rep(NA_character_, nrow(tbl))
  tbl$sender_node <- ifelse(!is.na(sm), sm, tbl$sender)
  tbl$receiver_node <- ifelse(!is.na(rm), rm, tbl$receiver)
  tbl
}

#' Directionnalité vélocité CCC : alignement T_ij et corrélation vitesse-récepteur
#'
#' Évalue de manière strictement descriptive (AUCUNE causalité) :
#' \itemize{
#'   \item l'alignement des probabilités de transition cellule-cellule
#'     \eqn{T_{ij}} entre cellules sender et receiver par rapport à la moyenne
#'     de fond hors-diagonale (si \code{transition_matrix} est fournie) ;
#'   \item la corrélation de Spearman intra-population (cellules receiver) entre
#'     la magnitude des vecteurs de vitesse et l'expression du/des récepteur(s)
#'     (si \code{receptor_expr} et les vecteurs de vitesse sont disponibles).
#' }
#'
#' @param communication_result Résultat canonique (\code{assert_communication_result}).
#' @param velocity_result Résultat velocity canonique (\code{assert_velocity_result},
#'   \code{view = "any"}).
#' @param cell_identities Vecteur character nommé (cellule -> identité).
#' @param transition_matrix Matrice numérique carrée optionnelle \eqn{T_{ij}}
#'   (rownames/colnames = identifiants de cellules). Si \code{NULL}, les
#'   colonnes de transition restent \code{NA} (\code{unavailable_no_transition_matrix}).
#' @param receptor_expr Vecteur numérique nommé (par cellule) ou matrice
#'   numérique gènes x cellules optionnelle. Si \code{NULL}, les colonnes de
#'   corrélation restent \code{NA}.
#' @param filtered_table Table canonique ou pré-filtrée optionnelle (défaut
#'   \code{NULL} = table canonique complète).
#' @return Liste structurée de type \code{"ccc_velocity_directionality"}.
#' @export
compute_communication_velocity_directionality <- function(communication_result,
                                                          velocity_result,
                                                          cell_identities,
                                                          transition_matrix = NULL,
                                                          receptor_expr = NULL,
                                                          filtered_table = NULL) {
  r <- assert_communication_result(communication_result,
                                   context = "directionnalité vélocité")
  vr <- assert_velocity_result(velocity_result, view = "any",
                               context = "directionnalité vélocité communication")
  .vdir_validate_identities(cell_identities)

  warnings_all <- character(0)

  # --- 1. Préparation de la matrice de transition T_ij (optionnelle) ---
  T_sub <- NULL
  bg_mean <- NA_real_
  t_cells <- character(0)
  if (is.null(transition_matrix)) {
    transition_status <- "unavailable_no_transition_matrix"
  } else {
    .vdir_validate_transition_matrix(transition_matrix)
    t_cells <- intersect(
      names(cell_identities),
      intersect(rownames(transition_matrix), colnames(transition_matrix))
    )
    if (length(t_cells) < 2L) {
      transition_status <- "insufficient_overlap"
      warnings_all <- c(warnings_all, paste0(
        "Moins de 2 cellules communes entre 'cell_identities' et ",
        "'transition_matrix' — alignement T_ij non calculé."
      ))
    } else {
      T_sub <- as.matrix(transition_matrix[t_cells, t_cells, drop = FALSE])
      off_diag <- T_sub[row(T_sub) != col(T_sub)]
      bg_mean <- if (any(is.finite(off_diag))) {
        mean(off_diag, na.rm = TRUE)
      } else {
        mean(T_sub, na.rm = TRUE)
      }
      transition_status <- "valid"
    }
  }

  # --- 2. Préparation des vecteurs de vitesse et expression récepteur ---
  rec_spec <- if (!is.null(receptor_expr)) .vdir_validate_receptor_expr(receptor_expr) else NULL
  vectors <- vr$velocity_vectors
  mag <- NULL
  if (is.null(rec_spec)) {
    cor_global_status <- "unavailable_no_receptor_expr"
  } else if (is.null(vectors)) {
    cor_global_status <- "unavailable_no_velocity_vectors"
    warnings_all <- c(warnings_all, paste0(
      "Résultat velocity sans vecteurs validés — corrélation vitesse-récepteur ",
      "indisponible (aucun vecteur fabriqué)."
    ))
  } else {
    mag <- setNames(sqrt(rowSums(as.matrix(vectors)^2)), rownames(vectors))
    v_cells <- intersect(names(cell_identities), names(mag))
    if (length(v_cells) < 3L) {
      cor_global_status <- "insufficient_overlap"
      warnings_all <- c(warnings_all, paste0(
        "Moins de 3 cellules communes entre 'cell_identities' et les vecteurs ",
        "velocity — corrélation Spearman non calculable."
      ))
    } else {
      cor_global_status <- "valid"
    }
  }

  # --- 3. Évaluation par paire (sender_node, receiver_node) ---
  tbl <- .vdir_resolve_nodes(
    if (is.null(filtered_table)) r$canonical_table else as.data.frame(filtered_table)
  )
  pairs <- unique(paste(tbl$sender_node, tbl$receiver_node, sep = "\r"))
  all_cells <- names(cell_identities)
  t_ids <- if (length(t_cells) > 0L) cell_identities[t_cells] else character(0)

  pair_rows <- list()
  unmapped <- character(0)

  for (k in pairs) {
    sp <- strsplit(k, "\r", fixed = TRUE)[[1L]]
    s_node <- sp[1L]
    r_node <- sp[2L]
    sub <- tbl[tbl$sender_node == s_node & tbl$receiver_node == r_node, , drop = FALSE]
    s_all <- all_cells[cell_identities == s_node]
    r_all <- all_cells[cell_identities == r_node]
    if (length(s_all) == 0L || length(r_all) == 0L) {
      unmapped <- c(unmapped, paste(sp, collapse = " -> "))
      next
    }

    # Transition alignment T_{sender -> receiver}
    t_mean <- NA_real_
    t_sum <- NA_real_
    t_bg <- if (identical(transition_status, "valid")) bg_mean else NA_real_
    t_enrich <- NA_real_
    t_aligned <- NA
    if (identical(transition_status, "valid")) {
      si <- t_cells[t_ids == s_node]
      ri <- t_cells[t_ids == r_node]
      if (length(si) > 0L && length(ri) > 0L) {
        sub_T <- T_sub[si, ri, drop = FALSE]
        t_mean <- mean(sub_T, na.rm = TRUE)
        t_sum <- sum(sub_T, na.rm = TRUE)
        if (is.finite(bg_mean) && bg_mean > 0) {
          t_enrich <- t_mean / bg_mean
        }
        if (is.finite(t_mean) && is.finite(bg_mean)) {
          t_aligned <- isTRUE(t_mean > bg_mean)
        }
      }
    }

    # Within-population Spearman correlation (velocity magnitude vs receptor expr)
    c_rho <- NA_real_
    c_n <- 0L
    c_st <- cor_global_status
    c_matched <- NA_character_
    if (identical(cor_global_status, "valid")) {
      ext <- .vdir_extract_receptor_vec(rec_spec, sub$receptor, r_all)
      if (length(ext$matched) == 0L) {
        c_st <- "receptor_not_found"
      } else {
        c_matched <- paste(ext$matched, collapse = ",")
        sp_res <- .vdir_spearman_within_pop(mag, ext$vec, r_all)
        c_rho <- sp_res$cor
        c_n <- sp_res$n
        c_st <- sp_res$status
      }
    }

    s_sc <- sub$score[!is.na(sub$score)]
    pair_rows[[length(pair_rows) + 1L]] <- data.frame(
      sender_node = s_node,
      receiver_node = r_node,
      n_interactions = as.integer(nrow(sub)),
      mean_imported_score = if (length(s_sc)) mean(s_sc) else NA_real_,
      n_sender_cells = as.integer(length(s_all)),
      n_receiver_cells = as.integer(length(r_all)),
      transition_mean = t_mean,
      transition_sum = t_sum,
      transition_bg_mean = t_bg,
      transition_enrichment = t_enrich,
      transition_aligned = t_aligned,
      cor_spearman = c_rho,
      cor_n_cells = as.integer(c_n),
      cor_status = c_st,
      receptors_matched = c_matched,
      stringsAsFactors = FALSE
    )
  }

  pair_table <- if (length(pair_rows) > 0L) {
    do.call(rbind, pair_rows)
  } else {
    data.frame(
      sender_node = character(0), receiver_node = character(0),
      n_interactions = integer(0), mean_imported_score = numeric(0),
      n_sender_cells = integer(0), n_receiver_cells = integer(0),
      transition_mean = numeric(0), transition_sum = numeric(0),
      transition_bg_mean = numeric(0), transition_enrichment = numeric(0),
      transition_aligned = logical(0), cor_spearman = numeric(0),
      cor_n_cells = integer(0), cor_status = character(0),
      receptors_matched = character(0), stringsAsFactors = FALSE
    )
  }

  overall_status <- if (nrow(pair_table) == 0L) {
    "insufficient_overlap"
  } else if (identical(transition_status, "valid") || identical(cor_global_status, "valid")) {
    "valid"
  } else {
    "unavailable"
  }

  qc <- list(
    n_cells_identities = as.integer(length(cell_identities)),
    n_cells_transition = as.integer(length(t_cells)),
    n_pairs_total = as.integer(length(pairs)),
    n_pairs_mapped = as.integer(nrow(pair_table)),
    n_pairs_unmapped = as.integer(length(unmapped)),
    unmapped_pairs = unmapped,
    n_pairs_transition_aligned = as.integer(sum(pair_table$transition_aligned %in% TRUE)),
    n_pairs_cor_valid = as.integer(sum(pair_table$cor_status == "valid"))
  )

  entry <- new_provenance_entry(
    analysis_id = "sc-communication-velocity-directional",
    method = "velocity_transition_alignment_spearman",
    parameters = list(
      transition_status = transition_status,
      cor_status = cor_global_status,
      n_pairs_mapped = qc$n_pairs_mapped,
      n_pairs_aligned = qc$n_pairs_transition_aligned,
      n_pairs_cor_valid = qc$n_pairs_cor_valid,
      parent_analysis_id = r$analysis_id
    ),
    warnings = warnings_all
  )
  entry$analysis_type <- "cell_cell_communication_velocity_directional"
  entry$import_only <- FALSE
  entry$parent_analysis_id <- r$analysis_id
  entry$status <- overall_status

  list(
    type = "ccc_velocity_directionality",
    status = overall_status,
    transition_status = transition_status,
    cor_status = cor_global_status,
    analysis_id = "sc-communication-velocity-directional",
    parent_analysis_id = r$analysis_id,
    source_method = r$source_method,
    velocity_status = vr$status,
    velocity_analysis_id = vr$analysis_id %||% NA_character_,
    pair_table = pair_table,
    qc = qc,
    warnings = warnings_all,
    epistemic_note = .vdir_epistemic_note(),
    provenance = entry,
    timestamp_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
  )
}

#' Visualisation directionnelle vélocité CCC (alignement T_ij / corrélation Spearman)
#'
#' Produit une figure \code{ggplot2} pure des paires (sender -> receiver) par
#' enrichissement de transition \eqn{T_{ij}} (ou corrélation de Spearman
#' vitesse-récepteur lorsque \eqn{T_{ij}} est absente), ou un graphe explicite
#' lorsque les métriques sont indisponibles.
#'
#' @param directionality_result Objet produit par
#'   \code{compute_communication_velocity_directionality()}.
#' @param top_n Nombre maximal de paires affichées (défaut \code{25L}).
#' @return Objet \code{ggplot}.
#' @export
plot_communication_velocity_directionality <- function(directionality_result,
                                                       top_n = 25L) {
  if (!is.list(directionality_result) ||
      !identical(directionality_result$type %||% NULL, "ccc_velocity_directionality")) {
    .vdir_stop(paste0(
      "plot_communication_velocity_directionality() : objet ",
      "'ccc_velocity_directionality' requis."
    ))
  }
  top_n <- max(1L, as.integer(top_n %||% 25L))
  note_wrap <- paste(
    strwrap(directionality_result$epistemic_note %||% .vdir_epistemic_note(), width = 76L),
    collapse = "\n"
  )

  pt <- directionality_result$pair_table
  has_t <- !is.null(pt) && nrow(pt) > 0L && any(is.finite(pt$transition_enrichment))
  has_c <- !is.null(pt) && nrow(pt) > 0L && any(is.finite(pt$cor_spearman))

  if (!identical(directionality_result$status, "valid") || (!has_t && !has_c)) {
    msg <- sprintf(
      "Directionnalité vélocité indisponible (statut : %s).",
      directionality_result$status %||% NA
    )
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                          size = 4.5, colour = "grey40") +
        ggplot2::theme_void() +
        ggplot2::labs(
          title = "Directionnalité vélocité CCC",
          subtitle = paste(msg, note_wrap, sep = "\n")
        )
    )
  }

  pt$pair_label <- paste(pt$sender_node, pt$receiver_node, sep = " -> ")

  if (has_t) {
    df <- pt[is.finite(pt$transition_enrichment), , drop = FALSE]
    df <- df[order(-df$transition_enrichment, -abs(df$cor_spearman %||% 0),
                   df$pair_label, na.last = TRUE), , drop = FALSE]
    if (nrow(df) > top_n) df <- df[seq_len(top_n), , drop = FALSE]
    df$pair_label <- factor(df$pair_label, levels = rev(df$pair_label))
    df$align_label <- ifelse(
      is.na(df$transition_aligned), "Indéterminé",
      ifelse(df$transition_aligned, "Aligné (> fond)", "Non aligné (<= fond)")
    )
    return(
      ggplot2::ggplot(df, ggplot2::aes(x = transition_enrichment, y = pair_label,
                                       colour = align_label)) +
        ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
        ggplot2::geom_segment(
          ggplot2::aes(x = 1, xend = transition_enrichment,
                       y = pair_label, yend = pair_label),
          linewidth = 0.8, alpha = 0.8
        ) +
        ggplot2::geom_point(size = 3.2) +
        ggplot2::scale_colour_manual(
          values = c("Aligné (> fond)" = "#2166AC",
                     "Non aligné (<= fond)" = "#B2182B",
                     "Indéterminé" = "grey60"),
          name = NULL
        ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          x = "Enrichissement de transition T_ij (ratio vs fond hors-diagonale)",
          y = NULL,
          title = sprintf("Alignement directionnel T_ij (top %d paires)", nrow(df)),
          subtitle = note_wrap
        )
    )
  }

  # Corrélation Spearman seule (sans matrice de transition T_ij)
  df <- pt[is.finite(pt$cor_spearman), , drop = FALSE]
  df <- df[order(-abs(df$cor_spearman), df$pair_label, na.last = TRUE), , drop = FALSE]
  if (nrow(df) > top_n) df <- df[seq_len(top_n), , drop = FALSE]
  df$pair_label <- factor(df$pair_label, levels = rev(df$pair_label))
  df$sign_label <- ifelse(df$cor_spearman >= 0, "Positive", "Négative")

  ggplot2::ggplot(df, ggplot2::aes(x = cor_spearman, y = pair_label,
                                   colour = sign_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = cor_spearman, y = pair_label, yend = pair_label),
      linewidth = 0.8, alpha = 0.8
    ) +
    ggplot2::geom_point(size = 3.2) +
    ggplot2::scale_colour_manual(
      values = c("Positive" = "#2166AC", "Négative" = "#B2182B"),
      name = NULL
    ) +
    ggplot2::coord_cartesian(xlim = c(-1, 1)) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Corrélation de Spearman (magnitude vitesse vs expression récepteur)",
      y = NULL,
      title = sprintf("Corrélation vitesse-récepteur intra-population (top %d paires)", nrow(df)),
      subtitle = note_wrap
    )
}
