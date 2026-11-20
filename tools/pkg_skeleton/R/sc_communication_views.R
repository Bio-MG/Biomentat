# =============================================================================
# R/sc/sc_communication_views.R — Communication 4D-2 : vues exploratoires des
# imports valides (Stage 12)
# =============================================================================
# Consommateurs PURES du resultat canonique (assert_communication_result) :
# filtres d'affichage, DotPlot sender-receiver, heatmap de pathways, reseau
# circulaire, centralite descriptive, provenance des filtres et exports.
#
# REGLES (Stage 12) :
#   - une vue ne cree, n'infere et ne repare JAMAIS un resultat ; elle lit la
#     table (fragmente par les filtres) et n'ecrit rien dans le resultat ;
#   - le type de score ET la methode source sont affiches sur CHAQUE figure
#     (sous-titre) : les scores de sources differentes ne sont PAS comparables ;
#   - les p-values sont etiquetees "importees" — aucune reinference ;
#   - une association ligand-receptor n'etablit aucune causalite (legende) ;
#   - la centralite est "derivede du reseau, descriptive" — PAS un controle
#     biologique ;
#   - les filtres sont des operations d'AFFICHAGE : la table canonique du
#     resultat n'est jamais modifiee ; les filtres actifs sont produit dans la
#     provenance (build_communication_filter_provenance) au moment de l'export.
#
# Nuds du reseau : identite harmonisee (sender_mapped/receiver_mapped) quand
# elle existe, sinon label brut — coalescence explicite (.communication_node_keys).
#
# Sourced in app.R AFTER R/sc/sc_communication.R. Pure ggplot2 (aucune
# dependance nouvelle : pas d'igraph/ggraph/circlize — courbes de Bezier
# quadratiques echantillonnees a la main).
#
# ── API PUBLIQUE FIGEE (Stage 12) ───────────────────────────────────────────
# communication_views_public_api() enumere la surface ; le test de freeze
# (test-communication-contract-freeze.R) la verifie. Helper interne :
#   .views_stop()            — erreur classee communication_import_error
#   .communication_node_keys() — coalescence mapped/raw
#   .communication_agg_pairs() — aggregation par paire de noeuds
#   .communication_empty_view_plot() — message explicite au lieu d'un graphe vide
# =============================================================================

#' Surface publique figee de R/sc/sc_communication_views.R (Stage 12)
#'
#' @return Vecteur character des noms de fonctions publiques.
#' @export
communication_views_public_api <- function() {
  c(
    "communication_views_public_api",
    "communication_apply_filters",
    "plot_communication_dotplot",
    "plot_communication_pathway_heatmap",
    "plot_communication_circle",
    "build_communication_centrality",
    "plot_communication_roles_dot",
    "plot_communication_centrality_quadrant",
    "build_communication_filter_provenance",
    "build_communication_filtered_export",
    # 5e pilule exploratoire « Synthèse multi-modale » (d183d11) — registre
    # reconcilié avec la surface du fichier (règle repo : pin au MÊME commit).
    "build_communication_multimodal_synthesis",
    "plot_communication_multimodal_synthesis",
    # CC-3 : graphe comparatif de flux par voie (@netP$prob).
    "plot_ccc_comparison_flow",
    # W2.2 §11.8 Gem #2 : consensus multi-méthodes par rangs intra-méthode.
    "build_communication_consensus",
    "plot_communication_consensus",
    # W2.2 §11.8 Gem #4 : vue alluviale 3 tiers Sender -> Pathway/Ligand -> Receiver.
    "plot_communication_alluvial",
    # W2.2 §11.8 Gem #3 : vérification descriptive d'expression L-R single-cell.
    "build_communication_lr_expression",
    "plot_communication_lr_expression",
    # W2.2 §11.8 Gem #1 : comparaison différentielle Condition A vs B.
    "build_communication_condition_diff",
    "plot_communication_condition_diff"
  )
}

#' Erreur d'une vue communication (classe du domaine, etat structure)
.views_stop <- function(message) {
  stop(errorCondition(
    message,
    state = "invalid_input",
    class = "communication_import_error"
  ))
}

#' Coalescence des noeuds : identite harmonisee sinon label brut
#'
#' Ajoute sender_node/receiver_node : la valeur harmonisee (exact match vers
#' la colonne choisie a l'import) quand elle existe, sinon le label brut —
#' la fusion de labels est donc TOUJOURS tracee par la colonne d'identites
#' citee dans les sous-titres.
.communication_node_keys <- function(table) {
  table$sender_node <- ifelse(
    !is.na(table$sender_mapped %||% rep(NA_character_, nrow(table))),
    table$sender_mapped, table$sender
  )
  table$receiver_node <- ifelse(
    !is.na(table$receiver_mapped %||% rep(NA_character_, nrow(table))),
    table$receiver_mapped, table$receiver
  )
  table
}

#' Aggregation descriptive par paire de noeuds (aucune reinference)
#'
#' @return data.frame(sender_node, receiver_node, n_interactions, mean_score,
#'   sum_score, n_with_score).
.communication_agg_pairs <- function(table) {
  keys <- unique(paste(table$sender_node, table$receiver_node, sep = "\r"))
  do.call(rbind, lapply(keys, function(k) {
    sp <- strsplit(k, "\r", fixed = TRUE)[[1L]]
    sub <- table[table$sender_node == sp[1L] & table$receiver_node == sp[2L], ]
    sc <- sub$score[!is.na(sub$score)]
    out <- data.frame(
      sender_node = sp[1L], receiver_node = sp[2L],
      n_interactions = nrow(sub),
      mean_score = if (length(sc)) mean(sc) else NA_real_,
      sum_score = if (length(sc)) sum(sc) else NA_real_,
      n_with_score = length(sc),
      stringsAsFactors = FALSE
    )
    if ("rank" %in% colnames(sub)) {
      rk <- suppressWarnings(as.numeric(sub$rank))
      rk <- rk[is.finite(rk)]
      out$mean_rank <- if (length(rk)) mean(rk) else NA_real_
      out$best_rank <- if (length(rk)) min(rk) else NA_real_
    }
    out
  }))
}

#' Garde-fou top-N generique sur un data.frame agrege (affichage uniquement)
.communication_top_n_rows <- function(df, n, primary, secondary = NULL) {
  if (is.null(n) || !is.numeric(n) || length(n) != 1L ||
      !is.finite(n) || n <= 0) {
    return(df)
  }
  n_keep <- as.integer(n)
  if (n_keep <= 0L || nrow(df) <= n_keep) {
    return(df)
  }
  ord <- if (is.null(secondary)) {
    order(primary, decreasing = TRUE, na.last = TRUE)
  } else {
    order(primary, secondary, decreasing = TRUE, na.last = TRUE)
  }
  df[ord[seq_len(n_keep)], , drop = FALSE]
}

#' Message explicite au lieu d'un graphique vide (aucun graphe trompeur)
#'
#' Le message est porte a la fois dans le sous-titre (visible et testable)
#' et en annotation centrale de la figure.
.communication_empty_view_plot <- function(message, title, subtitle) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.5, label = message,
                      size = 4.5, colour = "grey40") +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, subtitle = paste(message, subtitle, sep = " — "))
}

# Sous-titre commun : type de score + methode source (garde-fou d'echelle).
.communication_score_caption <- function(communication_result) {
  if ("rank" %in% colnames(communication_result$canonical_table)) {
    return(sprintf(
      paste0("Source : %s — rang importe, rang 1 = meilleur ",
             "(echelle de la source ; non comparable entre modes). ",
             "P-values importees, non re-inferrees. Association ",
             "ligand-receptor = aucune causalite."),
      communication_result$source_method %||% NA_character_
    ))
  }
  sprintf(
    paste0("Source : %s — score importe, echelle de la source ",
           "(non comparable entre sources). P-values importees, non ",
           "re-inferrees. Association ligand-receptor = aucune causalite."),
    communication_result$source_method %||% NA_character_
  )
}

#' Appliquer les filtres d'affichage (operation d'affichage, jamais destructive)
#'
#' Filtre la table canonique pour les vues. La table du RESULTAT n'est jamais
#' modifiee : une copie filtree est retournee. Semantique explicite :
#'   - score_min / p_value_max : une ligne sans valeur du champ filtre est
#'     retiree (un filtre non verifiable ne passe pas) — le retrait est
#'     comptabilise dans summary ;
#'   - pathways / senders / receivers : appartenance ; les filtres
#'     sender/receiver portent sur le NOEUD coalescent (identite harmonisee
#'     sinon label brut) ;
#'   - include_self = FALSE retire les auto-interactions (sender = receiver).
#'
#' @param communication_result Resultat canonique (assert_communication_result).
#' @param filters list(score_min = NA_real_, p_value_max = NA_real_,
#'   pathways = character(0), senders = character(0), receivers = character(0),
#'   include_self = TRUE) — champs absents pris par defaut.
#' @return list(table = copie filtree (+ sender_node/receiver_node),
#'   filters = filtres normalises, summary = list(n_before, n_after,
#'   dropped_score, dropped_p_value, dropped_pathway, dropped_sender,
#'   dropped_receiver, dropped_self), description = chaine lisible).
#' @export
communication_apply_filters <- function(communication_result,
                                        filters = list()) {
  r <- assert_communication_result(
    communication_result, context = "filtres communication"
  )
  rank_max_raw <- filters$rank_max %||% NA_real_
  if (length(rank_max_raw) != 1L) {
    .views_stop("communication_apply_filters() : rank_max doit etre un scalaire.")
  }
  rank_max <- suppressWarnings(as.numeric(rank_max_raw))
  if (is.na(rank_max) && !is.na(rank_max_raw)) {
    .views_stop("communication_apply_filters() : rank_max doit etre numerique.")
  }
  if (!is.na(rank_max) && rank_max < 0) {
    .views_stop("communication_apply_filters() : rank_max doit etre positif ou nul.")
  }
  f <- list(
    score_min = filters$score_min %||% NA_real_,
    rank_max = rank_max,
    p_value_max = filters$p_value_max %||% NA_real_,
    pathways = as.character(filters$pathways %||% character(0)),
    senders = as.character(filters$senders %||% character(0)),
    receivers = as.character(filters$receivers %||% character(0)),
    include_self = filters$include_self %||% TRUE
  )

  table <- .communication_node_keys(r$canonical_table)
  n_before <- nrow(table)
  keep <- rep(TRUE, n_before)
  score_values <- table$score
  score_available <- any(!is.na(score_values))
  rank_available <- "rank" %in% colnames(table)
  score_requested <- !is.na(f$score_min)
  rank_requested <- !is.na(f$rank_max)

  if (score_requested && rank_requested) {
    .views_stop(paste0(
      "communication_apply_filters() : score_min et rank_max sont incompatibles ",
      "sur une meme selection ; un rang n'est pas un score."
    ))
  }
  if (rank_requested && !rank_available) {
    .views_stop(paste0(
      "communication_apply_filters() : rank_max demande mais aucune colonne ",
      "rank n'est presente dans le resultat."
    ))
  }
  if (score_requested && rank_available && !score_available) {
    .views_stop(paste0(
      "communication_apply_filters() : score_min demande mais le resultat ",
      "ne contient aucune valeur de score ; un rang n'est pas un score."
    ))
  }

  if (score_requested) {
    keep <- keep & !is.na(score_values) & score_values >= f$score_min
  }
  dropped_score <- sum(!keep)
  dropped_rank <- 0L
  if (rank_requested) {
    rank_values <- suppressWarnings(as.numeric(table$rank))
    keep2 <- is.finite(rank_values) & rank_values <= f$rank_max
    dropped_rank <- sum(keep & !keep2)
    keep <- keep & keep2
  }
  if (!is.na(f$p_value_max)) {
    keep2 <- !is.na(table$p_value) & table$p_value <= f$p_value_max
    dropped_p <- sum(keep & !keep2)
    keep <- keep & keep2
  } else dropped_p <- 0L
  if (length(f$pathways)) {
    keep2 <- !is.na(table$pathway) & table$pathway %in% f$pathways
    dropped_pw <- sum(keep & !keep2)
    keep <- keep & keep2
  } else dropped_pw <- 0L
  if (length(f$senders)) {
    keep2 <- table$sender_node %in% f$senders
    dropped_s <- sum(keep & !keep2)
    keep <- keep & keep2
  } else dropped_s <- 0L
  if (length(f$receivers)) {
    keep2 <- table$receiver_node %in% f$receivers
    dropped_r <- sum(keep & !keep2)
    keep <- keep & keep2
  } else dropped_r <- 0L
  if (isFALSE(f$include_self)) {
    keep2 <- !(table$sender_node == table$receiver_node)
    dropped_self <- sum(keep & !keep2)
    keep <- keep & keep2
  } else dropped_self <- 0L

  table <- table[keep, , drop = FALSE]

  description <- sprintf(
    paste0("score_min=%s ; rank_max=%s ; p_value_max=%s ; pathways=[%s] ; ",
           "senders=[%s] ; receivers=[%s] ; include_self=%s ; dropped_rank=%s"),
    ifelse(is.na(f$score_min), "NA", format(f$score_min)),
    ifelse(is.na(f$rank_max), "NA", format(f$rank_max)),
    ifelse(is.na(f$p_value_max), "NA", format(f$p_value_max)),
    paste(f$pathways, collapse = ","), paste(f$senders, collapse = ","),
    paste(f$receivers, collapse = ","), as.character(f$include_self),
    as.character(dropped_rank)
  )

  list(
    table = table,
    filters = f,
    summary = list(
      n_before = as.integer(n_before),
      n_after = as.integer(nrow(table)),
      dropped_score = as.integer(dropped_score),
      dropped_rank = as.integer(dropped_rank),
      dropped_p_value = as.integer(dropped_p),
      dropped_pathway = as.integer(dropped_pw),
      dropped_sender = as.integer(dropped_s),
      dropped_receiver = as.integer(dropped_r),
      dropped_self = as.integer(dropped_self)
    ),
    description = description
  )
}

#' Vue DotPlot sender-receiver (aggregation descriptive)
#'
#' Taille = nombre d'interactions importees par paire de noeuds ; couleur =
#' score moyen importe (absent si aucun score — le point porte alors
#' uniquement l'effectif, sans substitution). Sous-titre : methode source et
#' type de score (garde-fou d'echelle).
#'
#' @param communication_result Resultat canonique.
#' @param filtered_table Table filtree (communication_apply_filters()$table)
#'   ou NULL pour utiliser la table complete du resultat.
#' @param seurat_obj Objet Seurat courant optionnel (peremption).
#' @param top_n_pairs Entier positif optionnel : limite aux N paires
#'   sender-receiver principales (NULL = toutes).
#' @return ggplot.
#' @export
plot_communication_dotplot <- function(communication_result, filtered_table = NULL,
                                       seurat_obj = NULL, top_n_pairs = NULL) {
  r <- assert_communication_result(
    communication_result, seurat_obj = seurat_obj,
    context = "DotPlot sender-receiver"
  )
  table <- .communication_node_keys(
    if (is.null(filtered_table)) r$canonical_table else filtered_table
  )
  subtitle_base <- .communication_score_caption(r)

  if (nrow(table) == 0L) {
    return(.communication_empty_view_plot(
      "Aucune interaction apres filtres : elargissez les filtres ou reimportez.",
      "DotPlot sender-receiver", subtitle_base
    ))
  }

  agg <- .communication_agg_pairs(table)
  has_rank <- "mean_rank" %in% colnames(agg) &&
    any(is.finite(agg$mean_rank))
  has_score <- any(!is.na(agg$mean_score))
  pair_weight <- if (has_rank) {
    agg$n_interactions
  } else if (has_score) {
    ifelse(is.na(agg$sum_score), agg$n_interactions, agg$sum_score)
  } else {
    agg$n_interactions
  }
  agg <- .communication_top_n_rows(agg, top_n_pairs, pair_weight, agg$n_interactions)

  p <- ggplot2::ggplot(
    agg,
    ggplot2::aes(x = receiver_node, y = sender_node,
                 size = n_interactions,
                 colour = if (has_rank) mean_rank else if (has_score) mean_score else n_interactions)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_size(range = c(2, 9), name = "Nb interactions") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      x = "Receiver (noeud harmonise)", y = "Sender (noeud harmonise)",
      title = "Interactions sender-receiver (resultats importes)",
      subtitle = sprintf(
        paste0("%s — %d interaction(s) affichee(s). %s"),
        subtitle_base, sum(agg$n_interactions),
        if (has_rank)
          "Couleur : rang moyen importe par paire (1 = meilleur)."
        else if (has_score)
          "Couleur : score moyen importe par paire."
        else
          "Aucun score dans la selection : l'effectif seule est affiche."
      )
    )
  if (has_rank) {
    p <- p + ggplot2::scale_colour_gradient(
      low = "#9ECAE1", high = "#084594",
      trans = scales::reverse_trans(),
      name = "rang moyen (1 = meilleur)"
    )
  } else if (has_score) {
    p <- p + ggplot2::scale_colour_gradient(
      low = "#9ECAE1", high = "#084594", name = "Score moyen importe"
    )
  } else {
    p <- p + ggplot2::scale_colour_gradient(
      low = "#9ECAE1", high = "#084594", name = "Nb interactions"
    )
  }
  p
}

#' Vue heatmap pathways x paires (aggregation descriptive)
#'
#' Une case = score moyen importe d'un pathway pour une paire de noeuds. Les
#' lignes sans pathway (NA) sont EXCLUES de la vue et comptabilisees dans le
#' sous-titre — jamais imputees a un pathway invente.
#'
#' @param communication_result Resultat canonique.
#' @param filtered_table Table filtree ou NULL (table complete).
#' @param seurat_obj Objet Seurat courant optionnel (peremption).
#' @return ggplot.
#' @export
plot_communication_pathway_heatmap <- function(communication_result,
                                               filtered_table = NULL,
                                               seurat_obj = NULL,
                                               top_n_pathways = NULL) {
  r <- assert_communication_result(
    communication_result, seurat_obj = seurat_obj,
    context = "heatmap pathways"
  )
  table <- .communication_node_keys(
    if (is.null(filtered_table)) r$canonical_table else filtered_table
  )
  subtitle_base <- .communication_score_caption(r)
  n_no_pathway <- sum(is.na(table$pathway))
  table <- table[!is.na(table$pathway), , drop = FALSE]

  if (nrow(table) == 0L) {
    return(.communication_empty_view_plot(
      paste0(
        if (n_no_pathway > 0L)
          paste0("Aucun pathway renseigne dans la source (", n_no_pathway,
                 " interaction(s) sans pathway). ")
        else "Aucune interaction apres filtres. ",
        "La vue heatmap necessite des pathways — utilisez une table ",
        "CellChat exportee avec la colonne pathway."
      ),
      "Heatmap pathways", subtitle_base
    ))
  }

  table$pair_label <- paste(table$sender_node, table$receiver_node, sep = " \u2192 ")
  agg <- do.call(rbind, lapply(
    split(table, paste(table$pathway, table$pair_label, sep = "\r")),
    function(sub) data.frame(
      pathway = sub$pathway[1L], pair_label = sub$pair_label[1L],
      mean_score = if (any(!is.na(sub$score))) mean(sub$score, na.rm = TRUE) else NA_real_,
      n_interactions = nrow(sub),
      stringsAsFactors = FALSE
    )
  ))

  if (!is.null(top_n_pathways) && is.finite(top_n_pathways) && top_n_pathways > 0L) {
    pw_rank <- tapply(
      ifelse(is.na(agg$mean_score), agg$n_interactions, agg$mean_score),
      agg$pathway, sum, na.rm = TRUE
    )
    if (length(pw_rank) > top_n_pathways) {
      keep_pw <- names(sort(pw_rank, decreasing = TRUE))[seq_len(as.integer(top_n_pathways))]
      agg <- agg[agg$pathway %in% keep_pw, , drop = FALSE]
    }
  }

  ggplot2::ggplot(agg, ggplot2::aes(x = pair_label, y = pathway,
                                    fill = mean_score)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient(low = "#F7F4F0", high = "#B2182B",
                                 na.value = "grey90",
                                 name = "Score moyen importe") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "Pathways par paire sender-receiver (resultats importes)",
      subtitle = sprintf(
        "%s — %d interaction(s) avec pathway ; %d sans pathway (exclues, jamais imputees).",
        subtitle_base, sum(agg$n_interactions), n_no_pathway
      )
    )
}

#' Vue reseau circulaire (ggplot pur, courbes de Bezier)
#'
#' Noeuds = identites harmonisees (sinon labels bruts) sur un cercle ; aretes =
#' paires sender-receiver aggregatees, largeur/transparence par poids importe,
#' fleche orientee sender -> receiver. Les auto-interactions ne sont pas
#' dessinables sur le cercle : elles sont comptabilisees dans le sous-titre.
#' VUE DESCRIPTIVE — aucune lecture causale.
#'
#' @param communication_result Resultat canonique.
#' @param filtered_table Table filtree ou NULL (table complete).
#' @param seurat_obj Objet Seurat courant optionnel (peremption).
#' @param top_n_edges Entier positif optionnel : limite aux N aretes les plus
#'   ponderees (NULL = toutes).
#' @return ggplot.
#' @export
plot_communication_circle <- function(communication_result, filtered_table = NULL,
                                      seurat_obj = NULL, top_n_edges = NULL) {
  r <- assert_communication_result(
    communication_result, seurat_obj = seurat_obj,
    context = "reseau circulaire"
  )
  table <- .communication_node_keys(
    if (is.null(filtered_table)) r$canonical_table else filtered_table
  )
  subtitle_base <- .communication_score_caption(r)
  n_self <- sum(table$sender_node == table$receiver_node)
  edges_table <- table[table$sender_node != table$receiver_node, , drop = FALSE]

  if (nrow(edges_table) == 0L) {
    return(.communication_empty_view_plot(
      sprintf(
        paste0("Aucune arete dessinable : %d interaction(s) dans la ",
               "selection dont %d auto-interaction(s) (sender = receiver, ",
               "non dessinables sur le cercle)."),
        nrow(table), n_self
      ),
      "Reseau sender-receiver", subtitle_base
    ))
  }

  agg <- .communication_agg_pairs(edges_table)
  has_rank <- "mean_rank" %in% colnames(agg) &&
    any(is.finite(agg$mean_rank))
  has_score <- any(!is.na(agg$sum_score))
  agg$weight <- if (has_rank) agg$n_interactions else if (has_score) agg$sum_score else agg$n_interactions
  agg <- .communication_top_n_rows(agg, top_n_edges, agg$weight, agg$n_interactions)

  nodes <- sort(unique(c(agg$sender_node, agg$receiver_node)))
  n_nodes <- length(nodes)
  angle <- seq(0, 2 * pi, length.out = n_nodes + 1L)[seq_len(n_nodes)]
  node_pos <- data.frame(
    node = nodes,
    x = sin(angle), y = cos(angle),
    stringsAsFactors = FALSE
  )

  # Bezier quadratique : point de controle tire vers le centre (0.35) —
  # courbure deterministe, echantillonnage fixe.
  n_pts <- 40L
  t_seq <- seq(0, 1, length.out = n_pts)
  edge_df <- do.call(rbind, lapply(seq_len(nrow(agg)), function(i) {
    p1 <- node_pos[node_pos$node == agg$sender_node[i], c("x", "y")]
    p2 <- node_pos[node_pos$node == agg$receiver_node[i], c("x", "y")]
    ctrl <- (p1 + p2) * 0.35
    data.frame(
      edge_id = i,
      t = t_seq,
      x = (1 - t_seq)^2 * p1$x + 2 * (1 - t_seq) * t_seq * ctrl$x + t_seq^2 * p2$x,
      y = (1 - t_seq)^2 * p1$y + 2 * (1 - t_seq) * t_seq * ctrl$y + t_seq^2 * p2$y,
      weight = agg$weight[i],
      stringsAsFactors = FALSE
    )
  }))

  node_tot <- do.call(rbind, lapply(node_pos$node, function(nd) {
    out_w <- sum(agg$weight[agg$sender_node == nd], na.rm = TRUE)
    in_w <- sum(agg$weight[agg$receiver_node == nd], na.rm = TRUE)
    data.frame(node = nd, total = out_w + in_w, stringsAsFactors = FALSE)
  }))
  node_pos$total <- node_tot$total[match(node_pos$node, node_tot$node)]

  ggplot2::ggplot() +
    ggplot2::geom_path(
      data = edge_df,
      ggplot2::aes(x = x, y = y, group = edge_id,
                   linewidth = weight, alpha = weight),
      lineend = "round",
      arrow = grid::arrow(length = grid::unit(7, "pt"), type = "closed")
    ) +
    ggplot2::geom_point(
      data = node_pos,
      ggplot2::aes(x = x, y = y, size = total),
      colour = "#2166AC", alpha = 0.9
    ) +
    ggplot2::geom_text(
      data = node_pos,
      ggplot2::aes(x = x * 1.18, y = y * 1.18, label = node),
      size = 3.4
    ) +
    ggplot2::scale_linewidth(range = c(0.2, 2.4), name = if (has_rank) "Nb interactions" else if (has_score) "Poids (score importe)" else "Nb interactions") +
    ggplot2::scale_alpha(range = c(0.25, 0.9), guide = "none") +
    ggplot2::scale_size(range = c(3, 9), guide = "none") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::labs(
      title = "Reseau sender-receiver (resultats importes)",
      subtitle = sprintf(
        paste0("%s — %d arete(s) ; %d auto-interaction(s) non dessinees. ",
               "Vue descriptive : aucune causalite.%s"),
        subtitle_base, nrow(agg), n_self,
        if (has_rank)
          " Poids : nombre d'interactions ; rang 1 = meilleur, jamais somme des rangs."
        else ""
      )
    )
}

#' Centralite reseau descriptive (NIVEAU RESEAU — pas un controle biologique)
#'
#' Par noeud : nombre d'interactions sortantes/entrantes, partenaires distincts,
#' somme des scores importes sortants/entrants. Ces quantites sont derivees du
#' graphe des resultats importes — elles ne mesurent AUCUNE activite
#' biologique.
#'
#' @param communication_result Resultat canonique.
#' @param filtered_table Table filtree ou NULL (table complete).
#' @return data.frame trie par poids total decroissant, trace par analysis_id.
#' @export
build_communication_centrality <- function(communication_result,
                                           filtered_table = NULL) {
  r <- assert_communication_result(
    communication_result, context = "centralite communication"
  )
  table <- .communication_node_keys(
    if (is.null(filtered_table)) r$canonical_table else filtered_table
  )
  nodes <- sort(unique(c(table$sender_node, table$receiver_node)))
  has_rank <- "rank" %in% colnames(table)
  if (length(nodes) == 0L) {
    out <- data.frame(
      node = character(),
      n_out_interactions = integer(),
      n_in_interactions = integer(),
      out_partners = integer(),
      in_partners = integer(),
      out_score_total = numeric(),
      in_score_total = numeric(),
      total_interactions = integer(),
      stringsAsFactors = FALSE
    )
    if (has_rank) {
      out$out_best_rank <- numeric()
      out$out_median_rank <- numeric()
      out$in_best_rank <- numeric()
      out$in_median_rank <- numeric()
    }
    out$analysis_id <- rep(r$analysis_id %||% NA_character_, nrow(out))
    out$source_method <- rep(r$source_method %||% NA_character_, nrow(out))
    out$identity_column <- rep(r$identity_column %||% NA_character_, nrow(out))
    return(out)
  }

  out <- do.call(rbind, lapply(nodes, function(nd) {
    s_rows <- table[table$sender_node == nd, ]
    r_rows <- table[table$receiver_node == nd, ]
    s_sc <- s_rows$score[!is.na(s_rows$score)]
    r_sc <- r_rows$score[!is.na(r_rows$score)]
    out_row <- data.frame(
      node = nd,
      n_out_interactions = nrow(s_rows),
      n_in_interactions = nrow(r_rows),
      out_partners = length(unique(s_rows$receiver_node)),
      in_partners = length(unique(r_rows$sender_node)),
      out_score_total = if (length(s_sc)) sum(s_sc) else NA_real_,
      in_score_total = if (length(r_sc)) sum(r_sc) else NA_real_,
      stringsAsFactors = FALSE
    )
    if (has_rank) {
      s_rk <- suppressWarnings(as.numeric(s_rows$rank))
      r_rk <- suppressWarnings(as.numeric(r_rows$rank))
      s_rk <- s_rk[is.finite(s_rk)]
      r_rk <- r_rk[is.finite(r_rk)]
      out_row$out_best_rank <- if (length(s_rk)) min(s_rk) else NA_real_
      out_row$out_median_rank <- if (length(s_rk)) stats::median(s_rk) else NA_real_
      out_row$in_best_rank <- if (length(r_rk)) min(r_rk) else NA_real_
      out_row$in_median_rank <- if (length(r_rk)) stats::median(r_rk) else NA_real_
    }
    out_row
  }))
  out$total_interactions <- out$n_out_interactions + out$n_in_interactions
  out <- out[order(-out$total_interactions), , drop = FALSE]
  rownames(out) <- NULL
  out$analysis_id <- r$analysis_id %||% NA_character_
  out$source_method <- r$source_method %||% NA_character_
  out$identity_column <- r$identity_column %||% NA_character_
  out
}

#' Entree de provenance pour une exploration filtree (PRODUITE a l'export)
#'
#' Les filtres sont des operations d'affichage ; leur trace provenance est
#' PRODUITE quand l'utilisateur exporte (filtres actifs figes dans les
#' parametres) — jamais reconstruite apres coup. L'appelant (module) append
#' l'entree a l'etat partage.
#'
#' @param communication_result Resultat canonique.
#' @param filter_summary Sortie de communication_apply_filters().
#' @return Entree new_provenance_entry() enrichie (analysis_type, status).
#' @export
build_communication_filter_provenance <- function(communication_result,
                                                  filter_summary) {
  r <- assert_communication_result(
    communication_result, context = "provenance filtres communication"
  )
  if (!is.list(filter_summary) || is.null(filter_summary$description)) {
    .views_stop(
      "build_communication_filter_provenance() : filter_summary doit provenir de communication_apply_filters()."
    )
  }
  entry <- new_provenance_entry(
    analysis_id = "sc-communication-explore",
    method = paste0("explore_", r$source_method %||% "unknown"),
    parameters = list(
      applied_filters = filter_summary$description,
      rank_max = filter_summary$filters$rank_max %||% NA_real_,
      n_rows_before = filter_summary$summary$n_before %||% NA_integer_,
      n_rows_after = filter_summary$summary$n_after %||% NA_integer_,
      dropped_score = filter_summary$summary$dropped_score %||% NA_integer_,
      dropped_rank = filter_summary$summary$dropped_rank %||% NA_integer_,
      dropped_p_value = filter_summary$summary$dropped_p_value %||% NA_integer_,
      dropped_pathway = filter_summary$summary$dropped_pathway %||% NA_integer_,
      dropped_sender = filter_summary$summary$dropped_sender %||% NA_integer_,
      dropped_receiver = filter_summary$summary$dropped_receiver %||% NA_integer_,
      dropped_self = filter_summary$summary$dropped_self %||% NA_integer_,
      identity_column = r$identity_column,
      base_analysis_id = r$analysis_id
    ),
    dataset = NULL,
    warnings = as.character(r$warnings %||% character(0))
  )
  entry$analysis_type <- "cell_cell_communication"
  entry$status <- r$status
  entry$import_only <- FALSE
  entry
}

#' Export CSV de la table filtree (filtres + analysis_id embarques)
#'
#' Chaque ligne porte l'identite de l'analyse de base, l'horodatage et la
#' description des filtres actifs — l'export reste reproductible.
#'
#' @param communication_result Resultat canonique.
#' @param filtered_table Table filtree (communication_apply_filters()$table).
#' @param filter_description Description lisible des filtres (champ description).
#' @return data.frame (n_after lignes).
#' @export
build_communication_filtered_export <- function(communication_result,
                                                filtered_table,
                                                filter_description) {
  r <- assert_communication_result(
    communication_result, context = "export table filtree"
  )
  if (!is.data.frame(filtered_table)) {
    .views_stop("build_communication_filtered_export() : table filtree requise (communication_apply_filters()$table).")
  }
  n <- nrow(filtered_table)
  data.frame(
    filtered_table,
    analysis_id = rep(r$analysis_id %||% NA_character_, n),
    base_timestamp_utc = rep(r$timestamp_utc %||% NA_character_, n),
    applied_filters = rep(as.character(filter_description)[1L], n),
    stringsAsFactors = FALSE
  )
}

#' Dot plot des rôles de signalisation (CC-1) — consommateur pur
#'
#' Visualise la table `cellchat_roles` produite par le moteur CellChat
#' (`run_cellchat()`, voir `.cellchat_engine_signaling_roles()`) : une ligne
#' par voie × groupe, taille = degré pondéré (`outdeg`/`indeg`), la couleur
#' porte le rôle classé (mediator/sender/receiver/minor). DESCRIPTIF : aucun
#' contrôle biologique, aucun seuil de signification — la classification est
#' re-calculable depuis les colonnes exportées.
#'
#' @param roles_table data.frame portant au moins (pathway, group, outdeg,
#'   indeg, role) — colonnes du contrat CC-1.
#' @return objet ggplot.
#' @export
plot_communication_roles_dot <- function(roles_table,
                                         hide_minor = FALSE,
                                         top_n_pathways = NULL) {
  if (!is.data.frame(roles_table) ||
      !all(c("pathway", "group", "outdeg", "indeg", "role") %in%
           colnames(roles_table))) {
    .views_stop("plot_communication_roles_dot() : table roles CC-1 requise (pathway, group, outdeg, indeg, role).")
  }
  df <- roles_table
  if (isTRUE(hide_minor)) {
    df <- df[as.character(df$role) != "minor" & (df$outdeg > 0 | df$indeg > 0), , drop = FALSE]
  }
  if (!is.null(top_n_pathways) && is.finite(top_n_pathways) && top_n_pathways > 0L && nrow(df) > 0L) {
    pw_deg <- tapply(pmax(df$outdeg, df$indeg, na.rm = TRUE), df$pathway, max, na.rm = TRUE)
    if (length(pw_deg) > top_n_pathways) {
      keep_pw <- names(sort(pw_deg, decreasing = TRUE))[seq_len(as.integer(top_n_pathways))]
      df <- df[df$pathway %in% keep_pw, , drop = FALSE]
    }
  }
  if (nrow(df) == 0L) {
    return(.communication_empty_view_plot(
      "Aucun role actif a afficher apres filtrage.",
      "Roles de signalisation par voie (CellChat)", ""
    ))
  }
  df$pathway <- factor(df$pathway, levels = rev(sort(unique(df$pathway))))
  df$role <- factor(df$role,
                    levels = c("sender", "receiver", "mediator", "minor"))
  ggplot2::ggplot(df, ggplot2::aes(x = .data$group, y = .data$pathway)) +
    ggplot2::geom_point(ggplot2::aes(size = pmax(.data$outdeg, .data$indeg),
                                    color = .data$role,
                                    alpha = .data$role)) +
    ggplot2::scale_color_manual(values = c(sender = "#2166AC",
                                           receiver = "#B2182B",
                                           mediator = "#4D4D4D",
                                           minor = "#BDBDBD"),
                                drop = FALSE) +
    ggplot2::scale_alpha_manual(values = c(sender = 0.9,
                                           receiver = 0.9,
                                           mediator = 0.9,
                                           minor = 0.25),
                                guide = "none",
                                drop = FALSE) +
    ggplot2::labs(x = NULL, y = NULL,
                  size = "degré pondéré",
                  color = "rôle",
                  title = "Rôles de signalisation par voie (CellChat)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right",
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

#' Quadrant de centralite sender/receiver (CellChat, ggplot pur)
#'
#' Trace `outdeg` (x) contre `indeg` (y) par groupe cellulaire (agrege sur
#' toutes les voies si `pathway = NULL`, ou pour une voie donnee), avec lignes
#' de reference medianes et etiquettes de groupe/role. VUE DESCRIPTIVE.
#'
#' @param roles_table data.frame portant au moins (pathway, group, outdeg,
#'   indeg, role) — colonnes du contrat CC-1.
#' @param pathway Nom de voie optionnel (NULL = agregation sur toutes voies).
#' @return objet ggplot.
#' @export
plot_communication_centrality_quadrant <- function(roles_table,
                                                   pathway = NULL) {
  if (!is.data.frame(roles_table) ||
      !all(c("pathway", "group", "outdeg", "indeg", "role") %in%
           colnames(roles_table))) {
    .views_stop(paste0(
      "plot_communication_centrality_quadrant() : table roles CC-1 requise ",
      "(pathway, group, outdeg, indeg, role)."
    ))
  }
  df <- roles_table
  pw_label <- NULL
  if (!is.null(pathway) && length(pathway) > 0L &&
      nzchar(as.character(pathway[1L]))) {
    pw_label <- as.character(pathway[1L])
    df <- df[as.character(df$pathway) == pw_label, , drop = FALSE]
  }
  if (nrow(df) == 0L) {
    return(.communication_empty_view_plot(
      if (is.null(pw_label)) "Aucun role actif a afficher."
      else sprintf("Aucun role actif pour la voie '%s'.", pw_label),
      "Quadrant de centralite sender/receiver (CellChat)",
      if (is.null(pw_label)) "" else paste0("Voie : ", pw_label)
    ))
  }

  agg <- do.call(rbind, lapply(
    split(df, as.character(df$group)),
    function(sub) {
      deg_max <- pmax(sub$outdeg, sub$indeg, na.rm = TRUE)
      idx <- if (any(is.finite(deg_max))) which.max(deg_max) else 1L
      data.frame(
        group = as.character(sub$group[1L]),
        outdeg = sum(sub$outdeg, na.rm = TRUE),
        indeg = sum(sub$indeg, na.rm = TRUE),
        role = as.character(sub$role[idx]),
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(agg) <- NULL
  agg$role <- factor(agg$role,
                     levels = c("sender", "receiver", "mediator", "minor"))

  med_out <- stats::median(agg$outdeg, na.rm = TRUE)
  med_in <- stats::median(agg$indeg, na.rm = TRUE)
  sub_txt <- if (is.null(pw_label)) {
    "Agregation sur toutes les voies — lignes pointillees : medianes outdeg / indeg."
  } else {
    sprintf("Voie : %s — lignes pointillees : medianes outdeg / indeg.", pw_label)
  }

  ggplot2::ggplot(agg, ggplot2::aes(x = .data$outdeg, y = .data$indeg)) +
    ggplot2::geom_vline(xintercept = med_out, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_hline(yintercept = med_in, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(size = pmax(.data$outdeg, .data$indeg),
                                     color = .data$role),
                        alpha = 0.9) +
    ggplot2::geom_text(ggplot2::aes(label = .data$group),
                       vjust = -0.7, size = 3.2, check_overlap = TRUE) +
    ggplot2::scale_color_manual(values = c(sender = "#2166AC",
                                           receiver = "#B2182B",
                                           mediator = "#4D4D4D",
                                           minor = "#BDBDBD"),
                                drop = FALSE) +
    ggplot2::labs(
      x = "Degre sortant pondere (outdeg)",
      y = "Degre entrant pondere (indeg)",
      size = "degre max",
      color = "role",
      title = "Quadrant de centralite sender/receiver (CellChat)",
      subtitle = sub_txt
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right",
                   panel.grid.minor = ggplot2::element_blank())
}

#' Synthèse multi-modale des contextes CCC (filtrage exploratoire)
#'
#' Consommateur pur qui triangule les 4 contextes calculés par l'utilisateur
#' (Spatial, Trajectoire, Vélocité, Perturbation) par paire `sender -> receiver`
#' sans recalcul de score ni agrégation statistique. STRICTEMENT EXPLORATOIRE :
#' `is_statistical_pvalue = FALSE`, aucune p-value combinée.
#'
#' @param communication_result Résultat canonique.
#' @param spatial_ctx Contexte spatial (`ccc_spatial_context`) ou NULL.
#' @param trajectory_ctx Contexte trajectoire (`ccc_trajectory_context`) ou NULL.
#' @param velocity_ctx Contexte vélocité (`ccc_velocity_context`) ou NULL.
#' @param perturbation_ctx Contexte perturbation (`ccc_perturbation`) ou NULL.
#' @param filtered_table Table filtrée optionnelle ou NULL (table complète).
#' @param min_spatial_frac Seuil exploratoire de fraction à portée (défaut 0.25).
#' @param min_traj_corr Seuil exploratoire de corrélation ligand/récepteur (défaut 0.3).
#' @return Liste `ccc_multimodal_synthesis` (`pair_table`, `qc`, `label`, `is_statistical_pvalue`).
#' @export
build_communication_multimodal_synthesis <- function(
    communication_result,
    spatial_ctx = NULL,
    trajectory_ctx = NULL,
    velocity_ctx = NULL,
    perturbation_ctx = NULL,
    filtered_table = NULL,
    min_spatial_frac = 0.25,
    min_traj_corr = 0.3) {
  r <- assert_communication_result(
    communication_result,
    context = "Synthèse multi-modale CCC"
  )
  table <- .communication_node_keys(
    if (is.null(filtered_table)) r$canonical_table else filtered_table
  )
  epistemic_label <- "Convergence des contextes (filtrage exploratoire)"

  if (nrow(table) == 0L) {
    return(list(
      type = "ccc_multimodal_synthesis",
      label = epistemic_label,
      is_statistical_pvalue = FALSE,
      pair_table = data.frame(
        sender_node = character(0), receiver_node = character(0),
        pair_label = character(0), n_interactions = integer(0),
        mean_imported_score = numeric(0),
        ctx_spatial = logical(0), ctx_trajectory = logical(0),
        ctx_velocity = logical(0), ctx_perturbation = logical(0),
        spatial_frac = numeric(0), traj_corr = numeric(0),
        velocity_recv_mean = numeric(0), pert_delta_score = numeric(0),
        n_contexts_evaluated = integer(0), n_contexts_convergent = integer(0),
        convergence_tier = character(0), context_badges = character(0),
        stringsAsFactors = FALSE
      ),
      qc = list(
        n_pairs = 0L,
        contexts_active = c(spatial = FALSE, trajectory = FALSE,
                            velocity = FALSE, perturbation = FALSE),
        n_contexts_active = 0L
      )
    ))
  }

  agg <- .communication_agg_pairs(table)
  keys <- paste(agg$sender_node, agg$receiver_node, sep = "\r")
  pair_labels <- paste(agg$sender_node, agg$receiver_node, sep = " -> ")

  has_spatial <- is.list(spatial_ctx) &&
    identical(spatial_ctx$type %||% NULL, "ccc_spatial_context") &&
    is.data.frame(spatial_ctx$pair_table) && nrow(spatial_ctx$pair_table) > 0L
  has_traj <- is.list(trajectory_ctx) &&
    identical(trajectory_ctx$type %||% NULL, "ccc_trajectory_context") &&
    is.data.frame(trajectory_ctx$pair_bin_table) && nrow(trajectory_ctx$pair_bin_table) > 0L
  has_vel <- is.list(velocity_ctx) &&
    identical(velocity_ctx$type %||% NULL, "ccc_velocity_context") &&
    identical(velocity_ctx$status %||% NULL, "valid") &&
    is.data.frame(velocity_ctx$pair_table) && nrow(velocity_ctx$pair_table) > 0L
  has_pert <- is.list(perturbation_ctx) &&
    identical(perturbation_ctx$type %||% NULL, "ccc_perturbation") &&
    is.data.frame(perturbation_ctx$delta_table) && nrow(perturbation_ctx$delta_table) > 0L

  # 1. Spatial proximity
  ctx_spatial <- rep(NA, nrow(agg))
  spatial_frac <- rep(NA_real_, nrow(agg))
  if (has_spatial) {
    sp_df <- spatial_ctx$pair_table
    sp_keys <- paste(sp_df$sender_node, sp_df$receiver_node, sep = "\r")
    m <- match(keys, sp_keys)
    for (i in seq_along(keys)) {
      if (is.na(m[i])) next
      row <- sp_df[m[i], , drop = FALSE]
      fr <- suppressWarnings(as.numeric(row$frac_receiver_within_radius[1L]))
      fs <- suppressWarnings(as.numeric(row$frac_sender_within_radius[1L]))
      pp <- suppressWarnings(as.numeric(row$p_perm[1L]))
      best_fr <- suppressWarnings(max(c(fr, fs), na.rm = TRUE))
      if (!is.finite(best_fr)) best_fr <- NA_real_
      spatial_frac[i] <- best_fr
      if (is.na(pp) && is.na(best_fr)) {
        ctx_spatial[i] <- NA
      } else {
        ctx_spatial[i] <- isTRUE(!is.na(pp) && pp <= 0.05) ||
          isTRUE(!is.na(best_fr) && best_fr >= min_spatial_frac)
      }
    }
  }

  # 2. Trajectory co-expression
  ctx_trajectory <- rep(NA, nrow(agg))
  traj_corr <- rep(NA_real_, nrow(agg))
  if (has_traj) {
    tb <- trajectory_ctx$pair_bin_table
    tb_keys <- paste(tb$sender_node, tb$receiver_node, sep = "\r")
    for (i in seq_along(keys)) {
      sub <- tb[tb_keys == keys[i], , drop = FALSE]
      if (nrow(sub) == 0L) next
      lx <- suppressWarnings(as.numeric(sub$mean_ligand_expression_senders))
      rx <- suppressWarnings(as.numeric(sub$mean_receptor_expression_receivers))
      ok <- is.finite(lx) & is.finite(rx)
      if (sum(ok) >= 3L && stats::sd(lx[ok]) > 0 && stats::sd(rx[ok]) > 0) {
        cc <- suppressWarnings(stats::cor(lx[ok], rx[ok]))
        traj_corr[i] <- cc
        ctx_trajectory[i] <- isTRUE(is.finite(cc) && cc >= min_traj_corr)
      } else if (any(ok)) {
        ctx_trajectory[i] <- any(lx[ok] > 0 & rx[ok] > 0)
      } else {
        ctx_trajectory[i] <- FALSE
      }
    }
  }

  # 3. Velocity-active receiver
  ctx_velocity <- rep(NA, nrow(agg))
  velocity_recv_mean <- rep(NA_real_, nrow(agg))
  if (has_vel) {
    vl_df <- velocity_ctx$pair_table
    vl_keys <- paste(vl_df$sender_node, vl_df$receiver_node, sep = "\r")
    med_recv <- stats::median(suppressWarnings(as.numeric(vl_df$receiver_velocity_mean)), na.rm = TRUE)
    m <- match(keys, vl_keys)
    for (i in seq_along(keys)) {
      if (is.na(m[i])) next
      rv <- suppressWarnings(as.numeric(vl_df$receiver_velocity_mean[m[i]]))
      velocity_recv_mean[i] <- rv
      ctx_velocity[i] <- isTRUE(is.finite(rv) && rv > 0 &&
                                  (!is.finite(med_recv) || rv >= med_recv))
    }
  }

  # 4. Perturbation-sensitive edge
  ctx_perturbation <- rep(NA, nrow(agg))
  pert_delta_score <- rep(NA_real_, nrow(agg))
  if (has_pert) {
    pt_df <- perturbation_ctx$delta_table
    pt_keys <- paste(pt_df$sender_node, pt_df$receiver_node, sep = "\r")
    m <- match(keys, pt_keys)
    for (i in seq_along(keys)) {
      if (is.na(m[i])) next
      pert_delta_score[i] <- suppressWarnings(as.numeric(pt_df$delta_score_total[m[i]]))
      ctx_perturbation[i] <- isTRUE(pt_df$affected[m[i]])
    }
  }

  mat <- cbind(ctx_spatial, ctx_trajectory, ctx_velocity, ctx_perturbation)
  n_eval <- apply(mat, 1L, function(x) sum(!is.na(x)))
  n_conv <- apply(mat, 1L, function(x) sum(x %in% TRUE))

  tier <- ifelse(
    n_eval == 0L, "Aucun contexte calculé",
    ifelse(n_conv >= 3L, "Exploratoire — fort (>=3)",
           ifelse(n_conv == 2L, "Exploratoire — modéré (2)",
                  ifelse(n_conv == 1L, "Exploratoire — partiel (1)",
                         "Aucun contexte convergent (0)")))
  )

  badge_names <- c("Spatial", "Trajectoire", "Vélocité", "Perturbation")
  badges <- apply(mat, 1L, function(x) {
    hit <- badge_names[which(x %in% TRUE)]
    if (length(hit)) paste(hit, collapse = " | ") else "—"
  })

  out <- data.frame(
    sender_node = agg$sender_node,
    receiver_node = agg$receiver_node,
    pair_label = pair_labels,
    n_interactions = as.integer(agg$n_interactions),
    mean_imported_score = as.numeric(agg$mean_score),
    ctx_spatial = ctx_spatial,
    ctx_trajectory = ctx_trajectory,
    ctx_velocity = ctx_velocity,
    ctx_perturbation = ctx_perturbation,
    spatial_frac = spatial_frac,
    traj_corr = traj_corr,
    velocity_recv_mean = velocity_recv_mean,
    pert_delta_score = pert_delta_score,
    n_contexts_evaluated = as.integer(n_eval),
    n_contexts_convergent = as.integer(n_conv),
    convergence_tier = tier,
    context_badges = badges,
    stringsAsFactors = FALSE
  )

  sc_ord <- ifelse(is.na(out$mean_imported_score), -Inf, out$mean_imported_score)
  ord <- order(-out$n_contexts_convergent, -out$n_contexts_evaluated, -sc_ord, out$pair_label)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  contexts_active <- c(
    spatial = has_spatial, trajectory = has_traj,
    velocity = has_vel, perturbation = has_pert
  )

  list(
    type = "ccc_multimodal_synthesis",
    label = epistemic_label,
    is_statistical_pvalue = FALSE,
    pair_table = out,
    qc = list(
      n_pairs = nrow(out),
      contexts_active = contexts_active,
      n_contexts_active = sum(contexts_active)
    )
  )
}

#' Matrice visuelle de convergence multi-modale CCC (filtrage exploratoire)
#'
#' @param synthesis Objet retourné par `build_communication_multimodal_synthesis()`.
#' @param top_n Nombre maximal de paires affichées (défaut 20).
#' @return objet ggplot.
#' @export
plot_communication_multimodal_synthesis <- function(synthesis, top_n = 20L) {
  if (!is.list(synthesis) ||
      !identical(synthesis$type %||% NULL, "ccc_multimodal_synthesis")) {
    .views_stop("plot_communication_multimodal_synthesis() : objet ccc_multimodal_synthesis requis.")
  }
  pt <- synthesis$pair_table
  if (!is.data.frame(pt) || nrow(pt) == 0L) {
    return(.communication_empty_view_plot(
      "Aucune paire disponible pour la synthèse multi-modale.",
      "Convergence des contextes (filtrage exploratoire)",
      "Aucune interaction filtrée"
    ))
  }
  n_show <- min(nrow(pt), max(1L, as.integer(top_n %||% 20L)))
  sub <- utils::head(pt, n_show)
  sub$pair_label <- factor(sub$pair_label, levels = rev(sub$pair_label))

  ctx_cols <- c("ctx_spatial", "ctx_trajectory", "ctx_velocity", "ctx_perturbation")
  ctx_labels <- c("Spatial", "Trajectoire", "Vélocité", "Perturbation")
  rows <- lapply(seq_along(ctx_cols), function(j) {
    val <- sub[[ctx_cols[j]]]
    status <- ifelse(is.na(val), "Non calculé",
                     ifelse(val %in% TRUE, "Convergent", "Non convergent"))
    data.frame(
      pair_label = sub$pair_label,
      context = ctx_labels[j],
      status = status,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  df$context <- factor(df$context, levels = ctx_labels)
  df$status <- factor(df$status, levels = c("Convergent", "Non convergent", "Non calculé"))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$context, y = .data$pair_label, fill = .data$status)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_manual(
      values = c("Convergent" = "#1B7837",
                 "Non convergent" = "#D9D9D9",
                 "Non calculé" = "#F5F5F5"),
      drop = FALSE,
      name = "Statut"
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "Convergence des contextes (filtrage exploratoire)",
      subtitle = sprintf(
        "Top %d paire(s) — %d/4 contexte(s) actif(s) — PAS une p-value statistique",
        n_show, as.integer(synthesis$qc$n_contexts_active %||% 0L)
      )
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold")
    )
}

#' Graphe comparatif des flux de signalisation par voie (CC-3)
#'
#' Visualise la table `flow_table` produite par `run_ccc_comparison()` :
#' voies les plus divergentes entre conditions (part relative `rel_<cond>` ou
#' poids brut `@netP$prob`), ordonnées par amplitude max-min décroissante.
#'
#' @param comparison_result Résultat retourné par `run_ccc_comparison()`.
#' @param top_n Nombre maximal de voies affichées (défaut 25).
#' @param metric `"relative"` (part du flux total) ou `"absolute"` (poids brut).
#' @return objet ggplot.
#' @export
plot_ccc_comparison_flow <- function(comparison_result,
                                     top_n = 25L,
                                     metric = c("relative", "absolute")) {
  metric <- match.arg(metric)
  if (!is.list(comparison_result) ||
      !is.data.frame(comparison_result$flow_table) ||
      !"pathway" %in% colnames(comparison_result$flow_table)) {
    .views_stop("plot_ccc_comparison_flow() : objet run_ccc_comparison()$flow_table requis.")
  }
  ft <- comparison_result$flow_table
  labels <- as.character(comparison_result$labels %||%
                           setdiff(colnames(ft), c("pathway", grep("^rel_", colnames(ft), value = TRUE))))
  if (length(labels) < 2L) {
    .views_stop("plot_ccc_comparison_flow() : au moins 2 conditions requises dans flow_table.")
  }
  val_cols <- if (identical(metric, "relative")) paste0("rel_", labels) else labels
  if (!all(val_cols %in% colnames(ft))) {
    .views_stop("plot_ccc_comparison_flow() : colonnes de métrique absentes de flow_table.")
  }
  if (nrow(ft) == 0L) {
    return(.communication_empty_view_plot(
      "Aucune voie commune à tracer entre les conditions.",
      "Comparaison des flux par voie (CC-3)",
      paste(labels, collapse = " vs ")
    ))
  }

  val_mat <- as.matrix(ft[, val_cols, drop = FALSE])
  spread <- apply(val_mat, 1L, function(x) {
    fin <- x[is.finite(x)]
    if (length(fin) < 2L) 0 else diff(range(fin))
  })
  ord <- order(-spread, ft$pathway)
  n_show <- min(nrow(ft), max(1L, as.integer(top_n %||% 25L)))
  sub <- ft[utils::head(ord, n_show), , drop = FALSE]
  sub$pathway <- factor(sub$pathway, levels = rev(sub$pathway))

  long_rows <- lapply(seq_along(labels), function(j) {
    data.frame(
      pathway = sub$pathway,
      condition = labels[j],
      value = as.numeric(sub[[val_cols[j]]]),
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, long_rows)
  df$condition <- factor(df$condition, levels = labels)

  x_lab <- if (identical(metric, "relative")) {
    "Part relative du flux d'information (@netP$prob)"
  } else {
    "Flux d'information cumulé (@netP$prob)"
  }
  sub_note <- sprintf(
    "%s — Top %d/%d voie(s) commune(s)%s",
    paste(labels, collapse = " vs "), n_show, nrow(ft),
    if (isTRUE(comparison_result$lifted)) " (groupes alignés via liftCellChat)" else ""
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$pathway, fill = .data$condition)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
    ggplot2::labs(
      x = x_lab, y = NULL, fill = "Condition",
      title = "Flux d'information par voie (comparaison inter-conditions)",
      subtitle = sub_note
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

# =============================================================================
# W2.2 §11.8 Gem #2 — Consensus multi-méthodes CCC par rangs intra-méthode
# =============================================================================

.consensus_stop <- function(message, state = "invalid_input") {
  stop(errorCondition(
    message,
    state = state,
    class = "communication_context_error"
  ))
}

.consensus_method_names <- function(results_list, validated_results) {
  raw_names <- names(results_list)
  if (is.null(raw_names)) raw_names <- rep("", length(validated_results))
  raw_names <- trimws(as.character(raw_names))
  raw_names[is.na(raw_names)] <- ""

  methods <- vapply(seq_along(validated_results), function(i) {
    if (nzchar(raw_names[i])) {
      raw_names[i]
    } else {
      sm <- as.character(validated_results[[i]]$source_method %||% "")
      if (length(sm) == 1L && !is.na(sm) && nzchar(trimws(sm))) {
        trimws(sm)
      } else {
        paste0("method_", i)
      }
    }
  }, character(1L))

  if (any(duplicated(methods))) {
    counts <- integer(0)
    for (i in seq_along(methods)) {
      m <- methods[i]
      k <- (counts[m] %||% 0L)
      if (is.na(k)) k <- 0L
      k <- k + 1L
      counts[m] <- k
    }
    seen <- integer(0)
    dup_set <- unique(methods[duplicated(methods)])
    for (i in seq_along(methods)) {
      m <- methods[i]
      if (m %in% dup_set) {
        idx <- (seen[m] %||% 0L)
        if (is.na(idx)) idx <- 0L
        idx <- idx + 1L
        seen[m] <- idx
        methods[i] <- paste0(m, "_", idx)
      }
    }
  }
  methods
}

.consensus_within_method_pct <- function(comm_res) {
  tbl <- .communication_node_keys(comm_res$canonical_table)
  if (nrow(tbl) == 0L) {
    return(data.frame(
      key = character(0),
      sender_node = character(0),
      receiver_node = character(0),
      ligand = character(0),
      receptor = character(0),
      within_pct = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  has_rank <- "rank" %in% colnames(tbl) &&
    any(is.finite(suppressWarnings(as.numeric(tbl$rank))))
  has_score <- "score" %in% colnames(tbl) &&
    any(is.finite(suppressWarnings(as.numeric(tbl$score))))

  if (has_rank) {
    val <- suppressWarnings(as.numeric(tbl$rank))
    keep <- is.finite(val)
    tbl <- tbl[keep, , drop = FALSE]
    val <- val[keep]
    # Direction : lower_is_better (rang 1 = meilleur)
    ord_metric <- val
  } else if (has_score) {
    val <- suppressWarnings(as.numeric(tbl$score))
    keep <- is.finite(val)
    tbl <- tbl[keep, , drop = FALSE]
    val <- val[keep]
    # Direction : higher_is_better (score élevé = meilleur -> tri sur -val)
    ord_metric <- -val
  } else {
    return(data.frame(
      key = character(0),
      sender_node = character(0),
      receiver_node = character(0),
      ligand = character(0),
      receptor = character(0),
      within_pct = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(tbl) == 0L) {
    return(data.frame(
      key = character(0),
      sender_node = character(0),
      receiver_node = character(0),
      ligand = character(0),
      receptor = character(0),
      within_pct = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  keys <- paste(
    as.character(tbl$sender_node),
    as.character(tbl$receiver_node),
    as.character(tbl$ligand),
    as.character(tbl$receptor),
    sep = "\r"
  )

  # Dédoublonnage par paire (sender_node, receiver_node, ligand, receptor) en
  # conservant la meilleure mesure intra-méthode avant calcul des rangs.
  if (any(duplicated(keys))) {
    ord_best <- order(ord_metric, seq_along(ord_metric))
    tbl <- tbl[ord_best, , drop = FALSE]
    ord_metric <- ord_metric[ord_best]
    keys <- keys[ord_best]
    keep_uniq <- !duplicated(keys)
    tbl <- tbl[keep_uniq, , drop = FALSE]
    ord_metric <- ord_metric[keep_uniq]
    keys <- keys[keep_uniq]
  }

  n_pairs <- nrow(tbl)
  # Rangs intra-méthode déterministes (ties.method = "min") dans (0, 1] :
  # meilleur = 1 / n_pairs, pire = 1.0.
  within_ranks <- rank(ord_metric, ties.method = "min")
  within_pct <- as.numeric(within_ranks) / as.numeric(n_pairs)

  data.frame(
    key = keys,
    sender_node = as.character(tbl$sender_node),
    receiver_node = as.character(tbl$receiver_node),
    ligand = as.character(tbl$ligand),
    receptor = as.character(tbl$receptor),
    within_pct = within_pct,
    stringsAsFactors = FALSE
  )
}

#' Consensus multi-méthodes CCC par rangs percentiles intra-méthode (W2.2 §11.8)
#'
#' Compare >= 2 résultats canoniques de communication intercellulaire (p. ex.
#' `CellChat × LIANA`) en calculant les rangs percentiles **intra-méthode** dans
#' `(0, 1]` (meilleur = plus petit percentile : outils par rang triés par
#' `rank` croissant, outils par score triés par `score` décroissant) et en
#' harmonisant les quadruplets `(sender_node, receiver_node, ligand, receptor)`
#' via `.communication_node_keys()`.
#'
#' **INVARIANT SCIENTIFIQUE** : les scores bruts ne sont JAMAIS moyennés entre
#' outils (échelles incompatibles : probabilité d'action de masse CellChat vs
#' rang agrégé RRA LIANA vs moyenne CellPhoneDB).
#'
#' @param results_list Liste (nommée ou non) d'au moins 2 résultats canoniques
#'   valides (`assert_communication_result()`).
#' @param top_pct Seuil de percentile intra-méthode dans `(0, 1]` définissant
#'   le top de chaque méthode (défaut `0.15`).
#' @param min_methods Nombre minimal de méthodes dans lesquelles une paire doit
#'   figurer dans le top `top_pct` pour être déclarée consensus (défaut `2L`,
#'   entier `>= 2L`).
#' @return Liste structurée avec `table`, `methods`, `params`, `qc` et
#'   `epistemic_note`.
#' @export
build_communication_consensus <- function(results_list,
                                          top_pct = 0.15,
                                          min_methods = 2L) {
  if (!is.list(results_list) || inherits(results_list, "data.frame") ||
      length(results_list) < 2L) {
    .consensus_stop(
      "build_communication_consensus() : 'results_list' doit être une liste d'au moins 2 résultats de communication."
    )
  }
  if (!is.numeric(top_pct) || length(top_pct) != 1L || !is.finite(top_pct) ||
      top_pct <= 0 || top_pct > 1) {
    .consensus_stop(
      "build_communication_consensus() : 'top_pct' doit être un numérique dans (0, 1]."
    )
  }
  if (!is.numeric(min_methods) || length(min_methods) != 1L ||
      !is.finite(min_methods) || as.integer(min_methods) < 2L) {
    .consensus_stop(
      "build_communication_consensus() : 'min_methods' doit être un entier >= 2."
    )
  }
  min_methods_int <- as.integer(min_methods)
  top_pct_num <- as.numeric(top_pct)

  validated <- lapply(seq_along(results_list), function(i) {
    tryCatch(
      assert_communication_result(
        results_list[[i]],
        context = sprintf("consensus communication (élément %d)", i)
      ),
      error = function(e) {
        .consensus_stop(
          sprintf(
            "build_communication_consensus() : élément %d invalide (%s).",
            i, conditionMessage(e)
          ),
          state = e$state %||% "invalid_input"
        )
      }
    )
  })

  methods <- .consensus_method_names(results_list, validated)
  n_methods <- length(methods)
  pct_dfs <- lapply(validated, .consensus_within_method_pct)
  epistemic_note <- "Consensus par rangs intra-méthode — scores bruts jamais moyennés entre outils"

  all_keys <- unique(unlist(lapply(pct_dfs, function(d) d$key), use.names = FALSE))
  if (length(all_keys) == 0L) {
    empty_tbl <- data.frame(
      sender_node = character(0),
      receiver_node = character(0),
      ligand = character(0),
      receptor = character(0),
      n_methods_present = integer(0),
      n_methods_top = integer(0),
      is_consensus = logical(0),
      mean_within_rank_pct = numeric(0),
      methods_present = character(0),
      stringsAsFactors = FALSE
    )
    for (m in methods) {
      empty_tbl[[paste0("pct_", m)]] <- numeric(0)
    }
    return(list(
      type = "ccc_consensus",
      table = empty_tbl,
      methods = methods,
      params = list(top_pct = top_pct_num, min_methods = min_methods_int),
      qc = list(n_total_pairs = 0L, n_consensus = 0L, n_methods = n_methods),
      epistemic_note = epistemic_note
    ))
  }

  key_parts <- do.call(rbind, strsplit(all_keys, "\r", fixed = TRUE))
  out <- data.frame(
    sender_node = key_parts[, 1L],
    receiver_node = key_parts[, 2L],
    ligand = key_parts[, 3L],
    receptor = key_parts[, 4L],
    stringsAsFactors = FALSE
  )

  pct_mat <- matrix(NA_real_, nrow = length(all_keys), ncol = n_methods)
  colnames(pct_mat) <- methods
  for (j in seq_along(methods)) {
    df_j <- pct_dfs[[j]]
    if (nrow(df_j) > 0L) {
      idx <- match(all_keys, df_j$key)
      pct_mat[, j] <- df_j$within_pct[idx]
    }
  }

  n_present <- as.integer(rowSums(!is.na(pct_mat)))
  n_top <- as.integer(rowSums(!is.na(pct_mat) & (pct_mat <= top_pct_num + 1e-12)))
  is_cons <- n_top >= min_methods_int
  mean_pct <- rowMeans(pct_mat, na.rm = TRUE)
  mean_pct[!is.finite(mean_pct)] <- NA_real_

  methods_present <- apply(!is.na(pct_mat), 1L, function(mask) {
    paste(methods[mask], collapse = ", ")
  })

  out$n_methods_present <- n_present
  out$n_methods_top <- n_top
  out$is_consensus <- is_cons
  out$mean_within_rank_pct <- mean_pct
  out$methods_present <- methods_present

  for (j in seq_along(methods)) {
    out[[paste0("pct_", methods[j])]] <- pct_mat[, j]
  }

  ord <- order(
    !out$is_consensus,
    -out$n_methods_top,
    out$mean_within_rank_pct,
    out$sender_node,
    out$receiver_node,
    out$ligand,
    out$receptor
  )
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  list(
    type = "ccc_consensus",
    table = out,
    methods = methods,
    params = list(top_pct = top_pct_num, min_methods = min_methods_int),
    qc = list(
      n_total_pairs = as.integer(nrow(out)),
      n_consensus = as.integer(sum(out$is_consensus)),
      n_methods = as.integer(n_methods)
    ),
    epistemic_note = epistemic_note
  )
}

#' Visualisation ggplot2 des interactions L-R consensus multi-méthodes (W2.2 §11.8)
#'
#' Affiche les meilleures interactions ligand-récepteur consensus ordonnées par
#' rang percentile intra-méthode moyen (`mean_within_rank_pct`) et nombre de
#' méthodes dans le top (`n_methods_top`). Si aucune paire n'atteint le seuil de
#' consensus, retourne un graphique d'état vide explicite portant la note
#' épistémique en sous-titre.
#'
#' @param consensus_result Résultat retourné par `build_communication_consensus()`.
#' @param top_n Nombre maximal d'interactions consensus affichées (défaut `25L`).
#' @return Objet `ggplot`.
#' @export
plot_communication_consensus <- function(consensus_result, top_n = 25L) {
  if (!is.list(consensus_result) ||
      !is.data.frame(consensus_result$table) ||
      !all(c("sender_node", "receiver_node", "ligand", "receptor",
             "n_methods_top", "is_consensus", "mean_within_rank_pct") %in%
           colnames(consensus_result$table))) {
    .consensus_stop(
      "plot_communication_consensus() : objet valide issu de build_communication_consensus() requis."
    )
  }

  epistemic_note <- as.character(
    consensus_result$epistemic_note %||%
      "Consensus par rangs intra-méthode — scores bruts jamais moyennés entre outils"
  )[1L]
  tbl <- consensus_result$table
  cons_rows <- tbl[tbl$is_consensus %in% TRUE, , drop = FALSE]

  if (nrow(cons_rows) == 0L) {
    return(.communication_empty_view_plot(
      "Aucune interaction consensus au seuil choisi.",
      "Consensus multi-méthodes CCC",
      epistemic_note
    ))
  }

  n_keep <- if (is.null(top_n) || !is.numeric(top_n) || length(top_n) != 1L ||
                !is.finite(top_n) || top_n <= 0) {
    nrow(cons_rows)
  } else {
    min(nrow(cons_rows), max(1L, as.integer(top_n)))
  }
  sub <- cons_rows[seq_len(n_keep), , drop = FALSE]
  sub$interaction_label <- sprintf(
    "%s -> %s : %s - %s",
    sub$sender_node, sub$receiver_node, sub$ligand, sub$receptor
  )
  # Déduplication d'étiquettes au cas où
  if (any(duplicated(sub$interaction_label))) {
    sub$interaction_label <- make.unique(sub$interaction_label, sep = " #")
  }
  sub$interaction_label <- factor(
    sub$interaction_label,
    levels = rev(sub$interaction_label)
  )
  sub$support_score <- 1 - sub$mean_within_rank_pct

  methods_str <- paste(as.character(consensus_result$methods %||% character(0)),
                       collapse = " \u00d7 ")
  top_pct_val <- as.numeric(consensus_result$params$top_pct %||% 0.15)
  min_m_val <- as.integer(consensus_result$params$min_methods %||% 2L)
  sub_full <- sprintf(
    "%s — %s (top %.0f%%, >= %d méthodes, %d/%d paires consensus)",
    epistemic_note,
    if (nzchar(methods_str)) methods_str else "Multi-méthodes",
    top_pct_val * 100,
    min_m_val,
    nrow(sub),
    nrow(cons_rows)
  )

  ggplot2::ggplot(
    sub,
    ggplot2::aes(
      x = .data$mean_within_rank_pct,
      y = .data$interaction_label,
      size = .data$n_methods_top,
      colour = .data$mean_within_rank_pct
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = .data$mean_within_rank_pct,
        y = .data$interaction_label,
        yend = .data$interaction_label
      ),
      linewidth = 0.6,
      colour = "grey70"
    ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_colour_gradient(
      low = "#084594",
      high = "#9ECAE1",
      name = "Rang pct moyen\n(plus bas = meilleur)"
    ) +
    ggplot2::scale_size_continuous(
      range = c(3.5, 7),
      breaks = sort(unique(sub$n_methods_top)),
      name = "Méthodes\n(top %)"
    ) +
    ggplot2::labs(
      x = "Rang percentile intra-méthode moyen (plus proche de 0 = meilleur)",
      y = NULL,
      title = "Consensus multi-méthodes CCC (interactions L-R)",
      subtitle = sub_full
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank()
    )
}

# =============================================================================
# W2.2 §11.8 Gem #4 : Vue alluviale 3 tiers (Sender -> Pathway/Ligand -> Receiver)
# =============================================================================

.alluvial_stack_tier <- function(agg, col_name, x_pos, gap_ratio = 0.04) {
  totals <- tapply(agg$weight, agg[[col_name]], sum)
  ord <- order(-as.numeric(totals), names(totals))
  nodes <- rev(names(totals)[ord])
  vals <- rev(as.numeric(totals)[ord])
  total_w <- sum(vals)
  gap <- if (length(nodes) > 1L) total_w * gap_ratio else 0

  ymins <- numeric(length(nodes))
  ymaxs <- numeric(length(nodes))
  cur <- 0
  for (i in seq_along(nodes)) {
    ymins[i] <- cur
    ymaxs[i] <- cur + vals[i]
    cur <- ymaxs[i] + gap
  }
  data.frame(
    tier = col_name,
    node = nodes,
    x = x_pos,
    ymin = ymins,
    ymax = ymaxs,
    ymid = (ymins + ymaxs) / 2,
    total = vals,
    stringsAsFactors = FALSE
  )
}

.alluvial_allocate_slices <- function(agg, tier_df, col_name, order_cols) {
  y0 <- numeric(nrow(agg))
  y1 <- numeric(nrow(agg))
  for (k in seq_len(nrow(tier_df))) {
    nd <- tier_df$node[k]
    idx <- which(agg[[col_name]] == nd)
    if (length(idx) > 1L) {
      ord_args <- lapply(order_cols, function(oc) agg[[oc]][idx])
      idx <- idx[do.call(order, ord_args)]
    }
    cur <- tier_df$ymin[k]
    for (j in idx) {
      y0[j] <- cur
      y1[j] <- cur + agg$weight[j]
      cur <- y1[j]
    }
  }
  list(y0 = y0, y1 = y1)
}

.alluvial_ribbon_poly <- function(poly_id, fill_group, x_start, x_end,
                                  y_start_bot, y_start_top,
                                  y_end_bot, y_end_top, n_pts = 25L) {
  t_seq <- seq(0, 1, length.out = n_pts)
  blend <- 3 * t_seq^2 - 2 * t_seq^3
  x_top <- x_start + (x_end - x_start) * t_seq
  y_top <- y_start_top + (y_end_top - y_start_top) * blend
  y_bot <- y_start_bot + (y_end_bot - y_start_bot) * blend
  data.frame(
    poly_id = poly_id,
    fill_group = fill_group,
    x = c(x_top, rev(x_top)),
    y = c(y_top, rev(y_bot)),
    stringsAsFactors = FALSE
  )
}

#' Vue alluviale 3 tiers Sender -> Pathway/Ligand -> Receiver (ggplot2 pur)
#'
#' Agrege les flux sur la table (filtree ou canonique) et trace un diagramme
#' alluvial pur ggplot2 sans dependance externe (rubans lisses par interpolation
#' cubique entre 3 colonnes). VUE DESCRIPTIVE — aucune lecture causale.
#'
#' @param communication_result Resultat canonique (assert_communication_result).
#' @param filtered_table Table filtree optionnelle ou NULL (table complete).
#' @param seurat_obj Objet Seurat courant optionnel (verification de peremption).
#' @param top_n_pathways Nombre maximal de noeuds intermediaires (pathway ou
#'   ligand) conserves par poids total decroissant (defaut `10L`).
#' @param middle_by Colonne du tier intermediaire : `"pathway"` (defaut) ou
#'   `"ligand"`.
#' @param top_n_flows Nombre maximal de chemins 3-tiers conserves par poids
#'   decroissant (defaut `30L`).
#' @return Objet `ggplot`.
#' @export
plot_communication_alluvial <- function(communication_result,
                                        filtered_table = NULL,
                                        seurat_obj = NULL,
                                        top_n_pathways = 10L,
                                        middle_by = c("pathway", "ligand"),
                                        top_n_flows = 30L) {
  r <- assert_communication_result(
    communication_result, seurat_obj = seurat_obj,
    context = "alluvial communication"
  )
  middle_by <- match.arg(middle_by)
  mid_col <- if (identical(middle_by, "pathway")) "pathway" else "ligand"
  mid_label <- if (identical(middle_by, "pathway")) "Pathway" else "Ligand"
  title_str <- sprintf("Flux Sender \u2192 %s \u2192 Receiver (resultats importes)", mid_label)
  subtitle_base <- .communication_score_caption(r)

  table <- .communication_node_keys(
    if (is.null(filtered_table)) r$canonical_table else filtered_table
  )

  if (nrow(table) == 0L) {
    return(.communication_empty_view_plot(
      "Aucune interaction apres filtres : elargissez les filtres ou reimportez.",
      title_str, subtitle_base
    ))
  }

  mid_vals <- as.character(table[[mid_col]])
  valid_nodes <- !is.na(table$sender_node) & nzchar(trimws(table$sender_node)) &
    !is.na(mid_vals) & nzchar(trimws(mid_vals)) &
    !is.na(table$receiver_node) & nzchar(trimws(table$receiver_node))
  table <- table[valid_nodes, , drop = FALSE]

  if (nrow(table) == 0L) {
    return(.communication_empty_view_plot(
      sprintf("Aucun flux exploitable : colonne %s ou noeuds vides/NA.", mid_col),
      title_str, subtitle_base
    ))
  }

  has_rank <- "rank" %in% colnames(table) &&
    any(is.finite(suppressWarnings(as.numeric(table$rank))))
  sc_num <- suppressWarnings(as.numeric(table$score))
  has_finite_score <- any(is.finite(sc_num))

  if (!has_rank && has_finite_score) {
    keep_pos <- is.finite(sc_num) & sc_num > 0
    table <- table[keep_pos, , drop = FALSE]
    row_weights <- sc_num[keep_pos]
  } else {
    row_weights <- rep(1, nrow(table))
  }

  if (nrow(table) == 0L) {
    return(.communication_empty_view_plot(
      "Aucun flux positif apres exclusion des scores nuls ou manquants.",
      title_str, subtitle_base
    ))
  }

  table$.mid_node <- as.character(table[[mid_col]])
  table$.row_weight <- row_weights

  keys <- paste(table$sender_node, table$.mid_node, table$receiver_node, sep = "\r")
  agg <- do.call(rbind, lapply(split(seq_len(nrow(table)), keys), function(idx) {
    data.frame(
      sender_node = table$sender_node[idx[1L]],
      mid_node = table$.mid_node[idx[1L]],
      receiver_node = table$receiver_node[idx[1L]],
      weight = sum(table$.row_weight[idx]),
      n_interactions = length(idx),
      stringsAsFactors = FALSE
    )
  }))
  agg <- agg[is.finite(agg$weight) & agg$weight > 0, , drop = FALSE]

  if (nrow(agg) == 0L) {
    return(.communication_empty_view_plot(
      "Aucun flux agrege positif disponible.",
      title_str, subtitle_base
    ))
  }

  if (!is.null(top_n_pathways) && is.numeric(top_n_pathways) &&
      length(top_n_pathways) == 1L && is.finite(top_n_pathways) &&
      top_n_pathways > 0) {
    mid_tot <- tapply(agg$weight, agg$mid_node, sum)
    n_mid_keep <- min(length(mid_tot), as.integer(top_n_pathways))
    keep_mid <- names(sort(mid_tot, decreasing = TRUE))[seq_len(n_mid_keep)]
    agg <- agg[agg$mid_node %in% keep_mid, , drop = FALSE]
  }

  agg <- .communication_top_n_rows(agg, top_n_flows, agg$weight, agg$n_interactions)
  if (nrow(agg) == 0L) {
    return(.communication_empty_view_plot(
      "Aucun flux apres application des limites top-N.",
      title_str, subtitle_base
    ))
  }

  t_sender <- .alluvial_stack_tier(agg, "sender_node", 1)
  t_mid <- .alluvial_stack_tier(agg, "mid_node", 2)
  t_receiver <- .alluvial_stack_tier(agg, "receiver_node", 3)

  agg$.s_ymid <- t_sender$ymid[match(agg$sender_node, t_sender$node)]
  agg$.m_ymid <- t_mid$ymid[match(agg$mid_node, t_mid$node)]
  agg$.r_ymid <- t_receiver$ymid[match(agg$receiver_node, t_receiver$node)]

  s_slices <- .alluvial_allocate_slices(agg, t_sender, "sender_node", c(".m_ymid", ".r_ymid"))
  m_slices <- .alluvial_allocate_slices(agg, t_mid, "mid_node", c(".s_ymid", ".r_ymid"))
  r_slices <- .alluvial_allocate_slices(agg, t_receiver, "receiver_node", c(".m_ymid", ".s_ymid"))

  node_w <- 0.18
  polys <- vector("list", nrow(agg) * 2L)
  for (i in seq_len(nrow(agg))) {
    polys[[2L * i - 1L]] <- .alluvial_ribbon_poly(
      poly_id = paste0("sm_", i),
      fill_group = agg$sender_node[i],
      x_start = 1 + node_w / 2,
      x_end = 2 - node_w / 2,
      y_start_bot = s_slices$y0[i],
      y_start_top = s_slices$y1[i],
      y_end_bot = m_slices$y0[i],
      y_end_top = m_slices$y1[i]
    )
    polys[[2L * i]] <- .alluvial_ribbon_poly(
      poly_id = paste0("mr_", i),
      fill_group = agg$sender_node[i],
      x_start = 2 + node_w / 2,
      x_end = 3 - node_w / 2,
      y_start_bot = m_slices$y0[i],
      y_start_top = m_slices$y1[i],
      y_end_bot = r_slices$y0[i],
      y_end_top = r_slices$y1[i]
    )
  }
  ribbon_df <- do.call(rbind, polys)

  node_df <- rbind(t_sender, t_mid, t_receiver)
  node_df$xmin <- node_df$x - node_w / 2
  node_df$xmax <- node_df$x + node_w / 2
  node_df$label_x <- ifelse(
    node_df$tier == "sender_node", node_df$xmin - 0.03,
    ifelse(node_df$tier == "receiver_node", node_df$xmax + 0.03, node_df$x)
  )
  node_df$hjust <- ifelse(
    node_df$tier == "sender_node", 1,
    ifelse(node_df$tier == "receiver_node", 0, 0.5)
  )

  max_y <- max(node_df$ymax, na.rm = TRUE)
  header_df <- data.frame(
    x = c(1, 2, 3),
    y = max_y * 1.04,
    label = c("Sender", mid_label, "Receiver"),
    stringsAsFactors = FALSE
  )

  weight_note <- if (has_rank || !has_finite_score) {
    "Epaisseur : nombre d'interactions (sans score scalaire)."
  } else {
    "Epaisseur : somme des scores importes."
  }

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = ribbon_df,
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$poly_id,
                   fill = .data$fill_group),
      alpha = 0.45,
      colour = NA
    ) +
    ggplot2::geom_rect(
      data = node_df,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax),
      fill = "#DCE6F1",
      colour = "#2C3E50",
      linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = node_df,
      ggplot2::aes(x = .data$label_x, y = .data$ymid, label = .data$node,
                   hjust = .data$hjust),
      size = 3.1,
      colour = "#1A1A1A"
    ) +
    ggplot2::geom_text(
      data = header_df,
      ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      fontface = "bold",
      size = 3.8,
      vjust = 0
    ) +
    ggplot2::coord_cartesian(xlim = c(0.45, 3.55), clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(10, 16, 10, 16)
    ) +
    ggplot2::labs(
      title = title_str,
      subtitle = sprintf(
        "%s — %d flux, %d %s(s). %s",
        subtitle_base, nrow(agg), nrow(t_mid), tolower(mid_label), weight_note
      )
    )
}

# =============================================================================
# W2.2 §11.8 Gem #3 : vérification descriptive d'expression L-R single-cell
# =============================================================================

.LR_EPISTEMIC_NOTE <- paste0(
  "V\u00e9rification descriptive d'expression L-R \u2014 la co-expression ",
  "ne prouve pas une interaction physique ni causale"
)

.lr_validate_nonempty_str <- function(val, arg_name) {
  if (!is.character(val) || length(val) != 1L || is.na(val) || !nzchar(trimws(val))) {
    .views_stop(sprintf(
      "build_communication_lr_expression() : '%s' doit etre une chaine non vide.",
      arg_name
    ))
  }
  trimws(val)
}

.lr_parse_subunits <- function(x, arg_name) {
  parts <- unlist(strsplit(as.character(x), "[_:+&,]"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  parts <- parts[!duplicated(parts)]
  if (length(parts) == 0L) {
    .views_stop(sprintf(
      "build_communication_lr_expression() : aucune sous-unite valide dans '%s'.",
      arg_name
    ))
  }
  parts
}

.lr_extract_seurat_data <- function(seurat_obj, assay, layer) {
  if (inherits(seurat_obj, "Seurat")) {
    meta <- seurat_obj@meta.data
    assays_avail <- names(seurat_obj@assays)
    if (length(assays_avail) == 0L) {
      .views_stop("build_communication_lr_expression() : aucun assay disponible dans l'objet Seurat.")
    }
    assay_use <- if (is.null(assay)) {
      da <- tryCatch(SeuratObject::DefaultAssay(seurat_obj), error = function(e) NULL)
      if (is.character(da) && length(da) == 1L && da %in% assays_avail) da else assays_avail[1L]
    } else {
      .lr_validate_nonempty_str(assay, "assay")
    }
    if (!assay_use %in% assays_avail) {
      .views_stop(sprintf(
        "build_communication_lr_expression() : assay '%s' introuvable (disponibles : %s).",
        assay_use, paste(assays_avail, collapse = ", ")
      ))
    }
    mat <- tryCatch(
      SeuratObject::LayerData(seurat_obj, assay = assay_use, layer = layer),
      error = function(e) NULL
    )
    if (is.null(mat)) {
      mat <- tryCatch(
        SeuratObject::GetAssayData(seurat_obj, assay = assay_use, layer = layer),
        error = function(e) NULL
      )
    }
    if (is.null(mat)) {
      mat <- tryCatch(
        SeuratObject::GetAssayData(seurat_obj, assay = assay_use, slot = layer),
        error = function(e) NULL
      )
    }
    if (is.null(mat) || (!is.matrix(mat) && !inherits(mat, "Matrix"))) {
      .views_stop(sprintf(
        "build_communication_lr_expression() : couche '%s' indisponible dans l'assay '%s'.",
        layer, assay_use
      ))
    }
    return(list(expr = mat, meta = meta, assay = assay_use))
  }

  if (is.list(seurat_obj) && !is.data.frame(seurat_obj)) {
    meta <- seurat_obj$meta.data %||% seurat_obj$metadata %||% seurat_obj$meta
    mat <- seurat_obj$expr %||% seurat_obj$expression %||%
      seurat_obj$matrix %||% seurat_obj[[layer]] %||%
      seurat_obj$data %||% seurat_obj$counts
    if (!is.data.frame(meta) || is.null(mat) ||
        (!is.matrix(mat) && !inherits(mat, "Matrix"))) {
      .views_stop(
        "build_communication_lr_expression() : liste compatible requise (matrice + meta.data)."
      )
    }
    assay_use <- if (is.null(assay)) NA_character_ else .lr_validate_nonempty_str(assay, "assay")
    return(list(expr = mat, meta = meta, assay = assay_use))
  }

  if (is.matrix(seurat_obj) || inherits(seurat_obj, "Matrix")) {
    meta <- attr(seurat_obj, "meta.data")
    if (!is.data.frame(meta)) {
      .views_stop(
        "build_communication_lr_expression() : attribut 'meta.data' requis sur une matrice."
      )
    }
    assay_use <- if (is.null(assay)) NA_character_ else .lr_validate_nonempty_str(assay, "assay")
    return(list(expr = seurat_obj, meta = meta, assay = assay_use))
  }

  .views_stop(
    "build_communication_lr_expression() : objet Seurat ou structure compatible requis."
  )
}

.lr_align_cells <- function(expr, meta) {
  if (!is.null(colnames(expr)) && !is.null(rownames(meta))) {
    common <- intersect(colnames(expr), rownames(meta))
    if (length(common) == 0L) {
      .views_stop(
        "build_communication_lr_expression() : aucune cellule commune entre expression et meta.data."
      )
    }
    return(list(
      expr = expr[, common, drop = FALSE],
      meta = meta[common, , drop = FALSE]
    ))
  }
  if (ncol(expr) == nrow(meta)) {
    return(list(expr = expr, meta = meta))
  }
  .views_stop(
    "build_communication_lr_expression() : dimensions incompatibles entre expression et meta.data."
  )
}

.lr_summarize_gene_pop <- function(expr, gene, role, population, cell_idx) {
  n_cells <- as.integer(length(cell_idx))
  gene_found <- gene %in% rownames(expr)
  if (!gene_found || n_cells == 0L) {
    return(data.frame(
      gene = gene,
      role = role,
      population = population,
      n_cells = n_cells,
      pct_expressing = NA_real_,
      mean_expr = NA_real_,
      mean_expr_nonzero = NA_real_,
      gene_found = gene_found,
      stringsAsFactors = FALSE
    ))
  }
  vals <- as.numeric(expr[gene, cell_idx])
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    return(data.frame(
      gene = gene,
      role = role,
      population = population,
      n_cells = n_cells,
      pct_expressing = NA_real_,
      mean_expr = NA_real_,
      mean_expr_nonzero = NA_real_,
      gene_found = TRUE,
      stringsAsFactors = FALSE
    ))
  }
  pos <- vals > 0
  data.frame(
    gene = gene,
    role = role,
    population = population,
    n_cells = n_cells,
    pct_expressing = 100 * sum(pos) / length(vals),
    mean_expr = mean(vals),
    mean_expr_nonzero = if (any(pos)) mean(vals[pos]) else NA_real_,
    gene_found = TRUE,
    stringsAsFactors = FALSE
  )
}

#' Verification descriptive d'expression single-cell Ligand-Recepteur
#'
#' Extrait la matrice d'expression pour les sous-unites du ligand et du
#' recepteur chez les cellules Sender et Receiver, sans jamais fabriquer de
#' valeurs pour les genes absents.
#'
#' @param seurat_obj Objet Seurat ou liste/matrice compatible avec `meta.data`.
#' @param sender Identite de la population emettrice.
#' @param receiver Identite de la population receptrice.
#' @param ligand Chaine ligand (mono-gene ou complexe ex. "ITGA4_ITGB1").
#' @param receptor Chaine recepteur (mono-gene ou complexe).
#' @param identity_column Colonne de `meta.data` identifiant les populations.
#' @param assay Assay Seurat optionnel (defaut = DefaultAssay).
#' @param layer Couche d'expression (defaut = "data").
#' @return Liste structuree : `table`, `qc`, `params`, `epistemic_note`.
#' @export
build_communication_lr_expression <- function(seurat_obj,
                                              sender,
                                              receiver,
                                              ligand,
                                              receptor,
                                              identity_column,
                                              assay = NULL,
                                              layer = "data") {
  sender_str <- .lr_validate_nonempty_str(sender, "sender")
  receiver_str <- .lr_validate_nonempty_str(receiver, "receiver")
  ligand_str <- .lr_validate_nonempty_str(ligand, "ligand")
  receptor_str <- .lr_validate_nonempty_str(receptor, "receptor")
  id_col <- .lr_validate_nonempty_str(identity_column, "identity_column")
  layer_str <- .lr_validate_nonempty_str(layer, "layer")

  extracted <- .lr_extract_seurat_data(seurat_obj, assay = assay, layer = layer_str)
  if (!id_col %in% colnames(extracted$meta)) {
    .views_stop(sprintf(
      "build_communication_lr_expression() : colonne '%s' absente de meta.data.",
      id_col
    ))
  }
  aligned <- .lr_align_cells(extracted$expr, extracted$meta)
  expr <- aligned$expr
  meta <- aligned$meta

  ligand_subunits <- .lr_parse_subunits(ligand_str, "ligand")
  receptor_subunits <- .lr_parse_subunits(receptor_str, "receptor")

  cell_ids <- as.character(meta[[id_col]])
  sender_idx <- which(!is.na(cell_ids) & cell_ids == sender_str)
  receiver_idx <- which(!is.na(cell_ids) & cell_ids == receiver_str)

  gene_specs <- rbind(
    data.frame(gene = ligand_subunits, role = "Ligand", stringsAsFactors = FALSE),
    data.frame(gene = receptor_subunits, role = "Receptor", stringsAsFactors = FALSE)
  )
  pop_specs <- list(
    list(name = "Sender", idx = sender_idx),
    list(name = "Receiver", idx = receiver_idx)
  )

  rows <- vector("list", nrow(gene_specs) * length(pop_specs))
  k <- 1L
  for (i in seq_len(nrow(gene_specs))) {
    for (p in pop_specs) {
      rows[[k]] <- .lr_summarize_gene_pop(
        expr = expr,
        gene = gene_specs$gene[i],
        role = gene_specs$role[i],
        population = p$name,
        cell_idx = p$idx
      )
      k <- k + 1L
    }
  }
  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL

  uniq_genes <- unique(c(ligand_subunits, receptor_subunits))
  found_mask <- uniq_genes %in% rownames(expr)
  rec_rows <- tbl[tbl$role == "Receptor" & tbl$population == "Receiver", , drop = FALSE]
  all_rec_expr <- isTRUE(
    nrow(rec_rows) > 0L &&
      all(rec_rows$gene_found %in% TRUE) &&
      all(!is.na(rec_rows$pct_expressing) & rec_rows$pct_expressing > 0)
  )

  list(
    table = tbl,
    qc = list(
      n_sender_cells = as.integer(length(sender_idx)),
      n_receiver_cells = as.integer(length(receiver_idx)),
      n_ligand_subunits = as.integer(length(ligand_subunits)),
      n_receptor_subunits = as.integer(length(receptor_subunits)),
      n_genes_found = as.integer(sum(found_mask)),
      n_genes_missing = as.integer(sum(!found_mask)),
      all_receptor_subunits_expressed = all_rec_expr
    ),
    params = list(
      sender = sender_str,
      receiver = receiver_str,
      ligand = ligand_str,
      receptor = receptor_str,
      ligand_subunits = ligand_subunits,
      receptor_subunits = receptor_subunits,
      identity_column = id_col,
      assay = extracted$assay,
      layer = layer_str,
      is_autocrine = identical(sender_str, receiver_str)
    ),
    epistemic_note = .LR_EPISTEMIC_NOTE
  )
}

#' DotPlot 2 lignes (Sender / Receiver) de verification d'expression L-R
#'
#' @param lr_expr_result Resultat produit par
#'   `build_communication_lr_expression()`.
#' @return Objet `ggplot2` (ou graphique d'etat vide explicite).
#' @export
plot_communication_lr_expression <- function(lr_expr_result) {
  req_cols <- c(
    "gene", "role", "population", "n_cells",
    "pct_expressing", "mean_expr", "gene_found"
  )
  if (!is.list(lr_expr_result) ||
      !is.data.frame(lr_expr_result$table) ||
      !all(req_cols %in% colnames(lr_expr_result$table))) {
    .views_stop(
      "plot_communication_lr_expression() : resultat valide issu de build_communication_lr_expression() requis."
    )
  }

  epistemic_note <- as.character(
    lr_expr_result$epistemic_note %||% .LR_EPISTEMIC_NOTE
  )[1L]
  sender_lbl <- as.character(lr_expr_result$params$sender %||% "Sender")[1L]
  receiver_lbl <- as.character(lr_expr_result$params$receiver %||% "Receiver")[1L]
  title_str <- sprintf("Expression L-R : %s -> %s", sender_lbl, receiver_lbl)

  tbl <- lr_expr_result$table
  if (nrow(tbl) == 0L ||
      !any(tbl$gene_found %in% TRUE) ||
      all(is.na(tbl$pct_expressing))) {
    return(.communication_empty_view_plot(
      "Aucun gene L-R detecte ou aucune cellule disponible.",
      title_str,
      epistemic_note
    ))
  }

  gene_levels <- unique(c(
    tbl$gene[tbl$role == "Ligand"],
    tbl$gene[tbl$role == "Receptor"]
  ))
  tbl$gene <- factor(tbl$gene, levels = gene_levels)
  tbl$role <- factor(tbl$role, levels = c("Ligand", "Receptor"))
  tbl$population <- factor(tbl$population, levels = c("Receiver", "Sender"))

  lig_str <- as.character(lr_expr_result$params$ligand %||% "")[1L]
  rec_str <- as.character(lr_expr_result$params$receptor %||% "")[1L]
  layer_str <- as.character(lr_expr_result$params$layer %||% "data")[1L]
  subtitle_str <- sprintf(
    "Ligand : %s | Recepteur : %s (couche : %s)",
    lig_str, rec_str, layer_str
  )

  ggplot2::ggplot(
    tbl,
    ggplot2::aes(x = .data$gene, y = .data$population)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$pct_expressing, colour = .data$mean_expr),
      na.rm = TRUE
    ) +
    ggplot2::facet_grid(. ~ role, scales = "free_x", space = "free_x") +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_size_continuous(
      range = c(1.5, 9),
      limits = c(0, 100),
      name = "% cellules > 0"
    ) +
    ggplot2::scale_colour_gradient(
      low = "#E0ECF4",
      high = "#08519C",
      name = "Expression moyenne",
      na.value = "grey85"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      plot.caption = ggplot2::element_text(size = 8.5, colour = "grey35")
    ) +
    ggplot2::labs(
      title = title_str,
      subtitle = subtitle_str,
      x = NULL,
      y = NULL,
      caption = epistemic_note
    )
}

# =============================================================================
# W2.2 §11.8 Gem #1 : comparaison différentielle Condition A vs Condition B
# =============================================================================

.DIFF_EPISTEMIC_NOTE <- paste0(
  "Comparaison différentielle descriptive intra-outil (delta_score et log2_fc ",
  "avec pseudocount) — aucune inférence causale ni comparaison inter-outils"
)

#' Valider un communication_result canonique pour la comparaison différentielle
.diff_validate_comm_result <- function(res, arg_name) {
  tryCatch(
    assert_communication_result(
      res,
      context = sprintf("comparaison différentielle (%s)", arg_name)
    ),
    error = function(e) {
      .views_stop(sprintf(
        "build_communication_condition_diff() : '%s' invalide (%s).",
        arg_name, conditionMessage(e)
      ))
    }
  )
  if (!is.data.frame(res$canonical_table)) {
    .views_stop(sprintf(
      "build_communication_condition_diff() : '%s$canonical_table' doit être un data.frame.",
      arg_name
    ))
  }
  req_cols <- c("sender", "receiver", "ligand", "receptor", "score")
  missing_cols <- setdiff(req_cols, colnames(res$canonical_table))
  if (length(missing_cols) > 0L) {
    .views_stop(sprintf(
      "build_communication_condition_diff() : colonnes absentes dans '%s$canonical_table' : %s.",
      arg_name, paste(missing_cols, collapse = ", ")
    ))
  }
  src <- res$source_method %||% res$source %||%
    unique(as.character(res$canonical_table$source_method %||% NA_character_))[1L]
  src <- as.character(src)[1L]
  if (is.na(src) || !nzchar(trimws(src))) {
    .views_stop(sprintf(
      "build_communication_condition_diff() : '%s' ne porte aucune source_method identifiable.",
      arg_name
    ))
  }
  list(result = res, source_method = trimws(src))
}

#' Extraire et dédoublonner les paires (sender, receiver, ligand, receptor)
.diff_extract_pairs <- function(res) {
  tbl <- res$canonical_table
  if (nrow(tbl) == 0L) {
    return(data.frame(
      key = character(0),
      sender = character(0),
      receiver = character(0),
      ligand = character(0),
      receptor = character(0),
      pathway = character(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    ))
  }
  tbl <- .communication_node_keys(tbl)
  sender_vec <- as.character(tbl$sender_node)
  receiver_vec <- as.character(tbl$receiver_node)
  ligand_vec <- as.character(tbl$ligand)
  receptor_vec <- as.character(tbl$receptor)
  pathway_vec <- if ("pathway" %in% colnames(tbl)) {
    as.character(tbl$pathway)
  } else {
    rep(NA_character_, nrow(tbl))
  }
  score_vec <- suppressWarnings(as.numeric(tbl$score))
  score_vec[!is.finite(score_vec)] <- NA_real_

  keys <- paste(sender_vec, receiver_vec, ligand_vec, receptor_vec, sep = "\r")
  # Dédoublonnage déterministe : conserver le score maximal par paire L-R
  if (any(duplicated(keys))) {
    ord_val <- ifelse(is.na(score_vec), -Inf, score_vec)
    ord <- order(-ord_val, seq_along(ord_val))
    keys <- keys[ord]
    sender_vec <- sender_vec[ord]
    receiver_vec <- receiver_vec[ord]
    ligand_vec <- ligand_vec[ord]
    receptor_vec <- receptor_vec[ord]
    pathway_vec <- pathway_vec[ord]
    score_vec <- score_vec[ord]
    keep <- !duplicated(keys)
    keys <- keys[keep]
    sender_vec <- sender_vec[keep]
    receiver_vec <- receiver_vec[keep]
    ligand_vec <- ligand_vec[keep]
    receptor_vec <- receptor_vec[keep]
    pathway_vec <- pathway_vec[keep]
    score_vec <- score_vec[keep]
  }

  data.frame(
    key = keys,
    sender = sender_vec,
    receiver = receiver_vec,
    ligand = ligand_vec,
    receptor = receptor_vec,
    pathway = pathway_vec,
    score = score_vec,
    stringsAsFactors = FALSE
  )
}

#' Comparaison différentielle CCC entre deux conditions (W2.2 §11.8 Gem #1)
#'
#' Compare deux `communication_result` canoniques issus du MÊME outil
#' (`source_method`), joint les interactions sur `(sender, receiver, ligand,
#' receptor)`, et calcule `delta_score = score_b - score_a` ainsi que
#' `log2_fc = log2((score_b + pseudocount) / (score_a + pseudocount))`.
#' Refuse toute comparaison inter-outils via `.views_stop()`.
#'
#' @param res_a Résultat canonique `communication_result` pour la condition A.
#' @param res_b Résultat canonique `communication_result` pour la condition B.
#' @param label_a Libellé de la condition A (défaut `"Condition A"`).
#' @param label_b Libellé de la condition B (défaut `"Condition B"`).
#' @param pseudocount Pseudocount strictement positif pour `log2_fc`
#'   (défaut `1e-6`).
#' @return `data.frame` (classe `c("ccc_condition_diff", "data.frame")`) trié
#'   par `abs(delta_score)` décroissant, avec attributs `source_method`,
#'   `label_a`, `label_b`, `pseudocount`, `epistemic_note` et `qc`.
#' @export
build_communication_condition_diff <- function(res_a,
                                               res_b,
                                               label_a = "Condition A",
                                               label_b = "Condition B",
                                               pseudocount = 1e-6) {
  if (!is.character(label_a) || length(label_a) != 1L ||
      is.na(label_a) || !nzchar(trimws(label_a))) {
    .views_stop("build_communication_condition_diff() : 'label_a' doit être une chaîne non vide.")
  }
  if (!is.character(label_b) || length(label_b) != 1L ||
      is.na(label_b) || !nzchar(trimws(label_b))) {
    .views_stop("build_communication_condition_diff() : 'label_b' doit être une chaîne non vide.")
  }
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      !is.finite(pseudocount) || pseudocount <= 0) {
    .views_stop("build_communication_condition_diff() : 'pseudocount' doit être un numérique fini > 0.")
  }
  label_a <- trimws(label_a)
  label_b <- trimws(label_b)
  pc <- as.numeric(pseudocount)

  val_a <- .diff_validate_comm_result(res_a, "res_a")
  val_b <- .diff_validate_comm_result(res_b, "res_b")

  if (!identical(tolower(val_a$source_method), tolower(val_b$source_method))) {
    .views_stop(sprintf(
      paste0(
        "build_communication_condition_diff() : comparaison inter-outils refusée ",
        "('%s' vs '%s') — les scores bruts ne sont pas comparables entre outils."
      ),
      val_a$source_method, val_b$source_method
    ))
  }

  df_a <- .diff_extract_pairs(val_a$result)
  df_b <- .diff_extract_pairs(val_b$result)
  all_keys <- unique(c(df_a$key, df_b$key))

  if (length(all_keys) == 0L) {
    empty_df <- data.frame(
      sender = character(0),
      receiver = character(0),
      ligand = character(0),
      receptor = character(0),
      pathway = character(0),
      score_a = numeric(0),
      score_b = numeric(0),
      delta_score = numeric(0),
      log2_fc = numeric(0),
      presence = character(0),
      stringsAsFactors = FALSE
    )
    attr(empty_df, "source_method") <- val_a$source_method
    attr(empty_df, "label_a") <- label_a
    attr(empty_df, "label_b") <- label_b
    attr(empty_df, "pseudocount") <- pc
    attr(empty_df, "epistemic_note") <- .DIFF_EPISTEMIC_NOTE
    attr(empty_df, "qc") <- list(n_total = 0L, n_both = 0L, n_only_a = 0L, n_only_b = 0L)
    class(empty_df) <- c("ccc_condition_diff", "data.frame")
    return(empty_df)
  }

  idx_a <- match(all_keys, df_a$key)
  idx_b <- match(all_keys, df_b$key)
  in_a <- !is.na(idx_a)
  in_b <- !is.na(idx_b)

  sender <- ifelse(in_a, df_a$sender[idx_a], df_b$sender[idx_b])
  receiver <- ifelse(in_a, df_a$receiver[idx_a], df_b$receiver[idx_b])
  ligand <- ifelse(in_a, df_a$ligand[idx_a], df_b$ligand[idx_b])
  receptor <- ifelse(in_a, df_a$receptor[idx_a], df_b$receptor[idx_b])

  pw_a <- ifelse(in_a, df_a$pathway[idx_a], NA_character_)
  pw_b <- ifelse(in_b, df_b$pathway[idx_b], NA_character_)
  pathway <- ifelse(!is.na(pw_a) & nzchar(pw_a), pw_a, pw_b)

  has_finite_a <- any(is.finite(df_a$score))
  has_finite_b <- any(is.finite(df_b$score))

  score_a <- ifelse(in_a, df_a$score[idx_a], if (has_finite_a || has_finite_b) 0 else NA_real_)
  score_b <- ifelse(in_b, df_b$score[idx_b], if (has_finite_a || has_finite_b) 0 else NA_real_)

  delta_score <- score_b - score_a
  valid_fc <- is.finite(score_a) & is.finite(score_b) &
    (score_a + pc > 0) & (score_b + pc > 0)
  log2_fc <- rep(NA_real_, length(all_keys))
  log2_fc[valid_fc] <- log2((score_b[valid_fc] + pc) / (score_a[valid_fc] + pc))

  presence <- ifelse(in_a & in_b, "both", ifelse(in_a, "only_a", "only_b"))

  out <- data.frame(
    sender = sender,
    receiver = receiver,
    ligand = ligand,
    receptor = receptor,
    pathway = pathway,
    score_a = as.numeric(score_a),
    score_b = as.numeric(score_b),
    delta_score = as.numeric(delta_score),
    log2_fc = as.numeric(log2_fc),
    presence = presence,
    stringsAsFactors = FALSE
  )

  abs_delta <- ifelse(is.finite(out$delta_score), abs(out$delta_score), -Inf)
  ord <- order(-abs_delta, out$sender, out$receiver, out$ligand, out$receptor)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "source_method") <- val_a$source_method
  attr(out, "label_a") <- label_a
  attr(out, "label_b") <- label_b
  attr(out, "pseudocount") <- pc
  attr(out, "epistemic_note") <- .DIFF_EPISTEMIC_NOTE
  attr(out, "qc") <- list(
    n_total = as.integer(nrow(out)),
    n_both = as.integer(sum(out$presence == "both")),
    n_only_a = as.integer(sum(out$presence == "only_a")),
    n_only_b = as.integer(sum(out$presence == "only_b"))
  )
  class(out) <- c("ccc_condition_diff", "data.frame")
  out
}

#' Visualisation ggplot2 des interactions CCC différentielles (W2.2 §11.8 Gem #1)
#'
#' Affiche les `top_n` interactions présentant la plus forte variation absolue
#' (`abs(delta_score)`) sous forme de graphique divergent (`delta_score` en
#' abscisse, `log2_fc` en taille/annotation, direction colorée par condition).
#'
#' @param diff_res `data.frame` retourné par `build_communication_condition_diff()`.
#' @param top_n Nombre maximal d'interactions affichées (défaut `25L`).
#' @return Objet `ggplot`.
#' @export
plot_communication_condition_diff <- function(diff_res, top_n = 25L) {
  req_cols <- c(
    "sender", "receiver", "ligand", "receptor",
    "score_a", "score_b", "delta_score", "log2_fc"
  )
  if (!is.data.frame(diff_res) || !all(req_cols %in% colnames(diff_res))) {
    .views_stop(
      "plot_communication_condition_diff() : data.frame valide issu de build_communication_condition_diff() requis."
    )
  }
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1) {
    .views_stop("plot_communication_condition_diff() : 'top_n' doit être un entier >= 1.")
  }

  label_a <- as.character(attr(diff_res, "label_a") %||% "Condition A")[1L]
  label_b <- as.character(attr(diff_res, "label_b") %||% "Condition B")[1L]
  src_method <- as.character(attr(diff_res, "source_method") %||% "communication")[1L]
  epistemic_note <- as.character(attr(diff_res, "epistemic_note") %||% .DIFF_EPISTEMIC_NOTE)[1L]

  title_str <- sprintf("CCC différentielle : %s vs %s", label_b, label_a)
  sub_str <- sprintf(
    "Outil : %s | %s -> %s (delta_score = %s - %s)",
    src_method, label_a, label_b, label_b, label_a
  )

  valid_rows <- diff_res[is.finite(diff_res$delta_score), , drop = FALSE]
  if (nrow(valid_rows) == 0L) {
    return(.communication_empty_view_plot(
      "Aucune interaction différentielle avec score fini disponible.",
      title_str,
      epistemic_note
    ))
  }

  ord <- order(-abs(valid_rows$delta_score), valid_rows$sender, valid_rows$receiver)
  valid_rows <- valid_rows[ord, , drop = FALSE]
  n_keep <- min(nrow(valid_rows), max(1L, as.integer(top_n)))
  sub <- valid_rows[seq_len(n_keep), , drop = FALSE]

  sub$interaction_label <- sprintf(
    "%s -> %s : %s - %s",
    sub$sender, sub$receiver, sub$ligand, sub$receptor
  )
  if (any(duplicated(sub$interaction_label))) {
    sub$interaction_label <- make.unique(sub$interaction_label, sep = " #")
  }
  sub <- sub[order(sub$delta_score), , drop = FALSE]
  sub$interaction_label <- factor(sub$interaction_label, levels = sub$interaction_label)

  dir_up <- paste0("Enrichi : ", label_b)
  dir_down <- paste0("Enrichi : ", label_a)
  dir_eq <- "Inchangé"
  sub$direction <- ifelse(
    sub$delta_score > 0, dir_up,
    ifelse(sub$delta_score < 0, dir_down, dir_eq)
  )
  sub$abs_log2_fc <- ifelse(is.finite(sub$log2_fc), abs(sub$log2_fc), 0)

  pal <- stats::setNames(
    c("#B2182B", "#2166AC", "grey55"),
    c(dir_up, dir_down, dir_eq)
  )

  ggplot2::ggplot(
    sub,
    ggplot2::aes(
      x = .data$delta_score,
      y = .data$interaction_label,
      colour = .data$direction
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = .data$delta_score,
        y = .data$interaction_label,
        yend = .data$interaction_label
      ),
      linewidth = 0.7
    ) +
    ggplot2::geom_point(ggplot2::aes(size = .data$abs_log2_fc), alpha = 0.9) +
    ggplot2::scale_colour_manual(values = pal, name = "Direction") +
    ggplot2::scale_size_continuous(range = c(2.5, 6.5), name = "|log2 FC|") +
    ggplot2::labs(
      title = title_str,
      subtitle = sub_str,
      x = sprintf("delta_score (%s - %s)", label_b, label_a),
      y = NULL,
      caption = epistemic_note
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(size = 8.5, colour = "grey35")
    )
}




