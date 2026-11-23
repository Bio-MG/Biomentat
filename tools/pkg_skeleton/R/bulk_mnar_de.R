# =============================================================================
# tools/pkg_skeleton/R/bulk_mnar_de.R — MNAR/MAR Missingness-Aware DE Engine
# =============================================================================
# Pure domain logic (zero Shiny symbols). Implements proteomics / bulk omics
# missingness regime classification (complete, MAR, MNAR), regime-aware
# imputation (Perseus left-censored normal draw for MNAR, row-mean for MAR),
# and PolySTest/MissTest-style combined intensity + missingness differential
# expression testing.
#
# Dependencies: limma (optional, checked via requireNamespace; falls back to
# Welch stats::t.test), stats.
# =============================================================================

#' Internal helper to construct a typed bulk_mnar_error condition
#'
#' @param msg Character scalar error message.
#' @return An error condition of class `"bulk_mnar_error"`.
#' @keywords internal
.bulk_mnar_stop <- function(msg) {
  errorCondition(msg, class = "bulk_mnar_error")
}

#' Validate expression/intensity matrix and 2-level group vector
#'
#' @param mat Numeric matrix (or numeric data.frame) with feature rownames and
#'   sample columns.
#' @param groups Factor or character vector of length `ncol(mat)` with 2 levels.
#' @return A list with `mat`, `groups_chr`, and `levels` (length-2 character).
#' @keywords internal
.validate_bulk_mnar_inputs <- function(mat, groups) {
  if (is.data.frame(mat)) {
    if (!all(vapply(mat, is.numeric, logical(1L)))) {
      stop(.bulk_mnar_stop("Input 'mat' data.frame must contain only numeric columns."))
    }
    mat <- as.matrix(mat)
  }
  if (is.null(mat) || !is.matrix(mat) || !is.numeric(mat)) {
    stop(.bulk_mnar_stop(sprintf(
      "Input 'mat' must be a numeric matrix or data.frame (got: %s).",
      if (is.null(mat)) "NULL" else paste(class(mat), collapse = "/")
    )))
  }
  if (nrow(mat) < 1L || ncol(mat) < 2L) {
    stop(.bulk_mnar_stop("Input 'mat' must have at least 1 feature (row) and 2 samples (columns)."))
  }
  feat_names <- rownames(mat)
  if (is.null(feat_names) || anyNA(feat_names) || any(!nzchar(feat_names))) {
    stop(.bulk_mnar_stop("Input 'mat' must have non-empty rownames corresponding to feature IDs."))
  }
  if (any(is.infinite(mat))) {
    stop(.bulk_mnar_stop("Input 'mat' must not contain infinite (Inf/-Inf) values."))
  }

  if (is.null(groups) || (!is.factor(groups) && !is.character(groups) && !is.numeric(groups))) {
    stop(.bulk_mnar_stop("Input 'groups' must be a factor or character vector of group labels."))
  }
  if (length(groups) != ncol(mat)) {
    stop(.bulk_mnar_stop(sprintf(
      "Length of 'groups' (%d) must equal ncol(mat) (%d).",
      length(groups), ncol(mat)
    )))
  }
  groups_chr <- as.character(groups)
  if (anyNA(groups_chr) || any(!nzchar(groups_chr))) {
    stop(.bulk_mnar_stop("Input 'groups' must not contain NA or empty group labels."))
  }

  lvls <- if (is.factor(groups)) {
    levels(droplevels(groups))
  } else {
    unique(groups_chr)
  }
  if (length(lvls) != 2L) {
    stop(.bulk_mnar_stop(sprintf(
      "Input 'groups' must have exactly 2 distinct levels (got %d).",
      length(lvls)
    )))
  }

  list(mat = mat, groups_chr = groups_chr, levels = lvls)
}

#' Classify per-feature missingness regime (complete, MNAR, or MAR)
#'
#' @param mat Numeric matrix (or coercible numeric data.frame) with non-empty
#'   rownames (features) and sample columns.
#' @param groups Factor or character vector of length `ncol(mat)` with 2 levels.
#' @param mnar_p_cutoff Numeric scalar in `(0, 1]`; Fisher's exact test p-value
#'   threshold below (or equal to) which group-skewed missingness is classified
#'   as `"MNAR"` (default `0.05`).
#' @return A `data.frame` with columns `feature`, `n_obs_g1`, `n_na_g1`,
#'   `pct_na_g1`, `n_obs_g2`, `n_na_g2`, `pct_na_g2`, `pvalue_missingness`,
#'   and `regime` (`"complete"`, `"MNAR"`, or `"MAR"`).
classify_missingness_regime <- function(mat, groups, mnar_p_cutoff = 0.05) {
  val <- .validate_bulk_mnar_inputs(mat, groups)
  mat <- val$mat
  groups_chr <- val$groups_chr
  lvls <- val$levels

  if (!is.numeric(mnar_p_cutoff) || length(mnar_p_cutoff) != 1L ||
      is.na(mnar_p_cutoff) || mnar_p_cutoff <= 0 || mnar_p_cutoff > 1) {
    stop(.bulk_mnar_stop("Parameter 'mnar_p_cutoff' must be a numeric scalar in (0, 1]."))
  }

  g1_idx <- which(groups_chr == lvls[1L])
  g2_idx <- which(groups_chr == lvls[2L])
  n_g1 <- length(g1_idx)
  n_g2 <- length(g2_idx)

  is_na_mat <- is.na(mat)
  n_na_g1 <- as.integer(rowSums(is_na_mat[, g1_idx, drop = FALSE]))
  n_na_g2 <- as.integer(rowSums(is_na_mat[, g2_idx, drop = FALSE]))
  n_obs_g1 <- n_g1 - n_na_g1
  n_obs_g2 <- n_g2 - n_na_g2

  pct_na_g1 <- n_na_g1 / n_g1
  pct_na_g2 <- n_na_g2 / n_g2
  pct_obs_g1 <- n_obs_g1 / n_g1
  pct_obs_g2 <- n_obs_g2 / n_g2

  n_feat <- nrow(mat)
  pval_miss <- rep(1.0, n_feat)
  regime <- rep("complete", n_feat)

  for (i in seq_len(n_feat)) {
    na_tot <- n_na_g1[i] + n_na_g2[i]
    if (na_tot == 0L) {
      pval_miss[i] <- 1.0
      regime[i] <- "complete"
    } else {
      obs_tot <- n_obs_g1[i] + n_obs_g2[i]
      if (obs_tot == 0L) {
        p_fish <- 1.0
      } else {
        tbl_2x2 <- matrix(
          c(n_obs_g1[i], n_na_g1[i], n_obs_g2[i], n_na_g2[i]),
          nrow = 2L,
          ncol = 2L
        )
        p_fish <- tryCatch(
          stats::fisher.test(tbl_2x2)$p.value,
          error = function(e) 1.0
        )
      }
      pval_miss[i] <- p_fish

      is_extreme_skew <- (pct_obs_g1[i] >= 0.5 && n_obs_g2[i] == 0L) ||
        (pct_obs_g2[i] >= 0.5 && n_obs_g1[i] == 0L)

      if ((!is.na(p_fish) && p_fish <= mnar_p_cutoff) || is_extreme_skew) {
        regime[i] <- "MNAR"
      } else {
        regime[i] <- "MAR"
      }
    }
  }

  data.frame(
    feature = rownames(mat),
    n_obs_g1 = n_obs_g1,
    n_na_g1 = n_na_g1,
    pct_na_g1 = pct_na_g1,
    n_obs_g2 = n_obs_g2,
    n_na_g2 = n_na_g2,
    pct_na_g2 = pct_na_g2,
    pvalue_missingness = pval_miss,
    regime = regime,
    row.names = rownames(mat),
    stringsAsFactors = FALSE
  )
}

#' Run PolySTest-style MNAR/MAR missingness-aware differential expression
#'
#' @param mat Numeric matrix (or coercible numeric data.frame) with non-empty
#'   rownames (features) and sample columns.
#' @param groups Factor or character vector of length `ncol(mat)` with 2 levels.
#' @param case_level Optional character scalar naming the case group level.
#' @param ctrl_level Optional character scalar naming the control group level.
#' @param mnar_shift_sd Numeric scalar >= 0; downward mean shift in SD units
#'   for Perseus-style MNAR left-censored normal imputation (default `1.8`).
#' @param mnar_width_sd Numeric scalar > 0; SD width multiplier for Perseus-style
#'   MNAR left-censored normal imputation (default `0.3`).
#' @param min_valid_total Integer scalar >= 1; minimum number of non-NA values
#'   across all samples required to retain a feature (default `2L`).
#' @param seed Optional integer seed for reproducible MNAR draws (default `42L`).
#'   Saves and restores `.Random.seed` on exit.
#' @param fdr_method Character scalar passed to `stats::p.adjust` (default `"BH"`).
#' @return A list with `table`, `n_features`, `n_mnar`, `n_mar`, `n_complete`,
#'   and `provenance`.
run_bulk_mnar_de <- function(mat,
                             groups,
                             case_level = NULL,
                             ctrl_level = NULL,
                             mnar_shift_sd = 1.8,
                             mnar_width_sd = 0.3,
                             min_valid_total = 2L,
                             seed = 42L,
                             fdr_method = "BH") {
  val <- .validate_bulk_mnar_inputs(mat, groups)
  mat <- val$mat
  groups_chr <- val$groups_chr
  lvls <- val$levels

  if (!is.numeric(mnar_shift_sd) || length(mnar_shift_sd) != 1L ||
      is.na(mnar_shift_sd) || mnar_shift_sd < 0) {
    stop(.bulk_mnar_stop("Parameter 'mnar_shift_sd' must be a non-negative numeric scalar."))
  }
  if (!is.numeric(mnar_width_sd) || length(mnar_width_sd) != 1L ||
      is.na(mnar_width_sd) || mnar_width_sd <= 0) {
    stop(.bulk_mnar_stop("Parameter 'mnar_width_sd' must be a positive numeric scalar."))
  }
  if (!is.numeric(min_valid_total) || length(min_valid_total) != 1L ||
      is.na(min_valid_total) || min_valid_total < 1) {
    stop(.bulk_mnar_stop("Parameter 'min_valid_total' must be a positive integer (>= 1)."))
  }
  min_valid_total <- as.integer(min_valid_total)

  if (!is.character(fdr_method) || length(fdr_method) != 1L ||
      !(fdr_method %in% stats::p.adjust.methods)) {
    stop(.bulk_mnar_stop(sprintf(
      "Parameter 'fdr_method' must be one of: %s.",
      paste(stats::p.adjust.methods, collapse = ", ")
    )))
  }

  # Resolve ctrl_level and case_level
  if (!is.null(case_level)) {
    case_level <- as.character(case_level)
    if (length(case_level) != 1L || !(case_level %in% lvls)) {
      stop(.bulk_mnar_stop(sprintf(
        "Specified 'case_level' ('%s') is not present in 'groups' levels (%s).",
        paste(case_level, collapse = ","), paste(lvls, collapse = ", ")
      )))
    }
  }
  if (!is.null(ctrl_level)) {
    ctrl_level <- as.character(ctrl_level)
    if (length(ctrl_level) != 1L || !(ctrl_level %in% lvls)) {
      stop(.bulk_mnar_stop(sprintf(
        "Specified 'ctrl_level' ('%s') is not present in 'groups' levels (%s).",
        paste(ctrl_level, collapse = ","), paste(lvls, collapse = ", ")
      )))
    }
  }

  if (is.null(ctrl_level) && is.null(case_level)) {
    ctrl_pat <- "^(ctrl|control|ref|reference|wt|untreated|normal|veh|vehicle)$"
    if (!is.factor(groups) && grepl(ctrl_pat, lvls[2L], ignore.case = TRUE) &&
        !grepl(ctrl_pat, lvls[1L], ignore.case = TRUE)) {
      ctrl_level <- lvls[2L]
      case_level <- lvls[1L]
    } else {
      ctrl_level <- lvls[1L]
      case_level <- lvls[2L]
    }
  } else if (is.null(ctrl_level)) {
    ctrl_level <- setdiff(lvls, case_level)[1L]
  } else if (is.null(case_level)) {
    case_level <- setdiff(lvls, ctrl_level)[1L]
  }

  if (identical(case_level, ctrl_level)) {
    stop(.bulk_mnar_stop("Parameters 'case_level' and 'ctrl_level' must be distinct group levels."))
  }

  # Save and restore .Random.seed on exit when seed is provided
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
      stop(.bulk_mnar_stop("Parameter 'seed' must be a single integer or NULL."))
    }
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

  # Filter features with < min_valid_total non-NA values
  n_valid <- rowSums(!is.na(mat))
  keep <- n_valid >= min_valid_total
  if (!any(keep)) {
    stop(.bulk_mnar_stop(sprintf(
      "No features have at least %d non-NA values across samples.",
      min_valid_total
    )))
  }
  mat_kept <- mat[keep, , drop = FALSE]

  # Order group factor as c(ctrl_level, case_level) so g1 = ctrl, g2 = case
  groups_fac <- factor(groups_chr, levels = c(ctrl_level, case_level))
  cls_df <- classify_missingness_regime(mat_kept, groups_fac, mnar_p_cutoff = 0.05)

  # Global fallback stats for MNAR imputation when a feature has < 2 observed values
  all_obs <- mat_kept[!is.na(mat_kept)]
  global_mean <- mean(all_obs)
  global_sd <- if (length(all_obs) >= 2L) stats::sd(all_obs) else 1.0
  if (is.na(global_sd) || global_sd <= 0) {
    global_sd <- 1.0
  }

  # Regime-specific imputation
  mat_imp <- mat_kept
  for (i in seq_len(nrow(mat_imp))) {
    reg_i <- cls_df$regime[i]
    if (reg_i == "complete") {
      next
    }
    na_idx <- which(is.na(mat_imp[i, ]))
    if (length(na_idx) == 0L) {
      next
    }
    obs_vals <- mat_imp[i, -na_idx]
    if (reg_i == "MAR") {
      mat_imp[i, na_idx] <- mean(obs_vals)
    } else if (reg_i == "MNAR") {
      mean_obs <- if (length(obs_vals) >= 1L) mean(obs_vals) else global_mean
      sd_obs <- if (length(obs_vals) >= 2L) stats::sd(obs_vals) else global_sd
      if (is.na(sd_obs) || sd_obs <= 0) {
        sd_obs <- global_sd
      }
      mat_imp[i, na_idx] <- stats::rnorm(
        length(na_idx),
        mean = mean_obs - mnar_shift_sd * sd_obs,
        sd = mnar_width_sd * sd_obs
      )
    }
  }

  ctrl_idx <- which(groups_chr == ctrl_level)
  case_idx <- which(groups_chr == case_level)

  mean_ctrl <- rowMeans(mat_imp[, ctrl_idx, drop = FALSE])
  mean_case <- rowMeans(mat_imp[, case_idx, drop = FALSE])
  log2fc <- as.numeric(mean_case - mean_ctrl)

  # Intensity p-values via limma (with Welch t.test fallback)
  pval_int <- rep(NA_real_, nrow(mat_imp))
  if (requireNamespace("limma", quietly = TRUE) && ncol(mat_imp) >= 3L) {
    design <- stats::model.matrix(~ groups_fac)
    fit <- tryCatch(
      limma::eBayes(limma::lmFit(mat_imp, design)),
      error = function(e) NULL
    )
    if (!is.null(fit) && !is.null(fit$p.value) && ncol(fit$p.value) >= 2L) {
      pval_int <- as.numeric(fit$p.value[, 2L])
    }
  }

  need_fallback <- which(is.na(pval_int))
  if (length(need_fallback) > 0L) {
    for (i in need_fallback) {
      pval_int[i] <- tryCatch(
        stats::t.test(mat_imp[i, case_idx], mat_imp[i, ctrl_idx])$p.value,
        error = function(e) 1.0
      )
      if (is.na(pval_int[i])) {
        pval_int[i] <- 1.0
      }
    }
  }

  # PolySTest / MissTest combination: Fisher's method (df = 4) for features with NAs
  pval_miss <- cls_df$pvalue_missingness
  pval_comb <- pval_int
  has_na <- cls_df$regime != "complete"
  if (any(has_na)) {
    p_i <- pmax(pval_int[has_na], .Machine$double.xmin)
    p_m <- pmax(pval_miss[has_na], .Machine$double.xmin)
    chi2 <- -2 * (log(p_i) + log(p_m))
    pval_comb[has_na] <- stats::pchisq(chi2, df = 4, lower.tail = FALSE)
  }

  fdr <- stats::p.adjust(pval_comb, method = fdr_method)

  res_tbl <- data.frame(
    feature = rownames(mat_imp),
    log2FC = log2fc,
    mean_case = as.numeric(mean_case),
    mean_ctrl = as.numeric(mean_ctrl),
    pct_na_ctrl = cls_df$pct_na_g1,
    pct_na_case = cls_df$pct_na_g2,
    pct_na_g1 = cls_df$pct_na_g1,
    pct_na_g2 = cls_df$pct_na_g2,
    regime = cls_df$regime,
    pvalue_intensity = pval_int,
    pvalue_missingness = pval_miss,
    pvalue_combined = pval_comb,
    FDR = fdr,
    row.names = rownames(mat_imp),
    stringsAsFactors = FALSE
  )

  ord <- order(res_tbl$FDR, res_tbl$pvalue_combined, -abs(res_tbl$log2FC))
  res_tbl <- res_tbl[ord, , drop = FALSE]

  list(
    table = res_tbl,
    n_features = as.integer(nrow(res_tbl)),
    n_mnar = as.integer(sum(res_tbl$regime == "MNAR")),
    n_mar = as.integer(sum(res_tbl$regime == "MAR")),
    n_complete = as.integer(sum(res_tbl$regime == "complete")),
    provenance = list(
      engine = "biomentat_mnar_de",
      case_level = case_level,
      ctrl_level = ctrl_level,
      mnar_shift_sd = mnar_shift_sd,
      mnar_width_sd = mnar_width_sd,
      min_valid_total = min_valid_total,
      fdr_method = fdr_method
    )
  )
}

#' Assign proteins to subcellular organelles from fractionation profiles (BK-14)
#'
#' @param profile_mat Numeric 2D matrix (proteins x fractions/samples) with
#'   non-empty `rownames(profile_mat)` and at least 2 columns.
#' @param markers Either a `data.frame` with columns `protein` and `organelle`,
#'   or a named character vector where names are protein IDs and values are
#'   organelle labels.
#' @param min_score Minimum cosine similarity (`0..1`) required to assign an
#'   organelle; profiles below `min_score` are labeled `"Unassigned"`.
#' @return Named list with `assignments`, `centroids`, `n_proteins`,
#'   `n_assigned`, and `n_organelles`.
#' @export
assign_subcellular_localization <- function(profile_mat, markers, min_score = 0.5) {
  if (is.null(profile_mat) || length(dim(profile_mat)) != 2L ||
      nrow(profile_mat) < 2L || ncol(profile_mat) < 2L || is.null(rownames(profile_mat))) {
    stop(.bulk_mnar_stop("'profile_mat' must be a 2D numeric matrix with rownames and >= 2 columns."))
  }
  mat <- as.matrix(profile_mat)
  storage.mode(mat) <- "double"
  mat[!is.finite(mat)] <- 0

  if (is.data.frame(markers)) {
    if (!all(c("protein", "organelle") %in% colnames(markers))) {
      stop(.bulk_mnar_stop("'markers' data.frame must contain columns 'protein' and 'organelle'."))
    }
    m_prot <- as.character(markers$protein)
    m_org  <- as.character(markers$organelle)
  } else if (is.character(markers) && !is.null(names(markers))) {
    m_prot <- names(markers)
    m_org  <- as.character(markers)
  } else {
    stop(.bulk_mnar_stop("'markers' must be a data.frame(protein, organelle) or named character vector."))
  }

  keep_m <- !is.na(m_prot) & nzchar(trimws(m_prot)) & !is.na(m_org) & nzchar(trimws(m_org))
  m_prot <- trimws(m_prot[keep_m])
  m_org  <- trimws(m_org[keep_m])
  in_mat <- m_prot %in% rownames(mat)
  m_prot <- m_prot[in_mat]
  m_org  <- m_org[in_mat]

  org_levels <- sort(unique(m_org))
  if (length(org_levels) < 2L) {
    stop(.bulk_mnar_stop("At least 2 distinct organelles with matched marker proteins are required."))
  }

  # L1 simplex normalization per protein
  row_sums <- rowSums(abs(mat))
  norm_mat <- mat
  pos <- row_sums > 0
  norm_mat[pos, ] <- mat[pos, , drop = FALSE] / row_sums[pos]
  norm_mat[!pos, ] <- 0

  # Organelle centroids
  centroids <- do.call(rbind, lapply(org_levels, function(org) {
    prots <- unique(m_prot[m_org == org])
    colMeans(norm_mat[prots, , drop = FALSE])
  }))
  rownames(centroids) <- org_levels
  colnames(centroids) <- colnames(mat)

  # Cosine similarity: proteins x organelles
  p_norm <- sqrt(rowSums(norm_mat^2))
  c_norm <- sqrt(rowSums(centroids^2))
  denom <- outer(p_norm, c_norm)
  sim <- norm_mat %*% t(centroids)
  sim <- ifelse(denom > 0, sim / denom, 0)

  top_idx <- max.col(sim, ties.method = "first")
  top_org <- org_levels[top_idx]
  top_sc  <- sim[cbind(seq_len(nrow(sim)), top_idx)]
  sec_sc  <- apply(sim, 1L, function(r) {
    sr <- sort(r, decreasing = TRUE)
    if (length(sr) >= 2L) sr[2L] else 0
  })

  marker_map <- stats::setNames(m_org, m_prot)
  is_m <- rownames(mat) %in% names(marker_map)
  assigned <- ifelse(top_sc >= as.numeric(min_score), top_org, "Unassigned")

  assign_df <- data.frame(
    protein = rownames(mat),
    top_organelle = top_org,
    score = round(as.numeric(top_sc), 4),
    second_score = round(as.numeric(sec_sc), 4),
    margin = round(as.numeric(top_sc - sec_sc), 4),
    assigned_organelle = assigned,
    is_marker = is_m,
    stringsAsFactors = FALSE
  )
  rownames(assign_df) <- rownames(mat)

  list(
    assignments = assign_df,
    centroids = centroids,
    norm_profiles = norm_mat,
    n_proteins = as.integer(nrow(assign_df)),
    n_assigned = as.integer(sum(assign_df$assigned_organelle != "Unassigned")),
    n_organelles = as.integer(length(org_levels))
  )
}

#' Detect subcellular translocation events between two conditions (BK-14)
#'
#' @param profile_ctrl Control condition protein-by-fraction matrix.
#' @param profile_treat Treated condition protein-by-fraction matrix.
#' @param markers Marker `data.frame(protein, organelle)` or named vector.
#' @param min_score Minimum cosine similarity threshold (default `0.5`).
#' @return Named list with `table`, `n_common`, and `n_translocated`.
#' @export
test_subcellular_translocation <- function(profile_ctrl, profile_treat, markers, min_score = 0.5) {
  common <- intersect(rownames(profile_ctrl), rownames(profile_treat))
  if (length(common) < 2L) {
    stop(.bulk_mnar_stop("At least 2 common proteins are required between control and treated profiles."))
  }
  res_c <- assign_subcellular_localization(profile_ctrl[common, , drop = FALSE], markers, min_score = min_score)
  res_t <- assign_subcellular_localization(profile_treat[common, , drop = FALSE], markers, min_score = min_score)

  ac <- res_c$assignments[common, , drop = FALSE]
  at <- res_t$assignments[common, , drop = FALSE]
  diff_mat <- res_t$norm_profiles[common, , drop = FALSE] - res_c$norm_profiles[common, , drop = FALSE]
  shift_dist <- sqrt(rowSums(diff_mat^2))
  transloc <- (ac$assigned_organelle != at$assigned_organelle) &
              (ac$assigned_organelle != "Unassigned") &
              (at$assigned_organelle != "Unassigned")

  tbl <- data.frame(
    protein = common,
    organelle_ctrl = ac$assigned_organelle,
    score_ctrl = ac$score,
    organelle_treat = at$assigned_organelle,
    score_treat = at$score,
    shift_distance = round(as.numeric(shift_dist), 4),
    translocated = transloc,
    stringsAsFactors = FALSE
  )
  tbl <- tbl[order(-as.integer(tbl$translocated), -tbl$shift_distance), , drop = FALSE]
  rownames(tbl) <- NULL

  list(
    table = tbl,
    n_common = as.integer(nrow(tbl)),
    n_translocated = as.integer(sum(tbl$translocated))
  )
}
