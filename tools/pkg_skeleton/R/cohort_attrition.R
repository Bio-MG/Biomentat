# =============================================================================
# tools/pkg_skeleton/R/cohort_attrition.R — OMOP / Cohort Attrition Engine (BK-11)
# =============================================================================
# Pure domain logic (zero Shiny symbols). Builds monotonic cohort/QC attrition
# tables, OMOP-style attrition funnels, Mermaid flowchart exports, and
# ggplot2 attrition diagrams for Bulk, Single-Cell, Spatial, and clinical
# cohorts.
#
# Dependencies: base R, stats, optional ggplot2 (checked via requireNamespace).
# =============================================================================

#' Internal typed error constructor for cohort attrition
#'
#' @param msg Character scalar error message.
#' @return An error condition of class `"cohort_attrition_error"`.
#' @keywords internal
.cohort_attrition_stop <- function(msg) {
  errorCondition(msg, class = "cohort_attrition_error")
}

#' Build a standardized cohort/QC attrition table from sequential steps
#'
#' @param steps A `data.frame` (or list of named lists) with columns/fields:
#'   `step` (character label), `n_remaining` (non-negative integer/numeric,
#'   monotonically non-increasing), and optional `reason` (character exclusion
#'   reason for that step).
#' @param unit Character scalar describing the counted unit (e.g., `"patients"`,
#'   `"samples"`, `"cells"`). Default `"samples"`.
#' @return A `data.frame` with columns `step_index`, `step`, `unit`, `n_in`,
#'   `n_excluded`, `n_remaining`, `pct_excluded`, `pct_retained`,
#'   `pct_cumulative_retained`, and `reason`.
#' @export
build_attrition_table <- function(steps, unit = "samples") {
  if (is.null(unit) || !is.character(unit) || length(unit) != 1L || is.na(unit) || !nzchar(trimws(unit))) {
    stop(.cohort_attrition_stop("'unit' must be a non-empty character scalar."))
  }
  unit <- trimws(unit)

  if (is.list(steps) && !is.data.frame(steps)) {
    if (length(steps) == 0L) {
      stop(.cohort_attrition_stop("'steps' must contain at least 1 step."))
    }
    steps <- do.call(rbind, lapply(steps, function(x) {
      if (!is.list(x) || is.null(x$step) || is.null(x$n_remaining)) {
        stop(.cohort_attrition_stop("Each element of 'steps' must have 'step' and 'n_remaining'."))
      }
      data.frame(
        step = as.character(x$step)[1L],
        n_remaining = as.numeric(x$n_remaining)[1L],
        reason = if (!is.null(x$reason)) as.character(x$reason)[1L] else "",
        stringsAsFactors = FALSE
      )
    }))
  }

  if (!is.data.frame(steps) || nrow(steps) < 1L) {
    stop(.cohort_attrition_stop("'steps' must be a non-empty data.frame or list of steps."))
  }
  if (!all(c("step", "n_remaining") %in% colnames(steps))) {
    stop(.cohort_attrition_stop("'steps' must contain columns 'step' and 'n_remaining'."))
  }

  step_labels <- as.character(steps$step)
  if (anyNA(step_labels) || any(!nzchar(trimws(step_labels)))) {
    stop(.cohort_attrition_stop("All 'step' labels must be non-empty strings."))
  }
  n_rem <- suppressWarnings(as.numeric(steps$n_remaining))
  if (anyNA(n_rem) || any(!is.finite(n_rem)) || any(n_rem < 0)) {
    stop(.cohort_attrition_stop("All 'n_remaining' values must be finite non-negative numbers."))
  }
  n_rem <- as.integer(round(n_rem))

  k <- length(n_rem)
  if (k > 1L && any(diff(n_rem) > 0L)) {
    stop(.cohort_attrition_stop("'n_remaining' must be monotonically non-increasing across steps."))
  }

  reasons <- if ("reason" %in% colnames(steps)) {
    r <- as.character(steps$reason)
    r[is.na(r)] <- ""
    r
  } else {
    rep("", k)
  }

  n_in <- c(n_rem[1L], if (k > 1L) n_rem[seq_len(k - 1L)] else integer(0L))
  n_excl <- n_in - n_rem
  n_start <- n_rem[1L]

  pct_excl <- if (k > 0L) ifelse(n_in > 0, round(100 * n_excl / n_in, 2), 0) else numeric(0L)
  pct_ret  <- if (k > 0L) ifelse(n_in > 0, round(100 * n_rem / n_in, 2), 100) else numeric(0L)
  pct_cum  <- if (n_start > 0L) round(100 * n_rem / n_start, 2) else rep(0, k)

  out <- data.frame(
    step_index = seq_len(k),
    step = trimws(step_labels),
    unit = rep(unit, k),
    n_in = as.integer(n_in),
    n_excluded = as.integer(n_excl),
    n_remaining = as.integer(n_rem),
    pct_excluded = pct_excl,
    pct_retained = pct_ret,
    pct_cumulative_retained = pct_cum,
    reason = reasons,
    stringsAsFactors = FALSE
  )
  class(out) <- c("biomentat_attrition", "data.frame")
  out
}

#' Build an OMOP-style clinical/transcriptomic cohort attrition table
#'
#' @param n_source Total records in source population / OMOP CDM.
#' @param n_phenotype Records matching index condition/phenotype.
#' @param n_observation Records passing observation window / clinical metadata.
#' @param n_qc Records passing RNA-seq / molecular QC thresholds.
#' @param n_final Final analytical cohort size (defaults to `n_qc`).
#' @param unit Unit label (default `"patients"`).
#' @return A `biomentat_attrition` data.frame.
#' @export
attrition_from_omop_cohort <- function(n_source,
                                       n_phenotype,
                                       n_observation,
                                       n_qc,
                                       n_final = n_qc,
                                       unit = "patients") {
  steps <- data.frame(
    step = c(
      "Source Population",
      "Index Phenotype / Condition",
      "Observation Period & Metadata",
      "Molecular QC & Completeness",
      "Final Analytical Cohort"
    ),
    n_remaining = c(n_source, n_phenotype, n_observation, n_qc, n_final),
    reason = c(
      "Initial cohort",
      "Missing index phenotype",
      "Incomplete observation / metadata",
      "Failed QC thresholds",
      "Excluded before modeling"
    ),
    stringsAsFactors = FALSE
  )
  build_attrition_table(steps, unit = unit)
}

#' Build a Single-Cell QC attrition table
#'
#' @param n_raw_cells Raw cells loaded.
#' @param n_post_qc_cells Cells remaining after mito/feature/count QC.
#' @param n_post_doublet_cells Cells remaining after doublet removal (defaults
#'   to `n_post_qc_cells`).
#' @param n_final_cells Cells retained in final clustered/annotated object
#'   (defaults to `n_post_doublet_cells`).
#' @return A `biomentat_attrition` data.frame.
#' @export
attrition_from_sc_qc <- function(n_raw_cells,
                                 n_post_qc_cells,
                                 n_post_doublet_cells = n_post_qc_cells,
                                 n_final_cells = n_post_doublet_cells) {
  steps <- data.frame(
    step = c(
      "Raw Cells",
      "Post-QC Filtering",
      "Doublet Removal",
      "Final Analyzed Cells"
    ),
    n_remaining = c(n_raw_cells, n_post_qc_cells, n_post_doublet_cells, n_final_cells),
    reason = c(
      "Initial matrix",
      "Low counts / high %MT / gene bounds",
      "Predicted doublets",
      "Unassigned / filtered"
    ),
    stringsAsFactors = FALSE
  )
  build_attrition_table(steps, unit = "cells")
}

#' Format a cohort attrition table as a Mermaid flowchart
#'
#' @param attrition_df Data.frame returned by `build_attrition_table()`.
#' @return Character scalar containing a Mermaid `flowchart TD` diagram.
#' @export
format_attrition_mermaid <- function(attrition_df) {
  if (!is.data.frame(attrition_df) || nrow(attrition_df) < 1L ||
      !all(c("step", "n_remaining", "n_excluded", "pct_cumulative_retained") %in% colnames(attrition_df))) {
    stop(.cohort_attrition_stop("'attrition_df' must be a valid table from build_attrition_table()."))
  }

  lines <- c("flowchart TD")
  k <- nrow(attrition_df)
  for (i in seq_len(k)) {
    node_id <- paste0("S", i)
    lbl <- sprintf("%s<br/>n = %d (%.1f%%)",
                   attrition_df$step[i],
                   attrition_df$n_remaining[i],
                   attrition_df$pct_cumulative_retained[i])
    lines <- c(lines, sprintf('  %s["%s"]', node_id, gsub('"', "'", lbl, fixed = TRUE)))
    if (i > 1L) {
      prev_id <- paste0("S", i - 1L)
      excl_n <- attrition_df$n_excluded[i]
      if (excl_n > 0L) {
        excl_id <- paste0("E", i)
        rsn <- if (nzchar(attrition_df$reason[i])) attrition_df$reason[i] else "Excluded"
        excl_lbl <- sprintf("-%d %s (%s)", excl_n, attrition_df$unit[i], rsn)
        lines <- c(
          lines,
          sprintf('  %s --> %s', prev_id, node_id),
          sprintf('  %s -.-> %s["%s"]', prev_id, excl_id, gsub('"', "'", excl_lbl, fixed = TRUE))
        )
      } else {
        lines <- c(lines, sprintf("  %s --> %s", prev_id, node_id))
      }
    }
  }
  paste(lines, collapse = "\n")
}

#' Plot a cohort attrition funnel diagram (ggplot2)
#'
#' @param attrition_df Data.frame returned by `build_attrition_table()`.
#' @param title Character plot title.
#' @param fill_color Hex color for retained bars (default `"#2C7BB6"`).
#' @return A `ggplot` object.
#' @export
plot_attrition_diagram <- function(attrition_df,
                                   title = "Cohort Attrition Diagram",
                                   fill_color = "#2C7BB6") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(.cohort_attrition_stop("Package 'ggplot2' is required for plot_attrition_diagram()."))
  }
  if (!is.data.frame(attrition_df) || nrow(attrition_df) < 1L ||
      !all(c("step", "n_remaining", "n_excluded", "pct_cumulative_retained") %in% colnames(attrition_df))) {
    stop(.cohort_attrition_stop("'attrition_df' must be a valid table from build_attrition_table()."))
  }

  df <- attrition_df
  df$step_f <- factor(df$step, levels = rev(unique(df$step)))
  df$bar_label <- ifelse(
    df$n_excluded > 0L,
    sprintf("%d (%.1f%%) | -%d", df$n_remaining, df$pct_cumulative_retained, df$n_excluded),
    sprintf("%d (%.1f%%)", df$n_remaining, df$pct_cumulative_retained)
  )

  x_max <- max(df$n_remaining, 1) * 1.28

  ggplot2::ggplot(df, ggplot2::aes(x = .data$n_remaining, y = .data$step_f)) +
    ggplot2::geom_col(fill = fill_color, width = 0.62, alpha = 0.9) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$bar_label),
      hjust = -0.05,
      size = 3.5
    ) +
    ggplot2::scale_x_continuous(limits = c(0, x_max), expand = c(0, 0)) +
    ggplot2::labs(
      title = title,
      x = sprintf("Remaining (%s)", df$unit[1L]),
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
}
