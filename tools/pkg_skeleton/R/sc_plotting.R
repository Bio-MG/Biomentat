# =============================================================================
# R/sc/sc_plotting.R — Single-Cell visualization builder
# =============================================================================
# Extracted from modules/sc/mod_sc_viz.R (Block 5 refactor).
# Pure plotting function: no Shiny reactivity, no module namespace.
# Called by mod_sc_viz_server() and sc_report_template.Rmd.
#
# Depends on: R/palettes.R (sc_discrete_scale, expression_continuous_scale,
#             bulk_diverging_ramp), Seurat, ggplot2, tidyr.
# Sourced in app.R AFTER R/palettes.R, BEFORE modules/sc/*.R.
# =============================================================================
#' Build a SC visualization ggplot from a config snapshot.
#'
#' Single source of truth used by BOTH the live render functions (which then
#' wrap with ggplotly or display static) AND the HTML/PDF report (via
#' params$saved_viz_list). Palette is applied via sc_discrete_scale() /
#' sc_continuous_scale() — both NULL-safe (NULL = Seurat/ggplot defaults,
#' which is what "default" palette produces).
#'
#' For the Volcano type, two extra attributes are attached to the returned
#' ggplot so the caller can construct a richer native plotly without re-running
#' FindMarkers:
#'   attr(p, "volcano_markers")  — the markers data.frame
#'   attr(p, "volcano_title")    — the plot title string
#'
#' @param obj     Seurat object.
#' @param cfg     Named list of input values (see .current_cfg() in server).
#' @param sc_palette Character palette name.
#' @param manual_colors Named character vector (level→hex), manual palette only.
#' @return ggplot or patchwork object.
build_sc_viz_plot <- function(obj, cfg, sc_palette = "default", manual_colors = NULL,
                              manual_gradient = NULL, manual_volcano_cols = NULL) {
  type    <- cfg$type     %||% "dim"
  pt_size <- as.numeric(cfg$pt_size %||% 0.5)
  grp     <- cfg$group_by %||% "seurat_clusters"

  # PLOT-S1 — shared resolver (R/plotting/theme.R). The default base_size (11)
  # is ggplot2's own default, i.e. EXACTLY what the bare theme_*() calls this
  # switch used before, so the rendering is unchanged. SC does not expose
  # base_size yet -> defensive fallback to the declared default.
  theme_fn <- ts_theme(cfg$plot_theme %||% TS_THEME_DEFAULT,
                       cfg$base_size %||% TS_BASE_SIZE_DEFAULT)

  pal_disc <- sc_discrete_scale(sc_palette, manual_colors, "color")
  pal_fill <- sc_discrete_scale(sc_palette, manual_colors, "fill")
  pal_cont <- sc_continuous_scale(sc_palette, "color", manual_gradient)

  # NULL-safe scale adder: suppressWarnings silences ggplot2's "replacing
  # existing scale" message (expected when overriding Seurat's built-in scale).
  .add <- function(p, s) if (is.null(s)) p else suppressWarnings(p + s)

  # Step-3.8B: past ~50k points, plain ggplot geom_point (what DimPlot/
  # FeaturePlot use by default) makes plotly conversion crawl or hang the
  # browser -- observed on the 1.3M-neurons dataset ("previews non
  # interactives / trop lentes" even though export PNG/PDF worked fine).
  # Seurat's raster=TRUE routes through scattermore (already a project
  # dependency, see global.R) -- point count no longer drives render cost.
  .raster_large <- ncol(obj) > 50000L

  # 1. DimPlot ----------------------------------------------------------------
  if (type == "dim") {
    red <- cfg$reduction %||% "umap"
    if (!red %in% names(obj@reductions)) stop(errorCondition(paste0("Réduction non calculée : ", red), class = "sc_plotting_error"))
    p <- DimPlot(obj, reduction = red, group.by = grp,
                 label = TRUE, pt.size = pt_size, raster = .raster_large) +
         theme_fn + ggtitle(paste(toupper(red), "\u2014", grp))
    return(.add(p, pal_disc))
  }

  # 2. FeaturePlot ------------------------------------------------------------
  if (type == "feature") {
    valid <- intersect(cfg$feat_sel %||% character(0), rownames(obj))
    if (!length(valid)) stop(errorCondition("Aucun gène valide", class = "sc_plotting_error"))
    # BK-04 opt-in 1st/99th percentile color clipping (robust to outlier cells)
    q_min <- if (isTRUE(cfg$clip_quantiles)) "q1" else NA
    q_max <- if (isTRUE(cfg$clip_quantiles)) "q99" else NA
    # multi-gene returns patchwork: scale only safe for single gene
    if (length(valid) == 1L) {
      p <- FeaturePlot(obj, features = valid[1], pt.size = pt_size,
                       min.cutoff = q_min, max.cutoff = q_max,
                       raster = .raster_large) + theme_fn
      return(.add(p, pal_cont))
    }
    return(FeaturePlot(obj, features = head(valid, 4), ncol = 2,
                       pt.size = pt_size, min.cutoff = q_min, max.cutoff = q_max,
                       raster = .raster_large) + theme_fn)
  }

  # 3. Scatter (enhanced) -----------------------------------------------------
  if (type == "scatter") {
    p <- plot_enhanced_scatter(
           obj, cfg$scatter_gene1, cfg$scatter_gene2, group.by = grp,
           method     = cfg$scatter_cor_method %||% "pearson",
           add_smooth = isTRUE(cfg$scatter_smooth),
           pt.size    = pt_size) + theme_fn
    return(.add(p, pal_disc))
  }

  # 4. Violin -----------------------------------------------------------------
  if (type == "violin") {
    valid <- intersect(cfg$feat_sel %||% character(0), rownames(obj))
    if (!length(valid)) stop(errorCondition("Aucun gène valide", class = "sc_plotting_error"))
    if (length(valid) == 1L) {
      p <- plot_violin_enhanced(obj, valid[1], grp,
                                add_boxplot = isTRUE(cfg$violin_boxplot)) + theme_fn
      return(.add(p, pal_fill))
    }
    plots <- lapply(head(valid, 4), function(g)
      .add(VlnPlot(obj, features = g, group.by = grp, pt.size = 0) +
           theme_fn + ggtitle(g), pal_fill))
    return(wrap_plots(plots, ncol = 2))
  }

  # 5. Stacked Violin ---------------------------------------------------------
  if (type == "stacked_violin") {
    valid <- intersect(cfg$feat_sel %||% character(0), rownames(obj))
    if (!length(valid)) stop(errorCondition("Aucun gène valide", class = "sc_plotting_error"))
    plots <- lapply(head(valid, 8), function(g)
      .add(VlnPlot(obj, features = g, group.by = grp, pt.size = 0) +
           theme(legend.position = "none",
                 axis.title.x   = element_blank(),
                 axis.text.x    = element_blank()) + ggtitle(g), pal_fill))
    return(wrap_plots(plots, ncol = 1) + theme_fn)
  }

  # 6. Ridge ------------------------------------------------------------------
  if (type == "ridge") {
    valid <- intersect(cfg$feat_sel %||% character(0), rownames(obj))
    if (!length(valid)) stop(errorCondition("Aucun gène valide", class = "sc_plotting_error"))
    p <- RidgePlot(obj, features = head(valid, 6), ncol = 2) + theme_fn
    return(.add(p, pal_fill))
  }

  # 7. DotPlot ----------------------------------------------------------------
  if (type == "dot") {
    valid <- intersect(cfg$feat_sel %||% character(0), rownames(obj))
    if (!length(valid)) stop(errorCondition("Aucun gène valide", class = "sc_plotting_error"))
    p <- DotPlot(obj, features = head(valid, 20), group.by = grp) + theme_fn +
         theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # DotPlot: color encodes avg expression (continuous), size encodes pct expressed.
    # Route through sc_continuous_scale() so i18n$t("Manuel") (sequential gradient) is
    # respected here too, same as FeaturePlot; size scale untouched.
    return(suppressWarnings(p + sc_continuous_scale(sc_palette, "color", manual_gradient)))
  }

  # 8. Heatmap ------------------------------------------------------------------
  # Step-3.7 FIX: DoHeatmap's group color bar is set via group.colors=, NOT via
  # a `+ scale_color_manual(...)` layer (unlike every other plot type here) —
  # so the selected palette (manual/Okabe-Ito/Set2/Viridis) used to be silently
  # ignored on the top annotation strip while the blue/red expression gradient
  # (unaffected, driven separately by DoHeatmap's internal fill scale) looked
  # correct. sc_discrete_colors() supplies the matching concrete color vector;
  # NULL (palette=="default") preserves the previous/original behaviour exactly.
  if (type == "heatmap") {
    valid <- intersect(cfg$feat_sel %||% character(0), rownames(obj))
    if (!length(valid)) stop(errorCondition("Aucun gène valide", class = "sc_plotting_error"))
    valid <- head(valid, 30)
    # Step-3.7A: smart_scale_data()/ScaleData() upstream now restricts scale.data
    # to VariableFeatures by default (RAM-safety) — marker genes picked here that
    # fall outside that set need scale.data extended on demand, without ever
    # re-scaling the whole assay.
    obj <- tryCatch(ensure_genes_scaled(obj, valid), error = function(e) obj)
    lvls     <- sort(unique(na.omit(as.character(obj@meta.data[[grp]]))))
    grp_cols <- tryCatch(sc_discrete_colors(lvls, sc_palette, manual_colors),
                         error = function(e) NULL)
    p <- if (!is.null(grp_cols))
      DoHeatmap(obj, features = valid, group.by = grp, group.colors = grp_cols)
    else
      DoHeatmap(obj, features = valid, group.by = grp)
    div_scale <- sc_diverging_scale(sc_palette, "fill", manual_gradient)
    if (!is.null(div_scale)) p <- suppressWarnings(p + div_scale)
    return(p + theme_fn)
  }

  # 9. Correlation Matrix -----------------------------------------------------
  if (type == "correlation_matrix") {
    valid <- intersect(cfg$feat_sel %||% character(0), rownames(obj))
    if (length(valid) < 2) stop(errorCondition("Au moins 2 gènes requis", class = "sc_plotting_error"))
    rc <- diverging_ramp_colors(sc_palette, manual_gradient)
    return(plot_correlation_matrix(obj, head(valid, 30), method = "pearson",
                                   low_color = rc$low, mid_color = rc$mid, high_color = rc$high) + theme_fn)
  }

  # 10. Multi-Sample ----------------------------------------------------------
  if (type == "multi_sample") {
    if (length(unique(obj$orig.ident)) < 2) stop(errorCondition("Au moins 2 échantillons requis", class = "sc_plotting_error"))
    p <- plot_multi_sample(obj, cfg$multi_gene %||% rownames(obj)[1],
                           cfg$multi_plot_type %||% "violin") + theme_fn
    return(.add(p, pal_fill))
  }

  # 11. Volcano (ggplot — also carries markers as attr for native plotly) ------
  if (type == "volcano") {
    grp_col <- grp
    if (!grp_col %in% colnames(obj@meta.data)) stop(errorCondition("Colonne de groupe introuvable", class = "sc_plotting_error"))
    Idents(obj) <- factor(as.character(obj@meta.data[[grp_col]]))
    ident1 <- cfg$volcano_group1 %||% levels(Idents(obj))[1]
    ident2 <- if (is.null(cfg$volcano_group2) ||
                  cfg$volcano_group2 == "rest") NULL else cfg$volcano_group2
    if (!ident1 %in% levels(Idents(obj))) stop(errorCondition(paste0("Groupe 1 invalide: ", ident1), class = "sc_plotting_error"))

    markers <- tryCatch(
      FindMarkers(obj, ident.1 = ident1, ident.2 = ident2,
                  only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0,
                  verbose = FALSE),
      error = function(e) stop(errorCondition(paste0("FindMarkers: ", e$message), class = "sc_plotting_error")))
    if (!nrow(markers)) stop(errorCondition("Aucun marqueur trouvé", class = "sc_plotting_error"))
    markers$gene <- rownames(markers)

    lfc  <- as.numeric(cfg$volcano_logfc %||% 0.25)
    pval <- as.numeric(cfg$volcano_pval  %||% 0.05)
    markers$status <- dplyr::case_when(
      markers$avg_log2FC >  lfc & markers$p_val_adj < pval ~ "Up",
      markers$avg_log2FC < -lfc & markers$p_val_adj < pval ~ "Down",
      TRUE ~ "NS")

    color_map <- bulk_role_colors(
      sc_palette,
      manual_colors = if (identical(sc_palette, "manual")) c(
        Up   = manual_volcano_cols[["up"]]   %||% "#E74C3C",
        Down = manual_volcano_cols[["down"]] %||% "#2980B9",
        NS   = manual_volcano_cols[["ns"]]   %||% "#BDC3C7"
      ) else NULL
    )
    show_lbl  <- isTRUE(cfg$volcano_show_labels %||% TRUE)
    top_genes <- c(head(markers$gene[markers$status == "Up"],   10),
                   head(markers$gene[markers$status == "Down"], 10))
    markers$label <- ifelse(show_lbl & markers$gene %in% top_genes, markers$gene, "")

    vtitle <- paste0("Volcano: ", ident1, " vs ", ident2 %||% "rest")
    p <- ggplot(markers, aes(x = avg_log2FC,
                              y = -log10(p_val_adj + 1e-300),
                              color = status)) +
         geom_point(alpha = 0.75, size = pt_size) +
         scale_color_manual(values = color_map) +
         geom_vline(xintercept = c(-lfc, lfc), linetype = "dotted", color = "grey40") +
         geom_hline(yintercept = -log10(pval),  linetype = "dotted", color = "grey40") +
         labs(title = vtitle, x = "Log2 Fold Change",
              y = "-log10(P.adj)", color = i18n$t("Statut")) + theme_fn
    if (show_lbl)
      # PLOT-Q1 — repel marker labels (same rationale as plot_volcano_bulk);
      # graceful fallback if ggrepel is unavailable.
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(aes(label = label), size = 2.8, na.rm = TRUE,
                                          max.overlaps = 20, box.padding = 0.4,
                                          segment.size = 0.3, show.legend = FALSE,
                                          min.segment.length = 0)
      } else {
        p <- p + geom_text(aes(label = label), size = 2.8, na.rm = TRUE,
                           vjust = -0.6, show.legend = FALSE)
      }

    # Attach markers + title as attributes so the caller can build native plotly
    # without re-running FindMarkers.
    attr(p, "volcano_markers") <- markers
    attr(p, "volcano_title")   <- vtitle
    return(p)
  }

  # 12. Heatmap Hierarchique -- returns a DRAWN ComplexHeatmap, NOT a ggplot.
  # Caller MUST route via the static render path, never ggsave()/ggplotly().
  if (type == "heatmap_hier") {
    valid <- intersect(cfg$feat_sel %||% character(0), rownames(obj))
    if (!length(valid)) stop(errorCondition("Aucun gène valide sélectionné pour la heatmap", class = "sc_plotting_error"))
    # PLOT-S4 — reglages de clustering. k < 2 (ou vide) = pas de decoupage,
    # ce qui reproduit exactement le rendu anterieur au jalon.
    k_row <- suppressWarnings(as.integer(cfg$hier_k_row %||% 0))
    if (length(k_row) != 1L || is.na(k_row) || k_row < 2L) k_row <- NULL
    return(build_sc_hierarchical_heatmap(
      obj, features = valid, group_by = grp, max_features = 50L,
      max_cells = as.integer(cfg$hier_max_cells %||% 5000L),
      palette = sc_palette, manual_colors = manual_colors, manual_gradient = manual_gradient,
      clustering_distance = cfg$hier_clust_distance %||% "euclidean",
      clustering_method   = cfg$hier_clust_method %||% "complete",
      k_row               = k_row
    ))
  }

  # 13. Densite d'Expression 2D (Nebulosa-like) --------------------------------
  if (type == "density_2d") {
    gene <- cfg$density_gene
    if (is.null(gene) || !nzchar(gene %||% "")) stop(errorCondition("Aucun gène sélectionné pour la densité", class = "sc_plotting_error"))
    red <- cfg$density_reduction %||% "umap"
    if (!red %in% names(obj@reductions)) stop(errorCondition(paste0("Réduction non calculée : ", red), class = "sc_plotting_error"))
    return(plot_sc_expression_density_2d(obj, feature = gene, reduction = red,
                                         max_cells = as.integer(cfg$density_max_cells %||% 50000L),
                                         palette = sc_palette, manual_gradient = manual_gradient) + theme_fn)
  }

  stop(errorCondition(paste0("Type de visualisation non supporté: ", type), class = "sc_plotting_error"))
}

# =============================================================================
# Adaptive UMAP Render Budget & Stratified Subsampling (W1.1 #13)
# =============================================================================

.validate_sc_budget_int <- function(val, name, min_val = 0L) {
  if (length(val) != 1L || !is.numeric(val) || is.na(val) || !is.finite(val) || val < min_val) {
    stop(errorCondition(paste0(name, " doit être un entier >= ", min_val),
                        class = "sc_plotting_error"))
  }
  as.integer(val)
}

#' Compute adaptive UMAP rendering budget for large single-cell datasets.
#'
#' Pure helper deciding whether interactive rendering should subsample
#' (`should_subsample`) or rasterize (`use_raster`), along with the target
#' point count and subsampling fraction.
#'
#' @param n_cells          Non-negative integer number of cells.
#' @param max_interactive  Maximum point count for interactive vector rendering
#'                         before subsampling is triggered (default 25000L).
#' @param raster_threshold Cell count above which rasterization (`raster = TRUE`)
#'                         is recommended (default 50000L).
#' @return Named list with `n_cells`, `max_interactive`, `raster_threshold`,
#'         `should_subsample`, `use_raster`, `target_points`,
#'         `subsample_fraction`, and `mode` (`"full"`, `"subsampled"`, or `"raster"`).
sc_umap_render_budget <- function(n_cells, max_interactive = 25000L, raster_threshold = 50000L) {
  n_cells          <- .validate_sc_budget_int(n_cells, "n_cells", 0L)
  max_interactive  <- .validate_sc_budget_int(max_interactive, "max_interactive", 1L)
  raster_threshold <- .validate_sc_budget_int(raster_threshold, "raster_threshold", 1L)

  should_subsample <- n_cells > max_interactive
  use_raster       <- n_cells > raster_threshold
  target_points    <- if (should_subsample) max_interactive else n_cells
  frac             <- if (n_cells == 0L) 1 else min(1, max_interactive / n_cells)
  mode             <- if (use_raster) "raster" else if (should_subsample) "subsampled" else "full"

  list(
    n_cells            = n_cells,
    max_interactive    = max_interactive,
    raster_threshold   = raster_threshold,
    should_subsample   = should_subsample,
    use_raster         = use_raster,
    target_points      = target_points,
    subsample_fraction = frac,
    mode               = mode
  )
}

#' Stratified deterministic subsampling of UMAP cell indices.
#'
#' When `n_cells <= max_points`, returns `seq_len(n_cells)`. Otherwise samples
#' `max_points` unique sorted indices. When `groups` is provided, guarantees
#' at least one cell per cluster (whenever `length(unique(groups)) <= max_points`)
#' and allocates the remaining budget proportionally via largest-remainder.
#'
#' @param n_cells    Non-negative integer number of cells.
#' @param max_points Positive integer target point budget (default 25000L).
#' @param groups     Optional vector/factor of length `n_cells` giving cluster
#'                   or group labels for stratified sampling.
#' @param seed       Integer RNG seed for deterministic sampling (default 42L).
#'                   Pass `NULL` to use the current RNG state.
#' @return Sorted integer vector of selected 1-based cell indices.
sc_subsample_umap_indices <- function(n_cells, max_points = 25000L, groups = NULL, seed = 42L) {
  n_cells    <- .validate_sc_budget_int(n_cells, "n_cells", 0L)
  max_points <- .validate_sc_budget_int(max_points, "max_points", 1L)
  if (!is.null(seed)) {
    seed <- .validate_sc_budget_int(seed, "seed", -.Machine$integer.max)
  }
  if (!is.null(groups) && length(groups) != n_cells) {
    stop(errorCondition(paste0("groups doit avoir une longueur égale à n_cells (", n_cells, ")"),
                        class = "sc_plotting_error"))
  }

  if (n_cells == 0L) return(integer(0))
  if (n_cells <= max_points) return(seq_len(n_cells))

  if (!is.null(seed)) {
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
    set.seed(seed)
  }

  if (is.null(groups)) {
    return(sort(sample.int(n_cells, size = max_points, replace = FALSE)))
  }

  grp_chr <- as.character(groups)
  grp_chr[is.na(grp_chr)] <- "<NA>"
  grp_split <- split(seq_len(n_cells), grp_chr)
  k <- length(grp_split)

  if (k > max_points) {
    chosen_groups <- sample.int(k, size = max_points, replace = FALSE)
    picked <- vapply(grp_split[chosen_groups], function(idx) {
      if (length(idx) == 1L) idx else idx[sample.int(length(idx), size = 1L)]
    }, integer(1))
    return(sort(unname(picked)))
  }

  # Base allocation of 1 cell per group preserves rare/singleton clusters.
  sizes <- lengths(grp_split)
  alloc <- rep(1L, k)
  rem   <- max_points - k
  cap   <- sizes - 1L
  total_cap <- sum(cap)

  if (rem > 0L && total_cap > 0L) {
    raw_extra   <- rem * (cap / total_cap)
    extra_floor <- pmin(cap, as.integer(floor(raw_extra)))
    alloc       <- alloc + extra_floor
    leftover    <- rem - sum(extra_floor)
    if (leftover > 0L) {
      eligible <- which((sizes - alloc) > 0L)
      ord      <- eligible[order(raw_extra[eligible] - extra_floor[eligible], decreasing = TRUE)]
      bump     <- ord[seq_len(min(leftover, length(ord)))]
      alloc[bump] <- alloc[bump] + 1L
    }
  }

  selected <- unlist(Map(function(idx, n_take) {
    if (n_take >= length(idx)) idx else idx[sample.int(length(idx), size = n_take, replace = FALSE)]
  }, grp_split, alloc), use.names = FALSE)

  sort(as.integer(selected))
}

