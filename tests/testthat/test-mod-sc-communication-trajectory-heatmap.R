# =============================================================================
# test-mod-sc-communication-trajectory-heatmap.R — W1.3 Lot 1 :
# Anti-collision verticale de la heatmap trajectoire (hauteur dynamique +
# plafond top-N paires + conteneur scrollable + taille des libelles Y)
# =============================================================================

suppressPackageStartupMessages(library(shiny))

source_project_file("tools/pkg_skeleton/R/sc_communication_trajectory.R")
source_project_file("modules/sc/mod_sc_communication_trajectory.R")

.make_large_trj_ctx <- function(n_pairs = 40L, n_bins = 5L) {
  senders <- rep(c("PopA", "PopB"), length.out = n_pairs)
  receivers <- sprintf("Target_%02d", seq_len(n_pairs))
  scores <- seq(10, 1, length.out = n_pairs)

  pair_summary <- data.frame(
    sender_node = senders,
    receiver_node = receivers,
    lineage = NA_character_,
    n_sender_cells = rep(12L, n_pairs),
    n_receiver_cells = rep(10L, n_pairs),
    ligand_genes = rep("LIG1", n_pairs),
    receptor_genes = rep("REC1", n_pairs),
    n_ligand_genes_missing = rep(0L, n_pairs),
    n_receptor_genes_missing = rep(0L, n_pairs),
    mean_imported_score = scores,
    stringsAsFactors = FALSE
  )

  rows <- vector("list", n_pairs * n_bins)
  idx <- 0L
  for (p in seq_len(n_pairs)) {
    for (b in seq_len(n_bins)) {
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        sender_node = senders[p],
        receiver_node = receivers[p],
        lineage = NA_character_,
        bin = b,
        bin_from = (b - 1) / n_bins,
        bin_to = b / n_bins,
        bin_mid = (b - 0.5) / n_bins,
        n_cells_in_bin = 20L,
        n_sender_cells = 4L,
        n_receiver_cells = 4L,
        frac_sender_of_bin = 0.2,
        frac_receiver_of_bin = 0.2,
        frac_sender_of_population = 0.25,
        frac_receiver_of_population = 0.25,
        mean_ligand_expression_senders = 1.5,
        mean_receptor_expression_receivers = 1.5,
        stringsAsFactors = FALSE
      )
    }
  }
  pair_bin_table <- do.call(rbind, rows)

  list(
    type = "ccc_trajectory_context",
    status = "valid",
    analysis_id = "sc-communication-trajectory",
    parent_analysis_id = "sc-communication-import",
    source_method = "cellchat",
    identity_column = "cell_type",
    params = list(n_bins_effective = n_bins, lineage_mode = "global"),
    pair_bin_table = pair_bin_table,
    pair_summary = pair_summary,
    qc = list(n_pairs_total = n_pairs, n_combinations_mapped = n_pairs, n_bins_effective = n_bins),
    warnings = character(0),
    provenance = list(),
    timestamp_utc = "2026-11-20T00:00:00Z"
  )
}

test_that("mod_sc_communication_trajectory_ui rend trj_heatmap en height=auto dans un conteneur scrollable", {
  e <- new.env(parent = globalenv())
  e$i18n <- list(t = function(x) x)
  e$.tr_plain <- function(x) x
  sys.source(file.path(ts_project_root(), "modules/sc/mod_sc_communication_trajectory.R"), envir = e)

  html <- as.character(e$mod_sc_communication_trajectory_ui("trj"))
  expect_true(grepl("overflow-y:\\s*auto", html, perl = TRUE))
  expect_false(grepl('id="trj-trj_heatmap"[^>]*height:\\s*340px', html, perl = TRUE))
})

test_that("heatmap trajectoire : plafond top-N, hauteur dynamique et taille Y adaptee pour >= 40 paires", {
  gd <- shiny::reactiveValues(language = "fr", i18n = NULL, sc_obj = NULL)
  shared <- shiny::reactiveValues()
  comm_state <- shiny::reactiveValues()

  shiny::testServer(mod_sc_communication_trajectory_server,
                    args = list(id = "trajectory", comm_state = comm_state,
                                global_data = gd, shared_rv = shared), {
    small_ctx <- .make_large_trj_ctx(n_pairs = 4L)
    prep_small <- .trj_heatmap_prepare_ctx(small_ctx)
    expect_identical(prep_small$n_pairs, 4L)
    expect_false(prep_small$capped)
    expect_identical(.trj_heatmap_height(prep_small$n_pairs), 340L)

    big_ctx <- .make_large_trj_ctx(n_pairs = 40L)
    prep_big <- .trj_heatmap_prepare_ctx(big_ctx)
    expect_true(prep_big$capped)
    expect_identical(prep_big$n_total, 40L)
    expect_lte(prep_big$n_pairs, 30L)
    expect_gt(.trj_heatmap_height(prep_big$n_pairs), 340L)

    # Top-N conserve les paires aux meilleurs scores (Target_01..Target_30)
    kept_targets <- unique(prep_big$ctx$pair_bin_table$receiver_node)
    expect_true("Target_01" %in% kept_targets)
    expect_false("Target_40" %in% kept_targets)

    p <- .trj_heatmap_plot(big_ctx)
    expect_s3_class(p, "ggplot")
    y_sz <- p$theme$axis.text.y$size
    expect_true(!is.null(y_sz) && y_sz <= 8.5)
  })
})
