# =============================================================================
# test-sc-communication-perturbation-multi.R — W2.2 §5.4 : perturbations
# multi-cibles, pathway et population (sender + receiver) IN SILICO
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_communication_perturbation.R")

.pert_multi_na_score_result <- function() {
  tab <- .ccx_cellchat_tab()
  tab$prob <- NA_real_
  parsed <- parse_cellchat_import(tab, source_file = "na_scores.csv")
  harm <- harmonize_communication_identities(
    parsed$table, unique(.ccx_identities()),
    "cell_type", context = "communication import"
  )
  qcr <- communication_import_qc(harm$table)
  finalize_communication_result(
    canonical_table = qcr$table, source_method = "cellchat",
    source_files = list(table = "na_scores.csv"), identity_column = "cell_type",
    identity_mapping = harm$mapping, identity_summary = harm$summary,
    column_mapping = parsed$column_mapping, qc = qcr$counts,
    n_input_rows = parsed$n_input_rows,
    seurat_obj = matrix(0, 2, 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))
  )
}

test_that("multi-ligand removal and attenuation produce exact deltas", {
  r <- .ccx_result()
  ctx_rm <- build_communication_perturbation(
    r, target = "ligand", value = c("L1", "L2"), mode = "remove"
  )
  expect_identical(ctx_rm$value, c("L1", "L2"))
  expect_identical(ctx_rm$baseline_summary$n_interactions, 5L)
  expect_identical(ctx_rm$perturbed_summary$n_interactions, 0L)
  expect_identical(ctx_rm$qc$n_removed, 5L)
  expect_identical(ctx_rm$qc$n_pairs_affected, 5L)
  expect_true(all(ctx_rm$delta_table$affected))
  expect_true(all(ctx_rm$delta_table$score_total_perturbed == 0))
  expect_true(all(abs(ctx_rm$delta_table$delta_fraction + 1) < 1e-12))
  expect_true(any(grepl("vide", ctx_rm$warnings)))
  expect_true(isTRUE(ctx_rm$params$multi_target))
  expect_identical(ctx_rm$params$n_targets, 2L)

  ctx_att <- build_communication_perturbation(
    r, target = "ligand", value = c("L1", "L2"), mode = "attenuate", factor = 0.25
  )
  expect_identical(ctx_att$perturbed_summary$n_interactions, 5L)
  expect_identical(ctx_att$qc$n_attenuated, 5L)
  expect_identical(ctx_att$qc$n_removed, 0L)
  expect_equal(
    ctx_att$delta_table$score_total_perturbed,
    ctx_att$delta_table$score_total_baseline * 0.25,
    tolerance = 1e-12
  )
  expect_equal(
    ctx_att$delta_table$delta_fraction,
    rep(-0.75, nrow(ctx_att$delta_table)),
    tolerance = 1e-12
  )
})

test_that("multi-receptor and multi-interaction perturbation work in both modes", {
  r <- .ccx_result()
  ctx_rec <- build_communication_perturbation(
    r, target = "receptor", value = c("R1", "R2"), mode = "attenuate", factor = 0.5
  )
  expect_identical(ctx_rec$qc$n_attenuated, 5L)
  expect_equal(ctx_rec$perturbed_summary$score_total,
               ctx_rec$baseline_summary$score_total * 0.5, tolerance = 1e-12)

  ctx_int <- build_communication_perturbation(
    r, target = "interaction", value = c("L1 -> R1", "L2 -> R2"), mode = "remove"
  )
  expect_identical(ctx_int$qc$n_removed, 5L)
  expect_identical(ctx_int$perturbed_summary$n_interactions, 0L)
})

test_that("whole-pathway perturbation removes or attenuates exact pathway rows", {
  r <- .ccx_result()
  # P1 = A->B (0.4) and C->B (0.2)
  ctx_p1 <- build_communication_perturbation(
    r, target = "pathway", value = "P1", mode = "remove"
  )
  expect_identical(ctx_p1$target, "pathway")
  expect_identical(ctx_p1$value, "P1")
  expect_identical(ctx_p1$qc$n_removed, 2L)
  expect_identical(ctx_p1$perturbed_summary$n_interactions, 3L)
  expect_setequal(ctx_p1$qc$affected_senders, c("A", "C"))
  expect_setequal(ctx_p1$qc$affected_receivers, "B")
  ab <- ctx_p1$delta_table[ctx_p1$delta_table$sender_node == "A" &
                             ctx_p1$delta_table$receiver_node == "B", ]
  expect_equal(ab$delta_score_total, -0.4, tolerance = 1e-12)
  cb <- ctx_p1$delta_table[ctx_p1$delta_table$sender_node == "C" &
                             ctx_p1$delta_table$receiver_node == "B", ]
  expect_equal(cb$delta_score_total, -0.2, tolerance = 1e-12)
  da <- ctx_p1$delta_table[ctx_p1$delta_table$sender_node == "D", ]
  expect_false(da$affected)

  # Multi-pathway P1 + P9 = A->B (0.4), C->B (0.2), D->A (0.9)
  ctx_multi <- build_communication_perturbation(
    r, target = "pathway", value = c("P1", "P9"), mode = "attenuate", factor = 0.5
  )
  expect_identical(ctx_multi$qc$n_attenuated, 3L)
  expect_equal(ctx_multi$perturbed_summary$score_total,
               ctx_multi$baseline_summary$score_total - 0.5 * (0.4 + 0.2 + 0.9),
               tolerance = 1e-12)
})

test_that("whole-population KO affects both sender and receiver roles without double-attenuating autocrine", {
  r <- .ccx_result()
  # Population A: sender in A->B (0.4), A->C (0.6); receiver in D->A (0.9)
  ctx_a <- build_communication_perturbation(
    r, target = "population", value = "A", mode = "remove"
  )
  expect_identical(ctx_a$target, "population")
  expect_identical(ctx_a$qc$n_removed, 3L)
  expect_identical(ctx_a$perturbed_summary$n_interactions, 2L)
  node_a <- ctx_a$node_table[ctx_a$node_table$node == "A", ]
  expect_identical(node_a$out_interactions_perturbed, 0L)
  expect_identical(node_a$in_interactions_perturbed, 0L)
  expect_equal(node_a$delta_out_score, -1.0, tolerance = 1e-12)
  expect_equal(node_a$delta_in_score, -0.9, tolerance = 1e-12)

  # Population C: sender in C->B (0.2), C->C (0.3); receiver in A->C (0.6), C->C (0.3)
  # Autocrine C->C must be attenuated ONCE (* 0.5 = 0.15), never squared (* 0.25).
  ctx_c <- build_communication_perturbation(
    r, target = "population", value = "C", mode = "attenuate", factor = 0.5
  )
  expect_identical(ctx_c$qc$n_attenuated, 3L)
  cc <- ctx_c$delta_table[ctx_c$delta_table$sender_node == "C" &
                            ctx_c$delta_table$receiver_node == "C", ]
  expect_equal(cc$score_total_perturbed, 0.15, tolerance = 1e-12)
  expect_equal(cc$delta_fraction, -0.5, tolerance = 1e-12)

  # Multi-population KO: A + B removes A->B, A->C, C->B, D->A (4 rows), leaving only C->C
  ctx_ab <- build_communication_perturbation(
    r, target = "population", value = c("A", "B"), mode = "remove"
  )
  expect_identical(ctx_ab$qc$n_removed, 4L)
  expect_identical(ctx_ab$perturbed_summary$n_interactions, 1L)

  # Unmapped population D still works via coalesced node keys
  ctx_d <- build_communication_perturbation(
    r, target = "population", value = "D", mode = "remove"
  )
  expect_identical(ctx_d$qc$n_removed, 1L)
  expect_identical(ctx_d$perturbed_summary$n_interactions, 4L)
})

test_that("parity with single-target build_communication_perturbation()", {
  r <- .ccx_result()
  # B only receives in .ccx_result(): population B == receiver B
  ctx_pop_b <- build_communication_perturbation(r, target = "population", value = "B", mode = "remove")
  ctx_rec_b <- build_communication_perturbation(r, target = "receiver", value = "B", mode = "remove")
  expect_identical(ctx_pop_b$delta_table, ctx_rec_b$delta_table)
  expect_identical(ctx_pop_b$node_table, ctx_rec_b$node_table)
  expect_identical(ctx_pop_b$baseline_summary, ctx_rec_b$baseline_summary)
  expect_identical(ctx_pop_b$perturbed_summary, ctx_rec_b$perturbed_summary)
  expect_identical(ctx_pop_b$qc, ctx_rec_b$qc)

  # Duplicate scalar value deduplicates to exact scalar parity
  ctx_dup <- build_communication_perturbation(r, target = "ligand", value = c("L1", "L1"), mode = "remove")
  ctx_one <- build_communication_perturbation(r, target = "ligand", value = "L1", mode = "remove")
  expect_identical(ctx_dup$value, ctx_one$value)
  expect_identical(ctx_dup$delta_table, ctx_one$delta_table)
  expect_identical(ctx_dup$node_table, ctx_one$node_table)
  expect_identical(ctx_dup$qc, ctx_one$qc)
})

test_that("invalid inputs throw communication_context_error", {
  r <- .ccx_result()
  expect_error(
    build_communication_perturbation(r, target = "gene", value = "L1"),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_perturbation(r, target = "ligand", value = character(0)),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_perturbation(r, target = "ligand", value = c("L1", NA_character_)),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_perturbation(r, target = "ligand", value = c("L1", "")),
    class = "communication_context_error"
  )
  expect_error(
    build_communication_perturbation(r, target = "ligand", value = 1),
    class = "communication_context_error"
  )
  # Partial missing in multi-target: strict exact match, never silent skip
  expect_error(
    build_communication_perturbation(r, target = "ligand", value = c("L1", "ZZZ")),
    "absente",
    class = "communication_context_error"
  )
  expect_error(
    build_communication_perturbation(r, target = "pathway", value = "MISSING_PW"),
    "absente",
    class = "communication_context_error"
  )
  expect_error(
    build_communication_perturbation(r, target = "population", value = c("A", "GHOST")),
    "absente",
    class = "communication_context_error"
  )
  # Pathway column absent from filtered_table
  no_pw <- r$canonical_table[, setdiff(names(r$canonical_table), "pathway"), drop = FALSE]
  expect_error(
    build_communication_perturbation(r, target = "pathway", value = "P1", filtered_table = no_pw),
    class = "communication_context_error"
  )
  # Invalid factor in attenuate mode
  expect_error(
    build_communication_perturbation(r, target = "pathway", value = "P1",
                                     mode = "attenuate", factor = 0),
    class = "communication_context_error"
  )
})

test_that("epistemic guard, immutability, NA preservation, plots and export hold for multi-target", {
  r <- .ccx_result()
  before <- r$canonical_table

  ctx <- build_communication_perturbation(
    r, target = "pathway", value = c("P1", "P9"), mode = "attenuate", factor = 0.5
  )
  expect_identical(r$canonical_table, before)
  expect_identical(ctx$params$label, "IN SILICO PERTURBATION")
  expect_identical(ctx$params$order, "first_order_no_propagation")
  expect_identical(ctx$provenance$parameters$label, "in_silico_not_ko")
  expect_true(any(grepl("IN SILICO", ctx$warnings)))
  expect_true(any(grepl("PREMIER ORDRE", ctx$warnings)))

  g1 <- plot_communication_perturbation_delta(ctx)
  g2 <- plot_communication_perturbation_nodes(ctx)
  expect_s3_class(g1, "ggplot")
  expect_s3_class(g2, "ggplot")
  expect_true(grepl("IN SILICO", g1$labels$subtitle %||% ""))
  expect_true(grepl("IN SILICO", g2$labels$subtitle %||% ""))
  expect_true(grepl("P1, P9", g1$labels$title %||% "", fixed = TRUE))
  expect_true(grepl("P1, P9", g2$labels$title %||% "", fixed = TRUE))

  ex <- build_communication_perturbation_export(ctx)
  expect_identical(nrow(ex), nrow(ctx$delta_table))
  expect_true(all(ex$target == "pathway"))
  expect_true(all(ex$value == "P1, P9"))
  expect_true(all(ex$label == "IN SILICO PERTURBATION"))

  # NA scores are never fabricated as 0
  r_na <- .pert_multi_na_score_result()
  ctx_na <- build_communication_perturbation(
    r_na, target = "population", value = c("A", "C"), mode = "remove"
  )
  expect_true(all(is.na(ctx_na$delta_table$score_total_baseline)))
  expect_true(all(is.na(ctx_na$delta_table$score_total_perturbed)))
  expect_true(all(is.na(ctx_na$delta_table$delta_score_total)))
  expect_true(all(is.na(ctx_na$node_table$delta_out_score)))
  expect_true(all(is.na(ctx_na$node_table$delta_in_score)))
})
