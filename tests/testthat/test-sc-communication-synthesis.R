# Tests pour la synthèse multi-modale CCC (5e pill exploratoire)

test_that("build_communication_multimodal_synthesis : cadrage épistémique strict et 0 contexte", {
  res <- .ccx_result()
  syn <- build_communication_multimodal_synthesis(res)
  expect_identical(syn$type, "ccc_multimodal_synthesis")
  expect_identical(syn$label, "Convergence des contextes (filtrage exploratoire)")
  expect_false(isTRUE(syn$is_statistical_pvalue))
  expect_identical(syn$qc$n_contexts_active, 0L)
  expect_true(nrow(syn$pair_table) >= 1L)
  expect_true(all(syn$pair_table$n_contexts_evaluated == 0L))
  expect_true(all(syn$pair_table$convergence_tier == "Aucun contexte calculé"))

  p <- plot_communication_multimodal_synthesis(syn)
  expect_s3_class(p, "ggplot")
})

test_that("build_communication_multimodal_synthesis : convergence sur 4 contextes simulés", {
  res <- .ccx_result()
  pt_pairs <- unique(.communication_node_keys(res$canonical_table)[, c("sender_node", "receiver_node")])
  p1 <- pt_pairs[1L, , drop = FALSE]
  p2 <- if (nrow(pt_pairs) >= 2L) pt_pairs[2L, , drop = FALSE] else p1

  spatial_ctx <- list(
    type = "ccc_spatial_context",
    pair_table = data.frame(
      sender_node = pt_pairs$sender_node,
      receiver_node = pt_pairs$receiver_node,
      frac_receiver_within_radius = c(0.60, rep(0.05, nrow(pt_pairs) - 1L)),
      frac_sender_within_radius = c(0.55, rep(0.05, nrow(pt_pairs) - 1L)),
      p_perm = c(0.01, rep(0.45, nrow(pt_pairs) - 1L)),
      stringsAsFactors = FALSE
    )
  )
  trajectory_ctx <- list(
    type = "ccc_trajectory_context",
    pair_bin_table = data.frame(
      sender_node = rep(p1$sender_node, 3L),
      receiver_node = rep(p1$receiver_node, 3L),
      bin = c(1L, 2L, 3L),
      mean_ligand_expression_senders = c(1, 2, 3),
      mean_receptor_expression_receivers = c(1, 2, 3),
      stringsAsFactors = FALSE
    )
  )
  velocity_ctx <- list(
    type = "ccc_velocity_context",
    status = "valid",
    pair_table = data.frame(
      sender_node = pt_pairs$sender_node,
      receiver_node = pt_pairs$receiver_node,
      receiver_velocity_mean = c(2.0, rep(0.5, nrow(pt_pairs) - 1L)),
      stringsAsFactors = FALSE
    )
  )
  perturbation_ctx <- list(
    type = "ccc_perturbation",
    delta_table = data.frame(
      sender_node = pt_pairs$sender_node,
      receiver_node = pt_pairs$receiver_node,
      delta_score_total = c(-0.9, rep(0, nrow(pt_pairs) - 1L)),
      affected = c(TRUE, rep(FALSE, nrow(pt_pairs) - 1L)),
      stringsAsFactors = FALSE
    )
  )

  syn <- build_communication_multimodal_synthesis(
    res,
    spatial_ctx = spatial_ctx,
    trajectory_ctx = trajectory_ctx,
    velocity_ctx = velocity_ctx,
    perturbation_ctx = perturbation_ctx
  )
  expect_identical(syn$qc$n_contexts_active, 4L)
  pt <- syn$pair_table
  expect_identical(pt$sender_node[1L], p1$sender_node)
  expect_identical(pt$receiver_node[1L], p1$receiver_node)
  expect_identical(pt$n_contexts_evaluated[1L], 4L)
  expect_identical(pt$n_contexts_convergent[1L], 4L)
  expect_identical(pt$convergence_tier[1L], "Exploratoire — fort (>=3)")
  expect_identical(pt$context_badges[1L], "Spatial | Trajectoire | Vélocité | Perturbation")

  p <- plot_communication_multimodal_synthesis(syn, top_n = 10L)
  expect_s3_class(p, "ggplot")
})
