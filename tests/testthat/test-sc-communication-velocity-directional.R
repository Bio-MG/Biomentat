# =============================================================================
# test-sc-communication-velocity-directional.R — W2.2 §4.3 : directionnalité
# vélocité CCC (alignement T_ij + corrélation Spearman vitesse-récepteur)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/sc_communication_velocity.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_velocity_directional.R")

# Fixture velocity canonique sur 9 cellules (3 par population A, B, C) afin de
# permettre des corrélations de Spearman avec n >= 3 cellules par population.
.vdir_cells <- paste0("cell", 1:9)

.vdir_identities <- function() {
  setNames(c("A", "A", "A", "B", "B", "B", "C", "C", "C"), .vdir_cells)
}

.vdir_velocity_result <- function(with_vectors = TRUE) {
  genes <- paste0("G", 1:4)
  m <- Matrix::Matrix(
    seq_len(36L), nrow = 4L, ncol = 9L, sparse = TRUE,
    dimnames = list(genes, .vdir_cells)
  )
  validated <- validate_velocity_matrices(
    spliced = m, unspliced = m * 2L,
    seurat_cells = .vdir_cells, seurat_genes = genes,
    orientation = "auto_strict"
  )
  validated$velocity_source <- "test"
  validated$velocity_method <- "precomputed"
  validated$input_orientation <- NA_character_
  validated$embedding_reduction <- "umap"
  vecs <- if (with_vectors) {
    # Magnitudes strictement croissantes dans A (1, 2, 3), B (1, 2, 3), C (3, 2, 1)
    dx <- c(1, 2, 3, 1, 2, 3, 3, 2, 1)
    dy <- rep(0, 9L)
    matrix(c(dx, dy), ncol = 2L, dimnames = list(.vdir_cells, c("dx", "dy")))
  } else {
    NULL
  }
  finalize_velocity_result(
    validated, input_mode = "rds", input_files = list(rds = "vdir.rds"),
    velocity_vectors = vecs, requested_reduction = "umap"
  )
}

# Matrice de transition T_ij carrée 9x9 déterministe :
# Forte masse A -> B (0.25 par paire de cellules), faible masse ailleurs (0.01).
.vdir_transition_matrix <- function() {
  T_mat <- matrix(0.01, nrow = 9L, ncol = 9L,
                  dimnames = list(.vdir_cells, .vdir_cells))
  diag(T_mat) <- 0.50
  # A (cell1..3) -> B (cell4..6) fortement aligné
  T_mat[1:3, 4:6] <- 0.25
  # C (cell7..9) -> B (cell4..6) très faible
  T_mat[7:9, 4:6] <- 0.002
  T_mat
}

# Matrice d'expression gènes x cellules (R1 corrélé +1 avec vitesse dans B,
# R2 corrélé -1 avec vitesse dans C).
.vdir_receptor_matrix <- function() {
  mat <- matrix(0.1, nrow = 2L, ncol = 9L,
                dimnames = list(c("R1", "R2"), .vdir_cells))
  # Dans B (cell4..6), vitesse = (1, 2, 3) -> R1 croissant (10, 20, 30) => rho = +1
  mat["R1", 4:6] <- c(10, 20, 30)
  # Dans C (cell7..9), vitesse = (3, 2, 1) -> R2 croissant (1, 2, 3) => rho = -1
  mat["R2", 7:9] <- c(1, 2, 3)
  # Dans A (cell1..3), constant
  mat["R1", 1:3] <- 5
  mat["R2", 1:3] <- 5
  mat
}

test_that("transition probability T_ij alignment computes exact mean, sum, enrichment and flag", {
  r <- .ccx_result(identities = .vdir_identities())
  vr <- .vdir_velocity_result()
  T_mat <- .vdir_transition_matrix()

  res <- compute_communication_velocity_directionality(
    communication_result = r,
    velocity_result = vr,
    cell_identities = .vdir_identities(),
    transition_matrix = T_mat
  )

  expect_identical(res$type, "ccc_velocity_directionality")
  expect_identical(res$status, "valid")
  expect_identical(res$transition_status, "valid")
  expect_identical(res$cor_status, "unavailable_no_receptor_expr")

  pt <- res$pair_table
  bg_expected <- mean(T_mat[row(T_mat) != col(T_mat)])

  ab <- pt[pt$sender_node == "A" & pt$receiver_node == "B", ]
  expect_equal(ab$transition_mean, 0.25, tolerance = 1e-12)
  expect_equal(ab$transition_sum, 0.25 * 9, tolerance = 1e-12)
  expect_equal(ab$transition_bg_mean, bg_expected, tolerance = 1e-12)
  expect_equal(ab$transition_enrichment, 0.25 / bg_expected, tolerance = 1e-12)
  expect_true(ab$transition_aligned)

  cb <- pt[pt$sender_node == "C" & pt$receiver_node == "B", ]
  expect_equal(cb$transition_mean, 0.002, tolerance = 1e-12)
  expect_false(cb$transition_aligned)

  # Sans receptor_expr, les colonnes de corrélation restent strictement NA
  expect_true(all(is.na(pt$cor_spearman)))
  expect_true(all(pt$cor_status == "unavailable_no_receptor_expr"))
})

test_that("within-population Spearman correlation works with gene x cell matrix and named numeric vector", {
  r <- .ccx_result(identities = .vdir_identities())
  vr <- .vdir_velocity_result()
  rmat <- .vdir_receptor_matrix()

  # Cas 1 : matrice gènes x cellules
  res_mat <- compute_communication_velocity_directionality(
    communication_result = r,
    velocity_result = vr,
    cell_identities = .vdir_identities(),
    receptor_expr = rmat
  )
  expect_identical(res_mat$status, "valid")
  expect_identical(res_mat$transition_status, "unavailable_no_transition_matrix")
  expect_identical(res_mat$cor_status, "valid")

  pt <- res_mat$pair_table
  expect_true(all(is.na(pt$transition_mean)))
  expect_true(all(is.na(pt$transition_aligned)))

  ab <- pt[pt$sender_node == "A" & pt$receiver_node == "B", ]
  expect_equal(ab$cor_spearman, 1.0, tolerance = 1e-12)
  expect_identical(ab$cor_n_cells, 3L)
  expect_identical(ab$cor_status, "valid")

  ac <- pt[pt$sender_node == "A" & pt$receiver_node == "C", ]
  expect_equal(ac$cor_spearman, -1.0, tolerance = 1e-12)
  expect_identical(ac$cor_n_cells, 3L)
  expect_identical(ac$cor_status, "valid")

  # Cas 2 : vecteur numérique nommé
  rvec <- setNames(c(1, 2, 3, 30, 20, 10, 5, 5, 5), .vdir_cells)
  res_vec <- compute_communication_velocity_directionality(
    communication_result = r,
    velocity_result = vr,
    cell_identities = .vdir_identities(),
    receptor_expr = rvec
  )
  pt_v <- res_vec$pair_table
  ab_v <- pt_v[pt_v$sender_node == "A" & pt_v$receiver_node == "B", ]
  expect_equal(ab_v$cor_spearman, -1.0, tolerance = 1e-12)
  expect_identical(ab_v$cor_status, "valid")

  # Dans C, rvec est constant (5, 5, 5) -> zero_variance, NA
  ac_v <- pt_v[pt_v$sender_node == "A" & pt_v$receiver_node == "C", ]
  expect_true(is.na(ac_v$cor_spearman))
  expect_identical(ac_v$cor_status, "zero_variance")
})

test_that("missing transition_matrix and receptor_expr degrade gracefully to NA without fabrication", {
  r <- .ccx_result(identities = .vdir_identities())
  vr <- .vdir_velocity_result()

  res <- compute_communication_velocity_directionality(
    communication_result = r,
    velocity_result = vr,
    cell_identities = .vdir_identities(),
    transition_matrix = NULL,
    receptor_expr = NULL
  )
  expect_identical(res$status, "unavailable")
  expect_identical(res$transition_status, "unavailable_no_transition_matrix")
  expect_identical(res$cor_status, "unavailable_no_receptor_expr")
  expect_true(all(is.na(res$pair_table$transition_mean)))
  expect_true(all(is.na(res$pair_table$transition_sum)))
  expect_true(all(is.na(res$pair_table$transition_aligned)))
  expect_true(all(is.na(res$pair_table$cor_spearman)))

  # Sans vecteurs velocity + avec receptor_expr -> cor_status explicite
  vr_novec <- .vdir_velocity_result(with_vectors = FALSE)
  res_nv <- compute_communication_velocity_directionality(
    communication_result = r,
    velocity_result = vr_novec,
    cell_identities = .vdir_identities(),
    transition_matrix = .vdir_transition_matrix(),
    receptor_expr = .vdir_receptor_matrix()
  )
  expect_identical(res_nv$status, "valid")
  expect_identical(res_nv$transition_status, "valid")
  expect_identical(res_nv$cor_status, "unavailable_no_velocity_vectors")
  expect_true(all(is.na(res_nv$pair_table$cor_spearman)))
})

test_that("filtered_table is honored and parent canonical_table is never mutated", {
  r <- .ccx_result(identities = .vdir_identities())
  before <- r$canonical_table
  fs <- communication_apply_filters(r, list(pathways = "P1"))

  res <- compute_communication_velocity_directionality(
    communication_result = r,
    velocity_result = .vdir_velocity_result(),
    cell_identities = .vdir_identities(),
    transition_matrix = .vdir_transition_matrix(),
    receptor_expr = .vdir_receptor_matrix(),
    filtered_table = fs$table
  )
  expect_identical(r$canonical_table, before)
  expect_identical(nrow(res$pair_table), 2L)
  expect_false(any(res$pair_table$sender_node == "A" & res$pair_table$receiver_node == "C"))
})

test_that("epistemic guardrail is present and invalid inputs raise communication_context_error", {
  r <- .ccx_result(identities = .vdir_identities())
  vr <- .vdir_velocity_result()
  ids <- .vdir_identities()

  res <- compute_communication_velocity_directionality(
    r, vr, ids, transition_matrix = .vdir_transition_matrix(),
    receptor_expr = .vdir_receptor_matrix()
  )
  expect_true(grepl("AUCUNE causalité", res$epistemic_note))

  # cell_identities non nommé
  expect_error(
    compute_communication_velocity_directionality(r, vr, c("A", "B")),
    class = "communication_context_error"
  )
  # transition_matrix non carrée
  bad_T <- matrix(0.1, nrow = 2L, ncol = 3L,
                  dimnames = list(c("cell1", "cell2"), c("cell1", "cell2", "cell3")))
  expect_error(
    compute_communication_velocity_directionality(r, vr, ids, transition_matrix = bad_T),
    "carrée", class = "communication_context_error"
  )
  # transition_matrix avec valeurs négatives
  neg_T <- .vdir_transition_matrix()
  neg_T[1L, 2L] <- -0.5
  expect_error(
    compute_communication_velocity_directionality(r, vr, ids, transition_matrix = neg_T),
    "négative", class = "communication_context_error"
  )
  # receptor_expr non nommé
  expect_error(
    compute_communication_velocity_directionality(r, vr, ids, receptor_expr = c(1, 2, 3)),
    class = "communication_context_error"
  )
})

test_that("plot_communication_velocity_directionality returns valid ggplot and stub when unavailable", {
  r <- .ccx_result(identities = .vdir_identities())
  vr <- .vdir_velocity_result()
  ids <- .vdir_identities()

  # Cas 1 : T_ij + corrélation
  res_full <- compute_communication_velocity_directionality(
    r, vr, ids,
    transition_matrix = .vdir_transition_matrix(),
    receptor_expr = .vdir_receptor_matrix()
  )
  p_full <- plot_communication_velocity_directionality(res_full, top_n = 2L)
  expect_s3_class(p_full, "ggplot")
  expect_true(grepl("AUCUNE causalité", p_full$labels$subtitle %||% ""))
  expect_lte(nrow(p_full$data), 2L)

  # Cas 2 : corrélation seule (sans T_ij)
  res_cor <- compute_communication_velocity_directionality(
    r, vr, ids,
    receptor_expr = .vdir_receptor_matrix()
  )
  p_cor <- plot_communication_velocity_directionality(res_cor)
  expect_s3_class(p_cor, "ggplot")
  expect_true(grepl("AUCUNE causalité", p_cor$labels$subtitle %||% ""))

  # Cas 3 : indisponible -> stub plot propre
  res_unavail <- compute_communication_velocity_directionality(r, vr, ids)
  p_stub <- plot_communication_velocity_directionality(res_unavail)
  expect_s3_class(p_stub, "ggplot")
  expect_true(grepl("indisponible", p_stub$labels$subtitle %||% ""))
  expect_true(grepl("AUCUNE causalité", p_stub$labels$subtitle %||% ""))

  # Cas 4 : objet invalide -> erreur classée
  expect_error(
    plot_communication_velocity_directionality(list(type = "wrong")),
    class = "communication_context_error"
  )
})
