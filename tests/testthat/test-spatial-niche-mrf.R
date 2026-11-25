# =============================================================================
# test-spatial-niche-mrf.R — BK-10 Bayesian MRF Spatial Niche Clustering Tests
# =============================================================================

source_project_file("tools/pkg_skeleton/R/spatial_niche.R")

.mrf_err <- function(expr) {
  tryCatch({ expr; NULL },
           error = function(e) list(msg = conditionMessage(e), class = class(e)))
}

.mrf_grid <- function(nx = 10L, ny = 10L) {
  g <- expand.grid(x = seq_len(nx), y = seq_len(ny))
  data.frame(
    id = sprintf("s_%02d_%02d", g$x, g$y),
    x = as.numeric(g$x),
    y = as.numeric(g$y),
    stringsAsFactors = FALSE
  )
}

test_that("compute_spatial_niches_mrf validates inputs with spatial_niche_error", {
  cd <- .mrf_grid(5L, 5L)
  labs <- stats::setNames(ifelse(cd$x <= 2.5, "A", "B"), cd$id)

  e_beta <- .mrf_err(compute_spatial_niches_mrf(cd, labs, beta = -0.5))
  expect_true("spatial_niche_error" %in% e_beta$class)

  e_iter <- .mrf_err(compute_spatial_niches_mrf(cd, labs, max_iter = 0L))
  expect_true("spatial_niche_error" %in% e_iter$class)

  e_tol <- .mrf_err(compute_spatial_niches_mrf(cd, labs, tol = 0))
  expect_true("spatial_niche_error" %in% e_tol$class)

  e_small <- .mrf_err(compute_spatial_niches_mrf(cd[1:4, ], labs[1:4]))
  expect_true("spatial_niche_error" %in% e_small$class)

  e_one <- .mrf_err(compute_spatial_niches_mrf(cd, stats::setNames(rep("A", nrow(cd)), cd$id)))
  expect_true("spatial_niche_error" %in% e_one$class)
})

test_that("compute_spatial_niches_mrf enforces RANN dependency cleanly", {
  e <- new.env(parent = globalenv())
  e$requireNamespace <- function(package, ...) FALSE
  f <- get("compute_spatial_niches_mrf", envir = globalenv())
  environment(f) <- e

  cd <- .mrf_grid(5L, 5L)
  labs <- stats::setNames(ifelse(cd$x <= 2.5, "A", "B"), cd$id)
  err <- .mrf_err(f(cd, labs))
  expect_true("spatial_niche_error" %in% err$class)
  expect_match(err$msg, "RANN", fixed = TRUE)
})

test_that("compute_spatial_niches_mrf converges and produces spatially coherent niches", {
  cd <- .mrf_grid(10L, 10L)
  # Two clean spatial domains with deterministic salt-and-pepper noise
  labs <- ifelse(cd$x <= 5, "L2/3 IT", "Tumor Core")
  flip_idx <- seq(5L, nrow(cd), by = 11L)
  labs[flip_idx] <- ifelse(labs[flip_idx] == "L2/3 IT", "Tumor Core", "L2/3 IT")
  names(labs) <- cd$id

  res_mrf <- compute_spatial_niches_mrf(
    coords = cd,
    group_labels = labs,
    k_neighbors = 8L,
    n_niches = 2L,
    beta = 2.0,
    max_iter = 30L
  )

  expect_identical(res_mrf$method, "bayesian_mrf")
  expect_true(isTRUE(res_mrf$converged))
  expect_true(res_mrf$iterations >= 1L && res_mrf$iterations <= 30L)
  expect_equal(nrow(res_mrf$assignments), nrow(cd))
  expect_true(all(res_mrf$assignments$niche %in% c("N1", "N2")))
  expect_true(all(c("L2/3 IT", "Tumor Core") %in% colnames(res_mrf$niche_composition)))
  expect_true(is.finite(res_mrf$spatial_coherence) && res_mrf$spatial_coherence >= 0.85)
  expect_true(is.finite(res_mrf$posterior_entropy) && res_mrf$posterior_entropy >= 0)
})

test_that("mod_spatial_niche.R wires MRF method selector and beta input", {
  src <- paste(readLines(file.path(ts_project_root(), "modules/spatial/mod_spatial_niche.R"),
                         encoding = "UTF-8", warn = FALSE), collapse = "\n")
  expect_true(grepl("niche_method", src, fixed = TRUE))
  expect_true(grepl("mrf_beta", src, fixed = TRUE))
  expect_true(grepl("compute_spatial_niches_mrf", src, fixed = TRUE))
})

test_that("mod_spatial_viz.R uses tagList (never paste) around i18n$t on btn_add_to_report", {
  viz_src <- paste(readLines(file.path(ts_project_root(), "modules/spatial/mod_spatial_viz.R"),
                             encoding = "UTF-8", warn = FALSE), collapse = "\n")
  expect_false(grepl("paste\\([^)]*i18n\\$t\\(\"Ajouter cette vue au rapport\"", viz_src))
  expect_true(grepl("tagList\\([^)]*i18n\\$t\\(\"Ajouter cette vue au rapport\"", viz_src))
})

test_that("mod_spatial.R updates .ts_wf_last_panel and arms .spatial_sync_guard in .activate_spatial_step", {
  sp_src <- paste(readLines(file.path(ts_project_root(), "modules/spatial/mod_spatial.R"),
                            encoding = "UTF-8", warn = FALSE), collapse = "\n")
  expect_true(grepl("\\.spatial_sync_guard", sp_src))
  expect_true(grepl("\\.ts_wf_last_panel\\(panel\\)", sp_src))
})
