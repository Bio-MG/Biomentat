# =============================================================================
# test-sc-ap-population-rarity.R — W1.1 #11 : opt-in Population Rarity in SC
# auto-pipeline (sc_ap_normalize_params + run_sc_auto_pipeline_core + mod_sc.R)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/provenance.R")
source_project_file("tools/pkg_skeleton/R/sc_helpers.R")
source_project_file("tools/pkg_skeleton/R/sc_bpcells.R")
source_project_file("tools/pkg_skeleton/R/sc_trajectory.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
source_project_file("tools/pkg_skeleton/R/sc_population_rarity.R")
source_project_file("tools/pkg_skeleton/R/memory.R")
source_project_file("tools/pkg_skeleton/R/sc_pipeline.R")
assign(".tr", function(key) key, envir = globalenv())

test_that("sc_ap_normalize_params: population rarity is opt-in FALSE with explicit rule/threshold", {
  p0 <- sc_ap_normalize_params(list())
  expect_true(isFALSE(p0$sc_ap_population_rarity))
  expect_identical(p0$sc_ap_rarity_rule, "absolute_n_cells")
  expect_true(is.na(p0$sc_ap_rarity_threshold))

  p1 <- sc_ap_normalize_params(list(
    sc_ap_population_rarity = TRUE,
    sc_ap_rarity_rule = "relative_fraction",
    sc_ap_rarity_threshold = 0.05
  ))
  expect_true(isTRUE(p1$sc_ap_population_rarity))
  expect_identical(p1$sc_ap_rarity_rule, "relative_fraction")
  expect_equal(p1$sc_ap_rarity_threshold, 0.05)
})

test_that("mod_sc.R modal exposes sc_ap_population_rarity, sc_ap_rarity_rule (C13), sc_ap_rarity_threshold, and completion wiring", {
  mod_path <- file.path(ts_project_root(), "modules/sc/mod_sc.R")
  ln <- readLines(mod_path, warn = FALSE, encoding = "UTF-8")
  txt <- paste(ln, collapse = "\n")

  expect_true(grepl('ns_m("sc_ap_population_rarity")', txt, fixed = TRUE))
  expect_true(grepl('ns_m("sc_ap_rarity_rule")', txt, fixed = TRUE))
  expect_true(grepl('ns_m("sc_ap_rarity_threshold")', txt, fixed = TRUE))
  expect_true(grepl('"absolute_n_cells"', txt, fixed = TRUE))
  expect_true(grepl('"relative_fraction"', txt, fixed = TRUE))
  expect_true(grepl('state_set(shared_rv, "population_rarity_result", res$population_rarity_result)',
                    txt, fixed = TRUE))

  pipe_path <- file.path(ts_project_root(), "tools/pkg_skeleton/R/sc_pipeline.R")
  pipe_txt <- paste(readLines(pipe_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_true(grepl('state_set(shared_rv, "population_rarity_result", res$population_rarity_result)',
                    pipe_txt, fixed = TRUE))
})

skip_if_not_installed("Seurat")
library(Seurat)

.make_rarity_ap_seurat <- function(n_cells = 80L, seed = 101) {
  set.seed(seed)
  n_genes <- 80L
  mat <- matrix(rpois(n_genes * n_cells, lambda = 2), nrow = n_genes,
                dimnames = list(paste0("G", seq_len(n_genes)),
                                paste0("C", seq_len(n_cells))))
  half <- floor(n_cells / 2L)
  mat[1:15, seq_len(half)] <- mat[1:15, seq_len(half)] + rpois(15L * half, 8)
  mat[16:30, (half + 1L):n_cells] <- mat[16:30, (half + 1L):n_cells] + rpois(15L * (n_cells - half), 8)
  suppressWarnings(CreateSeuratObject(counts = mat))
}

.base_ap_params <- list(
  sc_ap_min_gene     = 5,
  sc_ap_max_gene     = 8000,
  sc_ap_mt           = 50,
  sc_ap_pca_dim      = 5,
  sc_ap_res          = 0.5,
  sc_ap_compute_umap = FALSE,
  sc_ap_markers      = FALSE,
  sc_ap_pathway      = FALSE,
  sc_ap_correlation  = FALSE,
  sc_ap_trajectory   = FALSE,
  sc_ap_singler      = FALSE
)

test_that("run_sc_auto_pipeline_core skips population rarity cleanly when disabled", {
  obj <- .make_rarity_ap_seurat(n_cells = 80L)
  res <- suppressMessages(run_sc_auto_pipeline_core(obj, .base_ap_params))
  expect_null(res$population_rarity_result)
})

test_that("run_sc_auto_pipeline_core computes population_rarity_result when opted in", {
  obj <- .make_rarity_ap_seurat(n_cells = 80L)
  params <- c(.base_ap_params, list(
    sc_ap_population_rarity = TRUE,
    sc_ap_rarity_rule       = "absolute_n_cells",
    sc_ap_rarity_threshold  = 20
  ))
  res <- suppressMessages(run_sc_auto_pipeline_core(obj, params))
  expect_false(is.null(res$population_rarity_result))
  expect_identical(res$population_rarity_result$type, "sc_population_rarity")
  expect_identical(res$population_rarity_result$identity_column, "seurat_clusters")
  expect_identical(res$population_rarity_result$rarity_rule$rule_type, "absolute_n_cells")
  expect_equal(res$population_rarity_result$rarity_rule$threshold, 20L)
  expect_true(is.data.frame(res$population_rarity_result$population_table))
})

test_that("run_sc_auto_pipeline_core logs classified warning instead of aborting on invalid threshold or < 50 cells", {
  obj_small <- .make_rarity_ap_seurat(n_cells = 30L)
  params_small <- c(.base_ap_params, list(
    sc_ap_population_rarity = TRUE,
    sc_ap_rarity_rule       = "absolute_n_cells",
    sc_ap_rarity_threshold  = 5
  ))
  warned_small <- FALSE
  res_small <- withCallingHandlers(
    suppressMessages(run_sc_auto_pipeline_core(obj_small, params_small)),
    sc_ap_population_rarity_warning = function(w) {
      warned_small <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_true(warned_small)
  expect_null(res_small$population_rarity_result)
  expect_s4_class(res_small$obj, "Seurat")

  obj_ok <- .make_rarity_ap_seurat(n_cells = 80L)
  params_bad_thr <- c(.base_ap_params, list(
    sc_ap_population_rarity = TRUE,
    sc_ap_rarity_rule       = "absolute_n_cells",
    sc_ap_rarity_threshold  = NA_real_
  ))
  warned_thr <- FALSE
  res_bad_thr <- withCallingHandlers(
    suppressMessages(run_sc_auto_pipeline_core(obj_ok, params_bad_thr)),
    sc_ap_population_rarity_warning = function(w) {
      warned_thr <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_true(warned_thr)
  expect_null(res_bad_thr$population_rarity_result)
  expect_s4_class(res_bad_thr$obj, "Seurat")
})
