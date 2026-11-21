# =============================================================================
# test-mod-sc-viz-umap-budget.R -- W1.1 #13 : câblage du budget UMAP interactif
# dans modules/sc/mod_sc_viz.R (sc_umap_render_budget + sc_subsample_umap_indices)
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
load_project_pkg()
source_project_file("tools/pkg_skeleton/R/sc_plotting.R")

.UB_FILE <- "modules/sc/mod_sc_viz.R"

.ub_src <- function() {
  paste(readLines(file.path(ts_project_root(), .UB_FILE), warn = FALSE, encoding = "UTF-8"),
        collapse = "\n")
}

# Extrait le helper top-level sans sourcer tout le module (isolation de shard).
.ub_load_helper <- function() {
  exprs <- parse(file = file.path(ts_project_root(), .UB_FILE), encoding = "UTF-8")
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        identical(as.character(e[[2]]), ".sc_viz_umap_budget_idx")) {
      return(eval(e[[3]], envir = globalenv()))
    }
  }
  NULL
}

test_that("le module référence les deux helpers gardés par exists(mode='function')", {
  src <- .ub_src()
  expect_match(src, "sc_umap_render_budget", fixed = TRUE)
  expect_match(src, "sc_subsample_umap_indices", fixed = TRUE)
  expect_match(src, 'exists("sc_umap_render_budget", mode = "function")', fixed = TRUE)
  expect_match(src, 'exists("sc_subsample_umap_indices", mode = "function")', fixed = TRUE)
})

test_that(".preview_obj consomme le helper de budget", {
  expect_match(.ub_src(), ".sc_viz_umap_budget_idx(ncol(obj)", fixed = TRUE)
})

test_that("<= 25000 cellules : dataset complet conservé (NULL)", {
  f <- .ub_load_helper()
  expect_true(is.function(f))
  expect_null(f(100L))
  expect_null(f(25000L, groups = rep("a", 25000L)))
})

test_that("> 25000 cellules : 25000 indices, triés, déterministes", {
  f <- .ub_load_helper()
  n <- 60000L
  idx <- f(n)
  expect_length(idx, 25000L)
  expect_false(is.unsorted(idx))
  expect_identical(idx, f(n))
})

test_that("sous-échantillonnage stratifié préserve les clusters rares", {
  f <- .ub_load_helper()
  n <- 40000L
  g <- c(rep("big", n - 3L), "rare1", "rare2", "rare2")
  idx <- f(n, groups = g)
  expect_length(idx, 25000L)
  expect_true(all(c("rare1", "rare2") %in% g[idx]))
})

test_that("0 nouvelle clé i18n (.tr / i18n$t) dans le câblage", {
  src <- .ub_src()
  blk <- regmatches(src, regexpr("(?s)\\.sc_viz_umap_budget_idx <- function.*?\\r?\\n\\}\\r?\\n", src, perl = TRUE))
  expect_length(blk, 1L)
  expect_false(grepl("\\.tr\\(|i18n\\$t\\(", blk))
})
