# =============================================================================
# test-pk3-convergence-caveat.R — W2.3 PK-3 methodological caveat in UI
# =============================================================================

source_project_file("tools/pkg_skeleton/R/pathway_helpers.R")

test_that("packs_convergence_content displays the PK-3 global BH post-cutoff caveat", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("DT")
  suppressPackageStartupMessages(library(shiny))
  suppressPackageStartupMessages(library(DT))

  conv_rx <- function() list(
    ok = TRUE,
    tbl = data.frame(
      cluster_id = 1L,
      term = "GO:0001",
      library = "GO_BP",
      p_adj_global = 0.01,
      stringsAsFactors = FALSE
    ),
    msg = NULL
  )

  ui <- packs_convergence_content(NS("p"), conv_rx, function(x) x)
  html <- as.character(ui)
  expect_match(html, "PK-3", fixed = TRUE)
  expect_match(html, "post-cutoff", fixed = TRUE)
})
