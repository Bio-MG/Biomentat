# =============================================================================
# test-spatial-progress-ribbon.R — Compact Spatial pipeline progress ribbon
# (PROPOSAL_UX_SPREAD_SPATIAL_SLIDE_PANEL_2026-11-19.md Part B)
# =============================================================================

suppressPackageStartupMessages(library(shiny))

source_project_file("modules/spatial/mod_spatial.R")

test_that(".spatial_pipeline_progress_items handles NULL, empty, and empty-container states", {
  expected_codes <- c("qc", "moran", "cluster", "deconv", "niche", "stats")

  for (st in list(
    NULL,
    list(),
    list(
      qc_metrics        = data.frame(),
      moran_results     = data.frame(),
      cluster_labels    = character(0),
      deconv_props      = list(),
      niche_composition = data.frame(),
      hotspot_result    = data.frame(),
      ripley_result     = list()
    )
  )) {
    items <- .spatial_pipeline_progress_items(st)
    expect_type(items, "list")
    expect_length(items, 6L)

    codes <- vapply(items, `[[`, character(1), "code")
    expect_identical(codes, expected_codes)

    dones <- vapply(items, `[[`, logical(1), "done")
    expect_true(all(!dones))

    syms <- vapply(items, `[[`, character(1), "symbol")
    expect_true(all(syms == "\u26aa"))

    expect_true(all(nzchar(vapply(items, `[[`, character(1), "step"))))
    expect_true(all(nzchar(vapply(items, `[[`, character(1), "short"))))
    expect_true(all(nzchar(vapply(items, `[[`, character(1), "label"))))
  }
})

test_that(".spatial_pipeline_progress_items reflects partial progress accurately", {
  items <- .spatial_pipeline_progress_items(list(
    qc_metrics     = data.frame(nCount = c(100, 200)),
    cluster_labels = factor(c("1", "2")),
    deconv_props   = NULL
  ))
  by_code <- stats::setNames(items, vapply(items, `[[`, character(1), "code"))

  expect_true(by_code$qc$done)
  expect_identical(by_code$qc$symbol, "\u2705")
  expect_false(by_code$moran$done)
  expect_identical(by_code$moran$symbol, "\u26aa")
  expect_true(by_code$cluster$done)
  expect_identical(by_code$cluster$symbol, "\u2705")
  expect_false(by_code$deconv$done)
  expect_false(by_code$niche$done)
  expect_false(by_code$stats$done)

  # Stats step completes when EITHER hotspot_result OR ripley_result is present
  hs_items <- .spatial_pipeline_progress_items(list(
    hotspot_result = data.frame(id = "s1", gi_star = 2.1)
  ))
  expect_true(hs_items[[6]]$done)
  expect_identical(hs_items[[6]]$symbol, "\u2705")

  rp_items <- .spatial_pipeline_progress_items(list(
    ripley_result = list(curve = data.frame(r = 1, obs = 1.2))
  ))
  expect_true(rp_items[[6]]$done)
  expect_identical(rp_items[[6]]$symbol, "\u2705")
})

test_that(".spatial_pipeline_progress_items reflects full pipeline completion", {
  full_state <- list(
    qc_metrics        = data.frame(nCount = 100),
    moran_results     = data.frame(gene = "G1", I = 0.5),
    cluster_labels    = c("C1", "C2"),
    deconv_props      = matrix(c(0.4, 0.6), nrow = 1),
    niche_composition = data.frame(niche = 1, frac = 1),
    hotspot_result    = data.frame(id = "s1", hotspot = "hot"),
    ripley_result     = list(curve = data.frame(r = 1))
  )
  items <- .spatial_pipeline_progress_items(full_state)
  expect_true(all(vapply(items, `[[`, logical(1), "done")))
  expect_true(all(vapply(items, `[[`, character(1), "symbol") == "\u2705"))
})

test_that("mod_spatial_ui renders spatial_progress_ribbon_ui inside .ts-spatial-side-head", {
  suppressPackageStartupMessages(library(bslib))

  had_i18n <- exists("i18n", envir = globalenv(), inherits = FALSE)
  old_i18n <- if (had_i18n) get("i18n", envir = globalenv(), inherits = FALSE) else NULL
  had_trp  <- exists(".tr_plain", envir = globalenv(), inherits = FALSE)
  old_trp  <- if (had_trp) get(".tr_plain", envir = globalenv(), inherits = FALSE) else NULL
  on.exit({
    if (had_i18n) assign("i18n", old_i18n, envir = globalenv())
    else if (exists("i18n", envir = globalenv(), inherits = FALSE)) rm("i18n", envir = globalenv())
    if (had_trp) assign(".tr_plain", old_trp, envir = globalenv())
    else if (exists(".tr_plain", envir = globalenv(), inherits = FALSE)) rm(".tr_plain", envir = globalenv())
  }, add = TRUE)
  assign("i18n", list(t = function(x) x), envir = globalenv())
  assign(".tr_plain", function(x) x, envir = globalenv())

  load_project_pkg()
  source_project_file("tools/pkg_skeleton/R/io_helpers.R")
  source_project_file("tools/pkg_skeleton/R/spatial_deconv_prep.R")
  for (f in c(
    "modules/spatial/mod_spatial_pipeline.R",
    "modules/spatial/mod_spatial_cluster.R",
    "modules/spatial/deconv/mod_spatial_deconv_ui.R",
    "modules/spatial/deconv/mod_spatial_deconv.R",
    "modules/spatial/mod_spatial_viz.R",
    "modules/spatial/mod_spatial_multi.R",
    "modules/spatial/mod_spatial_niche.R",
    "modules/spatial/mod_spatial_export.R",
    "modules/spatial/mod_spatial_report.R",
    "modules/spatial/mod_spatial_qc.R",
    "modules/spatial/mod_spatial.R"
  )) {
    source_project_file(f)
  }

  html <- as.character(mod_spatial_ui("spatial"))
  expect_true(grepl('id="spatial-spatial_progress_ribbon_ui"', html, fixed = TRUE))

  # Verify the ribbon uiOutput sits inside .ts-spatial-side-head (before the
  # <hr/> separator and #spatial-steps accordion that follow it in the sidebar)
  pos_head <- regexpr('class="ts-spatial-side-head"', html, fixed = TRUE)[1]
  expect_gt(pos_head, 0L)
  after_head <- substr(html, pos_head, nchar(html))
  pos_ribbon <- regexpr('id="spatial-spatial_progress_ribbon_ui"', after_head, fixed = TRUE)[1]
  pos_hr     <- regexpr('<hr/>', after_head, fixed = TRUE)[1]
  pos_steps  <- regexpr('id="spatial-steps"', after_head, fixed = TRUE)[1]

  expect_gt(pos_ribbon, 0L)
  expect_gt(pos_hr, pos_ribbon)
  expect_gt(pos_steps, pos_hr)
})
