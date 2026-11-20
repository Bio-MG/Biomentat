# =============================================================================
# test-mod-sc-viz-clip-quantiles.R — W2.3 BK-04 UI wiring in mod_sc_viz.R
# =============================================================================
# Le moteur `sc_build_plot()` (sc_plotting.R:75-76) supporte déjà
# `cfg$clip_quantiles = TRUE` (q1-q99) pour FeaturePlot (BK-04), mais
# l'interface `mod_sc_viz.R` n'exposait encore aucun contrôle ni câblage cfg.
# =============================================================================

test_that("W2.3 BK-04 : mod_sc_viz.R expose clip_quantiles (UI + .current_cfg + plot_interactive)", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "modules/sc/mod_sc_viz.R"),
                         warn = FALSE), collapse = "\n")
  expect_match(src, 'checkboxInput\\(ns\\("clip_quantiles"\\)')
  expect_match(src, 'clip_quantiles\\s*=\\s*input\\$clip_quantiles')
  expect_match(src, '"clip_quantiles"')
})
