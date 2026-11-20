# =============================================================================
# test-sc-communication-spatial-self.R — W1.3 Lot 1 : honnêteté des
# auto-boucles (autocrine / self) dans le contexte spatial 2D
# =============================================================================

if (!exists("source_project_file", mode = "function")) {
  test_dir <- if (exists("test_path", mode = "function")) testthat::test_path() else "tests/testthat"
  helper_src <- file.path(test_dir, "helper-source.R")
  if (file.exists(helper_src)) source(helper_src)
}

for (.dep in c(
  "tools/pkg_skeleton/R/io_helpers.R",
  "tools/pkg_skeleton/R/rdata_io.R",
  "tools/pkg_skeleton/R/state.R",
  "tools/pkg_skeleton/R/provenance.R",
  "tools/pkg_skeleton/R/sc_velocity.R",
  "tools/pkg_skeleton/R/sc_communication.R",
  "tools/pkg_skeleton/R/sc_communication_views.R",
  "tools/pkg_skeleton/R/sc_communication_spatial.R"
)) {
  source_project_file(.dep)
}
rm(.dep)

test_that("helpers internes .spc_is_autocrine et .spc_distance_label identifient les auto-boucles", {
  expect_true(exists(".spc_is_autocrine", mode = "function"))
  expect_true(exists(".spc_distance_label", mode = "function"))

  expect_identical(
    .spc_is_autocrine(c("C", "A", "B"), c("C", "B", "D"), c(0, 10, 0)),
    c(TRUE, FALSE, TRUE)
  )
  expect_identical(
    .spc_distance_label(c(0, 10, 1.5, NA_real_), c("C", "A", "A", "A"), c("C", "B", "C", "D")),
    c("0 (autocrine / self)", "10", "1.5", NA_character_)
  )
  expect_identical(.spc_distance_label(0), "0 (autocrine / self)")
})

test_that("le contexte spatial préserve les distances numériques figées pour les auto-boucles", {
  ctx <- build_communication_spatial_context(
    .ccx_result(), .ccx_coords(), .ccx_identities(),
    identity_column = "cell_type", coordinate_source = "spatial_object"
  )
  pt <- ctx$pair_table
  cc <- pt[pt$sender_node == "C" & pt$receiver_node == "C", ]
  expect_equal(nrow(cc), 1L)
  expect_true(is.numeric(cc$centroid_distance))
  expect_equal(cc$centroid_distance, 0)
  expect_true(is.numeric(cc$median_recv_to_sender_nn))
  expect_equal(cc$median_recv_to_sender_nn, 0)
  expect_true(is.numeric(cc$median_send_to_recv_nn))
  expect_equal(cc$median_send_to_recv_nn, 0)
})

test_that("l'export spatial annote is_autocrine et distance_label = '0 (autocrine / self)'", {
  ctx <- build_communication_spatial_context(
    .ccx_result(), .ccx_coords(), .ccx_identities(),
    identity_column = "cell_type", coordinate_source = "spatial_object"
  )
  ex <- build_communication_spatial_export(ctx)
  expect_true(all(c("is_autocrine", "distance_label") %in% names(ex)))
  expect_type(ex$is_autocrine, "logical")
  expect_type(ex$distance_label, "character")

  cc <- ex[ex$sender_node == "C" & ex$receiver_node == "C", ]
  ab <- ex[ex$sender_node == "A" & ex$receiver_node == "B", ]
  expect_identical(cc$is_autocrine, TRUE)
  expect_identical(cc$distance_label, "0 (autocrine / self)")
  expect_identical(ab$is_autocrine, FALSE)
  expect_identical(ab$distance_label, "10")
})

test_that("plot_communication_spatial_edges étiquette honnêtement une auto-boucle", {
  ctx <- build_communication_spatial_context(
    .ccx_result(), .ccx_coords(), .ccx_identities(),
    identity_column = "cell_type", coordinate_source = "spatial_object"
  )
  p <- plot_communication_spatial_edges(ctx, "C", "C", .ccx_coords(), .ccx_identities())
  expect_s3_class(p, "ggplot")
  expect_true(grepl("0 (autocrine / self)", p$labels$subtitle %||% "", fixed = TRUE))
})

test_that("plot_communication_spatial_distance_summary annote les paires autocrines", {
  ctx <- build_communication_spatial_context(
    .ccx_result(), .ccx_coords(), .ccx_identities(),
    identity_column = "cell_type", coordinate_source = "spatial_object"
  )
  p <- plot_communication_spatial_distance_summary(ctx)
  expect_s3_class(p, "ggplot")
  expect_true(any(grepl("autocrine / self", p$data$pair_label, fixed = TRUE)))
})
