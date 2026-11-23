# ==============================================================================
# test-bk09-ome-ngff-webgl.R — BK-09 OME-NGFF/OpenSeadragon & WebGL scattergl
# ==============================================================================

.bk09_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("read_ome_ngff_metadata and build_openseadragon_tile_config parse 3-level pyramid", {
  source(.bk09_root("tools/pkg_skeleton/R/spatial_hd_ingest.R"), local = globalenv())

  zdir <- tempfile("ome_ngff_")
  dir.create(zdir, recursive = TRUE)
  on.exit(unlink(zdir, recursive = TRUE), add = TRUE)

  zattrs <- list(
    multiscales = list(list(
      version = "0.4",
      name = "histology_pyramid",
      axes = list(list(name = "c", type = "channel"),
                  list(name = "y", type = "space"),
                  list(name = "x", type = "space")),
      datasets = list(
        list(path = "0", coordinateTransformations = list(list(type = "scale", scale = list(1, 0.5, 0.5)))),
        list(path = "1", coordinateTransformations = list(list(type = "scale", scale = list(1, 1.0, 1.0)))),
        list(path = "2", coordinateTransformations = list(list(type = "scale", scale = list(1, 2.0, 2.0))))
      )
    ))
  )
  jsonlite::write_json(zattrs, file.path(zdir, ".zattrs"), auto_unbox = TRUE)

  meta <- read_ome_ngff_metadata(zdir)
  expect_s3_class(meta, "biomentat_ome_ngff")
  expect_equal(meta$n_levels, 3L)
  expect_equal(meta$levels$path, c("0", "1", "2"))

  osd <- build_openseadragon_tile_config(meta, base_url = "/zarr/sample1", tile_size = 256L,
                                         base_width = 4096L, base_height = 2048L)
  expect_equal(osd$type, "ome-ngff-pyramid")
  expect_equal(osd$tileSize, 256L)
  expect_equal(osd$maxLevel, 2L)
  expect_equal(osd$levels[[1L]]$width, 4096L)
  expect_equal(osd$levels[[2L]]$width, 2048L)
  expect_equal(osd$levels[[3L]]$width, 1024L)
  expect_equal(osd$levels[[1L]]$url, "/zarr/sample1/0")
})

test_that("read_ome_ngff_metadata errors cleanly on missing directory or .zattrs", {
  source(.bk09_root("tools/pkg_skeleton/R/spatial_hd_ingest.R"), local = globalenv())

  expect_error(read_ome_ngff_metadata(tempfile()), class = "spatial_hd_ingest_error")
  empty_d <- tempfile("empty_zarr_")
  dir.create(empty_d, recursive = TRUE)
  on.exit(unlink(empty_d, recursive = TRUE), add = TRUE)
  expect_error(read_ome_ngff_metadata(empty_d), class = "spatial_hd_ingest_error")
})

test_that("mod_sc_viz uses plotly::toWebGL for high-N interactive rendering", {
  src <- paste(readLines(.bk09_root("modules/sc/mod_sc_viz.R"), warn = FALSE), collapse = "\n")
  expect_true(grepl("plotly::toWebGL(", src, fixed = TRUE))
})
