# =============================================================================
# test-demo-loader.R — P3.3: 1-click offline demo dataset loader
# =============================================================================

source_project_file("tools/pkg_skeleton/R/demo_loader.R")

test_that("biomentat_demo_manifest lists bulk, sc, and spatial offline demos", {
  mf <- biomentat_demo_manifest()
  expect_true(is.data.frame(mf))
  expect_setequal(mf$modality, c("bulk", "sc", "spatial"))
  expect_true(all(mf$offline_only))
})

test_that("load_biomentat_demo rejects unknown modality", {
  err <- tryCatch(load_biomentat_demo("proteomics"), error = function(e) e)
  expect_s3_class(err, "demo_loader_error")
})

test_that("load_biomentat_demo generates deterministic offline fixtures for bulk, sc, and spatial", {
  b1 <- load_biomentat_demo("bulk", seed = 42L)
  b2 <- load_biomentat_demo("bulk", seed = 42L)
  expect_identical(b1, b2)
  expect_identical(dim(b1$counts), c(200L, 8L))
  expect_true(all(c("condition", "batch", "sample_id") %in% names(b1$metadata)))

  sc <- load_biomentat_demo("sc", seed = 42L)
  expect_identical(dim(sc$counts), c(100L, 60L))
  expect_true(inherits(sc$counts, "dgCMatrix"))
  expect_true(all(c("sample_id", "condition", "cluster") %in% names(sc$metadata)))

  sp <- load_biomentat_demo("spatial", seed = 42L)
  expect_identical(dim(sp$counts), c(80L, 50L))
  expect_true(all(c("id", "x", "y", "fov") %in% names(sp$coords)))
  expect_identical(sp$coords$id, colnames(sp$counts))
})
