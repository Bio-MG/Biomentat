# =============================================================================
# test-demo-loader.R — P3.3: 1-click offline demo dataset loader
# =============================================================================

source_project_file("tools/pkg_skeleton/R/demo_loader.R")
suppressPackageStartupMessages(library(Matrix))

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

test_that("load_biomentat_demo loads compact real offline fixtures for bulk, sc, and spatial", {
  b1 <- load_biomentat_demo("bulk")
  expect_identical(dim(b1$counts), c(2000L, 6L))
  expect_false(any(grepl("^Gene_", rownames(b1$counts))))
  expect_true(all(c("sample_id", "condition") %in% names(b1$metadata)))
  expect_true(isTRUE(b1$provenance$real_data))

  sc <- load_biomentat_demo("sc")
  expect_identical(dim(sc$counts), c(1500L, 400L))
  expect_true(inherits(sc$counts, "dgCMatrix"))
  expect_false(any(grepl("^Gene_", rownames(sc$counts))))
  expect_true(all(c("sample_id", "condition", "celltype", "cluster") %in% names(sc$metadata)))
  expect_true(isTRUE(sc$provenance$real_data))

  sp <- load_biomentat_demo("spatial")
  expect_identical(dim(sp$counts), c(1500L, 400L))
  expect_true(inherits(sp$counts, "dgCMatrix"))
  expect_true(all(c("id", "x", "y", "fov") %in% names(sp$coords)))
  expect_identical(sp$coords$id, colnames(sp$counts))
  expect_true(isTRUE(sp$provenance$real_data))
})
