# =============================================================================
# test-sc-h5ad-legacy-csc.R — Legacy CSC/CSR .h5ad fallback reader (W2.4 §2.7)
# =============================================================================

source_project_file("modules/import/mod_import_sc.R")

.write_test_h5ad <- function(path, mat_cxg, encoding = "csc_matrix",
                             obs_ids = rownames(mat_cxg),
                             var_ids = colnames(mat_cxg)) {
  f <- hdf5r::H5File$new(path, mode = "w")
  on.exit(try(f$close_all(), silent = TRUE), add = TRUE)

  x_grp <- f$create_group("X")
  if (identical(encoding, "csc_matrix")) {
    sp <- methods::as(Matrix::Matrix(mat_cxg, sparse = TRUE), "dgCMatrix")
    x_grp[["data"]] <- sp@x
    x_grp[["indices"]] <- as.integer(sp@i)
    x_grp[["indptr"]] <- as.integer(sp@p)
  } else {
    sp_r <- methods::as(Matrix::t(Matrix::Matrix(mat_cxg, sparse = TRUE)), "dgCMatrix")
    x_grp[["data"]] <- sp_r@x
    x_grp[["indices"]] <- as.integer(sp_r@i)
    x_grp[["indptr"]] <- as.integer(sp_r@p)
  }
  x_grp$create_attr("encoding-type", encoding)
  x_grp$create_attr("shape", as.integer(dim(mat_cxg)))

  obs_grp <- f$create_group("obs")
  obs_grp[["_index"]] <- obs_ids
  var_grp <- f$create_group("var")
  var_grp[["_index"]] <- var_ids
  invisible(path)
}

test_that(".read_legacy_h5ad_hdf5r reads legacy CSC .h5ad into genes x cells dgCMatrix", {
  skip_if_not_installed("hdf5r")
  mat_cxg <- matrix(c(1, 0, 2, 0, 3, 4, 0, 5, 0, 0, 6, 7), nrow = 4L, ncol = 3L,
                    dimnames = list(paste0("Cell", 1:4), paste0("Gene", 1:3)))
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  .write_test_h5ad(tmp, mat_cxg, encoding = "csc_matrix")

  res <- .read_legacy_h5ad_hdf5r(tmp)
  expect_s4_class(res, "dgCMatrix")
  expect_identical(dim(res), c(3L, 4L))
  expect_identical(rownames(res), paste0("Gene", 1:3))
  expect_identical(colnames(res), paste0("Cell", 1:4))
  expect_equal(as.matrix(res), t(mat_cxg))
})

test_that(".read_legacy_h5ad_hdf5r reads CSR .h5ad into genes x cells dgCMatrix", {
  skip_if_not_installed("hdf5r")
  mat_cxg <- matrix(c(0, 2, 0, 4, 5, 0, 7, 0, 9, 10, 0, 12), nrow = 4L, ncol = 3L,
                    dimnames = list(paste0("C", 1:4), paste0("G", 1:3)))
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp), add = TRUE)
  .write_test_h5ad(tmp, mat_cxg, encoding = "csr_matrix")

  res <- .read_legacy_h5ad_hdf5r(tmp)
  expect_s4_class(res, "dgCMatrix")
  expect_identical(dim(res), c(3L, 4L))
  expect_identical(rownames(res), paste0("G", 1:3))
  expect_identical(colnames(res), paste0("C", 1:4))
  expect_equal(as.matrix(res), t(mat_cxg))
})
