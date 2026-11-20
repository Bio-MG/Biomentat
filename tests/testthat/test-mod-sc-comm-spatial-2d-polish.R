# =============================================================================
# test-mod-sc-comm-spatial-2d-polish.R — W1.3 Lot 1 CCC Spatial 2D polish
# =============================================================================
# En mode Spatial 2D (`spc_source == "reduction"`), `spc_radius` et
# `spc_nperm` doivent être grisés (`shinyjs::disable`) et `spc_table` doit
# masquer les colonnes 100 % NA (tout en conservant `scrollY` / `scrollCollapse`).
# =============================================================================

test_that("W1.3 : mod_sc_communication_spatial.R grise spc_radius/spc_nperm en mode 2D et filtre les colonnes 100% NA", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "modules/sc/mod_sc_communication_spatial.R"),
                         warn = FALSE), collapse = "\n")
  expect_match(src, 'shinyjs::disable\\("spc_radius"\\)')
  expect_match(src, 'shinyjs::disable\\("spc_nperm"\\)')
  expect_match(src, 'shinyjs::enable\\("spc_radius"\\)')
  expect_match(src, 'shinyjs::enable\\("spc_nperm"\\)')
  expect_match(src, '\\.spc_drop_all_na_cols\\s*<-\\s*function')
  expect_match(src, '\\.spc_drop_all_na_cols\\(ctx\\$pair_table\\)')
})

test_that("W1.3 : .spc_drop_all_na_cols retire exactement les colonnes 100% NA", {
  src_lines <- readLines(file.path(ts_project_root(),
                                   "modules/sc/mod_sc_communication_spatial.R"),
                         warn = FALSE)
  idx <- grep("^\\s*\\.spc_drop_all_na_cols\\s*<-\\s*function", src_lines)[1]
  expect_false(is.na(idx))
  env <- new.env(parent = globalenv())
  eval(parse(text = src_lines[idx:(idx + 4L)]), envir = env)

  df <- data.frame(
    sender = c("A", "B"),
    dist_mean = c(1.2, 2.5),
    radius_prop = c(NA_real_, NA_real_),
    perm_pval = c(NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  out <- env$.spc_drop_all_na_cols(df)
  expect_identical(colnames(out), c("sender", "dist_mean"))
})
