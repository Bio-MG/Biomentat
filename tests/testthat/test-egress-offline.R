# =============================================================================
# test-egress-offline.R — P1.3c : registre d'egress central & mode Offline strict
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/pathway_helpers.R")
source_project_file("modules/import/mod_geo.R")
source_project_file("modules/sc/mod_sc_annotation.R")

test_that("P1.3c : ts_egress_registry() recense les 4 chemins egress canoniques", {
  reg <- ts_egress_registry()
  expect_true(is.data.frame(reg))
  expect_true(all(c("egress_id", "module", "host", "opt_in_only") %in% colnames(reg)))
  expect_true(all(c("geo_ncbi", "celldex_hub", "enrichr_cloud", "renv_cran") %in% reg$egress_id))
  expect_true(all(nzchar(reg$host)))
  expect_true(all(nzchar(reg$module)))
})

test_that("P1.3c : ts_is_offline_strict() repond a option et variable d'environnement", {
  withr::local_options(list(biomentat.offline_strict = NULL))
  withr::local_envvar(c(BIOMENTAT_OFFLINE_STRICT = ""))
  expect_false(ts_is_offline_strict())

  withr::with_options(list(biomentat.offline_strict = TRUE), {
    expect_true(ts_is_offline_strict())
  })
  withr::with_envvar(c(BIOMENTAT_OFFLINE_STRICT = "1"), {
    expect_true(ts_is_offline_strict())
  })
  expect_false(ts_is_offline_strict())
})

test_that("P1.3c : ts_require_online() et les appelants levent biomentat_offline_error avant tout reseau", {
  withr::local_options(list(biomentat.offline_strict = FALSE))
  withr::local_envvar(c(BIOMENTAT_OFFLINE_STRICT = ""))
  expect_true(isTRUE(ts_require_online("enrichr_cloud")))

  withr::local_options(list(biomentat.offline_strict = TRUE))
  expect_error(ts_require_online("enrichr_cloud"), class = "biomentat_offline_error")
  expect_error(enrichr_cloud_submit(c("TP53", "GAPDH"), consent = TRUE),
               class = "biomentat_offline_error")
  expect_error(enrichr_cloud_enrich("GO_Biological_Process_2023", cookie_jar = tempfile()),
               class = "biomentat_offline_error")
  expect_error(enrichr_cloud_libraries(), class = "biomentat_offline_error")
  expect_error(.geo_fetch_supp("GSE123456", destdir = tempdir()),
               class = "biomentat_offline_error")
  expect_error(.geo_fetch("GSE123456"), class = "biomentat_offline_error")
  expect_error(.load_ref("hpca"), class = "biomentat_offline_error")
})
