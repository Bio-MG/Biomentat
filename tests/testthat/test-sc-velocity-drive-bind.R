# =============================================================================
# test-sc-velocity-drive-bind.R
#   W1.1 #2: drive bindings for `velocity_validate` and `velocity_scvelo_run`
#   (allowlist + ts_drive_bind_button + ts_drive_publish_token + readiness).
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
source_project_file("modules/sc/mod_sc_velocity.R")

test_that("velocity_validate and velocity_scvelo_run are allowlisted for sc_velocity", {
  expect_true(exists("TS_DRIVE_SC_VELOCITY_VALIDATE_BUTTON"))
  expect_true(exists("TS_DRIVE_SC_VELOCITY_SCVELO_BUTTON"))
  expect_true(exists("TS_DRIVE_SC_VELOCITY_MODULE"))

  expect_identical(TS_DRIVE_SC_VELOCITY_VALIDATE_BUTTON, "sc-velocity-velocity_validate")
  expect_identical(TS_DRIVE_SC_VELOCITY_SCVELO_BUTTON, "sc-velocity-velocity_scvelo_run")
  expect_identical(TS_DRIVE_SC_VELOCITY_MODULE, "sc_velocity")

  expect_true(TS_DRIVE_SC_VELOCITY_MODULE %in% TS_DRIVE_MODULES)
  expect_true(TS_DRIVE_SC_VELOCITY_VALIDATE_BUTTON %in% TS_DRIVE_BUTTONS)
  expect_true(TS_DRIVE_SC_VELOCITY_SCVELO_BUTTON %in% TS_DRIVE_BUTTONS)

  for (btn in c(TS_DRIVE_SC_VELOCITY_VALIDATE_BUTTON, TS_DRIVE_SC_VELOCITY_SCVELO_BUTTON)) {
    expect_true(ts_drive_allowlisted(btn), info = btn)
    expect_identical(TS_DRIVE_ALLOWLIST[[btn]]$kind, "button", info = btn)
    expect_identical(ts_drive_module_of(btn), TS_DRIVE_SC_VELOCITY_MODULE, info = btn)
    expect_identical(ts_drive_button_module(btn), TS_DRIVE_SC_VELOCITY_MODULE, info = btn)
  }

  expect_true(isTRUE(ts_drive_allowlist_problems()))

  wrong_mod <- setNames(list(list(kind = "button", module = "sc_pipeline")),
                        TS_DRIVE_SC_VELOCITY_VALIDATE_BUTTON)
  expect_true(any(grepl("must declare module",
                        ts_drive_allowlist_problems(wrong_mod, character(0)))))
})

test_that("ts_drive_bind_button accepts both velocity buttons", {
  expect_identical(
    ts_drive_bind_button("sc-velocity-velocity_validate"),
    "tsdrive-token-sc-velocity-velocity_validate"
  )
  expect_identical(
    ts_drive_bind_button("sc-velocity-velocity_scvelo_run"),
    "tsdrive-token-sc-velocity-velocity_scvelo_run"
  )
})

test_that("mod_sc_velocity_server publishes tokens, readiness guards, and triggers on drive counter", {
  genes <- paste0("G", 1:4)
  cells <- paste0("Cell", 1:6)
  set.seed(42)
  sp <- matrix(rpois(24, 5), nrow = 4, ncol = 6, dimnames = list(genes, cells))
  un <- matrix(rpois(24, 2), nrow = 4, ncol = 6, dimnames = list(genes, cells))
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts = sp))
  umap_coords <- matrix(seq_len(12) * 0.1, nrow = 6, ncol = 2,
                        dimnames = list(cells, c("UMAP_1", "UMAP_2")))
  seu[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap_coords,
                                                key = "UMAP_", assay = "RNA")

  rds_path <- tempfile(fileext = ".rds")
  on.exit(unlink(rds_path), add = TRUE)
  saveRDS(list(spliced = sp, unspliced = un), rds_path)

  gd <- shiny::reactiveValues(
    sc_obj = NULL,
    language = "fr",
    i18n = NULL,
    drive_registry = new.env(parent = emptyenv())
  )
  srv <- shiny::reactiveValues()

  shiny::testServer(mod_sc_velocity_server, args = list(global_data = gd, shared_rv = srv), {
    reg <- ts_drive_registry(gd)
    val_entry <- reg[[TS_DRIVE_SC_VELOCITY_VALIDATE_BUTTON]]
    scv_entry <- reg[[TS_DRIVE_SC_VELOCITY_SCVELO_BUTTON]]

    expect_true(is.list(val_entry))
    expect_true(is.function(val_entry$counter))
    expect_true(is.function(val_entry$ready))
    expect_true(is.function(val_entry$state))

    expect_true(is.list(scv_entry))
    expect_true(is.function(scv_entry$counter))
    expect_true(is.function(scv_entry$ready))
    expect_true(is.function(scv_entry$state))

    # Not ready when sc_obj is NULL
    expect_type(val_entry$ready(), "character")
    expect_type(scv_entry$ready(), "character")

    # Provide sc_obj, still not ready without file input
    gd$sc_obj <- seu
    session$flushReact()
    expect_type(val_entry$ready(), "character")

    # Provide file input -> validate becomes ready
    session$setInputs(
      velocity_import_mode = "rds",
      velocity_rds_file = data.frame(
        name = "vel.rds",
        size = file.info(rds_path)$size,
        type = "application/octet-stream",
        datapath = rds_path,
        stringsAsFactors = FALSE
      )
    )
    expect_identical(val_entry$ready(), TRUE)
    expect_type(scv_entry$ready(), "character")

    # Trigger velocity_validate via drive counter
    val_entry$counter(val_entry$counter() + 1L)
    session$flushReact()

    expect_false(is.null(velocity_state$result))
    expect_identical(val_entry$state()$has_result, TRUE)

    # Now scvelo_run is ready (canonical result + umap reduction exist)
    expect_identical(scv_entry$ready(), TRUE)

    # Missing reduction makes scvelo_run not-ready; if forced, firing counter hits .vel_scvelo_fail
    session$setInputs(velocity_reduction = "tsne_missing")
    expect_type(scv_entry$ready(), "character")

    scv_entry$counter(scv_entry$counter() + 1L)
    session$flushReact()
    expect_identical(velocity_state$optin_status, "error")
  })
})
