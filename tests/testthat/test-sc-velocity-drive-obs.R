# =============================================================================
# test-sc-velocity-drive-obs.R
#   W1.1 #1 & #5: drive observability (`has_vectors`, `optin_status`) and
#   stale status/source label refresh on `import_velocity` drive completion.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
source_project_file("modules/sc/mod_sc_velocity.R")

test_that(".velocity_drive_probe_fields returns bounded has_vectors & optin_status", {
  p_null <- .velocity_drive_probe_fields(NULL)
  expect_identical(p_null$has_vectors, FALSE)
  expect_identical(p_null$optin_status, "unavailable")

  p_empty <- .velocity_drive_probe_fields(list(umap_vectors = NULL))
  expect_identical(p_empty$has_vectors, FALSE)
  expect_identical(p_empty$optin_status, "idle")

  vecs <- matrix(c(0.1, -0.1, 0.2, -0.2), nrow = 2,
                 dimnames = list(c("C1", "C2"), c("d1", "d2")))
  for (st in c("idle", "running", "done", "error", "unavailable")) {
    p <- .velocity_drive_probe_fields(list(umap_vectors = vecs), optin_status = st)
    expect_true(is.logical(p$has_vectors) && length(p$has_vectors) == 1L)
    expect_identical(p$has_vectors, TRUE)
    expect_identical(p$optin_status, st)
  }

  p_err <- .velocity_drive_probe_fields(list(umap_vectors = NULL,
                                             optin_error = list(message = "fail")))
  expect_identical(p_err$optin_status, "error")

  p_done <- .velocity_drive_probe_fields(list(
    umap_vectors = vecs,
    vector_validation = list(field = "scvelo_optin", ok = TRUE)
  ))
  expect_identical(p_done$has_vectors, TRUE)
  expect_identical(p_done$optin_status, "done")
})

test_that("mod_sc_velocity.R keeps its exact 3-function public surface", {
  exprs <- parse(file.path(ts_project_root(), "modules/sc/mod_sc_velocity.R"),
                 keep.source = FALSE)
  top_names <- character(0)
  for (e in exprs) {
    if (is.call(e) && as.character(e[[1L]]) %in% c("<-", "=")) {
      lhs <- e[[2L]]
      if (is.symbol(lhs)) top_names <- c(top_names, as.character(lhs))
    }
  }
  internal <- grep("^\\.", top_names, value = TRUE)
  public <- setdiff(top_names, internal)
  expect_setequal(public, c("mod_sc_velocity_ui",
                            "mod_sc_velocity_output_ui",
                            "mod_sc_velocity_server"))
  expect_true(".velocity_drive_probe_fields" %in% internal)
})

test_that("import_velocity drive transitions publish has_vectors/optin_status and refresh stale labels", {
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

  rds1 <- tempfile(fileext = ".rds")
  rds2 <- tempfile(fileext = ".rds")
  on.exit(unlink(c(rds1, rds2)), add = TRUE)

  saveRDS(list(spliced = sp, unspliced = un), rds1)
  vecs <- matrix(0.05, nrow = 6, ncol = 2,
                 dimnames = list(cells, c("d1", "d2")))
  saveRDS(list(spliced = sp, unspliced = un,
               umap_vectors = vecs, embedding_reduction = "umap"), rds2)

  gd <- shiny::reactiveValues(
    sc_obj = seu,
    language = "fr",
    i18n = NULL,
    drive_registry = new.env(parent = emptyenv())
  )
  srv <- shiny::reactiveValues()

  shiny::testServer(mod_sc_velocity_server, args = list(global_data = gd, shared_rv = srv), {
    reg_entry <- ts_drive_registry(gd)[[paste0(TS_DRIVE_IMPORTER_PREFIX,
                                               TS_DRIVE_SC_VELOCITY_IMPORT_MODULE)]]
    expect_true(is.list(reg_entry))
    expect_true(is.function(reg_entry$importer))
    expect_true(is.function(reg_entry$state))

    # Initial probe before any import
    st0 <- reg_entry$state()
    expect_identical(st0$has_result, FALSE)
    expect_identical(st0$has_vectors, FALSE)
    expect_identical(st0$optin_status, "idle")

    # 1. Drive import #1 without vectors
    res1 <- reg_entry$importer(list(path = rds1, name = "first_vel.rds"))
    session$flushReact()
    expect_true(isTRUE(res1$ok))
    expect_identical(res1$has_vectors, FALSE)
    expect_identical(res1$optin_status, "idle")
    expect_match(res1$source_label, "first_vel.rds (rds)", fixed = TRUE)
    expect_match(res1$status_msg, "first_vel.rds (rds)", fixed = TRUE)

    expect_identical(velocity_state$has_vectors, FALSE)
    expect_identical(velocity_state$optin_status, "idle")
    expect_identical(velocity_state$source_label, "first_vel.rds (rds)")
    expect_identical(velocity_state$status_msg, velocity_status_rv())
    expect_true(any(grepl("first_vel.rds (rds)", velocity_log_rv(), fixed = TRUE)))
    expect_identical(srv$velocity_panel_state$has_vectors, FALSE)
    expect_identical(srv$velocity_panel_state$optin_status, "idle")
    expect_identical(srv$velocity_panel_state$source_label, "first_vel.rds (rds)")

    st1 <- reg_entry$state()
    expect_identical(st1$has_vectors, FALSE)
    expect_identical(st1$optin_status, "idle")
    expect_identical(st1$source_label, "first_vel.rds (rds)")

    # 2. scvelo_run failure & success transitions via internal handlers
    .vel_scvelo_fail(simpleError("mock scvelo failure"))
    expect_identical(velocity_state$optin_status, "error")
    expect_identical(velocity_state$has_vectors, FALSE)
    expect_identical(srv$velocity_panel_state$optin_status, "error")
    expect_identical(reg_entry$state()$optin_status, "error")

    .apply_scvelo_result(
      list(
        embedding_vectors = vecs,
        banner = "scVelo opt-in",
        provenance = new_provenance_entry(
          analysis_id = "sc-velocity",
          method = "scvelo_optin",
          parameters = list(mode = "stochastic")
        )
      ),
      sel_red = "umap",
      canon = velocity_state$result
    )
    session$flushReact()
    expect_identical(velocity_state$optin_status, "done")
    expect_identical(velocity_state$has_vectors, TRUE)
    expect_identical(srv$velocity_panel_state$optin_status, "done")
    expect_identical(srv$velocity_panel_state$has_vectors, TRUE)
    expect_identical(reg_entry$state()$optin_status, "done")
    expect_identical(reg_entry$state()$has_vectors, TRUE)

    # 3. Drive import #2 with precomputed vectors: stale labels & optin_status refreshed
    res2 <- reg_entry$importer(list(path = rds2, name = "second_vel.rds"))
    session$flushReact()
    expect_true(isTRUE(res2$ok))
    expect_identical(res2$has_vectors, TRUE)
    expect_identical(res2$optin_status, "idle")
    expect_match(res2$source_label, "second_vel.rds (rds)", fixed = TRUE)
    expect_match(res2$status_msg, "second_vel.rds (rds)", fixed = TRUE)
    expect_false(grepl("first_vel.rds", res2$status_msg, fixed = TRUE))

    expect_identical(velocity_state$source_label, "second_vel.rds (rds)")
    expect_identical(velocity_state$status_msg, velocity_status_rv())
    expect_false(grepl("first_vel.rds", velocity_status_rv(), fixed = TRUE))
    expect_true(grepl("second_vel.rds (rds)", utils::tail(velocity_log_rv(), 1L), fixed = TRUE))
    expect_identical(srv$velocity_panel_state$has_vectors, TRUE)
    expect_identical(srv$velocity_panel_state$optin_status, "idle")
    expect_identical(srv$velocity_panel_state$source_label, "second_vel.rds (rds)")

    # 4. Reset transition on Seurat object change
    gd$sc_obj <- subset(seu, cells = cells[1:4])
    session$flushReact()
    expect_identical(velocity_state$has_vectors, FALSE)
    expect_identical(velocity_state$optin_status, "idle")
    expect_null(velocity_state$source_label)
    expect_null(srv$velocity_panel_state)
    st_reset <- reg_entry$state()
    expect_identical(st_reset$has_vectors, FALSE)
    expect_identical(st_reset$optin_status, "idle")
    expect_null(st_reset$source_label)
  })
})
