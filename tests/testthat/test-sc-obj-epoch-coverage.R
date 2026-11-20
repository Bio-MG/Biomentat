# =============================================================================
# test-sc-obj-epoch-coverage.R — W2.4 / ROADMAP_SC_FUNCTION_TEST.md M2
# =============================================================================
# Audit & lock of `global_data$sc_obj_epoch` coverage across all `modules/sc/`
# and `modules/import/mod_import_sc.R` writers of `global_data$sc_obj`.
#
# Architectural invariant (`tools/pkg_skeleton/R/state.R:80-110`,
# `modules/sc/mod_sc.R:1087-1099`):
#   - `observeEvent(global_data$sc_obj_epoch, ...)` in `mod_sc.R:1094` calls
#     `sc_purge_shared_results(shared_rv)` + `qc_applied_rv(NULL)` on the
#     reactive flush following an epoch bump.
#   - Tier 1 (Dataset replacement / pipeline-reset writers) MUST bump
#     `global_data$sc_obj_epoch` immediately when replacing `global_data$sc_obj`:
#       * `modules/import/mod_import_sc.R` (5 sites: single import, prepared
#         import, Option A merge, drive importer, Option B merge)
#       * `modules/sc/mod_sc_datasets.R` (1 site: `ds_activate`)
#       * `modules/sc/mod_sc_mapping.R` (1 site: `run_mapping`, where
#         `remap_seurat_ids_to_symbol()` builds a brand-new raw
#         `CreateSeuratObject` from counts, stripping all clusters/reductions
#         and resetting the pipeline).
#   - Tier 2 (In-lineage enrichment writers) MUST NOT bump `sc_obj_epoch`,
#     because bumping would fire `observeEvent(global_data$sc_obj_epoch)` on the
#     next reactive flush and wipe out `shared_rv` results (`markers_data`,
#     `pathway_results`, `qc_snapshot`, `traj_reduction`, `traj_method`,
#     `da_design_result`, `qc_applied_rv`) computed in the same lineage:
#       * `modules/sc/mod_sc.R` (2 sites: drive `run_annot`, auto-pipeline commit)
#       * `modules/sc/mod_sc_annotation.R` (2 sites: `run_annot`, `run_mapquery`)
#       * `modules/sc/mod_sc_metadata.R` (1 site: `apply_design`, also frozen by
#         `test-sc-metadata-contract-freeze.R:115`)
#       * `modules/sc/mod_sc_pipeline.R` (1 site: `.apply_manual_result`)
#       * `modules/sc/mod_sc_trajectory.R` (1 site: `run_traj`).
# =============================================================================

source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/state.R")
source_project_file("tools/pkg_skeleton/R/sc_helpers.R")
source_project_file("modules/sc/mod_sc_mapping.R")

suppressPackageStartupMessages({
  library(shiny)
  library(Matrix)
})

# ── AST helpers for `global_data$sc_obj <-` and `global_data$sc_obj_epoch <-` ──

.is_gd_slot_assign <- function(node, slot_name) {
  if (!is.call(node) || length(node) < 3L) return(FALSE)
  h <- tryCatch(node[[1L]], error = function(e) NULL)
  if (!is.symbol(h) || !as.character(h) %in% c("<-", "=")) return(FALSE)
  lhs <- tryCatch(node[[2L]], error = function(e) NULL)
  if (!is.call(lhs) || length(lhs) != 3L) return(FALSE)
  lhs_op <- tryCatch(lhs[[1L]], error = function(e) NULL)
  lhs_obj <- tryCatch(lhs[[2L]], error = function(e) NULL)
  lhs_slot <- tryCatch(lhs[[3L]], error = function(e) NULL)
  if (!is.symbol(lhs_op) || !as.character(lhs_op) %in% c("$", "[[")) return(FALSE)
  if (!is.symbol(lhs_obj) || !identical(as.character(lhs_obj), "global_data")) return(FALSE)
  slot_str <- if (is.symbol(lhs_slot) || is.character(lhs_slot)) as.character(lhs_slot) else ""
  length(slot_str) == 1L && identical(slot_str, slot_name)
}

.ast_count_slot_assigns <- function(relpath, slot_name) {
  p <- ts_ast_parse(relpath)
  hits <- ts_ast_find_all(p, function(x) .is_gd_slot_assign(x, slot_name))
  length(hits)
}

# Checks that every `{ ... }` block containing `global_data$sc_obj <- <rhs>`
# also contains `global_data$sc_obj_epoch <- <bump>` AFTER the `sc_obj` write.
.ast_blocks_with_sc_obj_assign <- function(relpath) {
  p <- ts_ast_parse(relpath)
  blocks <- ts_ast_find_all(p, function(x) {
    if (!ts_ast_is_call_to(x, "{")) return(FALSE)
    l <- as.list(x)
    for (i in seq_along(l)) {
      st <- tryCatch(l[[i]], error = function(e) NULL)
      if (.is_gd_slot_assign(st, "sc_obj")) return(TRUE)
    }
    FALSE
  })
  out <- list()
  for (b in blocks) {
    l <- as.list(b)
    sc_idx <- integer(0)
    ep_idx <- integer(0)
    for (i in seq_along(l)) {
      st <- tryCatch(l[[i]], error = function(e) NULL)
      if (.is_gd_slot_assign(st, "sc_obj")) sc_idx <- c(sc_idx, i)
      if (.is_gd_slot_assign(st, "sc_obj_epoch")) ep_idx <- c(ep_idx, i)
    }
    out[[length(out) + 1L]] <- list(
      sc_idx = sc_idx,
      ep_idx = ep_idx,
      has_bump_after = length(sc_idx) > 0L &&
        length(ep_idx) >= length(sc_idx) &&
        all(vapply(sc_idx, function(s) any(ep_idx > s), logical(1)))
    )
  }
  out
}

# ── 1. Exhaustive AST inventory across `modules/sc/` ────────────────────────

test_that("M2 : the set of global_data$sc_obj writers in modules/sc/ is closed and classified", {
  sc_dir <- file.path(ts_project_root(), "modules", "sc")
  sc_files <- sort(list.files(sc_dir, pattern = "\\.R$", full.names = FALSE))
  expect_true(length(sc_files) > 0L)

  writers <- character(0)
  counts <- integer(0)
  for (f in sc_files) {
    rel <- file.path("modules", "sc", f)
    n <- .ast_count_slot_assigns(rel, "sc_obj")
    if (n > 0L) {
      writers <- c(writers, f)
      counts <- c(counts, setNames(n, f))
    }
  }

  # Exact set of 7 files (9 assignment sites total) in modules/sc/:
  # - Tier 1 (bump required): mod_sc_datasets.R (1), mod_sc_mapping.R (1)
  # - Tier 2 (in-lineage enrichment, no bump): mod_sc.R (2),
  #   mod_sc_annotation.R (2), mod_sc_metadata.R (1), mod_sc_pipeline.R (1),
  #   mod_sc_trajectory.R (1)
  expected_counts <- c(
    "mod_sc.R"            = 2L,
    "mod_sc_annotation.R" = 2L,
    "mod_sc_datasets.R"   = 1L,
    "mod_sc_mapping.R"    = 1L,
    "mod_sc_metadata.R"   = 1L,
    "mod_sc_pipeline.R"   = 1L,
    "mod_sc_trajectory.R" = 1L
  )
  expect_identical(counts, expected_counts)
})

test_that("M2 Tier 1 : every replacement/reset writer bumps sc_obj_epoch in the same block", {
  tier1 <- c(
    "modules/import/mod_import_sc.R" = 5L,
    "modules/sc/mod_sc_datasets.R"   = 1L,
    "modules/sc/mod_sc_mapping.R"    = 1L
  )

  for (rel in names(tier1)) {
    expected_n <- tier1[[rel]]
    n_sc <- .ast_count_slot_assigns(rel, "sc_obj")
    n_ep <- .ast_count_slot_assigns(rel, "sc_obj_epoch")

    expect_identical(n_sc, expected_n, info = paste(rel, "sc_obj write count"))
    expect_identical(n_ep, expected_n, info = paste(rel, "sc_obj_epoch bump count"))

    blocks <- .ast_blocks_with_sc_obj_assign(rel)
    expect_identical(length(blocks), expected_n,
                     info = paste(rel, "AST block count"))
    for (k in seq_along(blocks)) {
      expect_true(blocks[[k]]$has_bump_after,
                  info = sprintf("%s block #%d missing sc_obj_epoch bump after sc_obj write",
                                 rel, k))
    }
  }
})

test_that("M2 Tier 2 : in-lineage enrichment writers in modules/sc/ NEVER assign sc_obj_epoch", {
  tier2 <- c(
    "modules/sc/mod_sc.R",
    "modules/sc/mod_sc_annotation.R",
    "modules/sc/mod_sc_metadata.R",
    "modules/sc/mod_sc_pipeline.R",
    "modules/sc/mod_sc_trajectory.R"
  )

  for (rel in tier2) {
    n_sc <- .ast_count_slot_assigns(rel, "sc_obj")
    n_ep <- .ast_count_slot_assigns(rel, "sc_obj_epoch")
    expect_gt(n_sc, 0L)
    expect_identical(n_ep, 0L,
                     info = paste(rel, "must not bump sc_obj_epoch (in-lineage enrichment)"))
  }

  # And no non-writer file in modules/sc/ may assign sc_obj_epoch either:
  sc_dir <- file.path(ts_project_root(), "modules", "sc")
  sc_files <- sort(list.files(sc_dir, pattern = "\\.R$", full.names = FALSE))
  for (f in setdiff(sc_files, c("mod_sc_datasets.R", "mod_sc_mapping.R"))) {
    rel <- file.path("modules", "sc", f)
    expect_identical(.ast_count_slot_assigns(rel, "sc_obj_epoch"), 0L,
                     info = paste(rel, "unexpected sc_obj_epoch write"))
  }
})

test_that("M2 listener contract : mod_sc.R observes global_data$sc_obj_epoch to purge shared_rv", {
  p <- ts_ast_parse("modules/sc/mod_sc.R")
  obs <- ts_ast_find_all(p, function(x) {
    if (!ts_ast_is_call_to(x, "observeEvent")) return(FALSE)
    l <- as.list(x)
    length(l) >= 3L && grepl("global_data\\$sc_obj_epoch", ts_ast_deparse(l[[2L]]))
  })
  expect_length(obs, 1L)
  body_dep <- ts_ast_deparse(as.list(obs[[1L]])[[3L]])
  expect_match(body_dep, "sc_purge_shared_results(shared_rv)", fixed = TRUE)
  expect_match(body_dep, "qc_applied_rv(NULL)", fixed = TRUE)
})

# ── 2. Targeted reactive verification (`testServer` + purge observer) ───────

.make_mapping_test_obj <- function() {
  skip_if_not_installed("SeuratObject")
  cnt <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 3L),
    j = c(1L, 1L, 2L, 2L),
    x = c(5, 3, 2, 8),
    dims = c(3L, 2L),
    dimnames = list(c("ENSG00000141510", "ENSG00000012048", "ENSG00000139618"),
                    c("cell_1", "cell_2"))
  )
  SeuratObject::CreateSeuratObject(counts = cnt)
}

# Wire the exact purge observer from mod_sc.R:1094 inside the active reactive domain
.install_sc_epoch_purge_observer <- function(gd, shared_rv, qc_applied_rv = NULL) {
  shiny::observeEvent(gd$sc_obj_epoch, {
    sc_purge_shared_results(shared_rv)
    if (is.function(qc_applied_rv)) qc_applied_rv(NULL)
  }, ignoreInit = TRUE)
}

test_that("M2 reactive : mod_sc_mapping_server bumps sc_obj_epoch and purges stale shared_rv results", {
  obj_before <- .make_mapping_test_obj()
  cnt_after <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L),
    j = c(1L, 1L, 2L),
    x = c(5, 3, 2),
    dims = c(2L, 2L),
    dimnames = list(c("TP53", "BRCA1"), c("cell_1", "cell_2"))
  )
  obj_after <- SeuratObject::CreateSeuratObject(counts = cnt_after)

  # Hermetic stub for remap_seurat_ids_to_symbol + showNotification
  old_remap <- get0("remap_seurat_ids_to_symbol", envir = globalenv())
  old_notif <- get0("showNotification", envir = globalenv())
  on.exit({
    if (!is.null(old_remap)) assign("remap_seurat_ids_to_symbol", old_remap, envir = globalenv())
    if (!is.null(old_notif)) assign("showNotification", old_notif, envir = globalenv())
  }, add = TRUE)

  assign("remap_seurat_ids_to_symbol",
         function(obj, from_type = "ensembl", organism = "human",
                  collapse_method = "sum", strip_version = TRUE) {
           list(object = obj_after, n_mapped = 2L, n_unmapped = 1L, n_collapsed = 0L)
         },
         envir = globalenv())
  assign("showNotification", function(...) invisible(NULL), envir = globalenv())

  gd <- shiny::reactiveValues(
    sc_obj = obj_before,
    sc_obj_epoch = 0L,
    i18n = NULL,
    language = "fr"
  )
  shared_rv <- create_sc_shared_state()
  shared_rv$markers_data <- data.frame(gene = "ENSG00000141510", p_val_adj = 0.001)
  shared_rv$pathway_results <- data.frame(ID = "GO:0001")
  shared_rv$qc_snapshot <- data.frame(min_genes = 100)
  shared_rv$traj_reduction <- "umap"
  shared_rv$da_design_result <- list(status = "valid")
  shared_rv$active_tab <- "tab_viz"
  shared_rv$selected_genes <- c("TP53")
  qc_applied_rv <- shiny::reactiveVal(list(source = "auto", min_genes = 100))

  testServer(mod_sc_mapping_server, args = list(id = "mapping", global_data = gd), {
    obs_purge <- .install_sc_epoch_purge_observer(gd, shared_rv, qc_applied_rv)
    on.exit(obs_purge$destroy(), add = TRUE)
    session$flushReact()  # initialize ignoreInit = TRUE observer at epoch 0L

    session$setInputs(
      map_organism = "human",
      map_from_type = "ensembl",
      collapse_method = "sum",
      strip_ensembl_version = TRUE,
      run_mapping = 1L
    )
    session$flushReact()
  })

  # 1. Active object replaced by mapped object and epoch incremented 0L -> 1L
  expect_identical(rownames(shiny::isolate(gd$sc_obj)), c("TP53", "BRCA1"))
  expect_identical(shiny::isolate(gd$sc_obj_epoch), 1L)

  # 2. All 13 shared result slots + qc_applied_rv were purged by the epoch bump
  for (fld in sc_shared_result_fields()) {
    expect_null(state_get(shared_rv, fld), info = paste("slot not purged:", fld))
  }
  expect_null(shiny::isolate(qc_applied_rv()))

  # 3. UI preferences survive the epoch bump
  expect_identical(state_get(shared_rv, "active_tab"), "tab_viz")
  expect_identical(state_get(shared_rv, "selected_genes"), c("TP53"))
})

test_that("M2 reactive : failed mapping does NOT bump sc_obj_epoch or purge shared_rv", {
  obj_before <- .make_mapping_test_obj()

  old_remap <- get0("remap_seurat_ids_to_symbol", envir = globalenv())
  old_notif <- get0("showNotification", envir = globalenv())
  on.exit({
    if (!is.null(old_remap)) assign("remap_seurat_ids_to_symbol", old_remap, envir = globalenv())
    if (!is.null(old_notif)) assign("showNotification", old_notif, envir = globalenv())
  }, add = TRUE)

  assign("remap_seurat_ids_to_symbol",
         function(...) {
           stop(errorCondition("Aucun gene converti", class = "sc_helpers_error"))
         },
         envir = globalenv())
  assign("showNotification", function(...) invisible(NULL), envir = globalenv())

  gd <- shiny::reactiveValues(
    sc_obj = obj_before,
    sc_obj_epoch = 3L,
    i18n = NULL,
    language = "fr"
  )
  shared_rv <- create_sc_shared_state()
  shared_rv$markers_data <- data.frame(gene = "ENSG00000141510", p_val_adj = 0.001)

  testServer(mod_sc_mapping_server, args = list(id = "mapping", global_data = gd), {
    obs_purge <- .install_sc_epoch_purge_observer(gd, shared_rv)
    on.exit(obs_purge$destroy(), add = TRUE)
    session$flushReact()

    session$setInputs(
      map_organism = "human",
      map_from_type = "ensembl",
      collapse_method = "sum",
      strip_ensembl_version = TRUE,
      run_mapping = 1L
    )
    session$flushReact()
  })

  expect_identical(shiny::isolate(gd$sc_obj_epoch), 3L)
  expect_identical(rownames(shiny::isolate(gd$sc_obj)), rownames(obj_before))
  expect_false(is.null(state_get(shared_rv, "markers_data")))
})

test_that("M2 reactive : Tier 2 in-lineage enrichment updates sc_obj while preserving shared_rv results", {
  obj_before <- .make_mapping_test_obj()
  suppressPackageStartupMessages(library(SeuratObject))
  source_project_file("tools/pkg_skeleton/R/sc_communication_input.R")
  source_project_file("tools/pkg_skeleton/R/sc_bpcells.R")
  source_project_file("modules/sc/mod_sc_annotation.R")
  source_project_file("modules/sc/mod_sc.R")

  old_singler <- get0(".run_singler_safe", envir = globalenv())
  on.exit({
    if (!is.null(old_singler)) assign(".run_singler_safe", old_singler, envir = globalenv())
  }, add = TRUE)

  assign(".run_singler_safe",
         function(obj, refcode, labellevel, maxcells = 50000L) {
           list(labels = factor(c("T_cell", "B_cell")), method = "fixture")
         },
         envir = globalenv())

  gd <- shiny::reactiveValues(
    sc_obj = obj_before,
    sc_obj_epoch = 2L,
    i18n = NULL,
    language = "fr"
  )
  shared_rv <- create_sc_shared_state()
  shared_rv$markers_data <- data.frame(gene = "ENSG00000141510", p_val_adj = 0.001)
  shared_rv$pathway_results <- data.frame(ID = "GO:0001")
  shared_rv$qc_snapshot <- data.frame(min_genes = 100)
  shared_rv$traj_reduction <- "umap"
  shared_rv$da_design_result <- list(status = "valid")
  qc_applied_rv <- shiny::reactiveVal(list(source = "auto", min_genes = 100))

  testServer(mod_sc_mapping_server, args = list(id = "mapping", global_data = gd), {
    obs_purge <- .install_sc_epoch_purge_observer(gd, shared_rv, qc_applied_rv)
    on.exit(obs_purge$destroy(), add = TRUE)
    session$flushReact()

    out <- run_annot(
      list(ref_singler = "hpca", label_level = "main", maxcells = 50000L),
      gd, shared_rv, session
    )
    expect_true(isTRUE(out$ok))
    session$flushReact()
  })

  # In-lineage enrichment updated sc_obj with SingleR column, left sc_obj_epoch at 2L,
  # and preserved all shared_rv results across reactive flush:
  expect_true("SingleR_hpca_main" %in% colnames(shiny::isolate(gd$sc_obj)@meta.data))
  expect_identical(shiny::isolate(gd$sc_obj_epoch), 2L)
  expect_false(is.null(state_get(shared_rv, "markers_data")))
  expect_false(is.null(state_get(shared_rv, "pathway_results")))
  expect_false(is.null(state_get(shared_rv, "qc_snapshot")))
  expect_identical(state_get(shared_rv, "traj_reduction"), "umap")
  expect_false(is.null(state_get(shared_rv, "da_design_result")))
  expect_false(is.null(shiny::isolate(qc_applied_rv())))
})

