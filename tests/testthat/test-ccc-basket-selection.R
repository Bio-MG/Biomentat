# =============================================================================
# test-ccc-basket-selection.R — W2.2 §11.8 Gem #5 : sélection click-to-basket
# et filtrage croisé sur communication_result canonique
# =============================================================================

if (!exists("communication_views_public_api", mode = "function") ||
    !exists("finalize_communication_result", mode = "function")) {
  if (exists("load_project_pkg", mode = "function")) {
    load_project_pkg()
  } else {
    repo_root <- if (file.exists("tools/pkg_skeleton/R/sc_communication_views.R")) "." else "../.."
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication.R"), local = FALSE)
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication_views.R"), local = FALSE)
  }
}

.basket_make_comm_result <- function(rows_df, source_method = "cellchat") {
  tbl <- data.frame(
    sender = as.character(rows_df$sender),
    receiver = as.character(rows_df$receiver),
    ligand = as.character(rows_df$ligand),
    receptor = as.character(rows_df$receptor),
    interaction = paste0(rows_df$ligand, " -> ", rows_df$receptor),
    pathway = if ("pathway" %in% colnames(rows_df)) as.character(rows_df$pathway)
              else paste0("PW_", rows_df$ligand),
    score = as.numeric(rows_df$score),
    p_value = if ("p_value" %in% colnames(rows_df)) as.numeric(rows_df$p_value) else 0.01,
    p_adjusted = if ("p_adjusted" %in% colnames(rows_df)) as.numeric(rows_df$p_adjusted) else 0.02,
    source_method = rep(source_method, nrow(rows_df)),
    source_file = rep(paste0("mock_", source_method, ".csv"), nrow(rows_df)),
    source_cell_identity_level = rep("celltype", nrow(rows_df)),
    stringsAsFactors = FALSE
  )
  if ("sender_mapped" %in% colnames(rows_df)) {
    tbl$sender_mapped <- as.character(rows_df$sender_mapped)
  }
  if ("receiver_mapped" %in% colnames(rows_df)) {
    tbl$receiver_mapped <- as.character(rows_df$receiver_mapped)
  }
  finalize_communication_result(
    canonical_table = tbl,
    source_method = source_method
  )
}

test_that("build_communication_basket_selection is in communication_views_public_api()", {
  api <- communication_views_public_api()
  expect_true("build_communication_basket_selection" %in% api)
})

test_that("build_communication_basket_selection validates res, selection, and max_items", {
  res <- .basket_make_comm_result(data.frame(
    sender = "CD4 T", receiver = "B", ligand = "LIG1", receptor = "REC1", score = 0.5,
    stringsAsFactors = FALSE
  ))

  expect_error(
    build_communication_basket_selection(NULL),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_basket_selection(list(type = "wrong")),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_basket_selection(res, selection = "CD4 T"),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_basket_selection(res, max_items = 0L),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_basket_selection(res, max_items = -5L),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_basket_selection(res, max_items = NA_integer_),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_basket_selection(res, max_items = "10"),
    class = "communication_import_error"
  )
})

test_that("build_communication_basket_selection deduplicates, sorts by score desc, and extracts unique genes", {
  res <- .basket_make_comm_result(data.frame(
    sender   = c("CD4 T", "CD4 T", "B",     "NK",    "Mono"),
    receiver = c("B",     "B",     "CD8 T", "CD4 T", "NK"),
    ligand   = c("LIG1",  "LIG1",  "LIG2",  "LIG1",  "LIG3"),
    receptor = c("REC1",  "REC1",  "REC2",  "REC2",  "REC3"),
    pathway  = c("PW1",   "PW1",   "PW2",   "PW1",   "PW3"),
    score    = c(0.20,    0.85,    0.60,    NA_real_, 0.95),
    stringsAsFactors = FALSE
  ))

  out <- build_communication_basket_selection(res)

  expect_true(is.list(out))
  expect_named(out, c("genes", "pairs", "n_selected"))
  expect_identical(out$n_selected, 4L)
  expect_true(is.data.frame(out$pairs))
  expect_identical(
    colnames(out$pairs),
    c("ligand", "receptor", "sender", "receiver", "pathway", "score")
  )
  expect_identical(nrow(out$pairs), 4L)

  # Ordered by score descending, NA last
  expect_equal(out$pairs$score, c(0.95, 0.85, 0.60, NA_real_))
  expect_identical(out$pairs$ligand, c("LIG3", "LIG1", "LIG2", "LIG1"))

  # Unique ligand + receptor symbols
  expect_setequal(out$genes, c("LIG1", "LIG2", "LIG3", "REC1", "REC2", "REC3"))
  expect_identical(length(out$genes), 6L)
})

test_that("build_communication_basket_selection cross-filters by selection and caps at max_items", {
  res <- .basket_make_comm_result(data.frame(
    sender   = c("CD4 T", "CD4 T", "CD4 T", "B",     "NK"),
    receiver = c("B",     "CD8 T", "NK",    "CD8 T", "CD4 T"),
    ligand   = c("LIG1",  "LIG2",  "LIG3",  "LIG4",  "LIG5"),
    receptor = c("REC1",  "REC2",  "REC1",  "REC4",  "REC5"),
    pathway  = c("PW1",   "PW1",   "PW2",   "PW1",   "PW2"),
    score    = c(0.40,    0.90,    0.70,    0.80,    0.50),
    stringsAsFactors = FALSE
  ))

  # Filter by sender + pathway (AND semantics), ignoring NULL/empty and extra click keys
  out_sp <- build_communication_basket_selection(
    res,
    selection = list(
      sender = "CD4 T",
      pathway = "PW1",
      receiver = NULL,
      ligand = "",
      curveNumber = 0L
    )
  )
  expect_identical(out_sp$n_selected, 2L)
  expect_identical(out_sp$pairs$ligand, c("LIG2", "LIG1"))
  expect_setequal(out_sp$genes, c("LIG1", "LIG2", "REC1", "REC2"))

  # Cap at max_items = 1L keeps only highest score and its genes
  out_cap <- build_communication_basket_selection(
    res,
    selection = list(sender = "CD4 T"),
    max_items = 1L
  )
  expect_identical(out_cap$n_selected, 1L)
  expect_identical(out_cap$pairs$ligand, "LIG2")
  expect_identical(out_cap$pairs$score, 0.90)
  expect_setequal(out_cap$genes, c("LIG2", "REC2"))

  # Multi-value filter on receptor
  out_rec <- build_communication_basket_selection(
    res,
    selection = list(receptor = c("REC1", "REC4"))
  )
  expect_identical(out_rec$n_selected, 3L)
  expect_equal(out_rec$pairs$score, c(0.80, 0.70, 0.40))
})

test_that("build_communication_basket_selection respects harmonized node keys and handles empty matches", {
  res <- .basket_make_comm_result(data.frame(
    sender          = c("CD4_raw", "B_raw"),
    receiver        = c("B_raw",   "NK_raw"),
    sender_mapped   = c("CD4 T",   "B"),
    receiver_mapped = c("B",       "NK"),
    ligand          = c("LIG1",    "LIG2"),
    receptor        = c("REC1",    "REC2"),
    pathway         = c("PW1",     "PW2"),
    score           = c(0.75,      0.25),
    stringsAsFactors = FALSE
  ))

  out_mapped <- build_communication_basket_selection(res, selection = list(sender = "CD4 T"))
  expect_identical(out_mapped$n_selected, 1L)
  expect_identical(out_mapped$pairs$sender, "CD4 T")
  expect_identical(out_mapped$pairs$receiver, "B")

  out_empty <- build_communication_basket_selection(res, selection = list(pathway = "NON_EXISTENT"))
  expect_identical(out_empty$n_selected, 0L)
  expect_identical(out_empty$genes, character(0))
  expect_identical(nrow(out_empty$pairs), 0L)
  expect_identical(
    colnames(out_empty$pairs),
    c("ligand", "receptor", "sender", "receiver", "pathway", "score")
  )
})
