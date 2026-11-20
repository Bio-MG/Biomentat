# =============================================================================
# test-ccc-diff-conditions.R — W2.2 §11.8 Gem #1 : comparaison différentielle
# Condition A vs Condition B sur des communication_result canoniques
# =============================================================================

if (!exists("communication_views_public_api", mode = "function") ||
    !exists("finalize_communication_result", mode = "function")) {
  load_project_pkg()
}

.diff_make_comm_result <- function(rows_df, source_method = "cellchat") {
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
  finalize_communication_result(
    canonical_table = tbl,
    source_method = source_method
  )
}

test_that("build_communication_condition_diff and plot_communication_condition_diff are in communication_views_public_api()", {
  api <- communication_views_public_api()
  expect_true("build_communication_condition_diff" %in% api)
  expect_true("plot_communication_condition_diff" %in% api)
})

test_that("build_communication_condition_diff validates inputs and refuses cross-tool comparisons", {
  res_cc <- .diff_make_comm_result(
    data.frame(sender = "CD4 T", receiver = "B", ligand = "LIG1", receptor = "REC1", score = 0.4),
    source_method = "cellchat"
  )
  res_cpdb <- .diff_make_comm_result(
    data.frame(sender = "CD4 T", receiver = "B", ligand = "LIG1", receptor = "REC1", score = 0.8),
    source_method = "cellphonedb"
  )

  expect_error(
    build_communication_condition_diff(NULL, res_cc),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_condition_diff(res_cc, list(type = "wrong")),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_condition_diff(res_cc, res_cpdb),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_condition_diff(res_cc, res_cc, label_a = ""),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_condition_diff(res_cc, res_cc, label_b = NA_character_),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_condition_diff(res_cc, res_cc, pseudocount = 0),
    class = "communication_import_error"
  )
  expect_error(
    build_communication_condition_diff(res_cc, res_cc, pseudocount = -1e-4),
    class = "communication_import_error"
  )
})

test_that("build_communication_condition_diff computes exact delta_score, log2_fc, join status, and metadata", {
  res_a <- .diff_make_comm_result(data.frame(
    sender   = c("CD4 T", "B",     "NK"),
    receiver = c("B",     "CD8 T", "CD4 T"),
    ligand   = c("LIG1",  "LIG2",  "LIG3"),
    receptor = c("REC1",  "REC2",  "REC3"),
    pathway  = c("PW1",   "PW2",   "PW3"),
    score    = c(0.20,    0.50,    0.30),
    stringsAsFactors = FALSE
  ))
  res_b <- .diff_make_comm_result(data.frame(
    sender   = c("CD4 T", "B",     "Mono"),
    receiver = c("B",     "CD8 T", "NK"),
    ligand   = c("LIG1",  "LIG2",  "LIG4"),
    receptor = c("REC1",  "REC2",  "REC4"),
    pathway  = c("PW1",   "PW2",   "PW4"),
    score    = c(0.80,    0.10,    0.60),
    stringsAsFactors = FALSE
  ))

  pc <- 1e-6
  diff_df <- build_communication_condition_diff(
    res_a, res_b,
    label_a = "Control",
    label_b = "Treated",
    pseudocount = pc
  )

  expect_s3_class(diff_df, "data.frame")
  expect_true(all(c(
    "sender", "receiver", "ligand", "receptor", "pathway",
    "score_a", "score_b", "delta_score", "log2_fc", "presence"
  ) %in% colnames(diff_df)))
  expect_identical(nrow(diff_df), 4L)

  # Attributes / metadata
  expect_identical(attr(diff_df, "source_method"), "cellchat")
  expect_identical(attr(diff_df, "label_a"), "Control")
  expect_identical(attr(diff_df, "label_b"), "Treated")
  expect_equal(attr(diff_df, "pseudocount"), pc)
  expect_true(nzchar(attr(diff_df, "epistemic_note")))

  # Shared up: LIG1 (0.2 -> 0.8)
  r1 <- diff_df[diff_df$ligand == "LIG1", ]
  expect_equal(r1$score_a, 0.20)
  expect_equal(r1$score_b, 0.80)
  expect_equal(r1$delta_score, 0.60, tolerance = 1e-12)
  expect_equal(r1$log2_fc, log2((0.80 + pc) / (0.20 + pc)), tolerance = 1e-12)
  expect_identical(r1$presence, "both")

  # Shared down: LIG2 (0.5 -> 0.1)
  r2 <- diff_df[diff_df$ligand == "LIG2", ]
  expect_equal(r2$delta_score, -0.40, tolerance = 1e-12)
  expect_equal(r2$log2_fc, log2((0.10 + pc) / (0.50 + pc)), tolerance = 1e-12)
  expect_identical(r2$presence, "both")

  # Only A: LIG3 (0.3 -> 0)
  r3 <- diff_df[diff_df$ligand == "LIG3", ]
  expect_equal(r3$score_a, 0.30)
  expect_equal(r3$score_b, 0.0)
  expect_equal(r3$delta_score, -0.30, tolerance = 1e-12)
  expect_equal(r3$log2_fc, log2((0 + pc) / (0.30 + pc)), tolerance = 1e-12)
  expect_identical(r3$presence, "only_a")

  # Only B: LIG4 (0 -> 0.6)
  r4 <- diff_df[diff_df$ligand == "LIG4", ]
  expect_equal(r4$score_a, 0.0)
  expect_equal(r4$score_b, 0.60)
  expect_equal(r4$delta_score, 0.60, tolerance = 1e-12)
  expect_equal(r4$log2_fc, log2((0.60 + pc) / (0 + pc)), tolerance = 1e-12)
  expect_identical(r4$presence, "only_b")

  # Sorted by abs(delta_score) descending
  expect_true(all(diff(abs(diff_df$delta_score)) <= 1e-12))
})

test_that("build_communication_condition_diff deduplicates repeated pairs by max score and handles NA scores", {
  res_a <- .diff_make_comm_result(data.frame(
    sender   = c("CD4 T", "CD4 T"),
    receiver = c("B",     "B"),
    ligand   = c("LIG1",  "LIG1"),
    receptor = c("REC1",  "REC1"),
    score    = c(0.15,    0.45),
    stringsAsFactors = FALSE
  ))
  res_b <- .diff_make_comm_result(data.frame(
    sender   = "CD4 T",
    receiver = "B",
    ligand   = "LIG1",
    receptor = "REC1",
    score    = 0.75,
    stringsAsFactors = FALSE
  ))

  diff_df <- build_communication_condition_diff(res_a, res_b)
  expect_identical(nrow(diff_df), 1L)
  expect_equal(diff_df$score_a, 0.45)
  expect_equal(diff_df$delta_score, 0.30, tolerance = 1e-12)

  # Rank-only (all NA scores) should not crash
  res_na_a <- .diff_make_comm_result(data.frame(
    sender = "CD4 T", receiver = "B", ligand = "LIG1", receptor = "REC1",
    score = NA_real_, stringsAsFactors = FALSE
  ))
  res_na_b <- .diff_make_comm_result(data.frame(
    sender = "CD4 T", receiver = "B", ligand = "LIG1", receptor = "REC1",
    score = NA_real_, stringsAsFactors = FALSE
  ))
  diff_na <- build_communication_condition_diff(res_na_a, res_na_b)
  expect_identical(nrow(diff_na), 1L)
  expect_true(is.na(diff_na$delta_score))
  expect_true(is.na(diff_na$log2_fc))
})

test_that("plot_communication_condition_diff validates input, enforces top_n, and renders ggplot2", {
  expect_error(
    plot_communication_condition_diff(NULL),
    class = "communication_import_error"
  )
  expect_error(
    plot_communication_condition_diff(data.frame(x = 1)),
    class = "communication_import_error"
  )

  res_a <- .diff_make_comm_result(data.frame(
    sender   = c("CD4 T", "B",     "NK",    "Mono"),
    receiver = c("B",     "CD8 T", "CD4 T", "NK"),
    ligand   = c("L1",    "L2",    "L3",    "L4"),
    receptor = c("R1",    "R2",    "R3",    "R4"),
    score    = c(0.10,    0.80,    0.30,    0.50),
    stringsAsFactors = FALSE
  ))
  res_b <- .diff_make_comm_result(data.frame(
    sender   = c("CD4 T", "B",     "NK",    "Mono"),
    receiver = c("B",     "CD8 T", "CD4 T", "NK"),
    ligand   = c("L1",    "L2",    "L3",    "L4"),
    receptor = c("R1",    "R2",    "R3",    "R4"),
    score    = c(0.90,    0.20,    0.35,    0.10),
    stringsAsFactors = FALSE
  ))

  diff_df <- build_communication_condition_diff(
    res_a, res_b,
    label_a = "Sain",
    label_b = "Malade"
  )
  p <- plot_communication_condition_diff(diff_df, top_n = 2L)
  expect_s3_class(p, "ggplot")
  expect_identical(nrow(p$data), 2L)
  expect_match(p$labels$subtitle, "Sain", fixed = TRUE)
  expect_match(p$labels$subtitle, "Malade", fixed = TRUE)
  expect_match(p$labels$subtitle, "cellchat", ignore.case = TRUE)

  # Empty / all-NA diff_res returns explicit empty ggplot
  diff_empty <- diff_df[0, , drop = FALSE]
  p_empty <- plot_communication_condition_diff(diff_empty)
  expect_s3_class(p_empty, "ggplot")
})
