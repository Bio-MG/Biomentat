# ==============================================================================
# test-biomentat-cli.R — BK-13 Headless Batch CLI Runner
# ==============================================================================

.cli_test_root <- function(rel = ".") {
  if (file.exists(file.path(rel, "tools/biomentat_cli.R"))) return(rel)
  file.path("..", "..")
}

test_that("biomentat_cli_main runs help, bulk_de (mnar), meta_de, and attrition headless", {
  root <- .cli_test_root()
  source(file.path(root, "tools/biomentat_cli.R"), local = globalenv())

  h <- biomentat_cli_main(c("--help"), root_dir = root)
  expect_equal(h$status, "help")
  expect_true(grepl("Biomentat Headless CLI", h$message, fixed = TRUE))

  # 1. attrition task
  tmp_steps <- tempfile(fileext = ".csv")
  tmp_att_out <- tempfile(fileext = ".csv")
  tmp_mmd <- tempfile(fileext = ".mmd")
  on.exit(unlink(c(tmp_steps, tmp_att_out, tmp_mmd)), add = TRUE)
  write.csv(data.frame(
    step = c("Source", "QC", "Final"),
    n_remaining = c(100, 80, 60),
    reason = c("Init", "Low QC", "Missing"),
    stringsAsFactors = FALSE
  ), tmp_steps, row.names = FALSE)

  att <- biomentat_cli_main(c(
    "--task", "attrition", "--steps", tmp_steps, "--unit", "patients",
    "--out", tmp_att_out, "--mermaid", tmp_mmd
  ), root_dir = root)
  expect_equal(att$status, "ok")
  expect_true(file.exists(tmp_att_out))
  expect_true(file.exists(tmp_mmd))
  expect_equal(nrow(att$table), 3L)

  # 2. bulk_de (mnar) task
  tmp_cnt <- tempfile(fileext = ".csv")
  tmp_meta <- tempfile(fileext = ".csv")
  tmp_de_out <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp_cnt, tmp_meta, tmp_de_out)), add = TRUE)

  set.seed(9)
  mat <- matrix(rnorm(48, mean = 10, sd = 0.5), nrow = 8, ncol = 6,
                dimnames = list(paste0("G", 1:8), paste0("S", 1:6)))
  mat[1, 1:3] <- NA_real_
  mat[1, 4:6] <- c(14, 14.2, 13.8)
  write.csv(cbind(feature = rownames(mat), as.data.frame(mat)), tmp_cnt, row.names = FALSE)
  write.csv(data.frame(
    sample = paste0("S", 1:6),
    group = rep(c("ctrl", "case"), each = 3),
    stringsAsFactors = FALSE
  ), tmp_meta, row.names = FALSE)

  de_res <- biomentat_cli_main(c(
    "--task", "bulk_de", "--counts", tmp_cnt, "--meta", tmp_meta,
    "--condition", "group", "--ref", "ctrl", "--target", "case",
    "--engine", "mnar", "--out", tmp_de_out
  ), root_dir = root)
  expect_equal(de_res$status, "ok")
  expect_equal(de_res$n_features, 8L)
  expect_true(file.exists(tmp_de_out))

  # 3. meta_de task
  tmp_s1 <- tempfile(fileext = ".csv")
  tmp_s2 <- tempfile(fileext = ".csv")
  tmp_meta_out <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp_s1, tmp_s2, tmp_meta_out)), add = TRUE)
  s_df <- data.frame(
    gene = paste0("G", 1:5),
    log2FoldChange = c(2.0, -1.5, 0.2, 1.8, -2.1),
    pvalue = c(1e-4, 2e-3, 0.7, 5e-4, 1e-5),
    padj = c(1e-3, 1e-2, 0.9, 2e-3, 1e-4)
  )
  write.csv(s_df, tmp_s1, row.names = FALSE)
  write.csv(s_df, tmp_s2, row.names = FALSE)

  m_res <- biomentat_cli_main(c(
    "--task", "meta_de", "--inputs", paste(tmp_s1, tmp_s2, sep = ","),
    "--method", "stouffer", "--out", tmp_meta_out
  ), root_dir = root)
  expect_equal(m_res$status, "ok")
  expect_equal(m_res$n_genes, 5L)
  expect_true(file.exists(tmp_meta_out))
})
