# =============================================================================
# test-bulk-meta-analysis.R — Horizontal multi-study DE meta-analysis (BK-07)
# =============================================================================
if (!exists("run_horizontal_meta_de", mode = "function")) {
  if (exists("source_project_file", mode = "function")) {
    source_project_file("tools/pkg_skeleton/R/bulk_meta_analysis.R")
  } else {
    source(file.path("..", "..", "tools", "pkg_skeleton", "R", "bulk_meta_analysis.R"))
  }
}

.make_study <- function(genes, lfc, pval) {
  data.frame(
    gene = genes,
    log2FC = lfc,
    pvalue = pval,
    stringsAsFactors = FALSE
  )
}

test_that("concordant up-regulated genes across 3 studies get low meta_pval, FDR < 0.05, concordance == 1", {
  genes <- paste0("G", 1:10)
  s1 <- .make_study(
    genes,
    lfc  = c(2.4, 1.9, -2.1, 0.05, -0.1, 0.2, -0.05, 0.1, -0.15, 0.02),
    pval = c(1e-4, 5e-4, 2e-4, 0.70, 0.80, 0.60, 0.90, 0.55, 0.75, 0.85)
  )
  s2 <- .make_study(
    genes,
    lfc  = c(2.1, 2.2, -1.8, -0.04, 0.12, -0.1, 0.08, -0.02, 0.11, -0.03),
    pval = c(2e-4, 3e-4, 4e-4, 0.65, 0.72, 0.81, 0.60, 0.92, 0.68, 0.77)
  )
  s3 <- .make_study(
    genes,
    lfc  = c(2.6, 1.7, -2.3, 0.08, -0.06, 0.15, -0.09, 0.04, -0.05, 0.01),
    pval = c(5e-5, 8e-4, 1e-4, 0.78, 0.64, 0.50, 0.84, 0.73, 0.88, 0.91)
  )

  studies <- list(CohortA = s1, CohortB = s2, CohortC = s3)

  for (m in c("stouffer", "fisher", "maxP", "minP")) {
    res <- run_horizontal_meta_de(studies, method = m)
    expect_type(res, "list")
    expect_identical(res$method, m)
    expect_true(res$direction_aware)
    expect_identical(res$n_studies, 3L)
    expect_identical(res$n_genes, 10L)
    expect_identical(res$study_names, c("CohortA", "CohortB", "CohortC"))
    expect_identical(res$provenance, list(
      engine = "biomentat_horizontal_meta_de",
      method = m
    ))

    tbl <- res$table
    expect_s3_class(tbl, "data.frame")
    expect_true(all(c(
      "gene", "meta_pval", "FDR", "mean_log2FC", "n_studies", "direction_concordance"
    ) %in% colnames(tbl)))

    # Sorted by FDR, then meta_pval
    expect_equal(tbl$FDR, sort(tbl$FDR))

    # Top concordant genes G1, G2 (up) and G3 (down)
    g1 <- tbl[tbl$gene == "G1", ]
    g2 <- tbl[tbl$gene == "G2", ]
    g3 <- tbl[tbl$gene == "G3", ]

    expect_lt(g1$meta_pval, 1e-3)
    expect_lt(g1$FDR, 0.05)
    expect_equal(g1$direction_concordance, 1)
    expect_gt(g1$mean_log2FC, 2)
    expect_identical(g1$n_studies, 3L)

    expect_lt(g2$meta_pval, 1e-3)
    expect_lt(g2$FDR, 0.05)
    expect_equal(g2$direction_concordance, 1)

    expect_lt(g3$meta_pval, 1e-3)
    expect_lt(g3$FDR, 0.05)
    expect_equal(g3$direction_concordance, 1)
    expect_lt(g3$mean_log2FC, -1.5)
  }
})

test_that("exact numerical formulas match reference calculations for all 4 combiners", {
  s1 <- .make_study(c("G1", "G2"), c(2.0, -1.5), c(0.01, 0.04))
  s2 <- .make_study(c("G1", "G2"), c(1.0, -0.5), c(0.02, 0.20))
  studies <- list(S1 = s1, S2 = s2)

  # Stouffer direction-aware
  res_st_dir <- run_horizontal_meta_de(studies, method = "stouffer", direction_aware = TRUE)
  z1 <- (stats::qnorm(1 - 0.01 / 2) + stats::qnorm(1 - 0.02 / 2)) / sqrt(2)
  p1_expected <- 2 * stats::pnorm(-abs(z1))
  expect_equal(
    res_st_dir$table$meta_pval[res_st_dir$table$gene == "G1"],
    p1_expected,
    tolerance = 1e-12
  )

  # Stouffer unsigned (direction_aware = FALSE)
  res_st_uns <- run_horizontal_meta_de(studies, method = "stouffer", direction_aware = FALSE)
  expect_false(res_st_uns$direction_aware)
  z1_u <- (stats::qnorm(1 - 0.01) + stats::qnorm(1 - 0.02)) / sqrt(2)
  p1_u_expected <- stats::pnorm(z1_u, lower.tail = FALSE)
  expect_equal(
    res_st_uns$table$meta_pval[res_st_uns$table$gene == "G1"],
    p1_u_expected,
    tolerance = 1e-12
  )

  # Fisher
  res_fi <- run_horizontal_meta_de(studies, method = "fisher")
  chi2_1 <- -2 * (log(0.01) + log(0.02))
  expect_equal(
    res_fi$table$meta_pval[res_fi$table$gene == "G1"],
    stats::pchisq(chi2_1, df = 4, lower.tail = FALSE),
    tolerance = 1e-12
  )

  # maxP (Wilkinson)
  res_max <- run_horizontal_meta_de(studies, method = "maxP")
  expect_equal(
    res_max$table$meta_pval[res_max$table$gene == "G1"],
    max(0.01, 0.02)^2,
    tolerance = 1e-12
  )

  # minP (Tippett)
  res_min <- run_horizontal_meta_de(studies, method = "minP")
  expect_equal(
    res_min$table$meta_pval[res_min$table$gene == "G1"],
    1 - (1 - min(0.01, 0.02))^2,
    tolerance = 1e-12
  )
})

test_that("discordant gene is penalized by direction-aware stouffer and maxP vs fisher", {
  # DISC_OPP: strong up in S1, equally strong down in S2
  # DISC_ONE: strong up in S1, null in S2
  # CONC_UP : strong up in both S1 and S2
  s1 <- .make_study(
    c("CONC_UP", "DISC_OPP", "DISC_ONE"),
    lfc  = c(2.5,  3.0, 2.8),
    pval = c(1e-5, 1e-6, 1e-6)
  )
  s2 <- .make_study(
    c("CONC_UP", "DISC_OPP", "DISC_ONE"),
    lfc  = c(2.4, -3.0, 0.05),
    pval = c(1e-5, 1e-6, 0.85)
  )
  studies <- list(S1 = s1, S2 = s2)

  res_st <- run_horizontal_meta_de(studies, method = "stouffer", direction_aware = TRUE)
  res_fi <- run_horizontal_meta_de(studies, method = "fisher")
  res_mx <- run_horizontal_meta_de(studies, method = "maxP")

  get_p <- function(res, g) res$table$meta_pval[res$table$gene == g]

  # Opposite-direction strong gene: direction-aware Stouffer cancels out (~1.0),
  # whereas Fisher ignores direction and returns a tiny p-value.
  expect_gt(get_p(res_st, "DISC_OPP"), 0.99)
  expect_lt(get_p(res_fi, "DISC_OPP"), 1e-9)
  expect_gt(get_p(res_st, "DISC_OPP"), get_p(res_fi, "DISC_OPP"))
  expect_equal(
    res_st$table$direction_concordance[res_st$table$gene == "DISC_OPP"],
    0
  )

  # Single-study hit (strong in S1, null in S2): maxP penalizes heavily
  # (0.85^2 = 0.7225) compared to Fisher (< 1e-4).
  expect_gt(get_p(res_mx, "DISC_ONE"), 0.7)
  expect_lt(get_p(res_fi, "DISC_ONE"), 1e-4)
  expect_gt(get_p(res_mx, "DISC_ONE"), get_p(res_fi, "DISC_ONE"))
})

test_that("min_studies filtering works across partial gene overlap", {
  s1 <- .make_study(c("ALL3", "IN2", "ONLY1"), c(2, 2, 2), c(1e-3, 1e-3, 1e-3))
  s2 <- .make_study(c("ALL3", "IN2"),          c(2, 2),    c(1e-3, 1e-3))
  s3 <- .make_study(c("ALL3"),                 c(2),       c(1e-3))
  studies <- list(S1 = s1, S2 = s2, S3 = s3)

  res2 <- run_horizontal_meta_de(studies, min_studies = 2L)
  expect_identical(res2$n_genes, 2L)
  expect_setequal(res2$table$gene, c("ALL3", "IN2"))
  expect_identical(res2$table$n_studies[res2$table$gene == "ALL3"], 3L)
  expect_identical(res2$table$n_studies[res2$table$gene == "IN2"], 2L)

  res3 <- run_horizontal_meta_de(studies, min_studies = 3L)
  expect_identical(res3$n_genes, 1L)
  expect_identical(res3$table$gene, "ALL3")

  # min_studies exceeding number of studies returns clean 0-row table
  res4 <- run_horizontal_meta_de(studies, min_studies = 4L)
  expect_identical(res4$n_genes, 0L)
  expect_identical(nrow(res4$table), 0L)
  expect_true(all(c(
    "gene", "meta_pval", "FDR", "mean_log2FC", "n_studies", "direction_concordance"
  ) %in% colnames(res4$table)))
})

test_that("weighted Stouffer reflects study weights and supports named weight vectors", {
  # Study 1 has strong up signal; Study 2 has modest down signal.
  s1 <- .make_study("G1", lfc =  2.0, pval = 1e-4)
  s2 <- .make_study("G1", lfc = -1.0, pval = 0.05)
  studies <- list(BigCohort = s1, SmallCohort = s2)

  res_eq <- run_horizontal_meta_de(studies, method = "stouffer")
  res_wt <- run_horizontal_meta_de(studies, method = "stouffer", weights = c(10, 1))
  res_wt_named <- run_horizontal_meta_de(
    studies, method = "stouffer",
    weights = c(SmallCohort = 1, BigCohort = 10)
  )

  # Giving BigCohort weight 10 makes combined p-value much smaller than equal weights
  expect_lt(res_wt$table$meta_pval, res_eq$table$meta_pval)
  expect_equal(res_wt$table$meta_pval, res_wt_named$table$meta_pval, tolerance = 1e-12)

  # Exact weighted Stouffer check
  z1 <- stats::qnorm(1e-4 / 2, lower.tail = FALSE)
  z2 <- -stats::qnorm(0.05 / 2, lower.tail = FALSE)
  z_expected <- (10 * z1 + 1 * z2) / sqrt(10^2 + 1^2)
  expect_equal(res_wt$table$meta_pval, 2 * stats::pnorm(-abs(z_expected)), tolerance = 1e-12)
})

test_that("column aliases, rownames, auto-naming, clamping, and zero-lfc edge cases work", {
  # Unnamed list with mixed column aliases and rownames
  df1 <- data.frame(
    gene_id = c("TP53", "BRCA1", "ZERO"),
    log2FoldChange = c(2.0, -1.5, 0.0),
    PValue = c(0, 0.01, 0.5),
    stringsAsFactors = FALSE
  )
  df2 <- data.frame(
    symbol = c("TP53", "BRCA1", "ZERO"),
    logFC = c(1.8, -1.2, 0.0),
    P.Value = c(1e-320, 0.02, 0.5),
    stringsAsFactors = FALSE
  )
  df3 <- data.frame(
    log2FC = c(2.2, -1.4, 0.0),
    pval = c(1e-10, 0.03, 0.5),
    row.names = c("TP53", "BRCA1", "ZERO"),
    stringsAsFactors = FALSE
  )

  res <- run_horizontal_meta_de(list(df1, df2, df3), method = "stouffer")
  expect_identical(res$study_names, c("study_1", "study_2", "study_3"))
  expect_identical(res$n_genes, 3L)
  expect_false(anyNA(res$table$meta_pval))
  expect_true(all(is.finite(res$table$meta_pval)))
  expect_equal(res$table$direction_concordance[res$table$gene == "ZERO"], 1)

  # Clamping p = 0 must also avoid Inf/NaN in Fisher, maxP, minP
  for (m in c("fisher", "maxP", "minP")) {
    rm <- run_horizontal_meta_de(list(df1, df2, df3), method = m)
    expect_false(anyNA(rm$table$meta_pval))
    expect_true(all(is.finite(rm$table$meta_pval)))
  }
})

test_that("input validation errors are raised cleanly with classed errorCondition", {
  valid_df <- .make_study(c("G1", "G2"), c(1, -1), c(0.01, 0.02))

  # Fewer than 2 studies or non-list
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df)),
    class = "bulk_meta_analysis_error"
  )
  expect_error(
    run_horizontal_meta_de(valid_df),
    class = "bulk_meta_analysis_error"
  )
  expect_error(
    run_horizontal_meta_de(NULL),
    class = "bulk_meta_analysis_error"
  )

  # Non-data.frame element
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = list(gene = "G1"))),
    class = "bulk_meta_analysis_error"
  )

  # Missing log2FC column
  no_lfc <- data.frame(gene = c("G1", "G2"), pvalue = c(0.01, 0.02))
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = no_lfc)),
    class = "bulk_meta_analysis_error"
  )

  # Missing pvalue column
  no_pval <- data.frame(gene = c("G1", "G2"), log2FC = c(1, -1))
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = no_pval)),
    class = "bulk_meta_analysis_error"
  )

  # Missing gene column and no explicit rownames
  no_gene <- data.frame(log2FC = c(1, -1), pvalue = c(0.01, 0.02))
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = no_gene)),
    class = "bulk_meta_analysis_error"
  )

  # Invalid method, weights, min_studies, fdr_method, direction_aware
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = valid_df), method = "invalid"),
    class = "bulk_meta_analysis_error"
  )
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = valid_df), weights = c(1, -2)),
    class = "bulk_meta_analysis_error"
  )
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = valid_df), weights = c(1, 2, 3)),
    class = "bulk_meta_analysis_error"
  )
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = valid_df), min_studies = 0L),
    class = "bulk_meta_analysis_error"
  )
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = valid_df), fdr_method = "bogus"),
    class = "bulk_meta_analysis_error"
  )
  expect_error(
    run_horizontal_meta_de(list(S1 = valid_df, S2 = valid_df), direction_aware = NA),
    class = "bulk_meta_analysis_error"
  )
})
