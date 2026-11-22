# =============================================================================
# test-bulk-provenance-bk08.R — SHA-256 dataset & anti-leakage CV folds (BK-08)
# =============================================================================
if (!exists("ts_dataset_sha256", mode = "function")) {
  source_project_file("tools/pkg_skeleton/R/bulk_provenance.R")
}

test_that("ts_dataset_sha256 : empreinte 64-hex deterministe et metadonnees fichier", {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(c("gene\ts1\ts2", "TP53\t10\t20", "EGFR\t5\t15"), tmp)

  res1 <- ts_dataset_sha256(tmp)
  res2 <- ts_dataset_sha256(tmp)

  expect_identical(res1$algo, "sha256")
  expect_match(res1$sha256, "^[0-9a-f]{64}$")
  expect_identical(res1$sha256, res2$sha256)
  expect_equal(res1$size_bytes, as.numeric(file.info(tmp)$size))
  expect_identical(res1$path, normalizePath(tmp, winslash = "/", mustWork = FALSE))

  writeLines(c("gene\ts1\ts2", "TP53\t11\t20", "EGFR\t5\t15"), tmp)
  res3 <- ts_dataset_sha256(tmp)
  expect_false(identical(res1$sha256, res3$sha256))
})

test_that("ts_dataset_sha256 : erreurs classees sur fichier absent, repertoire ou entree invalide", {
  missing_path <- tempfile(fileext = ".csv")
  expect_error(ts_dataset_sha256(missing_path), class = "biomentat_provenance_error")

  dir_path <- tempfile()
  dir.create(dir_path)
  on.exit(unlink(dir_path, recursive = TRUE), add = TRUE)
  expect_error(ts_dataset_sha256(dir_path), class = "biomentat_provenance_error")

  expect_error(ts_dataset_sha256(NULL), class = "biomentat_provenance_error")
  expect_error(ts_dataset_sha256(""), class = "biomentat_provenance_error")
  expect_error(ts_dataset_sha256(c("a", "b")), class = "biomentat_provenance_error")
  expect_error(ts_dataset_sha256(NA_character_), class = "biomentat_provenance_error")
})

test_that("ts_cv_leakage_free_folds : zero fuite de groupe entre les plis", {
  # 10 donneurs, 3 echantillons chacun = 30 echantillons
  groups  <- rep(paste0("donor_", seq_len(10)), each = 3)
  samples <- paste0("S_", seq_along(groups))

  cv <- ts_cv_leakage_free_folds(samples, groups, k = 5L, seed = 42L)

  expect_true(isTRUE(cv$leakage_verified))
  expect_identical(cv$k, 5L)
  expect_identical(cv$n_samples, 30L)
  expect_identical(cv$n_groups, 10L)
  expect_length(cv$fold_id, 30L)
  expect_length(cv$folds, 5L)
  expect_setequal(unlist(cv$folds, use.names = FALSE), samples)

  # Aucun donneur ne doit apparaitre dans plus d'un pli
  fold_groups <- lapply(seq_len(5L), function(f) unique(groups[cv$fold_id == f]))
  expect_identical(sum(lengths(fold_groups)), 10L)
  for (i in seq_len(4L)) {
    for (j in seq((i + 1L), 5L)) {
      expect_length(intersect(fold_groups[[i]], fold_groups[[j]]), 0L)
    }
  }
})

test_that("ts_cv_leakage_free_folds : stratification par strate et restauration de .Random.seed", {
  # 12 donneurs (6 cas, 6 temoins), 2 echantillons par donneur
  donor_ids <- paste0("D", seq_len(12))
  donor_stratum <- rep(c("case", "ctrl"), each = 6)
  groups  <- rep(donor_ids, each = 2)
  strata  <- rep(donor_stratum, each = 2)
  samples <- paste0("samp_", seq_along(groups))

  set.seed(123)
  seed_before <- .Random.seed
  cv1 <- ts_cv_leakage_free_folds(samples, groups, strata = strata, k = 3L, seed = 99L)
  expect_identical(.Random.seed, seed_before)

  cv2 <- ts_cv_leakage_free_folds(samples, groups, strata = strata, k = 3L, seed = 99L)
  expect_identical(cv1$fold_id, cv2$fold_id)
  expect_true(isTRUE(cv1$leakage_verified))

  # Chaque pli doit contenir exactement 2 donneurs 'case' et 2 donneurs 'ctrl' (4 echantillons chacun)
  for (f in seq_len(3L)) {
    idx <- which(cv1$fold_id == f)
    f_groups <- unique(groups[idx])
    f_strata <- vapply(f_groups, function(g) strata[match(g, groups)], character(1))
    expect_identical(sum(f_strata == "case"), 2L)
    expect_identical(sum(f_strata == "ctrl"), 2L)
  }
})

test_that("ts_cv_leakage_free_folds : erreurs de validation classees", {
  expect_error(ts_cv_leakage_free_folds(character(0), character(0)),
               class = "biomentat_provenance_error")
  expect_error(ts_cv_leakage_free_folds(c("s1", "s1"), c("g1", "g2"), k = 2L),
               class = "biomentat_provenance_error")
  expect_error(ts_cv_leakage_free_folds(c("s1", NA), c("g1", "g2"), k = 2L),
               class = "biomentat_provenance_error")
  expect_error(ts_cv_leakage_free_folds(c("s1", "s2"), c("g1"), k = 2L),
               class = "biomentat_provenance_error")
  expect_error(ts_cv_leakage_free_folds(c("s1", "s2"), c("g1", NA), k = 2L),
               class = "biomentat_provenance_error")
  expect_error(ts_cv_leakage_free_folds(c("s1", "s2"), c("g1", "g2"), strata = "A", k = 2L),
               class = "biomentat_provenance_error")
  # Un meme groupe dans 2 strates differentes -> erreur
  expect_error(ts_cv_leakage_free_folds(c("s1", "s2", "s3", "s4"),
                                        c("g1", "g1", "g2", "g2"),
                                        strata = c("A", "B", "A", "B"),
                                        k = 2L),
               class = "biomentat_provenance_error")
  # k invalide ou > n_groups
  expect_error(ts_cv_leakage_free_folds(c("s1", "s2"), c("g1", "g2"), k = 1L),
               class = "biomentat_provenance_error")
  expect_error(ts_cv_leakage_free_folds(c("s1", "s2", "s3"), c("g1", "g1", "g2"), k = 3L),
               class = "biomentat_provenance_error")
  expect_error(ts_cv_leakage_free_folds(c("s1", "s2"), c("g1", "g2"), k = 2L, seed = NA),
               class = "biomentat_provenance_error")
})
