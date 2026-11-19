# =============================================================================
# test-shard-merge.R — M-SH2 RED-first (unique écrivain du relevé mergé).
# Toutes les fixtures vivent dans tempdir() : le fichier canonique
# QC/full_suite_results.txt n'est JAMAIS ecrit par ce test.
# =============================================================================

library(testthat)

.tm_root <- local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "tools", "run_full_suite.R"))) break
    p <- dirname(d); if (identical(p, d)) stop("racine introuvable", call. = FALSE)
    d <- p
  }
  d
})
.tm_tool <- file.path(.tm_root, "tools", "merge_shard_results.R")
.tm_rscript <- file.path(R.home("bin"), "Rscript")

.tm_seam <- function() file.exists(.tm_tool)
.tm_skip <- function() skip_if_not(.tm_seam(),
  "merge_shard_results.R absent : RED attendu")

.tm_fixture <- function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  # test-check_duplication.R est le discriminant de collation de la fiche §4.0 :
  # en C il trie AVANT test-check-alpha.R, en francais APRES.
  pop <- c("test-check-alpha.R", "test-check_duplication.R",
           "test-zebra.R", "test-shinytest2-sc.R")
  popfile <- file.path(dir, "population.txt")
  writeLines(pop, popfile, useBytes = TRUE)
  s1 <- file.path(dir, "shard_1_of_2.txt")
  s2 <- file.path(dir, "shard_2_of_2.txt")
  writeLines(c(
    "start 00:00:00 | 2 files",
    sprintf("%-58s fail=0 pass=3 err=0 skip=0 t=1.0s", pop[1]),
    sprintf("%-58s fail=1 pass=2 err=0 skip=0 t=2.0s", pop[2]),
    "CENSUS: begin", "CENSUS: total=2",
    "BILAN: failed=1 passed=5 error=0 skipped=0 elapsed_s=3 | 00:00:03"), s1, useBytes = TRUE)
  writeLines(c(
    "start 00:00:05 | 2 files",
    sprintf("%-58s fail=0 pass=7 err=0 skip=1 t=4.0s", pop[3]),
    sprintf("%-58s fail=0 pass=1 err=0 skip=0 t=5.0s", pop[4]),
    "CENSUS: begin", "CENSUS: total=2",
    "BILAN: failed=0 passed=8 error=0 skipped=1 elapsed_s=5 | 00:00:10"), s2, useBytes = TRUE)
  list(pop = pop, popfile = popfile, s1 = s1, s2 = s2)
}
.tm_merge <- function(fx, out, shards = c(fx$s1, fx$s2)) {
  old <- getwd(); setwd(.tm_root); on.exit(setwd(old), add = TRUE)
  args <- c(shQuote(.tm_tool), shQuote(paste0("--out=", out)),
            shQuote(paste0("--population-file=", fx$popfile)),
            shQuote(paste0("--n=", length(shards))), shQuote(shards))
  log <- tempfile("merge-log-"); on.exit(unlink(log), add = TRUE)
  st <- suppressWarnings(system2(.tm_rscript, args, stdout = log, stderr = log,
                                 timeout = 120))
  list(status = st, log = paste(gsub("\r$", "", readLines(log, warn = FALSE)),
                                collapse = "\n"), out = out)
}

test_that("M-SH2 RED : l'outil de merge existe (sinon rien n'est lance)", {
  expect_true(.tm_seam(), info = "merge_shard_results.R non implemente — RED attendu")
})

test_that("merge : BILAN unique mode=sharded_2, sommes exactes, census complet", {
  .tm_skip()
  fx <- .tm_fixture(tempfile("merge-ok-"))
  out <- tempfile("merged-")
  r <- .tm_merge(fx, out)
  expect_identical(r$status, 0L, info = r$log)
  lines <- gsub("\r$", "", readLines(out, warn = FALSE))
  expect_identical(sum(grepl("^BILAN:", lines)), 1L)
  bl <- grep("^BILAN:", lines, value = TRUE)
  expect_match(bl, "mode=sharded_2", fixed = TRUE)
  expect_match(bl, "failed=1 passed=13 error=0 skipped=1", fixed = TRUE)
  expect_match(bl, "elapsed_s=10", fixed = TRUE)
  expect_identical(sum(grepl("^CENSUS: test-", lines)), 4L)
  expect_true(any(grepl("^CENSUS: total=4$", lines)))
})

test_that("merge : ordre per-file = ordre sort() du runner sous la locale imposee", {
  .tm_skip()
  fx <- .tm_fixture(tempfile("merge-order-"))
  out <- tempfile("merged-order-")
  r <- .tm_merge(fx, out)
  expect_identical(r$status, 0L, info = r$log)
  lines <- gsub("\r$", "", readLines(out, warn = FALSE))
  got <- sub(" .*$", "", grep("^test-", lines, value = TRUE))
  expect_identical(got, sort(fx$pop))
  expect_true(identical(got[1], "test-check-alpha.R"))
})

test_that("merge : CRLF purs et CRLF final (contrat byte du relevé canonique)", {
  .tm_skip()
  fx <- .tm_fixture(tempfile("merge-crlf-"))
  out <- tempfile("merged-crlf-")
  r <- .tm_merge(fx, out)
  expect_identical(r$status, 0L, info = r$log)
  raw <- readBin(out, "raw", n = file.info(out)$size)
  txt <- rawToChar(raw)
  expect_identical(nchar(gsub("[^\r]", "", txt)), nchar(gsub("[^\n]", "", txt)))
  expect_true(substr(txt, nchar(txt) - 1L, nchar(txt)) == "\r\n")
})

test_that("merge : shard manquant / sans BILAN / doublon => erreur dure, sortie non ecrite", {
  .tm_skip()
  fx <- .tm_fixture(tempfile("merge-errors-"))
  out <- tempfile("merged-err-")
  missing <- .tm_merge(fx, out, shards = c(fx$s1, file.path(dirname(fx$s1), "absent.txt")))
  expect_false(identical(missing$status, 0L))
  expect_false(file.exists(out))
  out2 <- tempfile("merged-err2-")
  nobilan <- file.path(dirname(fx$s1), "no_bilan.txt")
  writeLines(c("start 00:00:00 | 1 files",
               sprintf("%-58s fail=0 pass=1 err=0 skip=0 t=1.0s", fx$pop[3])), nobilan)
  bad <- .tm_merge(fx, out2, shards = c(nobilan, fx$s2))
  expect_false(identical(bad$status, 0L))
  expect_false(file.exists(out2))
})
