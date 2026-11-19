# =============================================================================
# test-shard-drive-root.R — M-SH1 : isolation drive au scope shard (fiche §3.1).
# Complementaire de test-shard-partition.R : ici on eprouve l'INVARIANT
# « aucun write host handshake », pas la géométrie de la partition.
# =============================================================================

library(testthat)

.tdr_root <- local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "tools", "run_full_suite.R"))) break
    p <- dirname(d); if (identical(p, d)) stop("racine introuvable", call. = FALSE)
    d <- p
  }
  d
})
.tdr_runner <- file.path(.tdr_root, "tools", "run_full_suite.R")
.tdr_rscript <- file.path(R.home("bin"), "Rscript")
.tdr_drive <- file.path(.tdr_root, "tools", "_drive")
.tdr_seam <- function() {
  src <- paste(readLines(.tdr_runner, warn = FALSE), collapse = "\n")
  grepl("TRANSCRIPTO_DRIVE_ROOT", src, fixed = TRUE) &&
    grepl(".shard_root", src, fixed = TRUE)
}
.tdr_skip <- function() skip_if_not(.tdr_seam(), "sandbox shard absent : RED attendu")
# `ready.json` est un battement de coeur VIVANT (mesuré 2026-11-19 02:44:19,
# réécrit toutes les quelques secondes par une session opérateur) : il ne peut
# PAS entrer dans un invariant byte-identique. Les fichiers de PROTOCOLE
# (state/arm/result), eux, doivent rester intacts — c'est la classe de
# contamination documentée (STATUS §0, réécriture 02:22:34).
.tdr_sha <- function() {
  files <- c("state.json", "arm.json", "result.json")
  vapply(files, function(f) {
    p <- file.path(.tdr_drive, f)
    if (file.exists(p)) tools::md5sum(p)[[1]] else "ABSENT"
  }, character(1))
}

test_that("M-SH1 RED : le runner pose un drive root par shard", {
  expect_true(.tdr_seam(), info = "sandbox shard non implementee — RED attendu")
})

test_that("plan shardié ne modifie aucun fichier d'état drive du dépôt", {
  .tdr_skip()
  before <- .tdr_sha()
  old <- getwd(); setwd(.tdr_root); on.exit(setwd(old), add = TRUE)
  log <- tempfile("drive-plan-"); on.exit(unlink(log), add = TRUE)
  st <- suppressWarnings(system2(.tdr_rscript,
    c(shQuote(.tdr_runner), "--plan", "--shard=1/2",
      shQuote(paste0("--timings=", file.path("QC", "full_suite_results.txt")))),
    stdout = log, stderr = log, timeout = 180))
  expect_identical(st, 0L)
  expect_identical(.tdr_sha(), before)
})

test_that("le sandbox est pose AVANT la boucle test_file (verrou source)", {
  .tdr_skip()
  src <- readLines(.tdr_runner, warn = FALSE)
  i_env <- grep("Sys\\.setenv\\(TRANSCRIPTO_DRIVE_ROOT", src)[1]
  i_loop <- grep("for \\(f in files\\)", src)[1]
  expect_true(!is.na(i_env) && !is.na(i_loop))
  expect_true(i_env < i_loop,
    info = "le root doit etre pose avant le premier test_file")
})
