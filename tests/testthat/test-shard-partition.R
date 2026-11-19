# =============================================================================
# test-shard-partition.R — M-SH1 RED-first (suite sharding, fiche
# docs/proposals/SUITE_SHARDING_DESIGN_2026-11-18.md §2.2).
#
# SECURITE BOOST : tant que le runner n'expose pas `--plan`, on ne lance
# JAMAIS tools/run_full_suite.R — sinon un test RED declencherait la suite
# COMPLETE (interdite sous MODE BOOST). Le premier bloc echoue par assertion,
# pas par erreur, et les suivants sont skippe's tant que le seam est absent.
# =============================================================================

.ts_shard_root <- local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "tools", "run_full_suite.R"))) break
    parent <- dirname(d)
    if (identical(parent, d)) stop("racine projet introuvable", call. = FALSE)
    d <- parent
  }
  d
})
.ts_shard_runner <- file.path(.ts_shard_root, "tools", "run_full_suite.R")
.ts_shard_rscript <- file.path(R.home("bin"), "Rscript")

.ts_shard_seam <- function() {
  src <- paste(readLines(.ts_shard_runner, warn = FALSE), collapse = "\n")
  all(vapply(c('"--plan"', '.shard_flag("shard")',
               '.shard_flag("timings")', '.shard_flag("out")'),
            function(tok) grepl(tok, src, fixed = TRUE), logical(1)))
}
.ts_shard_skip <- function() skip_if_not(.ts_shard_seam(),
  "runner sans seam --plan/--shard/--timings : RED attendu, suite non lancee")

.ts_shard_plan_lines <- function(timings, n = 2L) {
  out <- tempfile("shard-plan-")
  on.exit(unlink(out), add = TRUE)
  never <- tempfile("shard-must-not-exist-")
  on.exit(unlink(never), add = TRUE)
  cmd <- c(shQuote(.ts_shard_runner), "--plan", shQuote(paste0("--shard=1/", n)),
           shQuote(paste0("--timings=", timings)), shQuote(paste0("--out=", never)))
  # testthat::test_dir() tourne depuis tests/testthat ; le runner exige cwd =
  # racine projet (AGENTS.md "Critical environment rule"), sinon population vide.
  old_wd <- getwd()
  setwd(.ts_shard_root)
  on.exit(setwd(old_wd), add = TRUE)
  st <- suppressWarnings(system2(.ts_shard_rscript, cmd, stdout = out, stderr = out,
                                 timeout = 180))
  expect_identical(st, 0L, info = paste("plan runner exit:", st))
  lines <- gsub("\r$", "", readLines(out, warn = FALSE))
  expect_true(any(grepl("^PLAN: shard ", lines)), info = "aucune ligne PLAN: shard")
  expect_false(file.exists(never), info = "--plan ne doit ecrire aucun fichier resultat")
  lines
}
.ts_shard_of <- function(lines) {
  m <- regmatches(lines, regexpr("^PLAN: file (\\d+)/(\\d+) (.+\\.R)$", lines))
  v <- grep("^PLAN: file ", lines, value = TRUE)
  expect_true(length(v) > 0L,
    info = paste("plan sans lignes 'PLAN: file' :", length(lines), "lignes;",
                 "head:", paste(head(lines, 3L), collapse = " || ")))
  i <- as.integer(sub("^PLAN: file (\\d+)/.*$", "\\1", v))
  bn <- sub("^PLAN: file \\d+/\\d+ ", "", v)
  stats::setNames(i, bn)
}
.ts_shard_timings <- function(dir, weight_first) {
  files <- sort(basename(list.files(file.path(.ts_shard_root, "tests", "testthat"),
                                    pattern = "^test-.*\\.R$", full.names = TRUE)))
  stopifnot(length(files) >= 4L)
  lines <- c("start 00:00:00 | 3 files",
    sprintf("%-58s fail=0 pass=1 err=0 skip=0 t=%.1fs", files[1], weight_first),
    sprintf("%-58s fail=0 pass=1 err=0 skip=0 t=%.1fs", files[2], 0.1),
    sprintf("%-58s fail=0 pass=1 err=0 skip=0 t=%.1fs", files[3], 0.1),
    "BILAN: failed=0 passed=3 error=0 skipped=0 elapsed_s=1000 | 00:00:01")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  p <- file.path(dir, "timings.txt")
  writeLines(lines, p, useBytes = TRUE)
  list(path = p, files = files)
}

test_that("M-SH1 RED : le runner expose les seams de sharding (sinon on ne lance rien)", {
  expect_true(.ts_shard_seam(),
    info = paste("sharding non implemente dans tools/run_full_suite.R —",
                 "refus de lancer la suite complete sous BOOST (RED attendu)"))
})

test_that("partition = population exacte, chaque fichier une seule fois, deterministe", {
  .ts_shard_skip()
  fx <- .ts_shard_timings(tempdir(), 1000)
  l1 <- .ts_shard_plan_lines(fx$path)
  l2 <- .ts_shard_plan_lines(fx$path)
  # Les roots sont des tempdir() : determinisme porte sur la PARTITION, pas sur
  # les noms de repertoires volatils (fiche §2.2 point 4).
  part1 <- grep("^PLAN: (shard|file) ", l1, value = TRUE)
  part2 <- grep("^PLAN: (shard|file) ", l2, value = TRUE)
  expect_identical(part1, part2)
  got <- names(.ts_shard_of(l1))
  pop <- basename(list.files(file.path(.ts_shard_root, "tests", "testthat"),
                             pattern = "^test-.*\\.R$", full.names = TRUE))
  expect_setequal(got, pop)
  expect_identical(anyDuplicated(got), 0L)
})

test_that("LPT : le fichier le plus lourd va dans le shard 1 ; muter le BILAN deplaces la partition", {
  .ts_shard_skip()
  fx <- .ts_shard_timings(tempdir(), 1000)
  a <- .ts_shard_of(.ts_shard_plan_lines(fx$path))
  expect_identical(unname(a[[fx$files[1]]]), 1L,
    info = "le poids 1000s doit placer le fichier dans le shard 1")
  # meme population, poids inversés : le lourd doit changer de shard
  fx2 <- .ts_shard_timings(file.path(tempdir(), "mut"), 0.1)
  lines2 <- c(readLines(fx2$path, warn = FALSE)[-4],
              sprintf("%-58s fail=0 pass=1 err=0 skip=0 t=1000.0s", fx2$files[3]),
              "BILAN: failed=0 passed=3 error=0 skipped=0 elapsed_s=1000 | 00:00:01")
  writeLines(lines2, fx2$path, useBytes = TRUE)
  b <- .ts_shard_of(.ts_shard_plan_lines(fx2$path))
  expect_false(identical(a, b), info = "changer le BILAN doit changer la partition")
})

test_that("bundle same-shard : les 4 e2e shinytest2 restent dans le MEME shard", {
  .ts_shard_skip()
  fx <- .ts_shard_timings(tempdir(), 1000)
  s <- .ts_shard_of(.ts_shard_plan_lines(fx$path))
  e2e <- c("test-shinytest2-bulk.R", "test-shinytest2-import.R",
           "test-shinytest2-sc.R", "test-shinytest2-spatial.R")
  present <- intersect(e2e, names(s))
  expect_identical(length(present), 4L, info = "les 4 fichiers e2e doivent exister")
  expect_identical(length(unique(s[present])), 1L,
    info = "les 4 e2e doivent etre coloces (un seul spawn Chrome)")
})

test_that("chaque shard a un drive root distinct, sous tempdir, jamais dans le depot", {
  .ts_shard_skip()
  fx <- .ts_shard_timings(tempdir(), 1000)
  lines <- .ts_shard_plan_lines(fx$path)
  roots <- sub("^PLAN: drive_root [0-9]+/[0-9]+ ", "", grep("^PLAN: drive_root ", lines, value = TRUE))
  layouts <- grep("^PLAN: drive_layout ", lines, value = TRUE)
  expect_identical(length(roots), 2L)
  expect_false(identical(roots[1], roots[2]))
  expect_identical(length(layouts), 2L)
  expect_true(all(grepl(" yes$", layouts)),
    info = "chaque root doit mimer la layout tools/_drive (prouvee avant nettoyage)")
  expect_false(any(grepl("SHINYAPP test", roots, fixed = TRUE)),
    info = "aucun sandbox ne doit pointer dans le depot")
})
