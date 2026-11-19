# =============================================================================
# test-shard-runner-args.R — M-SH1 : validation des flags et NON-écriture du
# relevé canonique en mode --plan (garde anti-suite-sous-BOOST).
# =============================================================================

library(testthat)

.tra_root <- local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "tools", "run_full_suite.R"))) break
    p <- dirname(d); if (identical(p, d)) stop("racine introuvable", call. = FALSE)
    d <- p
  }
  d
})
.tra_runner <- file.path(.tra_root, "tools", "run_full_suite.R")
.tra_rscript <- file.path(R.home("bin"), "Rscript")
.tra_canonical <- file.path(.tra_root, "QC", "full_suite_results.txt")
.tra_seam <- function() {
  src <- paste(readLines(.tra_runner, warn = FALSE), collapse = "\n")
  grepl('"--plan"', src, fixed = TRUE) &&
    grepl("shard doit etre de la forme", src, fixed = TRUE)
}
.tra_skip <- function() skip_if_not(.tra_seam(), "seam runner absent : RED attendu")
.tra_call <- function(flag, extra = character(0)) {
  old <- getwd(); setwd(.tra_root); on.exit(setwd(old), add = TRUE)
  log <- tempfile("runner-args-"); on.exit(unlink(log), add = TRUE)
  st <- suppressWarnings(system2(.tra_rscript,
    c(shQuote(.tra_runner), flag, extra), stdout = log, stderr = log, timeout = 180))
  list(status = st, out = paste(gsub("\r$", "", readLines(log, warn = FALSE)),
                                collapse = "\n"))
}

test_that("M-SH1 RED : runner expose --plan et valide --shard", {
  expect_true(.tra_seam(), info = "sharding non implemente — RED attendu")
})

test_that("--shard invalide refuse avec remede, SANS lancer de suite", {
  .tra_skip()
  for (bad in c("--shard=3/2", "--shard=0/2", "--shard=1/0", "--shard=x/2")) {
    r <- .tra_call(bad)
    expect_false(identical(r$status, 0L), info = paste(bad, "aurait du etre refuse"))
    expect_match(r$out, "--shard", fixed = TRUE)
  }
})

test_that("--plan n'ecrit JAMAIS le relevé canonique (meme sans --out)", {
  .tra_skip()
  before <- if (file.exists(.tra_canonical)) {
    tools::md5sum(.tra_canonical)[[1]]
  } else NA_character_
  r <- .tra_call("--plan", c("--shard=2/2"))
  expect_identical(r$status, 0L, info = r$out)
  after <- if (file.exists(.tra_canonical)) {
    tools::md5sum(.tra_canonical)[[1]]
  } else NA_character_
  expect_identical(after, before)
  expect_match(r$out, "PLAN: file 2/2 ", fixed = TRUE)
})
