# =============================================================================
# test-launch_shards_defaults.R — O40 (Q7-bis) : les DEFAUTS du launcher
# tools/launch_shards.R sont le point de fonctionnement ADOPTÉ (N=6, poids du
# contexte shardié QC/full_suite_results_sharded.txt — levier L1), pas le
# brouillon M-SH1 (N=2, poids sériel canonique QC/full_suite_results.txt).
# Portée : --dry-run UNIQUEMENT — ce test ne lance jamais de suite.
# La ligne « poids / §14 » du dry-run est de la VISIBILITÉ, pas de
# l'automatisation : N reste 6 par défaut, le max-safe N de la table §14
# n'est jamais auto-appliqué (O40 : « N=7 réservé machine calme » ; la table
# §14 reste le pré-vol MANUEL de l'opérateur — runbook
# docs/audits/BASELINES.txt : « choisir N selon RAM libre, table §14 »).
# Chaînes asserties volontairement ASCII (piège locale Windows mesuré,
# AGENTS.md : toujours poser la locale avant parse/regex).
# =============================================================================

library(testthat)

.ls_root <- local({
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "tools", "launch_shards.R"))) break
    p <- dirname(d); if (identical(p, d)) stop("racine introuvable", call. = FALSE)
    d <- p
  }
  d
})
.ls_launcher <- file.path(.ls_root, "tools", "launch_shards.R")
.ls_rscript <- file.path(R.home("bin"), "Rscript")
.ls_dry <- function(extra = character(0)) {
  old <- getwd(); setwd(.ls_root); on.exit(setwd(old), add = TRUE)
  log <- tempfile("launch-dry-"); on.exit(unlink(log), add = TRUE)
  st <- suppressWarnings(system2(.ls_rscript,
    c(shQuote(.ls_launcher), "--dry-run", extra),
    stdout = log, stderr = log, timeout = 120))
  list(status = st, out = paste(gsub("\r$", "", readLines(log, warn = FALSE)),
                                collapse = "\n"))
}

test_that("O40 : le dry-run nu planifie shard 1/6 avec les poids du contexte shardié", {
  r <- .ls_dry()
  expect_match(r$out, "--shard=1/6", fixed = TRUE,
               info = "défaut du launcher != N=6 adopté (O40) — brouillon N=2 encore en place ?")
  expect_match(r$out, "--timings=QC/full_suite_results_sharded.txt", fixed = TRUE,
               info = "poids = relevé sériel canonique au lieu du contexte shardié (levier L1, O40)")
})

test_that("O40 : le dry-run imprime l'artefact de poids, sa population et le max-safe N §14 (informatif)", {
  skip_if_not(file.exists(file.path(.ls_root, "QC", "full_suite_results_sharded.txt")),
              "artefact de poids shardié absent sur ce poste (QC/ local-only)")
  r <- .ls_dry()
  expect_match(r$out, "poids : QC/full_suite_results_sharded.txt \\([0-9]+ fichiers",
               info = "artefact de poids chargé + sa population par fichier non imprimés au dry-run")
  expect_match(r$out, "max-safe N=[0-9]+",
               info = "max-safe N §14 pour la RAM libre mesurée non imprimé au dry-run")
})

test_that("les overrides explicites --n/--timings restent maitres", {
  r <- .ls_dry(c("--n=2", "--timings=QC/full_suite_results.txt"))
  expect_match(r$out, "--shard=1/2", fixed = TRUE,
               info = "un --n explicite n'est plus honoré")
  expect_match(r$out, "--timings=QC/full_suite_results.txt", fixed = TRUE,
               info = "un --timings explicite n'est plus honoré")
  expect_match(r$out, "poids : QC/full_suite_results.txt \\([0-9]+ fichiers",
               info = "la ligne poids ne reflète pas l'override --timings explicite")
})
