# tools/run_full_suite.R — full test suite, crash-resilient per-file bilan.
# Rationale: R exits with SIGSEGV (139) on this environment when some
# package combos unload (STATUS.md §2l) — with file-redirected stdout the
# test_dir summary is lost. Running file-by-file inside ONE process keeps
# the shared helper state AND flushes a result line per file, so a crash
# is isolated to its file and everything before it survives.
# Locale UTF-8 BEFORE any parse(): Git Bash exports LC_ALL=C.UTF-8, a name R
# does NOT recognise on Windows, so R silently falls back to "C", where it
# cannot parse some UTF-8 sources of the repo (STATUS.md, "Garde C13 +
# semantique testServer() MESUREE"). Not fatal if the locale is missing, but
# the warning must stay visible: a wrong locale yields FALSE failures.
.localectl <- Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8")
if (!nzchar(.localectl)) {
  warning("LC_CTYPE fr_FR.UTF-8 unavailable: parse() false failures are likely")
}

# --- M-SH1 (2026-11-19) : sharding seams, ALL opt-in -------------------------
# Sans AUCUN flag, le comportement reste byte-identique à la version sérielle
# (chemin canonique, population complète, format de ligne et BILAN inchangés).
# `--plan` est un dry-run : il imprime la partition sur stdout et quitte SANS
# écrire le moindre fichier résultats (garde anti-suite-under-BOOST des tests).
args <- commandArgs(trailingOnly = TRUE)
.shard_flag <- function(name) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^--", name, "="), "", a[1]) else NA_character_
}
shard_spec   <- .shard_flag("shard")
out_spec     <- .shard_flag("out")
timings_spec <- .shard_flag("timings")
plan_only    <- "--plan" %in% args
shard_i <- NA_integer_
shard_n <- NA_integer_
if (!is.na(shard_spec)) {
  m <- regmatches(shard_spec, regexpr("^([0-9]+)/([0-9]+)$", shard_spec))
  if (!length(m)) {
    stop("--shard doit etre de la forme i/N (ex. --shard=1/2)", call. = FALSE)
  }
  parts <- as.integer(strsplit(m, "/")[[1]])
  shard_n <- parts[2L]; shard_i <- parts[1L]
  if (is.na(shard_n) || shard_n < 1L || is.na(shard_i) ||
      shard_i < 1L || shard_i > shard_n) {
    stop("--shard : 1 <= i <= N et N >= 1 requis (recu : ", shard_spec, ")",
         call. = FALSE)
  }
}
timings_file <- if (is.na(timings_spec)) "QC/full_suite_results.txt" else timings_spec
# W_new : poids des fichiers sans timing (fiche §2.2). Ancrage mesuré :
# 205 des 305 fichiers du run #15 font <= 2,5 s ; 10 s ne agglutine pas la queue.
W_new <- 10.0
# Bundle same-shard : 4 e2e qui bootent l'app via Chrome/CDP — un seul shard
# spawn Chrome (fiche §2.3). Colocation = un seul item de poids sommé.
.e2e_bundle <- c("test-shinytest2-bulk.R", "test-shinytest2-import.R",
                 "test-shinytest2-sc.R", "test-shinytest2-spatial.R")

# Poids lus DEPUIS LE BILAN (jamais une liste cuite). Un artefact sans BILAN
# => tous les poids = W_new : equilibrage moindre, jamais de fausse exactitude.
.read_weights <- function(path, pop) {
  w <- stats::setNames(rep(W_new, length(pop)), pop)
  if (!file.exists(path)) return(w)
  txt <- readLines(path, warn = FALSE)
  if (!any(grepl("^BILAN:", txt))) return(w)
  idx <- grep("^(test-[^ ]+\\.R)[[:space:]]+fail=.*t=[0-9.]+s[[:space:]]*$", txt)
  for (line in txt[idx]) {
    bn <- sub("^([^[:space:]]+)[[:space:]].*$", "\\1", line)
    tm <- regmatches(line, regexpr("t=[0-9.]+s", line))
    if (!length(tm) || !(bn %in% names(w))) next
    secs <- suppressWarnings(as.numeric(sub("s$", "", sub("^t=", "", tm))))
    if (!is.na(secs)) w[[bn]] <- secs
  }
  w
}

# LPT (longest-processing-time-first) avec bundle de contrainte. Deterministe :
# memes entrees => meme partition ; which.min() departage les egalites par index.
.partition_shards <- function(files, weights, n) {
  bn <- basename(files)
  items <- list()
  used <- character(0)
  bundle <- intersect(.e2e_bundle, bn)
  if (length(bundle)) {
    items[[length(items) + 1L]] <- list(names = bundle,
                                        weight = sum(weights[bundle]))
    used <- c(used, bundle)
  }
  for (b in setdiff(bn, used)) {
    items[[length(items) + 1L]] <- list(names = b, weight = weights[[b]])
  }
  wts <- vapply(items, function(it) it$weight, numeric(1))
  keys <- vapply(items, function(it) sort(it$names)[1L], character(1))
  items <- items[order(-wts, keys)]
  loads <- rep(0, n)
  groups <- vector("list", n)
  for (it in items) {
    k <- which.min(loads)
    groups[[k]] <- c(groups[[k]], it$names)
    loads[k] <- loads[k] + it$weight
  }
  lapply(groups, function(g) files[match(g, bn)])
}

# Un drive root par shard, hors du depot, avec la layout tools/_drive)
# (fiche §3.1 — ts_drive_root() lit TRANSCRIPTO_DRIVE_ROOT).
.shard_root <- function(tag) {
  r <- normalizePath(tempfile(paste0("ts_drive_shard_", tag, "_")),
                     winslash = "/", mustWork = FALSE)
  dir.create(file.path(r, "tools", "_drive"), recursive = TRUE, showWarnings = FALSE)
  r
}

if (plan_only) {
  all_files <- sort(list.files("tests/testthat", pattern = "^test-.*\\.R$",
                               full.names = TRUE))
  n <- if (is.na(shard_n)) 1L else shard_n
  w <- .read_weights(timings_file, basename(all_files))
  groups <- .partition_shards(all_files, w, n)
  roots <- vapply(seq_len(n), function(i) .shard_root(i), character(1))
  on.exit(unlink(roots, recursive = TRUE, force = TRUE), add = TRUE)
  for (i in seq_len(n)) {
    cat(sprintf("PLAN: shard %d/%d files=%d weight=%.1f\n",
                i, n, length(groups[[i]]), sum(w[basename(groups[[i]])])))
    cat(sprintf("PLAN: drive_root %d/%d %s\n", i, n, roots[i]))
    cat(sprintf("PLAN: drive_layout %d/%d %s\n", i, n,
                if (dir.exists(file.path(roots[i], "tools", "_drive"))) "yes" else "no"))
  }
  for (i in seq_len(n)) {
    for (f in groups[[i]]) {
      cat(sprintf("PLAN: file %d/%d %s\n", i, n, basename(f)))
    }
  }
  quit(save = "no", status = 0L)
}
# ----------------------------------------------------------------------------

files <- sort(list.files("tests/testthat", pattern = "^test-.*\\.R$",
                         full.names = TRUE))
weights <- NULL
if (!is.na(shard_i)) {
  weights <- .read_weights(timings_file, basename(files))
  groups <- .partition_shards(files, weights, shard_n)
  files <- groups[[shard_i]]
  root <- .shard_root(shard_i)
  old_drive_root <- Sys.getenv("TRANSCRIPTO_DRIVE_ROOT", unset = NA_character_)
  on.exit({
    if (is.na(old_drive_root)) Sys.unsetenv("TRANSCRIPTO_DRIVE_ROOT")
    else Sys.setenv(TRANSCRIPTO_DRIVE_ROOT = old_drive_root)
    unlink(root, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(TRANSCRIPTO_DRIVE_ROOT = root)
}
cat(sprintf("files: %d\n", length(files)))
if (!is.na(shard_i)) {
  cat(sprintf("shard %d/%d: %d files, weight %.1fs\n", shard_i, shard_n,
              length(files), sum(weights[basename(files)])))
}
# ROOT_HYGIENE item 2 (2026-10-28) : le bilan vit dans QC/ (plus la racine).
dir.create("QC", showWarnings = FALSE)
results_file <- if (is.na(out_spec)) "QC/full_suite_results.txt" else out_spec
con <- file(results_file, open = "wt")
writeLines(sprintf("start %s | %d files", format(Sys.time(), "%H:%M:%S"), length(files)), con)
flush(con)

total_failed <- 0L; total_passed <- 0L; total_error <- 0L; total_skipped <- 0L
.suite_t0 <- proc.time()[["elapsed"]]
crashed_after <- NA_character_
for (f in files) {
  bn <- basename(f)
  # Timing par fichier (H2#9, mesurer avant d'optimiser) : sans lui, « où vont
  # les 45 min » reste une devinette — le CENSUS reste lui déterministe.
  .t0 <- proc.time()[["elapsed"]]
  r <- tryCatch({
    res <- testthat::test_file(f, reporter = "silent", stop_on_failure = FALSE,
                               env = testthat::test_env(),
                               load_package = "none")
    df <- as.data.frame(res)
    sprintf("%-58s fail=%d pass=%d err=%d skip=%d t=%.1fs",
            bn, sum(df$failed), sum(df$passed), sum(df$error), sum(df$skipped),
            proc.time()[["elapsed"]] - .t0)
  }, error = function(e) sprintf("%-58s RUNNER-ERROR: %.120s t=%.1fs", bn,
                                 conditionMessage(e),
                                 proc.time()[["elapsed"]] - .t0))
  writeLines(r, con); flush(con)
  if (grepl("RUNNER-ERROR", r)) next
  # Tally anchored on the labelled fields (a naive digit grab breaks on
  # filenames containing digits, e.g. test-shinytest2-bulk.R).
  getn <- function(tag) {
    m <- regmatches(r, regexpr(paste0(tag, "=(-?\\d+)"), r))
    if (!length(m)) return(0L)
    as.integer(sub(paste0(tag, "="), "", m))
  }
  total_failed  <- total_failed  + getn("fail")
  total_passed  <- total_passed  + getn("pass")
  total_error   <- total_error   + getn("err")
  total_skipped <- total_skipped + getn("skip")
  rm(r)
}
# Census (audit 2026-10-01, R2): emit the exact list of test files that were
# candidates for this run, next to the BILAN verdict, so the bilan can never be
# read as covering a different population than it actually ran. Deterministic:
# sorted, no timestamps, no locale-dependent values. NOTE: the pattern below is
# DUPLICATED in tools/build_census.R on purpose — if you change one, change the
# other (a shared helper was rejected to keep this crash-resilient script
# free of source() dependencies).
writeLines("CENSUS: begin", con)
for (bn in basename(files)) writeLines(paste0("CENSUS: ", bn), con)
writeLines(sprintf("CENSUS: total=%d", length(files)), con)
writeLines(sprintf("BILAN: failed=%d passed=%d error=%d skipped=%d elapsed_s=%.0f | %s",
                   total_failed, total_passed, total_error, total_skipped,
                   proc.time()[["elapsed"]] - .suite_t0,
                   format(Sys.time(), "%H:%M:%S")), con)
close(con)
