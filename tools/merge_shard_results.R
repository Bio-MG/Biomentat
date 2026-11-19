# tools/merge_shard_results.R — M-SH2 : unique écrivain du relevé mergé.
# Reconstruit UN relevé byte-compatible de la suite sérielle à partir de N
# fichiers de shard. Aucun vert partiel : shard manquant, shard tronqué (sans
# BILAN), doublon inter-shards, population incomplète ou somme incohérente =>
# erreur dure, aucun fichier de sortie n'est créé.
#
# Usage : LC_ALL=French_France.65001 Rscript tools/merge_shard_results.R \
#           --out=QC/full_suite_results_sharded.txt \
#           [shard_1.txt shard_2.txt ...]
#   --population-file=<path> : override de la population (FIXTURES DE TEST
#     uniquement — en production la population est re-dérivée de tests/testthat,
#     exactement comme le runner, avec la même grammaire).
#   --n=<N> : nombre de shards attendu (défaut = nombre de fichiers passés).
# Le BILAN mergé porte mode=sharded_<N> et elapsed_s = MAKESPAN (horloge murale
# du start global au BILAN du dernier shard), jamais la somme des elapsed_s.
.localectl <- Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8")
if (!nzchar(.localectl)) {
  warning("LC_CTYPE fr_FR.UTF-8 unavailable: parse() false failures are likely")
}
args <- commandArgs(trailingOnly = TRUE)
.flag <- function(name) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^--", name, "="), "", a[1]) else NA_character_
}
out_spec   <- .flag("out")
pop_spec   <- .flag("population-file")
n_spec     <- .flag("n")
shard_files <- setdiff(args, grep("^--(out|population-file|n)=", args, value = TRUE))
if (!length(shard_files)) {
  cat("ERREUR : aucun fichier de shard fourni\n", file = stderr())
  quit(save = "no", status = 1L)
}
out_file <- if (is.na(out_spec)) "QC/full_suite_results_sharded.txt" else out_spec
die <- function(msg) {
  cat(paste0("ERREUR merge : ", msg, "\n"), file = stderr())
  quit(save = "no", status = 1L)
}
getn <- function(txt, tag) {
  m <- regmatches(txt, regexpr(paste0(tag, "=(-?\\d+)"), txt))
  if (!length(m)) return(NA_integer_)
  as.integer(sub(paste0(tag, "="), "", m))
}
secs <- function(hms) {
  p <- as.integer(strsplit(hms, ":")[[1]])
  if (length(p) != 3L || any(is.na(p))) return(NA_real_)
  p[1L] * 3600 + p[2L] * 60 + p[3L]
}

# --- LECTURE DES SHARDS -----------------------------------------------------
per_file <- list()
bilans <- character(0)
starts <- character(0)
ends <- character(0)
for (sf in shard_files) {
  if (!file.exists(sf)) die(paste("shard introuvable :", sf))
  lines <- gsub("\r$", "", readLines(sf, warn = FALSE))
  bl <- grep("^BILAN:", lines, value = TRUE)
  if (length(bl) != 1L) {
    die(sprintf("shard %s : %d lignes BILAN (exactement 1 requis — shard tronqué ?)",
                basename(sf), length(bl)))
  }
  bilans <- c(bilans, bl)
  hd <- grep("^start [0-9:]+ \\| ", lines, value = TRUE)
  if (length(hd) != 1L) die(sprintf("shard %s : header start manquant", basename(sf)))
  starts <- c(starts, sub("^start ([0-9:]+) .*$", "\\1", hd))
  ends <- c(ends, sub("^.*\\| ([0-9:]+)$", "\\1", bl))
  fl <- grep("^(test-[^[:space:]]+\\.R)[[:space:]]+fail=", lines, value = TRUE)
  bns <- sub("^([^[:space:]]+)[[:space:]].*$", "\\1", fl)
  if (!length(fl)) die(sprintf("shard %s : aucune ligne per-file", basename(sf)))
  per_file[[sf]] <- stats::setNames(fl, bns)
}
all_names <- unlist(lapply(per_file, names), use.names = FALSE)
if (anyDuplicated(all_names)) {
  die(paste("fichier dupliqué entre shards :",
            paste(unique(all_names[duplicated(all_names)]), collapse = ", ")))
}

# --- POPULATION -------------------------------------------------------------
pop <- if (!is.na(pop_spec)) {
  gsub("\r$", "", readLines(pop_spec, warn = FALSE))
} else {
  basename(sort(list.files("tests/testthat", pattern = "^test-.*\\.R$",
                           full.names = TRUE)))
}
pop <- sort(pop[nzchar(pop)])
if (length(all_names) != length(pop) || !setequal(all_names, pop)) {
  die(sprintf("union des shards (%d) != population (%d) — manquants : %s ; externes : %s",
              length(all_names), length(pop),
              paste(setdiff(pop, all_names), collapse = ","),
              paste(setdiff(all_names, pop), collapse = ",")))
}
# Somme des tallies per-file == somme des BILAN de shards (garde anti-vert-partiel).
tags <- c("failed", "passed", "error", "skipped")
file_tags <- c("fail", "pass", "err", "skip")
sums <- vapply(seq_along(tags), function(j) {
  vals <- unlist(lapply(per_file, function(v) getn(v, file_tags[j])), use.names = FALSE)
  sum(vals, na.rm = TRUE)
}, numeric(1))
bilan_sums <- vapply(tags, function(t) sum(vapply(bilans, getn, integer(1), t)), numeric(1))
if (!identical(as.integer(sums), as.integer(bilan_sums))) {
  die(sprintf("somme per-file (%s) != somme BILAN (%s)",
              paste(sums, collapse = "/"), paste(bilan_sums, collapse = "/")))
}

# --- SORTIE (écrite UNIQUEMENT après toutes les validations) ----------------
n <- if (is.na(n_spec)) length(shard_files) else as.integer(n_spec)
if (is.na(n) || n < 1L) die("--n invalide")
t0 <- min(vapply(starts, secs, numeric(1)), na.rm = TRUE)
t1 <- max(vapply(ends, secs, numeric(1)), na.rm = TRUE)
makespan <- if (t1 < t0) t1 + 86400 - t0 else t1 - t0
lines <- c(
  sprintf("start %s | %d files", starts[which.min(vapply(starts, secs, numeric(1)))], length(pop)),
  unlist(lapply(pop, function(b) per_file[[names(which(vapply(per_file, function(v) b %in% names(v), logical(1))))[1L]]][[b]]), use.names = FALSE),
  "CENSUS: begin", paste0("CENSUS: ", pop), sprintf("CENSUS: total=%d", length(pop)),
  sprintf("BILAN: mode=sharded_%d failed=%d passed=%d error=%d skipped=%d elapsed_s=%.0f | %s",
          n, sums[1L], sums[2L], sums[3L], sums[4L], makespan,
          ends[which.max(vapply(ends, secs, numeric(1)))])
)
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
con <- file(out_file, open = "wt")
writeLines(lines, con)
close(con)
cat(sprintf("merge OK : %d shards -> %d fichiers, failed=%d passed=%d error=%d skipped=%d elapsed_s=%.0f | %s\n",
            n, length(pop), sums[1L], sums[2L], sums[3L], sums[4L], makespan, out_file))
