# Per-file equivalence protocol (fiche SUITE_SHARDING_DESIGN §6).
# Usage: Rscript compare_equiv.R <serial.txt> <sharded.txt> <report.txt>
# Verdict = EXACT equality of (fail,pass,err,skip) per basename. `t=` excluded
# by design (that is exactly what we parallelise). Totals alone are NOT enough.
Sys.setlocale("LC_CTYPE", "French_France.65001")
a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) == 3L)
read_releve <- function(p) {
  lines <- gsub("\r$", "", readLines(p, warn = FALSE))
  per <- grep("^test-.*\\.R[[:space:]]+fail=", lines, value = TRUE)
  bn <- sub("^([^[:space:]]+)[[:space:]].*$", "\\1", per)
  num <- function(tag) as.integer(sub(paste0("^.*", tag, "=(-?\\d+).*$"), "\\1", per))
  df <- data.frame(file = bn, fail = num("fail"), pass = num("pass"),
                   err = num("err"), skip = num("skip"), stringsAsFactors = FALSE)
  bl <- grep("^BILAN:", lines, value = TRUE)
  stopifnot(length(bl) == 1L)
  cen <- sub("^CENSUS: ", "", grep("^CENSUS: test-", lines, value = TRUE))
  list(df = df[order(df$file), ], bilan = bl, census = sort(cen), lines = lines)
}
s <- read_releve(a[1]); h <- read_releve(a[2])
problems <- character(0)
if (!setequal(s$df$file, h$df$file)) {
  problems <- c(problems, sprintf("populations differ: serial=%d sharded=%d, missing=%s extra=%s",
    nrow(s$df), nrow(h$df),
    paste(setdiff(s$df$file, h$df$file), collapse=","),
    paste(setdiff(h$df$file, s$df$file), collapse=",")))
}
m <- merge(s$df, h$df, by = "file", suffixes = c("_s", "_h"))
for (tag in c("fail", "pass", "err", "skip")) {
  bad <- m[[paste0(tag, "_s")]] != m[[paste0(tag, "_h")]]
  if (any(bad)) {
    problems <- c(problems, sprintf("%s differs on %d files: %s", tag, sum(bad),
      paste(head(m$file[bad], 12L), collapse=",")))
  }
}
if (!identical(s$census, h$census)) problems <- c(problems, "CENSUS sets differ")
for (tag in c("failed", "passed", "error", "skipped")) {
  gs <- as.integer(sub(paste0("^.*", tag, "=(\\d+).*$"), "\\1", s$bilan))
  gh <- as.integer(sub(paste0("^.*", tag, "=(\\d+).*$"), "\\1", h$bilan))
  if (!identical(gs, gh)) problems <- c(problems, sprintf("BILAN %s differs: %d vs %d", tag, gs, gh))
}
el_s <- as.numeric(sub("^.*elapsed_s=([0-9.]+).*$", "\\1", s$bilan))
el_h <- as.numeric(sub("^.*elapsed_s=([0-9.]+).*$", "\\1", h$bilan))
report <- c(
  sprintf("EQUIVALENCE REPORT %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("serial  : %s", a[1]), sprintf("sharded : %s", a[2]),
  sprintf("files   : serial=%d sharded=%d (compared=%d)", nrow(s$df), nrow(h$df), nrow(m)),
  sprintf("BILAN s : %s", s$bilan), sprintf("BILAN h : %s", h$bilan),
  sprintf("wall    : serial=%.0fs sharded=%.0fs speedup=%.2fx", el_s, el_h, el_s / el_h),
  "",
  if (length(problems)) c("VERDICT: FALSIFIED — per-file equality FAILED", paste0("  - ", problems)) else
    c("VERDICT: EQUIVALENT — per-file (fail,pass,err,skip) identical on every file")
)
writeLines(report, a[3], useBytes = TRUE)
cat(paste(report, collapse="\n"), "\n")
quit(save = "no", status = if (length(problems)) 1L else 0L)
