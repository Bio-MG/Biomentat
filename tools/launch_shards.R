# tools/launch_shards.R — M-SH1 : orchestrateur des shards de la suite complète.
# Porteur du garde RAM, de l'ordre de lancement, du préfixe locale, du polling
# des fichiers shard ET de l'interdiction de juger sur les codes de sortie
# (fiche docs/proposals/SUITE_SHARDING_DESIGN_2026-11-18.md §3.2/§10.4).
#
# Usage (GO opérateur requis — BOOST) :
#   LC_ALL=French_France.65001 Rscript tools/launch_shards.R --n=2
#   --dry-run : imprime les commandes et le verdict RAM sans rien lancer.
# Le verdict d'un shard = présence d'une ligne ^BILAN: dans SON fichier, jamais
# son exit code (SIGSEGV 139 possible même en succès, AGENTS.md:21).
.localectl <- Sys.setlocale("LC_CTYPE", "fr_FR.UTF-8")
if (!nzchar(.localectl)) {
  warning("LC_CTYPE fr_FR.UTF-8 unavailable: parse() false failures are likely")
}
args <- commandArgs(trailingOnly = TRUE)
.flag <- function(name) {
  a <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^--", name, "="), "", a[1]) else NA_character_
}
n_spec <- .flag("n")
out_dir <- if (!is.na(.flag("out-dir"))) .flag("out-dir") else "QC"
timings <- if (!is.na(.flag("timings"))) .flag("timings") else "QC/full_suite_results.txt"
threshold <- suppressWarnings(as.numeric(if (!is.na(.flag("ram-gb"))) .flag("ram-gb") else "8"))
dry <- "--dry-run" %in% args
n <- suppressWarnings(as.integer(if (is.na(n_spec)) "2" else n_spec))
die <- function(msg) {
  cat(paste0("ERREUR launch_shards : ", msg, "\n"), file = stderr())
  quit(save = "no", status = 1L)
}
if (is.na(n) || n < 1L) die("--n doit etre un entier >= 1")
if (is.na(threshold) || threshold <= 0) die("--ram-gb doit etre un nombre > 0")

# --- GARDE RAM (mesure réelle, pas une constante) ---------------------------
ram_lines <- suppressWarnings(system2(
  "powershell",
  c("-NoProfile", "-Command",
    "Get-CimInstance Win32_OperatingSystem | ForEach-Object { 'FreeRAM_GB={0:N2}' -f ($_.FreePhysicalMemory/1MB) }"),
  stdout = TRUE))
free_txt <- grep("FreeRAM_GB=", ram_lines, value = TRUE)
free <- if (length(free_txt)) {
  suppressWarnings(as.numeric(gsub(",", ".", sub("FreeRAM_GB=", "", free_txt[1]))))
} else NA_real_
if (is.na(free)) die("impossible de mesurer la RAM libre")
cat(sprintf("RAM libre = %.2f Go (seuil %.2f Go) — %s\n", free, threshold,
            if (free >= threshold) "OK" else "REFUS"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
runner <- "tools/run_full_suite.R"
cmds <- lapply(seq_len(n), function(i) {
  out_i <- file.path(out_dir, sprintf("full_suite_results_shard_%d_of_%d.txt", i, n))
  log_i <- file.path(out_dir, sprintf("shard_%d_of_%d.log", i, n))
  c(shQuote(rscript), shQuote(runner), shQuote(paste0("--shard=", i, "/", n)),
    shQuote(paste0("--out=", out_i)), shQuote(paste0("--timings=", timings)))
})
if (dry) {
  for (i in seq_len(n)) cat(sprintf("[dry] shard %d/%d : %s\n", i, n,
                                    paste(cmds[[i]], collapse = " ")))
  if (free < threshold) quit(save = "no", status = 1L)
  quit(save = "no", status = 0L)
}
if (free < threshold) {
  die(sprintf("%.2f Go libres < %.2f Go : refus de lancer %d shards", free, threshold, n))
}

# Préfixe locale : l'enfant hérite des variables d'environnement (propriété déjà
# utilisée par ts_e2e_with_child_locale(), helper-app-driver.R:22-34).
Sys.setenv(LC_ALL = "French_France.65001")
Sys.setenv(LC_CTYPE = "French_France.65001")
t0 <- Sys.time()
handles <- vector("list", n)
for (i in seq_len(n)) {
  out_i <- file.path(out_dir, sprintf("full_suite_results_shard_%d_of_%d.txt", i, n))
  log_i <- file.path(out_dir, sprintf("shard_%d_of_%d.log", i, n))
  cat(sprintf("lancement shard %d/%d -> %s (log %s)\n", i, n, out_i, log_i))
  handles[[i]] <- processx::process$new(rscript, cmds[[i]],
                                        stdout = log_i, stderr = log_i)
  Sys.sleep(1)  # léger étalement : réduit le pic RAM au démarrage
}
# Polling sur LES FICHIERS, jamais sur les exit codes.
repeat {
  ready <- vapply(seq_len(n), function(i) {
    p <- file.path(out_dir, sprintf("full_suite_results_shard_%d_of_%d.txt", i, n))
    file.exists(p) && any(grepl("^BILAN:", readLines(p, warn = FALSE)))
  }, logical(1))
  if (all(ready)) break
  if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > 9000) {
    die(paste("timeout 9000s — shards sans BILAN :",
              paste(which(!ready), collapse = ",")))
  }
  Sys.sleep(5)
}
for (i in seq_len(n)) {
  try(handles[[i]]$wait(timeout = 1000), silent = TRUE)  # récolte, PAS verdict
}
cat(sprintf("Tous les shards ont écrit leur BILAN (%d/%d), makespan=%.0fs.\n",
            n, n, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat("Prochaine étape (hors de cet outil) : tools/merge_shard_results.R puis protocole d'équivalence §6.\n")
