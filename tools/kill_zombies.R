#!/usr/bin/env Rscript
# =============================================================================
# tools/kill_zombies.R — Inventaire et purge contrôlée des processus zombies R/Chrome
# =============================================================================
# Référence : docs/ideation/IDEATION_ZOMBIE_PROCESS_LIFECYCLE_2026-11-18.md (W1.4 étape 1)
#
# Par défaut : DRY-RUN uniquement (aucun processus terminé).
# Usage :
#   Rscript tools/kill_zombies.R
#   Rscript tools/kill_zombies.R --protect-ports=6573,8100
#   Rscript tools/kill_zombies.R --kill --protect-ports=6573,8100
# =============================================================================

parse_zombie_cli_args <- function(args = character(0)) {
  kill <- "--kill" %in% args
  protect_ports <- integer(0)
  protect_pids <- integer(0)
  help <- any(args %in% c("-h", "--help"))

  for (a in args) {
    if (startsWith(a, "--protect-ports=")) {
      raw <- sub("^--protect-ports=", "", a)
      vals <- suppressWarnings(as.integer(trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])))
      protect_ports <- unique(vals[!is.na(vals) & vals > 0L])
    } else if (startsWith(a, "--protect-pids=")) {
      raw <- sub("^--protect-pids=", "", a)
      vals <- suppressWarnings(as.integer(trimws(strsplit(raw, ",", fixed = TRUE)[[1L]])))
      protect_pids <- unique(vals[!is.na(vals) & vals > 0L])
    }
  }

  list(
    dry_run = !kill,
    kill = kill,
    protect_ports = protect_ports,
    protect_pids = protect_pids,
    help = help
  )
}

classify_zombie_processes <- function(procs,
                                      live_pids,
                                      listening_pids = integer(0),
                                      protect_pids = integer(0),
                                      self_pid = Sys.getpid()) {
  if (!is.data.frame(procs) || nrow(procs) == 0L) {
    return(data.frame(
      pid = integer(0),
      ppid = integer(0),
      name = character(0),
      status = character(0),
      reason = character(0),
      rss_mb = numeric(0),
      cmdline = character(0),
      stringsAsFactors = FALSE
    ))
  }

  pids <- as.integer(procs$pid)
  ppids <- as.integer(procs$ppid)
  names_vec <- as.character(procs$name)
  cmds <- as.character(procs$cmdline)
  cmds[is.na(cmds)] <- ""
  rss_bytes <- as.numeric(procs$rss_bytes)
  rss_bytes[is.na(rss_bytes)] <- 0

  # Protéger self_pid et sa chaîne d'ancêtres (anti-suicide)
  ancestor_pids <- integer(0)
  cur <- as.integer(self_pid)
  for (step in seq_len(32L)) {
    if (is.na(cur) || cur <= 0L || cur %in% ancestor_pids) break
    ancestor_pids <- c(ancestor_pids, cur)
    idx <- match(cur, pids)
    if (is.na(idx)) break
    cur <- ppids[[idx]]
  }

  protected_set <- unique(c(as.integer(protect_pids), as.integer(listening_pids), ancestor_pids))
  live_set <- unique(as.integer(live_pids))

  n <- nrow(procs)
  status <- character(n)
  reason <- character(n)

  for (i in seq_len(n)) {
    pid_i <- pids[[i]]
    ppid_i <- ppids[[i]]
    nm <- names_vec[[i]]
    cmd <- cmds[[i]]

    if (pid_i %in% protected_set) {
      status[[i]] <- "PROTECTED"
      reason[[i]] <- if (pid_i %in% ancestor_pids) {
        "self_or_ancestor"
      } else if (pid_i %in% listening_pids) {
        "listening_port"
      } else {
        "explicit_protect"
      }
      next
    }

    if (!is.na(ppid_i) && ppid_i %in% live_set) {
      status[[i]] <- "ALIVE"
      reason[[i]] <- "live_parent"
      next
    }

    is_r <- grepl("^R(script)?(\\.exe)?$", nm, ignore.case = TRUE)
    is_browser <- grepl("chrome|chromedriver", nm, ignore.case = TRUE) &&
      grepl("headless|chromedriver|remote-debugging-port", paste(nm, cmd), ignore.case = TRUE)

    if (is_r && grepl("mirai", cmd, ignore.case = TRUE)) {
      status[[i]] <- "ZOMBIE"
      reason[[i]] <- "mirai_orphan"
    } else if (is_r && grepl("Sys\\.getpid\\(|port", cmd, ignore.case = TRUE)) {
      status[[i]] <- "ZOMBIE"
      reason[[i]] <- "probe_orphan"
    } else if (is_browser) {
      status[[i]] <- "ZOMBIE"
      reason[[i]] <- "headless_browser_orphan"
    } else if (is_r) {
      status[[i]] <- "ZOMBIE"
      reason[[i]] <- "r_orphan"
    } else {
      status[[i]] <- "ALIVE"
      reason[[i]] <- "non_target"
    }
  }

  data.frame(
    pid = pids,
    ppid = ppids,
    name = names_vec,
    status = status,
    reason = reason,
    rss_mb = round(rss_bytes / (1024 * 1024), 1),
    cmdline = cmds,
    stringsAsFactors = FALSE
  )
}

summarize_zombie_candidates <- function(classified_df) {
  if (!is.data.frame(classified_df) || nrow(classified_df) == 0L) {
    return(list(n_total = 0L, n_zombies = 0L, zombie_mb = 0, n_protected = 0L, n_alive = 0L))
  }
  is_z <- classified_df$status == "ZOMBIE"
  list(
    n_total = as.integer(nrow(classified_df)),
    n_zombies = as.integer(sum(is_z)),
    zombie_mb = round(sum(classified_df$rss_mb[is_z], na.rm = TRUE), 1),
    n_protected = as.integer(sum(classified_df$status == "PROTECTED")),
    n_alive = as.integer(sum(classified_df$status == "ALIVE"))
  )
}

collect_windows_processes <- function() {
  ps_cmd <- paste(
    "Get-CimInstance Win32_Process |",
    "Select-Object ProcessId,ParentProcessId,Name,WorkingSetSize,CommandLine |",
    "ConvertTo-Csv -NoTypeInformation"
  )
  raw <- tryCatch(
    system2("powershell", c("-NoProfile", "-Command", ps_cmd), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  )
  if (length(raw) < 2L) {
    return(list(procs = data.frame(), live_pids = integer(0)))
  }
  df <- tryCatch(
    utils::read.csv(text = paste(raw, collapse = "\n"), stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (nrow(df) == 0L) {
    return(list(procs = data.frame(), live_pids = integer(0)))
  }
  all_pids <- unique(as.integer(df$ProcessId))
  all_pids <- all_pids[!is.na(all_pids)]

  keep <- grepl("^(R|Rscript|chrome|chromedriver)(\\.exe)?$", df$Name, ignore.case = TRUE)
  sub_df <- df[keep, , drop = FALSE]
  procs <- data.frame(
    pid = as.integer(sub_df$ProcessId),
    ppid = as.integer(sub_df$ParentProcessId),
    name = as.character(sub_df$Name),
    cmdline = as.character(sub_df$CommandLine),
    rss_bytes = as.numeric(sub_df$WorkingSetSize),
    stringsAsFactors = FALSE
  )
  list(procs = procs, live_pids = all_pids)
}

collect_listening_pids <- function(protect_ports = integer(0)) {
  ps_cmd <- "Get-NetTCPConnection -State Listen -ErrorAction SilentlyContinue | Select-Object LocalPort,OwningProcess | ConvertTo-Csv -NoTypeInformation"
  raw <- tryCatch(
    system2("powershell", c("-NoProfile", "-Command", ps_cmd), stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  )
  if (length(raw) < 2L) return(integer(0))
  df <- tryCatch(
    utils::read.csv(text = paste(raw, collapse = "\n"), stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (nrow(df) == 0L) return(integer(0))
  if (length(protect_ports) > 0L) {
    df <- df[as.integer(df$LocalPort) %in% protect_ports, , drop = FALSE]
  }
  unique(as.integer(df$OwningProcess[!is.na(df$OwningProcess)]))
}

main_kill_zombies <- function(args = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_zombie_cli_args(args)
  if (isTRUE(opts$help)) {
    cat("Usage: Rscript tools/kill_zombies.R [--kill] [--protect-ports=6573,8100] [--protect-pids=123]\n")
    return(invisible(0L))
  }

  if (!identical(.Platform$OS.type, "windows")) {
    stop("tools/kill_zombies.R ne supporte actuellement que Windows (Win32_Process).", call. = FALSE)
  }

  snap <- collect_windows_processes()
  listening_pids <- collect_listening_pids(opts$protect_ports)
  classified <- classify_zombie_processes(
    procs = snap$procs,
    live_pids = snap$live_pids,
    listening_pids = listening_pids,
    protect_pids = opts$protect_pids,
    self_pid = Sys.getpid()
  )
  summary <- summarize_zombie_candidates(classified)

  mode_lbl <- if (opts$dry_run) "DRY-RUN" else "KILL"
  cat(sprintf("=== kill_zombies.R [%s] ===\n", mode_lbl))
  cat(sprintf(
    "Inspectes: %d | ZOMBIE: %d (%.1f MB) | PROTECTED: %d | ALIVE: %d\n",
    summary$n_total, summary$n_zombies, summary$zombie_mb,
    summary$n_protected, summary$n_alive
  ))

  z <- classified[classified$status == "ZOMBIE", , drop = FALSE]
  if (nrow(z) > 0L) {
    for (i in seq_len(nrow(z))) {
      cat(sprintf(
        "  - PID=%d PPID=%d (%s) %.1f MB [%s] %s\n",
        z$pid[[i]], z$ppid[[i]], z$name[[i]], z$rss_mb[[i]],
        z$reason[[i]], substr(z$cmdline[[i]], 1L, 100L)
      ))
      if (isTRUE(opts$kill)) {
        suppressWarnings(system2("taskkill", c("/F", "/PID", as.character(z$pid[[i]])),
                                 stdout = FALSE, stderr = FALSE))
      }
    }
  }
  invisible(0L)
}

if (sys.nframe() == 0L) {
  quit(status = main_kill_zombies(), save = "no")
}
