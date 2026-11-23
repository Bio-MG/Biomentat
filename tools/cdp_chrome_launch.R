#!/usr/bin/env Rscript
# =============================================================================
# tools/cdp_chrome_launch.R — Loopback Chrome/Edge CDP Isolated-Profile Helper
# =============================================================================
# W2.4 (OBSCURA_CDP_BROWSER_CANDIDATE_2026-11-13.md O1): Chromium >=136 blocks
# --remote-debugging-port on the default OS user profile. This utility builds
# and launches an isolated temporary profile strictly bound to 127.0.0.1 for
# reproducible local E2E runs without touching the user's browser profile.
# =============================================================================

#' Locate an installed Chromium-based browser executable
#'
#' @param candidates Optional character vector of executable paths to check.
#'   When `NULL`, checks `Sys.getenv("CHROME_PATH")` and standard Windows
#'   Chrome/Edge/Brave install locations.
#' @return First existing file path, or `NA_character_` if none found.
#' @export
cdp_find_browser <- function(candidates = NULL) {
  if (is.null(candidates)) {
    env_p <- Sys.getenv("CHROME_PATH", unset = "")
    local_app <- Sys.getenv("LOCALAPPDATA", unset = "")
    pf <- Sys.getenv("ProgramFiles", unset = "C:/Program Files")
    pf86 <- Sys.getenv("ProgramFiles(x86)", unset = "C:/Program Files (x86)")
    candidates <- c(
      if (nzchar(env_p)) env_p,
      file.path(pf, "Google/Chrome/Application/chrome.exe"),
      file.path(pf86, "Google/Chrome/Application/chrome.exe"),
      file.path(local_app, "Google/Chrome/Application/chrome.exe"),
      file.path(pf86, "Microsoft/Edge/Application/msedge.exe"),
      file.path(pf, "Microsoft/Edge/Application/msedge.exe"),
      file.path(pf, "BraveSoftware/Brave-Browser/Application/brave.exe")
    )
  }
  candidates <- as.character(candidates)
  candidates <- candidates[!is.na(candidates) & nzchar(trimws(candidates))]
  for (p in candidates) {
    if (file.exists(p)) {
      return(normalizePath(p, winslash = "/", mustWork = FALSE))
    }
  }
  NA_character_
}

#' Build validated Chrome/Edge CLI flags for isolated loopback CDP
#'
#' @param port Integer TCP port in `[1024, 65535]` (default `9333L`).
#' @param user_data_dir Isolated temporary profile directory.
#' @param url Initial URL to open (default `"about:blank"`).
#' @param headless Logical; include `--headless=new` when `TRUE`.
#' @return Character vector of command-line arguments.
#' @export
cdp_build_chrome_args <- function(port = 9333L,
                                  user_data_dir = tempfile("biomentat_chrome_prof_"),
                                  url = "about:blank",
                                  headless = TRUE) {
  p_int <- suppressWarnings(as.integer(port))
  if (length(p_int) != 1L || is.na(p_int) || p_int < 1024L || p_int > 65535L) {
    stop("cdp_build_chrome_args(): 'port' must be an integer in [1024, 65535].", call. = FALSE)
  }
  if (is.null(user_data_dir) || !is.character(user_data_dir) ||
      length(user_data_dir) != 1L || is.na(user_data_dir) || !nzchar(trimws(user_data_dir))) {
    stop("cdp_build_chrome_args(): 'user_data_dir' must be a non-empty string.", call. = FALSE)
  }
  if (is.null(url) || !is.character(url) || length(url) != 1L ||
      is.na(url) || !nzchar(trimws(url))) {
    stop("cdp_build_chrome_args(): 'url' must be a non-empty string.", call. = FALSE)
  }

  prof_norm <- normalizePath(trimws(user_data_dir), winslash = "/", mustWork = FALSE)
  args <- c(
    "--remote-debugging-address=127.0.0.1",
    sprintf("--remote-debugging-port=%d", p_int),
    sprintf("--user-data-dir=%s", prof_norm),
    "--no-first-run",
    "--no-default-browser-check",
    "--disable-gpu",
    "--disable-extensions",
    "--disable-background-networking"
  )
  if (isTRUE(headless)) {
    args <- c(args, "--headless=new")
  }
  c(args, trimws(url))
}

#' Probe a local loopback CDP endpoint (/json/version)
#'
#' @param port Integer TCP port on `127.0.0.1`.
#' @param timeout_sec Numeric timeout in seconds (default `2`).
#' @return Named list with `ok` (logical), `port`, and `raw` (character).
#' @export
cdp_probe_endpoint <- function(port = 9333L, timeout_sec = 2) {
  p_int <- suppressWarnings(as.integer(port))
  if (length(p_int) != 1L || is.na(p_int) || p_int < 1L || p_int > 65535L) {
    return(list(ok = FALSE, port = NA_integer_, raw = character(0)))
  }
  target <- sprintf("http://127.0.0.1:%d/json/version", p_int)
  old_to <- getOption("timeout")
  on.exit(options(timeout = old_to), add = TRUE)
  options(timeout = timeout_sec)

  res <- tryCatch({
    con <- suppressWarnings(base::url(target, open = "rb"))
    on.exit(close(con), add = TRUE)
    lines <- readLines(con, warn = FALSE)
    list(ok = length(lines) > 0L, port = p_int, raw = lines)
  }, error = function(e) {
    list(ok = FALSE, port = p_int, raw = character(0))
  })
  res
}

#' Launch Chrome/Edge with an isolated temporary profile on 127.0.0.1
#'
#' @param port Integer TCP port in `[1024, 65535]`.
#' @param url Initial URL to open.
#' @param headless Logical; run headless when `TRUE`.
#' @param browser Path to browser executable.
#' @param user_data_dir Isolated temporary profile directory.
#' @return Handle list with `process`, `port`, `user_data_dir`, `browser`, `args`.
#' @export
cdp_launch_chrome <- function(port = 9333L,
                              url = "about:blank",
                              headless = TRUE,
                              browser = cdp_find_browser(),
                              user_data_dir = tempfile("biomentat_chrome_prof_")) {
  if (is.null(browser) || is.na(browser) || !nzchar(browser) || !file.exists(browser)) {
    stop("cdp_launch_chrome(): no Chromium/Chrome/Edge executable found.", call. = FALSE)
  }
  args <- cdp_build_chrome_args(port = port, user_data_dir = user_data_dir,
                                url = url, headless = headless)
  if (!dir.exists(user_data_dir)) {
    dir.create(user_data_dir, recursive = TRUE, showWarnings = FALSE)
  }

  proc <- if (requireNamespace("processx", quietly = TRUE)) {
    processx::process$new(command = browser, args = args, cleanup = TRUE)
  } else {
    system2(browser, args = args, wait = FALSE)
    NULL
  }

  list(
    process = proc,
    port = as.integer(port),
    user_data_dir = normalizePath(user_data_dir, winslash = "/", mustWork = FALSE),
    browser = browser,
    args = args
  )
}

#' Stop a launched CDP browser handle and optionally clean its temp profile
#'
#' @param handle Handle list returned by `cdp_launch_chrome()`.
#' @param clean_profile Logical; removes `user_data_dir` only when it resides
#'   inside `tempdir()`.
#' @return Logical `TRUE` invisibly.
#' @export
cdp_stop_chrome <- function(handle, clean_profile = TRUE) {
  if (is.list(handle) && !is.null(handle$process)) {
    tryCatch({
      if (is.function(handle$process$is_alive) && isTRUE(handle$process$is_alive())) {
        handle$process$kill()
      }
    }, error = function(e) NULL)
  }

  if (isTRUE(clean_profile) && is.list(handle) && !is.null(handle$user_data_dir) &&
      is.character(handle$user_data_dir) && length(handle$user_data_dir) == 1L &&
      nzchar(handle$user_data_dir) && dir.exists(handle$user_data_dir)) {
    tmp_root <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
    prof_norm <- normalizePath(handle$user_data_dir, winslash = "/", mustWork = FALSE)
    if (startsWith(tolower(prof_norm), tolower(tmp_root)) && !identical(tolower(prof_norm), tolower(tmp_root))) {
      unlink(prof_norm, recursive = TRUE, force = TRUE)
    }
  }
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  cat("cdp_chrome_launch.R loaded. Browser detected:", cdp_find_browser(), "\n")
}
