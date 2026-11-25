# =============================================================================
# test-kill-zombies.R — W1.4 controlled zombie process inventory & classifier
# =============================================================================

source_project_file("tools/kill_zombies.R")

test_that("classify_zombie_processes keeps self, ancestors, protected PIDs, and live-parent workers", {
  procs <- data.frame(
    pid = c(100L, 101L, 102L, 200L, 201L, 300L, 301L, 400L),
    ppid = c(1L, 100L, 9999L, 9998L, 9997L, 9996L, 100L, 9995L),
    name = c("Rscript.exe", "Rscript.exe", "Rscript.exe", "Rscript.exe",
             "Rscript.exe", "chrome.exe", "chromedriver.exe", "Rscript.exe"),
    cmdline = c(
      "Rscript app.R",
      "Rscript -e mirai::daemon('ipc://abc')",
      "Rscript -e mirai::daemon('ipc://dead')",
      "Rscript -e \"cat(Sys.getpid(), file='probe.pid')\"",
      "Rscript run_shiny.R",
      "chrome.exe --headless --remote-debugging-port=9222",
      "chromedriver.exe --port=9515",
      "Rscript tools/kill_zombies.R"
    ),
    rss_bytes = c(500e6, 600e6, 634e6, 2e6, 1800e6, 300e6, 40e6, 50e6),
    stringsAsFactors = FALSE
  )
  live_pids <- c(1L, procs$pid)

  res <- classify_zombie_processes(
    procs = procs,
    live_pids = live_pids,
    listening_pids = c(100L, 201L),
    protect_pids = c(201L),
    self_pid = 400L
  )

  expect_true(is.data.frame(res))
  expect_true(all(c("pid", "status", "reason", "rss_mb") %in% names(res)))

  by_pid <- setNames(res$status, as.character(res$pid))
  expect_identical(unname(by_pid["100"]), "PROTECTED") # listening & live parent
  expect_identical(unname(by_pid["101"]), "ALIVE")     # mirai with live parent 100
  expect_identical(unname(by_pid["102"]), "ZOMBIE")    # orphan mirai daemon
  expect_identical(unname(by_pid["200"]), "ZOMBIE")    # orphan probe worker
  expect_identical(unname(by_pid["201"]), "PROTECTED") # explicitly protected / listening
  expect_identical(unname(by_pid["300"]), "ZOMBIE")    # orphan headless chrome
  expect_identical(unname(by_pid["301"]), "ALIVE")     # chromedriver with live parent 100
  expect_identical(unname(by_pid["400"]), "PROTECTED") # self_pid always protected
})

test_that("summarize_zombie_candidates computes count and total MB accurately", {
  res <- data.frame(
    pid = c(102L, 200L, 100L),
    status = c("ZOMBIE", "ZOMBIE", "PROTECTED"),
    rss_mb = c(634.0, 2.0, 500.0),
    stringsAsFactors = FALSE
  )
  s <- summarize_zombie_candidates(res)
  expect_identical(s$n_zombies, 2L)
  expect_equal(s$zombie_mb, 636.0)
  expect_identical(s$n_protected, 1L)
})

test_that("parse_zombie_cli_args defaults to dry_run = TRUE and parses --kill and --protect-ports", {
  d <- parse_zombie_cli_args(character(0))
  expect_true(d$dry_run)
  expect_false(d$kill)
  expect_identical(d$protect_ports, integer(0))

  k <- parse_zombie_cli_args(c("--kill", "--protect-ports=6573,8100"))
  expect_false(k$dry_run)
  expect_true(k$kill)
  expect_false(k$include_workrsock)
  expect_identical(k$protect_ports, c(6573L, 8100L))

  w <- parse_zombie_cli_args(c("--include-workrsock"))
  expect_true(w$include_workrsock)
})

test_that("classify_zombie_processes protects workRSOCK workers by default", {
  procs <- data.frame(
    pid = c(500L),
    ppid = c(9999L),
    name = c("Rscript.exe"),
    cmdline = c("Rscript -e parallel:::.slaveRSOCK() workRSOCK"),
    rss_bytes = c(2e6),
    stringsAsFactors = FALSE
  )
  res_def <- classify_zombie_processes(procs, live_pids = integer(0), self_pid = 1L)
  expect_identical(res_def$status[[1L]], "PROTECTED")
  expect_identical(res_def$reason[[1L]], "workrsock_protected")

  res_inc <- classify_zombie_processes(procs, live_pids = integer(0), self_pid = 1L,
                                       include_workrsock = TRUE)
  expect_identical(res_inc$status[[1L]], "ZOMBIE")
})
