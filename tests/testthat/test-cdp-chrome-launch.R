# ==============================================================================
# test-cdp-chrome-launch.R — W2.4 Isolated-Profile CDP Launcher
# ==============================================================================

.cdp_test_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("cdp_build_chrome_args builds isolated loopback flags and validates inputs", {
  source(.cdp_test_root("tools/cdp_chrome_launch.R"), local = globalenv())

  prof <- file.path(tempdir(), "prof_test_1")
  args_hl <- cdp_build_chrome_args(port = 9333L, user_data_dir = prof,
                                   url = "http://127.0.0.1:3838", headless = TRUE)
  expect_true("--remote-debugging-address=127.0.0.1" %in% args_hl)
  expect_true("--remote-debugging-port=9333" %in% args_hl)
  expect_true(any(startsWith(args_hl, "--user-data-dir=")))
  expect_true("--headless=new" %in% args_hl)
  expect_equal(utils::tail(args_hl, 1L), "http://127.0.0.1:3838")

  args_vis <- cdp_build_chrome_args(port = 9334L, user_data_dir = prof, headless = FALSE)
  expect_false("--headless=new" %in% args_vis)

  expect_error(cdp_build_chrome_args(port = 80L))
  expect_error(cdp_build_chrome_args(port = 99999L))
  expect_error(cdp_build_chrome_args(port = NA_integer_))
  expect_error(cdp_build_chrome_args(port = 9333L, user_data_dir = ""))
})

test_that("cdp_find_browser, cdp_probe_endpoint, and cdp_stop_chrome behave safely offline", {
  source(.cdp_test_root("tools/cdp_chrome_launch.R"), local = globalenv())

  fake_exe <- tempfile(fileext = ".exe")
  file.create(fake_exe)
  on.exit(unlink(fake_exe), add = TRUE)

  expect_equal(
    cdp_find_browser(candidates = c(file.path(tempdir(), "missing_browser.exe"), fake_exe)),
    normalizePath(fake_exe, winslash = "/", mustWork = FALSE)
  )
  expect_true(is.na(cdp_find_browser(candidates = file.path(tempdir(), "nonexistent_999.exe"))))

  pr <- cdp_probe_endpoint(port = 59999L, timeout_sec = 0.2)
  expect_false(pr$ok)

  tmp_prof <- tempfile("cdp_prof_cleanup_")
  dir.create(tmp_prof, recursive = TRUE)
  expect_true(dir.exists(tmp_prof))
  cdp_stop_chrome(list(process = NULL, user_data_dir = tmp_prof), clean_profile = TRUE)
  expect_false(dir.exists(tmp_prof))

  # Refuses to delete paths outside tempdir()
  outside_dir <- normalizePath(".", winslash = "/", mustWork = FALSE)
  cdp_stop_chrome(list(process = NULL, user_data_dir = outside_dir), clean_profile = TRUE)
  expect_true(dir.exists(outside_dir))
})
