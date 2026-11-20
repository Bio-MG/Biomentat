# =============================================================================
# test-status-ledger.R — Tests RED-first du generateur et des gardes C6-C10
#                        (tools/status_ledger.py + docs/agent/status_ledger.json)
# =============================================================================

.sl_root <- local({
  d <- normalizePath(getwd(), mustWork = FALSE)
  repeat {
    if (file.exists(file.path(d, "AGENTS.md"))) break
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  d
})

.sl_script <- file.path(.sl_root, "tools", "status_ledger.py")
.sl_profile <- file.path(
  .sl_root, "docs", "agent", "skills", "scan-inventory",
  "references", "profile.example.json"
)
.sl_make_fixture <- file.path(
  .sl_root, "docs", "agent", "skills", "scan-inventory",
  "scripts", "make_fixture.py"
)

.sl_run_py <- function(args) {
  out <- tempfile("status-ledger-probe-")
  on.exit(unlink(out), add = TRUE)
  st <- suppressWarnings(system2(
    "python",
    c("-B", args),
    stdout = out, stderr = out, timeout = 120
  ))
  txt <- if (file.exists(out)) {
    paste(readLines(out, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  } else {
    ""
  }
  list(status = as.integer(st), out = txt)
}

.sl_run_ledger <- function(mode = "--check", root = .sl_root,
                           card = "docs/ROADMAP.md",
                           ledger = "docs/agent/status_ledger.json",
                           profile = .sl_profile) {
  if (!file.exists(.sl_script)) {
    return(list(status = 127L, out = "MISSING_SCRIPT: tools/status_ledger.py"))
  }
  .sl_run_py(c(
    shQuote(.sl_script),
    mode,
    "--root", shQuote(root),
    "--card", shQuote(card),
    "--ledger", shQuote(ledger),
    "--profile", shQuote(profile)
  ))
}

.sl_setup_fixtures <- function() {
  base <- tempfile("sl_fixtures_")
  dir.create(base, recursive = TRUE, showWarnings = FALSE)
  res <- .sl_run_py(c(shQuote(.sl_make_fixture), "--out", shQuote(base)))
  list(base = base, init_status = res$status)
}

test_that("tools/status_ledger.py existe et --build puis --check passent sur fixture clean + idempotence", {
  expect_true(file.exists(.sl_script))
  fx <- .sl_setup_fixtures()
  on.exit(unlink(fx$base, recursive = TRUE, force = TRUE), add = TRUE)
  expect_identical(fx$init_status, 0L)

  clean_root <- file.path(fx$base, "clean")
  b1 <- .sl_run_ledger("--build", root = clean_root)
  expect_identical(b1$status, 0L)

  ledger_p <- file.path(clean_root, "docs", "agent", "status_ledger.json")
  card_p <- file.path(clean_root, "docs", "ROADMAP.md")
  expect_true(file.exists(ledger_p))
  expect_true(file.exists(card_p))

  if (file.exists(ledger_p) && file.exists(card_p)) {
    hash_ledger_1 <- unname(tools::md5sum(ledger_p))
    hash_card_1 <- unname(tools::md5sum(card_p))
    card_lines <- readLines(card_p, warn = FALSE, encoding = "UTF-8")
    expect_true(any(grepl("<!-- LEDGER:BEGIN -->", card_lines, fixed = TRUE)))
    expect_true(any(grepl("<!-- LEDGER:END -->", card_lines, fixed = TRUE)))

    # Deuxieme --build : doit etre 100% octet-identique (idempotence)
    b2 <- .sl_run_ledger("--build", root = clean_root)
    expect_identical(b2$status, 0L)
    expect_identical(unname(tools::md5sum(ledger_p)), hash_ledger_1)
    expect_identical(unname(tools::md5sum(card_p)), hash_card_1)
  }

  chk <- .sl_run_ledger("--check", root = clean_root)
  expect_identical(chk$status, 0L)
})

test_that("C6 (vocabulaire ferme) : un etat hors enum (STATE:BOGUS) echoue en exit 1 avec C6", {
  expect_true(file.exists(.sl_script))
  fx <- .sl_setup_fixtures()
  on.exit(unlink(fx$base, recursive = TRUE, force = TRUE), add = TRUE)

  clean_root <- file.path(fx$base, "clean")
  .sl_run_ledger("--build", root = clean_root)

  card_p <- file.path(clean_root, "docs", "ROADMAP.md")
  if (file.exists(card_p)) {
    txt <- readLines(card_p, warn = FALSE, encoding = "UTF-8")
    txt <- sub("STATE:ACTIVE", "STATE:BOGUS", txt, fixed = TRUE)
    writeLines(txt, card_p, useBytes = TRUE)
  }

  chk <- .sl_run_ledger("--check", root = clean_root)
  expect_identical(chk$status, 1L)
  expect_match(chk$out, "C6", fixed = TRUE)
})

test_that("C7 (TERMINATED uniquement dans archive/) : STATE:TERMINATED dans proposals/ echoue en exit 1 avec C7", {
  expect_true(file.exists(.sl_script))
  fx <- .sl_setup_fixtures()
  on.exit(unlink(fx$base, recursive = TRUE, force = TRUE), add = TRUE)

  clean_root <- file.path(fx$base, "clean")
  .sl_run_ledger("--build", root = clean_root)

  card_p <- file.path(clean_root, "docs", "ROADMAP.md")
  if (file.exists(card_p)) {
    txt <- readLines(card_p, warn = FALSE, encoding = "UTF-8")
    txt <- sub("STATE:ACTIVE", "STATE:TERMINATED", txt, fixed = TRUE)
    writeLines(txt, card_p, useBytes = TRUE)
  }

  chk <- .sl_run_ledger("--check", root = clean_root)
  expect_identical(chk$status, 1L)
  expect_match(chk$out, "C7", fixed = TRUE)
})

test_that("C8 (ancrage et preuves obligatoires) : DELIVERED sans preuve et PARTIAL avec 0 jalon clos echouent en exit 1 avec C8", {
  expect_true(file.exists(.sl_script))
  fx <- .sl_setup_fixtures()
  on.exit(unlink(fx$base, recursive = TRUE, force = TRUE), add = TRUE)

  dirty_root <- file.path(fx$base, "dirty")
  .sl_run_ledger("--build", root = dirty_root)
  chk_dirty <- .sl_run_ledger("--check", root = dirty_root)
  expect_identical(chk_dirty$status, 1L)
  expect_match(chk_dirty$out, "C8", fixed = TRUE)

  # Cas isole : DELIVERED sans preuve sur fixture clean
  clean_root <- file.path(fx$base, "clean")
  .sl_run_ledger("--build", root = clean_root)
  card_p <- file.path(clean_root, "docs", "ROADMAP.md")
  if (file.exists(card_p)) {
    txt <- readLines(card_p, warn = FALSE, encoding = "UTF-8")
    txt <- sub(
      "STATE:ACTIVE \u00b7 JALONS:0/2 \u00b7 NEXT:GO M1 \u00b7 FLAGS:needs_go",
      "STATE:DELIVERED \u00b7 JALONS:2/2",
      txt,
      fixed = TRUE
    )
    writeLines(txt, card_p, useBytes = TRUE)
  }
  chk_deliv <- .sl_run_ledger("--check", root = clean_root)
  expect_identical(chk_deliv$status, 1L)
  expect_match(chk_deliv$out, "C8", fixed = TRUE)

  # Cas isole : PARTIAL avec 0/3 jalon clos sur fixture clean
  if (file.exists(card_p)) {
    txt <- readLines(card_p, warn = FALSE, encoding = "UTF-8")
    txt <- sub(
      "STATE:DELIVERED \u00b7 JALONS:2/2",
      "STATE:PARTIAL \u00b7 JALONS:0/3",
      txt,
      fixed = TRUE
    )
    writeLines(txt, card_p, useBytes = TRUE)
  }
  chk_partial <- .sl_run_ledger("--check", root = clean_root)
  expect_identical(chk_partial$status, 1L)
  expect_match(chk_partial$out, "C8", fixed = TRUE)
})

test_that("C9 (population > 0) : famille requise vide ou registre vide echoue en exit 1 avec C9", {
  expect_true(file.exists(.sl_script))
  fx <- .sl_setup_fixtures()
  on.exit(unlink(fx$base, recursive = TRUE, force = TRUE), add = TRUE)

  empty_root <- file.path(fx$base, "empty")
  .sl_run_ledger("--build", root = empty_root)
  chk_empty <- .sl_run_ledger("--check", root = empty_root)
  expect_identical(chk_empty$status, 1L)
  expect_match(chk_empty$out, "C9", fixed = TRUE)

  # Registre JSON vide {"items": []} sur fixture clean
  clean_root <- file.path(fx$base, "clean")
  .sl_run_ledger("--build", root = clean_root)
  ledger_p <- file.path(clean_root, "docs", "agent", "status_ledger.json")
  if (file.exists(ledger_p)) {
    writeLines("{\"items\": []}", ledger_p, useBytes = TRUE)
  }
  chk_zero <- .sl_run_ledger("--check", root = clean_root)
  expect_identical(chk_zero$status, 1L)
  expect_match(chk_zero$out, "C9", fixed = TRUE)
})

test_that("C10 (accord ROADMAP.md <-> status_ledger.json) : desaccord d'etat ou de bloc LEDGER echoue en exit 1 avec C10", {
  expect_true(file.exists(.sl_script))
  fx <- .sl_setup_fixtures()
  on.exit(unlink(fx$base, recursive = TRUE, force = TRUE), add = TRUE)

  clean_root <- file.path(fx$base, "clean")
  .sl_run_ledger("--build", root = clean_root)

  ledger_p <- file.path(clean_root, "docs", "agent", "status_ledger.json")
  if (file.exists(ledger_p)) {
    txt <- readLines(ledger_p, warn = FALSE, encoding = "UTF-8")
    txt <- sub("\"state\": \"ACTIVE\"", "\"state\": \"FROZEN\"", txt, fixed = TRUE)
    writeLines(txt, ledger_p, useBytes = TRUE)
  }

  chk <- .sl_run_ledger("--check", root = clean_root)
  expect_identical(chk$status, 1L)
  expect_match(chk$out, "C10", fixed = TRUE)
})

test_that("Depot reel : tools/status_ledger.py --check passe en exit 0 et le registre JSON est valide", {
  expect_true(file.exists(.sl_script))
  ledger_p <- file.path(.sl_root, "docs", "agent", "status_ledger.json")
  expect_true(file.exists(ledger_p))

  chk <- .sl_run_ledger("--check", root = .sl_root)
  expect_identical(chk$status, 0L)
})
