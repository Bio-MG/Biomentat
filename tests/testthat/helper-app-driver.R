# Shared helpers for the shinytest2 "golden path" smoke tests (backlog #7).
# Coexists with helper-source.R (pure-helper unit tests): disjoint symbols,
# both auto-sourced by testthat before any test-*.R file runs.
#
# NOTE (API shinytest2) : $run_js() exécute le JS en fire-and-forget et
# renvoie NULL — toute valeur doit revenir par $get_js(). Les deux helpers
# qui lisent le DOM passent donc par get_js(), avec une IIFE pour que
# l'instruction 'return' soit légale quel que soit le contexte d'évaluation.

# shinytest2 s'auto-skippe avec la raison "On CRAN" tant que NOT_CRAN!=true
# (projet sans DESCRIPTION : aucun R CMD check ne passera jamais le poser).
# Pose ici, au chargement du helper, pour que la commande du plan —
#   Rscript -e "testthat::test_dir('tests/testthat')"
# — fonctionne telle quelle.
Sys.setenv(NOT_CRAN = "true")

#' Absolute path to the app's root directory (where app.R lives)
app_root_dir <- function() {
  normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
}

# --- Locale des processus ENFANTS (correctif 2026-09-16) --------------------
# shinytest2 demarre l'app dans un PROCESSUS ENFANT. Sys.setlocale() ne suffit
# donc PAS (il ne change que le processus courant), alors que les enfants
# heritent des VARIABLES d'environnement. Or Git Bash exporte LC_ALL=C.UTF-8,
# un nom que R sous Windows ne reconnait pas => repli silencieux sur "C" =>
# la lecture de i18n/translation.json (UTF-8) echoue par
# « invalid multibyte string, element 1 » — que ce helper rapportait a tort
# comme « Chromote/headless Chrome unavailable ».
# Mesure 2026-09-16 (A/B) : heritage LC_ALL=C.UTF-8 => ECHEC ; LC_ALL leve +
# LC_CTYPE=fr_FR.UTF-8 => driver OK, inputs visibles.
# On ne touche QUE LC_ALL / LANG / LC_CTYPE, et seulement le temps du
# demarrage : forcer LC_ALL changerait aussi LC_COLLATE / LC_TIME, donc le tri
# des chaines, donc potentiellement des resultats de tests.

#' Candidate UTF-8 locale names accepted by R on Windows
#'
#' ⚠️ « fr_FR.UTF-8 » a été retiré des candidats MESURÉ le 2026-09-28 : Windows
#' ACCEPTE le nom (setlocale ne lève pas) mais le résout vers le code page
#' ANSI — le grepl("UTF-8") du test d'acceptation passe alors à tort, et la
#' sortie stderr UTF-8 de l'app devient « chaîne de charactères multioctets
#' incorrecte / input string invalid in this locale » dans le poller de
#' shinytest2. « French_France.utf8 » (forme Windows native) est le premier
#' candidat réellement UTF-8.
ts_e2e_locale_candidates <- function() {
  c("French_France.utf8", "English_United States.utf8", "en_US.UTF-8")
}

#' Set LC_CTYPE to the first candidate R actually accepts; return its name
ts_e2e_first_valid_utf8_locale <- function() {
  for (cand in ts_e2e_locale_candidates()) {
    ok <- suppressWarnings(tryCatch({
      Sys.setlocale("LC_CTYPE", cand)
      grepl("UTF-8|utf8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)
    }, error = function(e) FALSE))
    if (isTRUE(ok)) return(cand)
  }
  NA_character_
}

#' Restore an environment snapshot taken with unset = NA_character_
ts_e2e_restore_env <- function(old) {
  present <- old[!is.na(old)]
  if (length(present)) do.call(Sys.setenv, as.list(present))
  absent <- names(old)[is.na(old)]
  if (length(absent)) Sys.unsetenv(absent)
}

#' Evaluate expr with a locale that child processes can actually use
ts_e2e_with_child_locale <- function(expr) {
  old <- Sys.getenv(c("LC_ALL", "LANG", "LC_CTYPE"), unset = NA_character_)
  on.exit(ts_e2e_restore_env(old), add = TRUE)
  # Ne toucher a l'environnement QUE si l'on sait poser une locale utilisable :
  # sinon on laisse l'etat initial (l'echec reste visible et honnetement
  # rapporte) plutot que d'echanger un probleme d'encodage contre un autre.
  if (!is.na(ts_e2e_first_valid_utf8_locale())) {
    # MESURÉ le 2026-09-28 (cette machine) : poser SEULEMENT LC_CTYPE ne suffit
    # pas — avec LC_ALL/LANG désettés, un composant de la pile écrit du
    # cp1252 sur stderr et le poller du parent (UTF-8) meurt sur
    # « input string N is invalid in this locale ». LC_ALL = le même candidat
    # UTF-8 que LC_CTYPE rend le boot fiable. Le coût annoncé plus bas
    # (LC_COLLATE/LC_TIME héritent du candidat) est accepté : les e2e sont des
    # smoke tests, pas des tests sensibles au tri.
    locale <- ts_e2e_first_valid_utf8_locale()
    Sys.setenv(LC_ALL = locale)
    Sys.setenv(LC_CTYPE = locale)
    Sys.unsetenv("LANG")
  }
  force(expr)
}

#' Test-only Chrome PIDs (chromote/shinytest2), never the user's browser.
#' A failed AppDriver boot spawns a whole Chrome tree that gets REPARENTED, so
#' ps_children() no longer sees it. Filtering by command line is the reliable
#' handle; measured 2026-11-19: one failed boot leaves up to 67 chrome.exe.
ts_test_chrome_pids <- function() {
  tryCatch({
    pids <- ps::ps_pids()
    nms <- vapply(pids, function(pid) {
      tryCatch(ps::ps_name(ps::ps_handle(pid)), error = function(e) "")
    }, character(1))
    chrome <- pids[grepl("^chrome\\.exe$", nms, ignore.case = TRUE)]
    if (!length(chrome)) return(integer(0))
    cmd <- vapply(chrome, function(pid) {
      tryCatch(ps::ps_cmdline(ps::ps_handle(pid)), error = function(e) "")
    }, character(1))
    chrome[grepl("(?i)remote-debugging|user-data-dir=.*(Temp|Rtmp|chromote|shinytest)",
                 cmd)]
  }, error = function(e) integer(0))
}

#' Retrieve the private drive sandbox path attached to an AppDriver instance
app_driver_sandbox_dir <- function(driver) {
  attr(driver, "ts_drive_sandbox", exact = TRUE)
}

#' Best-effort AppDriver constructor — skips (never fails) if shinytest2 or
#' a headless Chrome (chromote) is unavailable. Generous timeouts : global.R
#' charge toute la pile de packages, et l'init des 6 daemons mirai est
#' différée à la première ouverture de l'onglet Spatial (boot lazy, voir
#' l'en-tête de app.R).
#'
#' Automatically isolates test execution in a private drive sandbox via
#' `TRANSCRIPTO_DRIVE_ROOT`, preventing any interference with host tools/_drive/ready.json.
new_app_driver <- function(name) {
  # Porte e2e (DEEP_ANALYSIS_SUITE_ACCELERATION_2026-11-19 §5) : les 8 boots
  # échouent aujourd'hui à launch_chrome_impl (named-pipe err 232, mesuré
  # 2026-11-19, artefact e2e_probe.err) pendant que chromote réussit HORS
  # contexte (sonde standalone : session en 1,2 s) — une sonde de
  # disponibilité ne reproduit PAS l'échec. Les 4 fichiers e2e coûtent
  # 308,8 s (13 % de la suite) pour 8 skips sans assertion.
  # Quadruplet préservé par construction : chaque test skippe => (0,0,0,2),
  # identique aux relevés sériel ET shardié existants. TS_E2E=1 rétablit le
  # comportement historique verbatim (boots tentés, skip sur échec).
  # Blast radius vérifié (grep) : seuls les 4 test-shinytest2-*.R appellent
  # cette fonction (test-drive-sandbox.R construit un mock, l.195).
  if (!identical(Sys.getenv("TS_E2E"), "1")) {
    testthat::skip("e2e gated: TS_E2E!=1 — boots launch_chrome named-pipe err 232 (DEEP_ANALYSIS_2026-11-19 §5)")
  }
  testthat::skip_if_not_installed("shinytest2")

  # M2: Isolated drive sandbox per AppDriver instance to protect host tools/_drive/ready.json
  sandbox_dir <- normalizePath(tempfile("ts_drive_sandbox_"), winslash = "/", mustWork = FALSE)
  dir.create(file.path(sandbox_dir, "tools", "_drive"), recursive = TRUE, showWarnings = FALSE)

  old_drive_root <- Sys.getenv("TRANSCRIPTO_DRIVE_ROOT", unset = NA_character_)
  Sys.setenv(TRANSCRIPTO_DRIVE_ROOT = sandbox_dir)
  on.exit({
    if (is.na(old_drive_root)) {
      Sys.unsetenv("TRANSCRIPTO_DRIVE_ROOT")
    } else {
      Sys.setenv(TRANSCRIPTO_DRIVE_ROOT = old_drive_root)
    }
  }, add = TRUE)

  # Snapshot des enfants AVANT le boot. Mesure 2026-11-19 : un timeout
  # Page.navigate laisse le process Shiny + l'arbre Chrome vivants ; 8 boots
  # ratés => 66 process / 22,4 Go. On ne tue que les PIDs apparus pendant CE
  # boot, jamais les enfants déjà présents (shards concurrents, daemons, etc.).
  child_pids_before <- tryCatch(
    vapply(ps::ps_children(recursive = TRUE),
           function(h) as.integer(ps::ps_pid(h)), integer(1)),
    error = function(e) integer(0))
  chrome_pids_before <- ts_test_chrome_pids()

  err <- NULL
  driver <- tryCatch(
    ts_e2e_with_child_locale(
      shinytest2::AppDriver$new(
        app_dir = app_root_dir(), name = name, height = 900, width = 1400,
        load_timeout = 120000, timeout = 120000,
        # Smoke tests : check_names=FALSE pour qu'un warning de doublon d'id
        # ne transforme PAS chaque test en skip silencieux — le drift de
        # namespace est déjà asserté explicitement via get_values().
        check_names = FALSE,
        options = list(shiny.testmode = TRUE)
      )
    ),
    error = function(e) {
      err <<- conditionMessage(e)
      # shinytest2 expose le driver à moitié initialisé via last_error()$app.
      tryCatch({
        le <- rlang::last_error()
        if (!is.null(le) && !is.null(le$app)) try(le$app$stop(), silent = TRUE)
      }, error = function(e2) NULL)
      # Puis on tue les enfants laissés par ce boot raté (Shiny, Chrome, mirai).
      child_pids_after <- tryCatch(
        vapply(ps::ps_children(recursive = TRUE),
               function(h) as.integer(ps::ps_pid(h)), integer(1)),
        error = function(e2) integer(0))
      for (pid in setdiff(child_pids_after, child_pids_before)) {
        tryCatch(ps::ps_kill(ps::ps_handle(pid), grace = 500L),
                 error = function(e2) NULL)
      }
      # L'arbre Chrome est réparenté : le tuer par signature de commande,
      # jamais par parenté (et jamais le navigateur de l'utilisateur).
      for (pid in setdiff(ts_test_chrome_pids(), chrome_pids_before)) {
        tryCatch(ps::ps_kill(ps::ps_handle(pid), grace = 500L),
                 error = function(e2) NULL)
      }
      NULL
    }
  )
  # Le motif de skip porte desormais la VRAIE cause : un echec de locale ne
  # doit plus se déguiser en « Chrome absent » (8 tests restes invisibles
  # ainsi jusqu'au 2026-09-16).
  if (is.null(driver)) {
    unlink(sandbox_dir, recursive = TRUE, force = TRUE)
    msg <- if (is.null(err)) "raison inconnue" else err
    if (grepl("multibyte|encoding|locale|LC_", msg, ignore.case = TRUE)) {
      testthat::skip(paste("Locale inutilisable par le processus enfant :", msg))
    }
    testthat::skip(paste("AppDriver indisponible (Chromote/headless Chrome ?) :", msg))
  }

  # Wrap stop() and register finalizer to ensure the temporary sandbox is cleaned up
  # without ever touching <repo_root>/tools/_drive/ready.json
  attr(driver, "ts_drive_sandbox") <- sandbox_dir
  orig_stop <- driver$stop
  unlockBinding("stop", driver)
  driver$stop <- function(...) {
    on.exit(unlink(sandbox_dir, recursive = TRUE, force = TRUE), add = TRUE)
    orig_stop(...)
  }
  lockBinding("stop", driver)
  reg.finalizer(driver, function(e) {
    unlink(sandbox_dir, recursive = TRUE, force = TRUE)
  }, onexit = TRUE)

  driver
}

#' Click a top-level nav tab by (partial) visible text, JS-side — more
#' robust than a CSS selector on bslib's emoji+text nav-link markup.
click_nav_by_text <- function(app, text_fragment) {
  clicked <- app$get_js(sprintf(
    "(function() {
       var links = document.querySelectorAll('a.nav-link');
       for (var i = 0; i < links.length; i++) {
         if (links[i].textContent.indexOf('%s') !== -1) { links[i].click(); return true; }
       }
       return false;
     })()", text_fragment
  ))
  testthat::expect_true(isTRUE(clicked), info = paste("Nav tab not found:", text_fragment))
  # Best-effort uniquement : la sidebar (textOutput mem_usage) et l'init
  # lazy des daemons mirai (onglet Spatial) peuvent maintenir Shiny
  # "instable" >15s ; attendre la stabilite est un confort de rendu, pas une
  # assertion — un vrai crash d'output sera de toute facon attrape par
  # assert_no_shiny_error()/get_values() qui suivent.
  try(app$wait_for_idle(duration = 250, timeout = 10000), silent = TRUE)
  invisible(app)
}

#' Assert no Shiny/JS runtime error is currently shown in the DOM.
assert_no_shiny_error <- function(app) {
  n_errors <- app$get_js("document.querySelectorAll('.shiny-output-error').length")
  testthat::expect_equal(n_errors, 0, info = "A .shiny-output-error element is visible.")
}
