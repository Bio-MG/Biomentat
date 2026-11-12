# =============================================================================
# test-mod-spatial-deconv.R — P-lane deconv (2026-11-12) : garde pool mirai +
# journal honnête pour modules/spatial/deconv/mod_spatial_deconv.R
# =============================================================================
# Diagnostic LIVE (instance opérateur 6561, soirée 2026-11-12) : les 4
# dernières sessions de déconvolution ont un journal à 0 octet —
# `reset_log()` tronque au clic, puis le worker mirai (mort, RAM : référence
# 73260x10237 + pipeline complet + 4 workers multisession + daemons mirai)
# n'écrit JAMAIS la 1re ligne. Le module orchestrator dispatchait sans AUCUN
# appel de `ensure_app_daemons()` : le pool mort n'était jamais détecté ni
# respawné (le watchdog d'init_spatial_daemons, changelog v7, ne court que si
# ON l'appelle) — timeout muet de 20-45 min puis « voir le log » sur un
# fichier vide.
#
# Contrats (ROUGE d'abord, règle 5) :
#   (1) pool indisponible => invoke NON appelé + notification type "error"
#       honnête, et le journal N'EST PAS tronqué par un clic qui ne part pas ;
#   (2) statut "error" => le journal reçoit le VRAI texte d'erreur de la
#       tâche, horodaté par le `write_mirai_log` RÉEL.
#
# Idiome : le VRAI bloc `observeEvent` est extrait par AST (helper-ast.R) et
# levé dans un environnement enfant — jamais re-copié (cf. la série D4(a) de
# test-mod-bulk-de-run.R). `mirai::daemons_set()` est un appel de NAMESPACE
# (non doublable) : le chemin « pool OK » arme un VRAI daemon.
# =============================================================================

.MDC_FILE <- "modules/spatial/deconv/mod_spatial_deconv.R"

# write_mirai_log RÉEL : le test de journal ne prouverait rien avec une copie.
source_project_file("tools/pkg_skeleton/R/spatial_async.R")

# Les constantes du bloc (définies AU SOMMET du module ou de son prep) :
# importées par AST depuis le dépôt, jamais re-copiées (cf. .mdr_env).
e_consts <- new.env(parent = globalenv())
e_consts$DECONV_DEFAULT_N_HVG <- ts_ast_assignment(
  "tools/pkg_skeleton/R/spatial_deconv_prep.R", "DECONV_DEFAULT_N_HVG", e_consts)
.mdc_const <- function(name) ts_ast_assignment(.MDC_FILE, name, e_consts)

.mdc_lift <- function(envir, needle) {
  f <- function() NULL
  body(f) <- ts_ast_observe_block(.MDC_FILE, needle)
  environment(f) <- envir
  f
}

.mdc_run <- function(e, needle) {
  f <- .mdc_lift(e, needle)
  shiny::withReactiveDomain(shiny::MockShinySession$new(), f())
  invisible(e)
}

.mdc_notif_recorder <- function(e) {
  function(message, type = "message", duration = NULL, id = NULL, ...) {
    e$notifs[[length(e$notifs) + 1L]] <- list(message = message, type = type)
    paste0("id", length(e$notifs))
  }
}

# req() avec la sémantique RÉELLE de Shiny (avortement silencieux) —
# cf. test-mod-bulk-de-run.R:92-97.
.mdc_req <- function(...) {
  vals <- list(...)
  dead <- vapply(vals, function(v) is.null(v) || length(v) == 0L || isFALSE(v),
                 logical(1))
  if (any(dead)) {
    stop(structure(class = c("shiny.silent.error", "error", "condition"),
                   list(message = "", call = NULL)))
  }
  if (length(vals)) vals[[1]] else invisible(NULL)
}

.mdc_env <- function(pool_ok) {
  e <- new.env(parent = e_consts)
  e$`%||%` <- function(a, b) if (is.null(a)) b else a
  e$RCTD_CELL_MIN_INSTANCE <- .mdc_const("RCTD_CELL_MIN_INSTANCE")
  e$LABEL_TRANSFER_MIN_SHARED_GENES <- .mdc_const("LABEL_TRANSFER_MIN_SHARED_GENES")
  e$notifs <- list()
  e$invoked <- list()
  e$showNotification <- .mdc_notif_recorder(e)
  e$.tr <- function(key) key
  e$.t_fmt <- function(fmt, ...) fmt
  e$req <- .mdc_req
  e$reset_log <- function(file) {
    e$reset_called <- TRUE
    cat("", file = file)
    invisible(file)
  }
  # Le point d'injection du garde : pool_ok = FALSE simule « mirai absent ou
  # init/vérification en échec » — le tryCatch du bloc rend le chemin entier.
  e$ensure_app_daemons <- function() isTRUE(pool_ok)
  e$log_file <- tempfile(fileext = ".log")
  # reactiveVal doublé : setter ignore, getter rend NULL (chemin « ? » du
  # plafond de délai dans le bloc de statut).
  e$last_deconv_mode <- function(x = NULL) invisible(NULL)
  e$deconv_task <- list(
    invoke = function(...) {
      e$invoked[[length(e$invoked) + 1L]] <- list(...)
      invisible(NULL)
    },
    status = function() NULL,
    result = function() NULL)
  # `shared_rv` a une sémantique de RÉFÉRENCE dans l'app (reactiveValues) :
  # un environnement la reproduit, une liste serait COPIÉE (leçon bulk).
  e$shared_rv <- list2env(list(qc_pass_idx = 1L, deconv_params = NULL,
                               deconv_props = NULL), parent = emptyenv())
  e$global_data <- list(
    i18n = NULL,
    session_checkpoint_publish = function(...) invisible(NULL),
    spatial_obj = list(bpcells_dir = "bp_dir",
                       coords = data.frame(id = "c1", x = 1, y = 2)))
  e$ref <- list(reference_is_ready = function() TRUE,
                using_shared_ref = function() FALSE,
                effective_ref_path = function() "ref.rds",
                ref_state = list(status = "loaded"))
  # mode "stdeconvolve" : le chemin court SANS référence — le garde doit
  # tirer AVANT reset_log/invoke quel que soit le mode.
  e$input <- list(btn_deconv = 1L, mode = "stdeconvolve", n_topics = 5L,
                  n_top_od = 100L, lt_npcs = 30L, lt_norm_method = "lognorm",
                  lt_ncells = 3000L, rctd_n_hvg = 2000L)
  e
}

# ---------------------------------------------------------------------------
# (1) POOL MORT => invoke NON appelé + notification honnête, journal intact
# ---------------------------------------------------------------------------
test_that("P-lane deconv : pool mirai mort => invoke NON appelé + notification error", {
  e <- .mdc_env(pool_ok = FALSE)
  .mdc_run(e, "btn_deconv")

  expect_length(e$invoked, 0L)
  expect_length(e$notifs, 1L)
  expect_identical(e$notifs[[1]]$type, "error",
                   info = "un pool mort est une ERREUR, pas un refus de référence")
  expect_true(grepl("Pool mirai indisponible", e$notifs[[1]]$message, fixed = TRUE),
              info = "le message honnête nomme la panne")
  expect_false(isTRUE(e$reset_called),
               info = "reset_log() AVANT la garde recrée le journal vide")
  expect_null(e$shared_rv$deconv_params,
              info = "aucun paramètre publié pour un run qui ne part pas")
})

# ---------------------------------------------------------------------------
# (1bis) POOL OK => la garde n'entrave pas le dispatch
# ---------------------------------------------------------------------------
test_that("P-lane deconv : pool OK => invoke appelé avec le log_file du module", {
  skip_if_not_installed("mirai")
  on.exit(tryCatch(mirai::daemons(0), error = function(e) NULL), add = TRUE)
  if (!mirai::daemons_set()) suppressMessages(mirai::daemons(1))
  e <- .mdc_env(pool_ok = TRUE)
  .mdc_run(e, "btn_deconv")

  expect_length(e$invoked, 1L)
  expect_identical(e$invoked[[1]]$mode, "stdeconvolve")
  expect_identical(e$invoked[[1]]$log_file, e$log_file)
})

# ---------------------------------------------------------------------------
# (2) STATUT "error" => le journal reçoit le VRAI texte d'erreur de la tâche
# ---------------------------------------------------------------------------
test_that("P-lane deconv : statut error => journal horodaté avec l'erreur RÉELLE", {
  log_file <- tempfile(fileext = ".log")
  cat("", file = log_file)  # reset_log a déjà couru : le fichier existe, vide
  on.exit(unlink(log_file), add = TRUE)
  e <- .mdc_env(pool_ok = TRUE)
  e$log_file <- log_file
  e$deconv_task <- list(
    status = function() "error",
    result = function() stop("mirai: timed out after 1200000 ms"))
  .mdc_run(e, "deconv_task$status")

  lines <- readLines(log_file, warn = FALSE)
  expect_true(any(grepl("ECHEC :", lines, fixed = TRUE)),
              info = paste("journal lu :", paste(lines, collapse = " || ")))
  expect_true(any(grepl("mirai: timed out after 1200000 ms", lines, fixed = TRUE)),
              info = "le TEXTE RÉEL de la tâche, pas un message générique")
  expect_true(any(grepl("^[0-9]{2}:[0-9]{2}:[0-9]{2} \\| ECHEC :", lines)),
              info = "ligne horodatée par write_mirai_log (même format que le daemon)")
  expect_true(length(e$notifs) >= 1L && e$notifs[[1]]$type == "error",
              info = "la notification « voir le log » reste en place")
  expect_null(e$shared_rv$deconv_props,
              info = "un échec ne publie aucun résultat de déconvolution")
})
