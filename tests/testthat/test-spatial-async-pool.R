# =============================================================================
# test-spatial-async-pool.R — PC (2026-11-09) : bornes de vie du pool mirai +
# watchdog « pool épuisé » (tools/pkg_skeleton/R/spatial_async.R, changelog v7)
# =============================================================================
# Contrats testés :
#  - ts_mirai_pool_lifetime_args() : OPT-IN strict — sans réglage TS_MIRAI_*,
#    liste vide => l'appel de daemons() reste mirai::daemons(n) (workers
#    permanents, zéro changement).
#  - .ts_mirai_pool_exhausted() : drapeau ready levé mais 0 worker vivant —
#    MESURE 2026-11-09 : après mort par idletime, status() reste synchrone et
#    renvoie connections = 0, alors qu'un mirai() émis dans cet état attend
#    indéfiniment (le dispatcher ne relance PAS les daemons locaux).
#  - init_spatial_daemons() : le pool épuisé est respawn (preload +
#    vérification complètes) au lieu de l'early-return aveugle.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/spatial_async.R")

skip_if_not_installed("mirai")

# --- ts_mirai_pool_lifetime_args ---------------------------------------------

test_that("PC : bornes de vie ABSENTES par defaut (workers permanents, zero changement)", {
  withr::local_envvar(c(TS_MIRAI_IDLETIME_MS = NA, TS_MIRAI_WALLTIME_MS = NA))
  withr::local_options(list(ts.mirai.idletime = NULL, ts.mirai.walltime = NULL))
  expect_identical(ts_mirai_pool_lifetime_args(), list())
})

test_that("PC : bornes de vie lues env d'abord, option ensuite, invalides ignorees", {
  withr::local_options(list(ts.mirai.idletime = NULL, ts.mirai.walltime = NULL))
  withr::with_envvar(c(TS_MIRAI_IDLETIME_MS = "1500", TS_MIRAI_WALLTIME_MS = "900000"),
                     expect_identical(ts_mirai_pool_lifetime_args(),
                                      list(idletime = 1500, walltime = 900000)))
  # option en repli quand l'env est absent
  withr::with_envvar(c(TS_MIRAI_IDLETIME_MS = NA), {
    withr::with_options(list(ts.mirai.idletime = 2000), {
      expect_identical(ts_mirai_pool_lifetime_args(), list(idletime = 2000))
    })
  })
  # invalides (0, negatif, non numerique, NA) => ignores, pas d'argument pose
  withr::with_envvar(c(TS_MIRAI_IDLETIME_MS = "0", TS_MIRAI_WALLTIME_MS = "-5"),
                     expect_identical(ts_mirai_pool_lifetime_args(), list()))
  withr::with_envvar(c(TS_MIRAI_IDLETIME_MS = "abc", TS_MIRAI_WALLTIME_MS = NA),
                     expect_identical(ts_mirai_pool_lifetime_args(), list()))
})

# --- .ts_mirai_pool_exhausted ------------------------------------------------

test_that("PC : pool inexistant ou vivant => PAS epuise (faux positif interdit)", {
  stop_spatial_daemons()
  tryCatch(if (mirai::daemons_set()) mirai::daemons(0), error = function(e) NULL)
  expect_false(.ts_mirai_pool_exhausted())
  mirai::daemons(2)
  on.exit(stop_spatial_daemons(), add = TRUE)
  expect_false(.ts_mirai_pool_exhausted())
})

test_that("PC : pool mort par idletime => EPUISE (status() synchrone, connections 0)", {
  mirai::daemons(2, idletime = 500)
  on.exit(stop_spatial_daemons(), add = TRUE)
  Sys.sleep(2)   # > idletime : les deux workers meurent
  expect_true(mirai::daemons_set())
  expect_identical(as.integer(mirai::status()$connections), 0L)
  expect_true(.ts_mirai_pool_exhausted())
})

# --- watchdog d'init ---------------------------------------------------------

test_that("PC : init_spatial_daemons respawn un pool epuise (watchdog, pas early-return)", {
  # le cableage PC-1 pose le formal `future` de DT::renderDT : restaure en
  # sortie pour ne pas fuiter dans la suite du process de test.
  on.exit(try(ts_dt_future_restore(), silent = TRUE), add = TRUE)

  mirai::daemons(2, idletime = 500)
  Sys.sleep(2)   # mort par idle
  # Etat "pret mais epuise" : ce que l'app voit apres une mort par idletime.
  .spatial_async_env$daemons_ready <- TRUE
  .spatial_async_env$n_daemons <- 2L
  on.exit(stop_spatial_daemons(), add = TRUE)

  # preload vide : la verification post-init echouera (warning attendu, les
  # fonctions projet ne sont pas preloadees) — le watchdog teste ICI est le
  # RESPAWN, pas la verification (couverte par les tests spatiaux existants).
  # PAS de tryCatch(warning=) : il avorterait init SUR l'avertissement du
  # watchdog lui-meme — testthat muffie les warnings, on appelle direct.
  init_spatial_daemons(n_daemons = 2L, source_files = character(0))
  expect_true(isTRUE(.spatial_async_env$daemons_ready))
  expect_true(mirai::daemons_set())
  expect_gt(as.integer(mirai::status()$connections), 0L)

  # un deuxieme init (pool vivant) reste un no-op idempotent :
  init_spatial_daemons(n_daemons = 2L, source_files = character(0))
  expect_gt(as.integer(mirai::status()$connections), 0L)
})
