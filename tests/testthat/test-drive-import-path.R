# =============================================================================
# test-drive-import-path.R — §H.6 (2026-11-14) : MCP dataset-ingest
# =============================================================================
# Demande opérateur (« MCP should be able to add data/dataset », P0.4 note fix)
# : un agent ne peut pas remplir un fileInput — le navigateur agent (IAB)
# refuse le file-chooser, et `import_file` ne couvre que les importeurs
# 10x/Spatial. Continuation P2.3 (FLAW V-FI, `b6f703e`) en trois couches :
#
#   (a) `ts_resolve_data_path()` — moteur PUR de résolution/confinement
#       (repli env, allowlist d'extensions PAR CIBLE, erreur classée) ;
#   (b) action drive `import_path` — la demande route vers les importeurs
#       PATH-driven (TS_DRIVE_PATH_MODULES) par la MÊME seam que `import_file`
#       (ts_drive_publish_importer), avec la règle UN ROUTE PAR MODULE ;
#   (c) `ts_velocity_resolve_loom_path()` DÉMÉNAGÉ de sc_velocity.R vers
#       drive_allowlist.R — wrapper du moteur (a). Le déménagement réconcilie
#       le pin « public API surface of sc_velocity.R is frozen »
#       (test-velocity-contract-freeze.R) : le nom était public au top-level
#       du moteur sans être déclaré dans `velocity_public_api()`.
#   (d) mod_sc_velocity.R publie l'importeur `import_velocity` : ingestion par
#       `.velocity_ingest()` — le pipeline EXACT du bouton humain (zéro voie
#       parallèle, motif G3 import_bulk) — après re-résolution CONFINÉE du
#       chemin par le moteur (a).
# =============================================================================

Sys.setenv(OMNIPATHR_LOGFILE = "none")

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
source_project_file("tools/pkg_skeleton/R/defaults.R")
source_project_file("tools/pkg_skeleton/R/drive_allowlist.R")
source_project_file("tools/pkg_skeleton/R/drive_watcher.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")

.read_src <- function(rel) {
  paste(readLines(file.path(ts_project_root(), rel),
                  warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

# ── (a) le moteur pur : résolution + confinement + allowlist d'extensions ───
test_that("ts_resolve_data_path resout et confine un chemin de donnees", {
  root <- withr::local_tempdir()
  loom_f <- file.path(root, "sample.loom")
  h5_f   <- file.path(root, "sample.h5ad")
  rds_f  <- file.path(root, "sample.rds")
  txt_f  <- file.path(root, "bad.exe")
  writeLines("dummy", loom_f)
  writeLines("dummy", h5_f)
  writeLines("dummy", rds_f)
  writeLines("dummy", txt_f)

  outside_dir <- withr::local_tempdir()
  outside_f <- file.path(outside_dir, "outside.loom")
  writeLines("dummy", outside_f)

  expect_identical(
    normalizePath(ts_resolve_data_path(loom_f, roots = root), winslash = "/"),
    normalizePath(loom_f, winslash = "/")
  )
  # repli env : path NULL + variable d'environnement
  withr::with_envvar(c(TS_H6_DATA_PATH = h5_f), {
    expect_identical(
      normalizePath(ts_resolve_data_path(env_var = "TS_H6_DATA_PATH",
                                         roots = root), winslash = "/"),
      normalizePath(h5_f, winslash = "/")
    )
  })

  # allowlist d'extensions PAR APPEL : un .rds refuse quand la cible ne veut
  # que loom/h5ad — le message est actionnable (liste les attendues)
  expect_error(
    ts_resolve_data_path(rds_f, roots = root, allowed_exts = c("loom", "h5ad")),
    class = "data_path_validation_error"
  )
  # allowlist PAR DEFAUT : loom/h5ad/rds (voie velocity)
  expect_identical(
    normalizePath(ts_resolve_data_path(rds_f, roots = root), winslash = "/"),
    normalizePath(rds_f, winslash = "/")
  )

  expect_error(
    ts_resolve_data_path(txt_f, roots = root),
    class = "data_path_validation_error"
  )
  expect_error(
    ts_resolve_data_path(file.path(root, "..", basename(outside_f)), roots = root),
    class = "data_path_validation_error"
  )
  expect_error(
    ts_resolve_data_path(outside_f, roots = root),
    class = "data_path_validation_error"
  )
  # ni path ni env : refusal classé (jamais un retour silencieux)
  withr::with_envvar(c(TS_H6_DATA_PATH = ""), {
    expect_error(
      ts_resolve_data_path(env_var = "TS_H6_DATA_PATH", roots = root),
      class = "data_path_validation_error"
    )
  })
})

# ── (c) le wrapper velocity : même contrat V-FI, déménagé ───────────────────
test_that("ts_velocity_resolve_loom_path delegue au moteur (meme classe, meme env)", {
  root <- withr::local_tempdir()
  h5_f <- file.path(root, "sample.h5ad")
  writeLines("dummy", h5_f)
  withr::with_envvar(c(TS_VELOCITY_LOOM_PATH = h5_f), {
    expect_identical(
      normalizePath(ts_velocity_resolve_loom_path(roots = root), winslash = "/"),
      normalizePath(h5_f, winslash = "/")
    )
  })
  txt_f <- file.path(root, "bad.exe")
  writeLines("dummy", txt_f)
  expect_error(
    ts_velocity_resolve_loom_path(txt_f, roots = root),
    class = "velocity_validation_error"
  )
  # le nom VIT dans la couche drive, pas dans le moteur velocity — la
  # réconciliation du pin freeze est un DÉMÉNAGEMENT, jamais une duplication
  expect_false(grepl("ts_velocity_resolve_loom_path <- function",
                     .read_src("tools/pkg_skeleton/R/sc_velocity.R"),
                     fixed = TRUE),
               info = "l'ancienne definition doit avoir QUITTE sc_velocity.R")
  expect_true(grepl("ts_velocity_resolve_loom_path <- function",
                    .read_src("tools/pkg_skeleton/R/drive_allowlist.R"),
                    fixed = TRUE),
              info = "le wrapper delegant vit dans drive_allowlist.R")
})

# ── (b) le protocole : action + registres + un route par module ─────────────
test_that("import_path est dans l'enum et epingle au token", {
  expect_true("import_path" %in% TS_DRIVE_ACTIONS)
  expect_true("import_path" %in% TS_DRIVE_TOKEN_PINNED_ACTIONS,
              info = "import_path REMPLACE l'etat d'analyse velocity — meme classe S3 que import_file")
})

test_that("le schema declare l'importeur path-driven velocity", {
  expect_true("import_velocity" %in% names(TS_DRIVE_IMPORT_SCHEMA))
  sch <- TS_DRIVE_IMPORT_SCHEMA$import_velocity
  expect_identical(sch$required, "path")
  expect_identical(sch$optional, "name")
  expect_identical(sch$exts, c("loom", "h5ad", "rds"))
  # `check` MANDATORY : sans lui le validateur partage retombe dans la branche
  # Bulk (counts_path) — mesure avant fix sur un scenario bien forme
  expect_identical(sch$check, "ts_drive_validate_velocity_path")
  expect_true(is.function(get0("ts_drive_validate_velocity_path", envir = globalenv())),
              info = "un check non appelable est un REFUS, jamais un skip silencieux (regle validate_import)")
  # TS_DRIVE_PATH_MODULES est le set DATA des importeurs path-driven, et
  # chaque membre est un importeur declare (le schema reste la source unique)
  expect_true(all(TS_DRIVE_PATH_MODULES %in% names(TS_DRIVE_IMPORT_SCHEMA)))
  expect_true(all(TS_DRIVE_PATH_MODULES %in% TS_DRIVE_MODULES),
              info = "un importeur hors TS_DRIVE_MODULES serait inatteignable (regle allowlist_problems)")
})

test_that("ts_drive_validate_velocity_path confine et trime (seam)", {
  root <- withr::local_tempdir()
  loom_f <- file.path(root, "wb_sample.loom")
  writeLines("dummy", loom_f)
  r <- ts_drive_validate_velocity_path(list(path = loom_f, name = " WB "),
                                       roots = root)
  expect_true(r$ok)
  expect_identical(normalizePath(r$import$path, winslash = "/"),
                   normalizePath(loom_f, winslash = "/"))
  expect_identical(r$import$name, "WB")
  # ext hors allowlist refusee au seam (actionable : liste les attendues)
  bad <- file.path(root, "bad.exe")
  writeLines("dummy", bad)
  r <- ts_drive_validate_velocity_path(list(path = bad), roots = root)
  expect_false(r$ok)
  expect_match(r$errors[1], "loom", fixed = TRUE)
  # `name` avec separateur : refuse (label d'affichage, jamais un chemin)
  r <- ts_drive_validate_velocity_path(list(path = loom_f, name = "a/b"),
                                       roots = root)
  expect_false(r$ok)
})

test_that("UN ROUTE PAR MODULE : import_file refuse import_velocity, import_path refuse les 10x", {
  # import_file -> import_velocity : refuse (l'importeur lit `path`, pas les
  # cles 10x ; le message nomme l'action correcte)
  scn_ff <- list(protocol = TS_DRIVE_PROTOCOL, seq = 1, module = "import_velocity",
                 action = "import_file", session_token = "tok",
                 import = list(path = "/tmp/x.loom"))
  v <- ts_drive_validate_scenario(scn_ff, "tok", last_seq = 0)
  expect_false(v$ok)
  expect_identical(v$status, "invalid")
  expect_match(v$errors[1], "import_path", fixed = TRUE)

  # import_path -> import_bulk : refuse (le 10x reste import_file ; le
  # message liste les importeurs path-driven)
  scn_pb <- list(protocol = TS_DRIVE_PROTOCOL, seq = 2, module = "import_bulk",
                 action = "import_path", session_token = "tok",
                 import = list(path = "/tmp/x.csv"))
  v <- ts_drive_validate_scenario(scn_pb, "tok", last_seq = 0)
  expect_false(v$ok)
  expect_identical(v$status, "invalid")
  expect_match(v$errors[1], "import_velocity", fixed = TRUE)

  # import_path -> import_velocity : la charge utile route par le schema
  # (validation generee : `path` requis, `name` optionnel). tempdir() est
  # dans les racines allowlistées (ts_drive_import_roots) — le fichier y est
  # créé pour que la validation d'EXISTENCE passe comme en live.
  ok_f <- file.path(tempdir(), "h6_sample.loom")
  writeLines("dummy", ok_f)
  on.exit(unlink(ok_f), add = TRUE)
  scn_ok <- list(protocol = TS_DRIVE_PROTOCOL, seq = 3, module = "import_velocity",
                 action = "import_path", session_token = "tok",
                 import = list(path = ok_f, name = "WB"))
  v <- ts_drive_validate_scenario(scn_ok, "tok", last_seq = 0)
  expect_true(v$ok)
  expect_identical(v$status, "applied-candidate")
  expect_identical(normalizePath(v$scenario$import$path, winslash = "/"),
                   normalizePath(ok_f, winslash = "/"))
  expect_identical(v$scenario$import$name, "WB")

  # cle inconnue refusee (whitelist schema, motif import_sc sample_name)
  scn_bad <- list(protocol = TS_DRIVE_PROTOCOL, seq = 4, module = "import_velocity",
                  action = "import_path", session_token = "tok",
                  import = list(path = ok_f, dir_path = "data_sc/"))
  v <- ts_drive_validate_scenario(scn_bad, "tok", last_seq = 0)
  expect_false(v$ok)
})

# ── (d) le module velocity publie l'importeur (contrat source) ──────────────
test_that("mod_sc_velocity publie l'importeur import_velocity sur la seam G3", {
  src <- .read_src("modules/sc/mod_sc_velocity.R")
  # garde défensive mesurée : la couche drive peut être absente des
  # environnements de test minimalistes — l'importeur n'est publié que si
  # ts_drive_publish_importer existe (motif ExtendedTask)
  expect_match(src, "if (.vel_drive_ok) ts_drive_publish_importer(global_data,",
               fixed = TRUE,
               info = "l'importeur doit etre publie sur la seam standard (garde défensive)")
  expect_match(src, '.vel_drive_ok <- exists("ts_drive_publish_importer"',
               fixed = TRUE,
               info = "la publication doit être gardée sur la présence de la couche drive")
  # le chemin arrive VALIDE du seam (motif import_sc) : l'importeur consomme,
  # il ne re-dérive pas — la validation confinement/exts vit dans
  # ts_drive_validate_velocity_path (testée ci-dessus)
  expect_match(src, "path <- request$path", fixed = TRUE,
               info = "l'importeur consomme le chemin validé, sans re-dérivation")
  # zéro voie parallèle : l'ingestion passe par le pipeline partagé
  n <- lengths(regmatches(src, gregexpr(".velocity_ingest(", src, fixed = TRUE)))[1]
  expect_true(n >= 3L,
              info = sprintf("definition + appel bouton + appel importeur attendus (trouve : %d)", n))
})

test_that("le plan de navigation import_velocity ouvre les accordions velocity", {
  nav <- ts_drive_nav_plan(TS_DRIVE_SC_VELOCITY_IMPORT_MODULE)
  expect_identical(nav$top, TS_DRIVE_SC_TOP_TAB)
  expect_identical(nav$panel, "8_velocity")
  expect_identical(nav$accordion_id, c("sc-acc_workflow", "sc-acc_dynamique"))
  expect_null(nav$tab)
})
