# =============================================================================
# test-packs-manifest.R — PK-1 : pilier 1 « Système de packs »
# =============================================================================
# Fiche (fait foi) :
#   * docs/proposals/PROPOSAL_PACKS_ECOSYSTEM_2026-11-09.md §3 (jalon PK-1 :
#     lecteur manifest.json PUR, validation stricte, refus honnêtes, scan au
#     démarrage, registre session), §4/D2 (emplacement data_sc/packs/),
#     §4/D6 (zéro nouvelle dép — jsonlite déjà au lock), §5 (non-goals :
#     PAS d'outil MCP, PAS d'API cloud, PAS de palettes ni d'enricher ici),
#     §6 (RED→GREEN : fixture pack valide/invalide, refus honnêtes) ;
#   * docs/ideation/PACKS_ECOSYSTEM_ENRICHR_PALETTES.md §3 (schéma manifest
#     « libraries » : type/id/name/version/source/downloaded/attribution/files),
#     §7.1 (schéma manifest « palette » : palettes hex embarquées + attribution
#     + license), §8 (data_sc/packs/<pack>/manifest.json), §9 (invalid manifest
#     → pack SKIPPÉ avec message visible, JAMAIS un crash ; attribution/license
#     obligatoires pour les palettes).
#
# Verrouillé ici :
#   * surface publique gelée + vocabulaire d'états + types de pack ;
#   * pack_read_manifest() : lectures valides (libraries + palette) et REFUS
#     HONNÊTES classés (manifest absent, JSON malformé, champ requis absent,
#     champ malformé, type inconnu, fichier GMT déclaré absent) — chaque refus
#     porte un `state` gelé ET un message qui dit POURQUOI (fiche §6) ;
#   * pack_scan_directory() : scan passif STRICTEMENT LECTURE-SEULE d'UNE
#     profondeur (motif ts_dataset_dir_scan), racine absente => résultat vide
#     jamais une erreur, pack invalide = refus ISOLÉ (jamais un crash, §9),
#     doublon d'id = refus (le registre ne peux pas mentir) ;
#   * registre session PUR (motif sc_multi.R) : register/remove/get/summary
#     retournent la NOUVELLE liste nommée invisiblement, ne mutent JAMAIS un
#     global — l'appelant assigne global_data$packs lui-même ;
#   * pureté Shiny du moteur (garde du freeze test sc_multi) ;
#   * ancres app.R : source() + champs registre + appel du scan.
#
# FIXTURES : construites dans tempdir() — data_sc/ est la zone DONNÉES
# UTILISATEUR (gitignorée, .gitignore:121) : AUCUN fixture n'y est écrit ni
# commité (discipline PK-1). Le scan par défaut (data_sc/packs) n'est pas
# touché : tous les appels passent une racine explicite tempdir().
# =============================================================================

source_project_file("tools/pkg_skeleton/R/session_checkpoint.R")  # .sc_project_root/.sc_scalar_text (réutilisés, zéro duplication)
source_project_file("tools/pkg_skeleton/R/packs.R")

# ── Fixtures tempdir() ───────────────────────────────────────────────────────

.pfx_root <- file.path(tempdir(), "pk1_packs_fixtures")
if (dir.exists(.pfx_root)) unlink(.pfx_root, recursive = TRUE)
dir.create(.pfx_root, recursive = TRUE)

.pfx_write_manifest <- function(pack_dir, lst, pretty = TRUE) {
  dir.create(pack_dir, recursive = TRUE, showWarnings = FALSE)
  txt <- jsonlite::toJSON(lst, auto_unbox = TRUE, pretty = pretty, null = "null")
  writeLines(txt, file.path(pack_dir, "manifest.json"), useBytes = TRUE)
  file.path(pack_dir, "manifest.json")
}

.pfx_write_gmt <- function(pack_dir, fname) {
  dir.create(pack_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("KEGG_2021_Human\t\tTP53\tBRCA1", "ChEA_2022\t\tMYC\tMAX"),
             file.path(pack_dir, fname), useBytes = TRUE)
  file.path(pack_dir, fname)
}

# Pack « libraries » VALIDE — schéma exact de la fiche idéation §3.
.pfx_valid_libraries <- function(pack_dir) {
  .pfx_write_gmt(pack_dir, "KEGG_2021_Human.gmt")
  .pfx_write_gmt(pack_dir, "ChEA_2022.gmt")
  .pfx_write_manifest(pack_dir, list(
    type = "libraries",
    id = "enrichr_core",
    name = "Enrichr core libraries",
    version = "1.0.0",
    source = "https://maayanlab.cloud/Enrichr/#libraries",
    downloaded = "2026-10-31",
    attribution = "Enrichr (Ma'ayan Lab, Icahn School of Medicine at Mount Sinai)",
    files = c("KEGG_2021_Human.gmt", "ChEA_2022.gmt")
  ))
}

# Pack « palette » VALIDE — schéma exact de la fiche idéation §7.1.
.pfx_valid_palette <- function(pack_dir) {
  .pfx_write_manifest(pack_dir, list(
    type = "palette",
    id = "wes_anderson",
    name = "Wes Anderson",
    version = "1.0.0",
    attribution = "wesanderson R package (Karthik Ram), GPL-3",
    license = "hexes redistribues avec attribution — verifier la licence de chaque source",
    palettes = list(
      Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"),
      GrandBudapest1 = c("#F1BB7B", "#FD6467", "#5B7218", "#E9E3D4")
    )
  ))
}

# Capture d'un refus classé (motif tryCatch du freeze test sc_multi).
.pfx_err <- function(code) {
  tryCatch(force(code), error = function(e) e)
}
.pfx_is_packs_error <- function(e, state) {
  inherits(e, "packs_error") && identical(e$state, state)
}

# ── Surface publique gelée ───────────────────────────────────────────────────

test_that("public API surface is frozen (names + signatures)", {
  expect_setequal(
    packs_public_api(),
    c("packs_public_api", "packs_error_states", "packs_types",
      "packs_dir", "pack_read_manifest", "pack_scan_directory",
      "packs_register", "packs_remove", "packs_get", "packs_summary")
  )
  for (fn_name in packs_public_api()) {
    expect_true(exists(fn_name, where = globalenv(), inherits = FALSE),
                info = paste("fonction absente :", fn_name))
  }
  expected_args <- list(
    packs_public_api     = list(),
    packs_error_states   = list(),
    packs_types          = list(),
    packs_dir            = list(dir = NULL),
    pack_read_manifest   = list(manifest_path = NULL),
    pack_scan_directory  = list(packs_root = NULL),
    packs_register       = list(registry = NULL, pack = NULL, overwrite = FALSE),
    packs_remove         = list(registry = NULL, id = NULL),
    packs_get            = list(registry = NULL, id = NULL),
    packs_summary        = list(registry = NULL)
  )
  for (fn_name in names(expected_args)) {
    fmls <- names(formals(get(fn_name, envir = globalenv())))
    expect_identical(fmls, names(expected_args[[fn_name]]),
                     info = paste("signature de", fn_name))
  }
})

test_that("error states and pack types are frozen", {
  expect_setequal(
    packs_error_states(),
    c("invalid_input", "missing_manifest", "malformed_json", "missing_field",
      "invalid_field", "unknown_type", "missing_file", "duplicate_id",
      "unknown_id", "scan_failed")
  )
  # Fiche idéation §3 : type = "libraries" | "palette" | <futur>. Le lecteur
  # ne connaît que les types IMPLÉMENTÉS — un type futur = refus honnête
  # (unknown_type), pas un chargement silencieux.
  expect_setequal(packs_types(), c("libraries", "palette"))
  # Toute erreur émise appartient à l'ensemble gelé (invariant).
  bad <- .pfx_err(pack_read_manifest(file.path(tempdir(), "n_existe_pas.json")))
  expect_true(bad$state %in% packs_error_states())
})

# ── pack_read_manifest : lectures valides ────────────────────────────────────

test_that("pack_read_manifest reads a valid 'libraries' pack (fiche ideation §3)", {
  d <- file.path(.pfx_root, "ok_libraries")
  man <- .pfx_valid_libraries(d)
  p <- pack_read_manifest(man)
  expect_s3_class(p, "biomentat_pack")
  expect_identical(p$id, "enrichr_core")
  expect_identical(p$type, "libraries")
  expect_identical(p$name, "Enrichr core libraries")
  expect_identical(p$version, "1.0.0")
  expect_match(p$attribution, "Ma'ayan Lab", fixed = TRUE)
  expect_identical(p$source, "https://maayanlab.cloud/Enrichr/#libraries")
  expect_identical(p$downloaded, "2026-10-31")
  expect_identical(p$files, c("KEGG_2021_Human.gmt", "ChEA_2022.gmt"))
  # Chemins résolus : le pack sait où il vit (PK-2 consommera files+dir).
  expect_identical(p$manifest, normalizePath(man, winslash = "/", mustWork = TRUE))
  expect_identical(p$dir, normalizePath(d, winslash = "/", mustWork = TRUE))
})

test_that("pack_read_manifest reads a valid 'palette' pack (fiche ideation §7.1)", {
  d <- file.path(.pfx_root, "ok_palette")
  man <- .pfx_valid_palette(d)
  p <- pack_read_manifest(man)
  expect_s3_class(p, "biomentat_pack")
  expect_identical(p$id, "wes_anderson")
  expect_identical(p$type, "palette")
  expect_identical(p$palettes$Zissou1,
                   c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"))
  expect_match(p$license, "attribution", fixed = TRUE)
})

# ── pack_read_manifest : refus honnêtes (chaque refus dit POURQUOI) ─────────

test_that("refusal: missing manifest states the expected path", {
  d <- file.path(.pfx_root, "sans_manifest")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  e <- .pfx_err(pack_read_manifest(file.path(d, "manifest.json")))
  expect_true(.pfx_is_packs_error(e, "missing_manifest"),
              info = paste(class(e), conditionMessage(e)))
  expect_match(conditionMessage(e), "manifest.json", fixed = TRUE)
})

test_that("refusal: malformed JSON names the parse failure", {
  d <- file.path(.pfx_root, "json_casse")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines("{ type: pas du json ]", file.path(d, "manifest.json"),
             useBytes = TRUE)
  e <- .pfx_err(pack_read_manifest(file.path(d, "manifest.json")))
  expect_true(.pfx_is_packs_error(e, "malformed_json"),
              info = conditionMessage(e))
  # un manifeste qui n'est PAS un objet JSON est également malformé
  d2 <- file.path(.pfx_root, "json_tableau")
  dir.create(d2, recursive = TRUE, showWarnings = FALSE)
  writeLines('["type"]', file.path(d2, "manifest.json"), useBytes = TRUE)
  e2 <- .pfx_err(pack_read_manifest(file.path(d2, "manifest.json")))
  expect_true(.pfx_is_packs_error(e2, "malformed_json"),
              info = conditionMessage(e2))
})

test_that("refusal: unknown type lists the allowed types", {
  d <- file.path(.pfx_root, "type_inconnu")
  man <- .pfx_write_manifest(d, list(
    type = "rapport_template", id = "x", name = "X", version = "1.0.0",
    attribution = "A"))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "unknown_type"),
              info = conditionMessage(e))
  expect_match(conditionMessage(e), "rapport_template", fixed = TRUE)
  expect_match(conditionMessage(e), "libraries", fixed = TRUE)
  expect_match(conditionMessage(e), "palette", fixed = TRUE)
})

test_that("refusal: each missing required field is named", {
  base <- list(type = "libraries", id = "x", name = "X", version = "1.0.0",
               attribution = "A", files = "f.gmt")
  for (fld in c("type", "id", "name", "version", "attribution", "files")) {
    d <- file.path(.pfx_root, paste0("champ_manquant_", fld))
    lst <- base[setdiff(names(base), fld)]
    man <- .pfx_write_manifest(d, lst)
    e <- .pfx_err(pack_read_manifest(man))
    expect_true(.pfx_is_packs_error(e, "missing_field"),
                info = sprintf("champ %s : %s", fld, conditionMessage(e)))
    expect_match(conditionMessage(e), fld, fixed = TRUE,
                 info = sprintf("le refus doit nommer le champ %s", fld))
  }
  # files est requis pour libraries même quand présent mais vide
  d <- file.path(.pfx_root, "files_vide")
  man <- .pfx_write_manifest(d, list(
    type = "libraries", id = "x", name = "X", version = "1.0.0",
    attribution = "A", files = list()))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "invalid_field"),
              info = conditionMessage(e))
})

test_that("refusal: malformed fields state expected vs got", {
  # id non scalaire / vide / trop long / caractères de contrôle (motif sc_multi)
  d <- file.path(.pfx_root, "id_vide")
  man <- .pfx_write_manifest(d, list(
    type = "libraries", id = "   ", name = "X", version = "1.0.0",
    attribution = "A", files = "f.gmt"))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "invalid_field"), info = conditionMessage(e))

  # version non textuelle
  d <- file.path(.pfx_root, "version_num")
  man <- .pfx_write_manifest(d, list(
    type = "libraries", id = "x", name = "X", version = 1, attribution = "A",
    files = "f.gmt"))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "invalid_field"), info = conditionMessage(e))

  # downloaded présent mais pas une date ISO
  d <- file.path(.pfx_root, "date_branlee")
  man <- .pfx_write_manifest(d, list(
    type = "libraries", id = "x", name = "X", version = "1.0.0",
    attribution = "A", downloaded = "31/10/2026", files = "f.gmt"))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "invalid_field"),
              info = conditionMessage(e))
  expect_match(conditionMessage(e), "downloaded", fixed = TRUE)

  # files non textuel (objet au lieu de tableau)
  d <- file.path(.pfx_root, "files_objet")
  man <- .pfx_write_manifest(d, list(
    type = "libraries", id = "x", name = "X", version = "1.0.0",
    attribution = "A", files = list(a = "f.gmt")))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "invalid_field"), info = conditionMessage(e))

  # palette : hex invalide nommé
  d <- file.path(.pfx_root, "hex_branle")
  man <- .pfx_write_manifest(d, list(
    type = "palette", id = "p", name = "P", version = "1.0.0",
    attribution = "A", license = "L",
    palettes = list(Zissou1 = c("#3B9AB2", "bleu"))))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "invalid_field"), info = conditionMessage(e))
  expect_match(conditionMessage(e), "Zissou1", fixed = TRUE)
  expect_match(conditionMessage(e), "bleu", fixed = TRUE)

  # palette : license absente (fiche §9 — attribution/license obligatoires)
  d <- file.path(.pfx_root, "palette_sans_license")
  man <- .pfx_write_manifest(d, list(
    type = "palette", id = "p", name = "P", version = "1.0.0",
    attribution = "A",
    palettes = list(Zissou1 = c("#3B9AB2"))))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "missing_field"), info = conditionMessage(e))
  expect_match(conditionMessage(e), "license", fixed = TRUE)

  # palettes : entrées sans nom
  d <- file.path(.pfx_root, "palettes_sans_noms")
  man <- .pfx_write_manifest(d, list(
    type = "palette", id = "p", name = "P", version = "1.0.0",
    attribution = "A", license = "L",
    palettes = list(c("#3B9AB2"))))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "invalid_field"), info = conditionMessage(e))
})

test_that("refusal: missing referenced GMT file names the file and the dir", {
  d <- file.path(.pfx_root, "gmt_absent")
  .pfx_write_gmt(d, "KEGG_2021_Human.gmt")   # ChEA_2022.gmt volontairement absent
  man <- .pfx_write_manifest(d, list(
    type = "libraries", id = "x", name = "X", version = "1.0.0",
    attribution = "A", files = c("KEGG_2021_Human.gmt", "ChEA_2022.gmt")))
  e <- .pfx_err(pack_read_manifest(man))
  expect_true(.pfx_is_packs_error(e, "missing_file"), info = conditionMessage(e))
  expect_match(conditionMessage(e), "ChEA_2022.gmt", fixed = TRUE)
  expect_match(conditionMessage(e), basename(d), fixed = TRUE)
})

test_that("refusal: bad arguments are honest, not cryptic", {
  e <- .pfx_err(pack_read_manifest(character(0)))
  expect_true(.pfx_is_packs_error(e, "invalid_input"), info = conditionMessage(e))
  # un RÉPERTOIRE passé n'est pas un manifest.json
  d <- file.path(.pfx_root, "ok_libraries")
  e2 <- .pfx_err(pack_read_manifest(d))
  expect_true(.pfx_is_packs_error(e2, "invalid_input"), info = conditionMessage(e2))
  expect_match(conditionMessage(e2), "manifest.json", fixed = TRUE)
})

# ── pack_scan_directory : scan passif lecture-seule (§8/§9) ──────────────────

test_that("scan of a missing root returns an empty result and CREATES NOTHING", {
  absent <- file.path(.pfx_root, "racine_absente")
  res <- pack_scan_directory(absent)
  expect_identical(res$packs, list())
  expect_identical(res$refusals, list())
  expect_false(dir.exists(absent),
               info = "le scan ne doit JAMAIS créer la racine (lecture seule)")
})

test_that("scan collects valid packs AND isolated refusals, never crashes", {
  root <- file.path(.pfx_root, "scan_mixte")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  .pfx_valid_libraries(file.path(root, "a_enrichr_core"))
  .pfx_valid_palette(file.path(root, "b_wes_anderson"))
  dir.create(file.path(root, "c_sans_manifest"), recursive = TRUE)
  dir.create(file.path(root, "d_json_casse"), recursive = TRUE)
  writeLines("{ casse ]", file.path(root, "d_json_casse", "manifest.json"),
             useBytes = TRUE)
  .pfx_write_manifest(file.path(root, "e_type_inconnu"), list(
    type = "futur", id = "t", name = "T", version = "1.0.0", attribution = "A"))
  .pfx_write_gmt(file.path(root, "f_gmt_absent"), "present.gmt")
  .pfx_write_manifest(file.path(root, "f_gmt_absent"), list(
    type = "libraries", id = "g", name = "G", version = "1.0.0",
    attribution = "A", files = "absent.gmt"))
  # doublon d'id : un registre ne peut pas porter deux packs de même clé
  .pfx_valid_libraries(file.path(root, "g_doublon_id"))

  res <- pack_scan_directory(root)
  expect_setequal(names(res$packs), c("enrichr_core", "wes_anderson"))
  expect_s3_class(res$packs$enrichr_core, "biomentat_pack")
  # 5 refus isolés : sans manifest / JSON cassé / type inconnu / GMT absent / doublon
  expect_length(res$refusals, 5L)
  states <- vapply(res$refusals, function(r) r$state, character(1))
  expect_setequal(states, c("missing_manifest", "malformed_json",
                            "unknown_type", "missing_file", "duplicate_id"))
  for (r in res$refusals) {
    expect_true(nzchar(r$message) && !is.na(r$message),
                info = "chaque refus porte un message qui dit pourquoi")
    expect_true(nzchar(r$dir), info = "chaque refus nomme son dossier")
  }
  # le message du doublon nomme l'id déjà pris
  dup <- res$refusals[states == "duplicate_id"][[1L]]
  expect_match(dup$message, "enrichr_core", fixed = TRUE)
})

test_that("scan is STRICTLY read-only (AC2 motif dataset_directory)", {
  root <- file.path(.pfx_root, "scan_lecture_seule")
  .pfx_valid_libraries(file.path(root, "a_pack"))
  .pfx_write_manifest(file.path(root, "b_ko"), list(type = "futur", id = "k"))
  before <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  mtimes <- file.info(list.files(root, full.names = TRUE, recursive = TRUE))$mtime
  res <- pack_scan_directory(root)
  after <- list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  mtimes_after <- file.info(list.files(root, full.names = TRUE, recursive = TRUE))$mtime
  expect_identical(before, after)
  expect_identical(unname(mtimes), unname(mtimes_after))
  expect_length(res$refusals, 1L)   # b_ko : type inconnu, isolé
})

# ── Registre session PUR (motif sc_multi.R : jamais de mutation de global) ──

test_that("packs_register/get/remove follow the pure-registry pattern", {
  d1 <- file.path(.pfx_root, "reg_libraries")
  man1 <- .pfx_valid_libraries(d1)
  p1 <- pack_read_manifest(man1)
  d2 <- file.path(.pfx_root, "reg_palette")
  p2 <- pack_read_manifest(.pfx_valid_palette(d2))

  reg0 <- packs_register(NULL, p1)
  expect_identical(names(reg0), "enrichr_core")
  expect_s3_class(reg0$enrichr_core, "biomentat_pack")
  expect_true(nzchar(reg0$enrichr_core$registered_at))

  # registre pur : reg0 n'est PAS muté par l'ajout de p2 (l'appelant assigne)
  reg1 <- packs_register(reg0, p2)
  expect_length(reg0, 1L)
  expect_length(reg1, 2L)

  # doublon d'id refusé, sauf overwrite explicite
  e <- .pfx_err(packs_register(reg1, p1))
  expect_true(.pfx_is_packs_error(e, "duplicate_id"), info = conditionMessage(e))
  reg2 <- packs_register(reg1, p1, overwrite = TRUE)
  expect_length(reg2, 2L)

  expect_identical(packs_get(reg2, "wes_anderson")$id, "wes_anderson")
  e2 <- .pfx_err(packs_get(reg2, "inconnu"))
  expect_true(.pfx_is_packs_error(e2, "unknown_id"), info = conditionMessage(e2))

  reg3 <- packs_remove(reg2, "wes_anderson")
  expect_identical(names(reg3), "enrichr_core")
  e3 <- .pfx_err(packs_remove(reg3, "wes_anderson"))
  expect_true(.pfx_is_packs_error(e3, "unknown_id"), info = conditionMessage(e3))

  # mauvais arguments : refus classé, pas un crash cryptique
  e4 <- .pfx_err(packs_register(NULL, list(id = "pas_un_pack")))
  expect_true(.pfx_is_packs_error(e4, "invalid_input"), info = conditionMessage(e4))
  e5 <- .pfx_err(packs_register(c(a = 1), p1))
  expect_true(.pfx_is_packs_error(e5, "invalid_input"), info = conditionMessage(e5))
})

test_that("packs_summary has frozen columns and a 0-row empty frame", {
  cols <- c("id", "type", "name", "version", "n_files", "n_palettes",
            "attribution", "registered_at")
  empty <- packs_summary(NULL)
  expect_s3_class(empty, "data.frame")
  expect_identical(names(empty), cols)
  expect_identical(nrow(empty), 0L)

  d <- file.path(.pfx_root, "sum_libraries")
  p1 <- pack_read_manifest(.pfx_valid_libraries(d))
  reg <- packs_register(NULL, p1)
  s <- packs_summary(reg)
  expect_identical(names(s), cols)
  expect_identical(nrow(s), 1L)
  expect_identical(s$id[1], "enrichr_core")
  expect_identical(s$type[1], "libraries")
  expect_identical(s$n_files[1], 2L)
  expect_identical(s$n_palettes[1], 0L)
})

# ── Pureté Shiny (garde du freeze test sc_multi) ─────────────────────────────

test_that("tools/pkg_skeleton/R/packs.R is pure (no Shiny symbols)", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "tools/pkg_skeleton/R/packs.R"),
                         warn = FALSE), collapse = "\n")
  for (sym in c("reactiveVal", "reactiveValues", "observeEvent", "output\\$",
                "input\\$", "moduleServer", "showNotification",
                "isolate\\(", "reactive\\(")) {
    expect_false(grepl(sym, src, perl = TRUE),
                 info = paste("symbole Shiny interdit trouvé :", sym))
  }
  # lecture seule : aucune écriture disque dans le moteur (motif AC2)
  for (sym in c("writeLines", "write.table", "file.copy", "dir.create",
                "unlink")) {
    expect_false(grepl(sym, src, fixed = TRUE),
                 info = paste("écriture disque interdite dans le moteur :", sym))
  }
})

# ── Ancres app.R (câblage boot : source + registre session + scan) ───────────

test_that("app.R anchors: source, registry fields, scan call", {
  app <- paste(readLines(file.path(ts_project_root(), "app.R"), warn = FALSE),
               collapse = "\n")
  expect_match(app, 'source("tools/pkg_skeleton/R/packs.R")', fixed = TRUE)
  expect_match(app, "packs = list(),", fixed = TRUE)
  expect_match(app, "packs_refusals = list(),", fixed = TRUE)
  expect_match(app, "global_data$packs <- .packs_scan$packs", fixed = TRUE)
  expect_match(app, "global_data$packs_refusals <- .packs_scan$refusals",
               fixed = TRUE)
  # EXACTEMENT un point d'appel du scan (le boot de session) — aucun autre
  # site ne doit rescanner en douce (motif comptage du freeze test sc_multi).
  calls <- length(gregexpr("pack_scan_directory(", app, fixed = TRUE)[[1L]])
  expect_identical(calls, 1L)
})
