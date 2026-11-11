# =============================================================================
# test-cqf4-ui-cosmetics.R — CQF-4 : passe cosmétique spatial/import
# =============================================================================
# Constats mesurés (triage UIUX Spatial 2026-11-08, idée I4 — tous [M] au
# source ; GO opérateur 2026-11-11 « go for the UI/UX task ») :
#   RC6 — boutons PRIMAIRES en btn-danger (sémantique Bootstrap inverse :
#         « Lancer le pipeline complet », « Appliquer les seuils »)
#   RC4 — textes explicatifs de sidebar sous le plancher de lisibilité
#         0.8rem (0.7rem / 0.68rem, text-muted)
#   RC4 — console d'import spatial sans white-space:pre-wrap (reco
#         UIUX_REVIEW_2026-10-29 §3.4 partiellement appliquée — pre-wrap
#         présent dans la console pipeline, absent de la console import)
#   RC4 — littéraux UI ASCII du domaine QC spatial (§3.4 : « Apercu du jeu
#         de donnees », « Genes spatialement variables » — accents absents)
#   RC4 — boutons fileInput en anglais « Browse... » : la clé
#         « Parcourir... » existe déjà (mod_bulk_pathways la consomme) ;
#         périmètre CQF-4 = sites HORS modules/sc/ — la session concurrente
#         CCC-DRAW-OPT possède modules/sc/ (protocole essaim, jamais -A)
# Périmètre volontairement limité : mod_spatial_pipeline.R, mod_spatial_qc.R,
# mod_import_spatial.R + fileInput hors modules/sc/. « Deconvolution »
# (9 fichiers + i18n + gelés) = CQF-4bis, hors lot.
# =============================================================================
.ts_cos_read <- function(rel) {
  paste(readLines(file.path(ts_project_root(), rel), warn = FALSE,
                  encoding = "UTF-8"), collapse = "\n")
}

# --- 1. Sémantique des boutons primaires (RC6) --------------------------------

test_that("CQF-4 : les 2 CTA primaires spatial ne sont plus btn-danger", {
  # « Lancer le pipeline complet » = action principale, jamais destructive
  pipe <- .ts_cos_read("modules/spatial/mod_spatial_pipeline.R")
  expect_false(grepl("btn-danger", pipe, fixed = TRUE),
               info = "mod_spatial_pipeline.R : btn-danger sur CTA primaire (RC6)")
  expect_match(pipe, 'btn_run_all.*\\n.*btn-primary w-100', fixed = FALSE)
  # « Appliquer les seuils » = application de filtres QC, pas une destruction
  qc <- .ts_cos_read("modules/spatial/mod_spatial_qc.R")
  expect_false(grepl("btn-danger", qc, fixed = TRUE),
               info = "mod_spatial_qc.R : btn-danger sur CTA primaire (RC6)")
  expect_match(qc, "btn-primary w-100 mt-2", fixed = TRUE)
})

# --- 2. Plancher de lisibilité des fontes de sidebar (RC4) ---------------------

test_that("CQF-4 : plus de fonte < 0.8rem dans la sidebar pipeline spatial", {
  pipe <- .ts_cos_read("modules/spatial/mod_spatial_pipeline.R")
  expect_false(grepl("font-size:0\\.7rem", pipe),
               info = "texte explicatif à 0.7rem sous le plancher 0.8rem")
  expect_false(grepl("font-size:0\\.68rem", pipe),
               info = "texte explicatif à 0.68rem sous le plancher 0.8rem")
  expect_match(pipe, "font-size:0.8rem", fixed = TRUE)
})

# --- 3. Console d'import spatial : pre-wrap (RC4, reco 2026-10-29 §3.4) --------

test_that("CQF-4 : la console d'import spatial passe en white-space:pre-wrap", {
  imp <- .ts_cos_read("modules/import/mod_import_spatial.R")
  expect_match(imp,
               'overflow-y:auto; white-space:pre-wrap;',
               fixed = TRUE,
               info = "console_log : les retours à la ligne du log doivent être préservés")
})

# --- 4. Accents des littéraux UI du domaine QC spatial (RC4, §3.4) -------------

test_that("CQF-4 : les labels nav QC spatial sont accentés (et les ASCII absents)", {
  qc <- .ts_cos_read("modules/spatial/mod_spatial_qc.R")
  # code seul (hors commentaires) — même convention que les greps de contrat
  qc_code <- paste(grep("^\\s*#", readLines(file.path(ts_project_root(),
                     "modules/spatial/mod_spatial_qc.R"), warn = FALSE,
                     encoding = "UTF-8"), value = TRUE, invert = TRUE),
                   collapse = "\n")
  expect_match(qc_code, "Aperçu du jeu de données", fixed = TRUE)
  expect_match(qc_code, "Gènes spatialement variables (Moran's I)", fixed = TRUE)
  expect_false(grepl("Apercu du jeu de donnees", qc_code, fixed = TRUE),
               info = "littéral ASCII subsistant")
  expect_false(grepl("Genes spatialement variables", qc_code, fixed = TRUE),
               info = "littéral ASCII subsistant")
})

# --- 5. buttonLabel FR sur les fileInput hors modules/sc/ (RC4) -----------------

test_that("CQF-4 : tout fileInput hors modules/sc/ porte buttonLabel Parcourir...", {
  # mod_bulk_pathways.R est le site de référence (buttonLabel déjà présent).
  files <- c(
    "modules/import/mod_geo.R",
    "modules/import/mod_import_sc.R",
    "modules/import/mod_import_bulk.R",
    "modules/import/mod_import_spatial.R",
    "modules/spatial/deconv/mod_spatial_deconv_ui.R",
    "modules/bulk/mod_bulk_signatures.R",
    "modules/bulk/mod_bulk_pathways.R"
  )
  offenders <- character(0)
  for (rel in files) {
    src <- .ts_cos_read(rel)
    n_in <- length(gregexpr("fileInput\\(", src, perl = TRUE)[[1]])
    n_lbl <- length(gregexpr("buttonLabel\\s*=", src, perl = TRUE)[[1]])
    if (n_in != n_lbl) offenders <- c(offenders, rel)
  }
  expect_length(offenders, 0L)
  info <- if (length(offenders))
    paste0("fileInput sans buttonLabel (Browse... anglais) : ",
           paste(offenders, collapse = ", ")) else ""
  expect_identical(info, "")
})
