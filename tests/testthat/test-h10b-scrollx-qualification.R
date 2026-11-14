# =============================================================================
# test-h10b-scrollx-qualification.R — §H.10b : qualification de l'opt-in
# scroll_x par table (suite du lot 2H/O28, arbitrage 2026-10-30).
# =============================================================================
# WHY THIS FILE EXISTS. Le 2H/O28 a fixé scroll_x = FALSE PAR DÉFAUT dans
# ts_datatable parce que la COMBINAISON `scrollX = TRUE` + `filter = "top"`
# produit le bug DataTables du double header / double filtre (le clonage du
# `thead` par le scroll horizontal duplique les filtres du header flottant).
# Le §H.10b (P0.4) demande de QUALIFIER l'opt-in pour les tables réellement
# larges — celles à colonnes DYNAMIQUES (une colonne par type cellulaire /
# cluster / condition), qui débordent sinon et clippent leur dernière colonne.
#
# INVARIANT QUALIFIÉ (mesuré par lecture des sources, pas par déduction) :
# tout opt-in `scroll_x = TRUE` dans les modules Spatial DOIT voyager avec
# `filter = "none"` — le combo documenté sûr. Un `scroll_x = TRUE` avec le
# filtre "top" (défaut du wrapper) RECREE le glitch et doit être refusé à
# la relecture. (La vérification visuelle du RENDU par table reste du
# ressort de la session opérateur avec les données correspondantes — même
# limitation honnête que le §H.9 : voir handoff §S.)
# =============================================================================

spatial_module_sources <- function() {
  root <- file.path("..", "..", "modules", "spatial")
  files <- list.files(root, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  stats::setNames(lapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n")),
                  files)
}

test_that("H.10b: every scroll_x = TRUE opt-in in Spatial travels with filter = 'none'", {
  srcs <- spatial_module_sources()
  expect_true(length(srcs) > 0L, info = "spatial module sources must be readable")

  findings <- list()
  for (i in seq_along(srcs)) {
    src <- srcs[[i]]
    # Every ts_datatable( call — approximate split on the call boundaries:
    # scan each occurrence of "ts_datatable(" and take up to the closing "))"
    # is fragile; instead require that WHEREVER scroll_x = TRUE appears, the
    # same call fragment also names filter = "none". We approximate a call
    # fragment as the 400 chars AFTER "ts_datatable(" up to the first
    # unbalanced ")" heuristic-free window; every measured opt-in in this
    # repo sits well inside that window.
    hits <- gregexpr("scroll_x\\s*=\\s*TRUE", src)[[1]]
    if (hits[1] == -1L) next
    for (h in hits) {
      # Call-fragment window: 200 chars BEFORE + 400 AFTER the opt-in — the
      # named args of a ts_datatable() call sit on either side of scroll_x.
      from <- max(1L, h - 200L)
      to <- min(h + 400L, nchar(src))
      window <- substr(src, from, to)
      findings[[length(findings) + 1L]] <- list(
        file = names(srcs)[i],
        has_filter_none = grepl('filter\\s*=\\s*"none"', window)
      )
    }
  }

  expect_true(length(findings) >= 2L,
              info = paste("at least the two H.10b opt-ins (deconv_outputs,",
                           "niche_composition) must exist; found", length(findings)))
  for (f in findings) {
    expect_true(f$has_filter_none,
                info = paste("scroll_x = TRUE in", f$file,
                             "must carry filter = 'none' (2H/O28 glitch guard)"))
  }
})

test_that("H.10b: the two widest dynamic-column tables are opted in", {
  # 1. deconv_outputs — deconv_props: UNE colonne PAR type cellulaire.
  src <- paste(readLines(file.path("..", "..", "modules", "spatial", "deconv",
                                   "mod_spatial_deconv_outputs.R"), warn = FALSE),
               collapse = "\n")
  expect_true(grepl("filename_base = \"deconv_outputs_table\", filter = \"none\",\n                 scroll_x = TRUE",
                    src, fixed = TRUE),
              info = "the deconv proportions table must opt into scroll_x = TRUE")

  # 2. niche_composition — UNE colonne PAR cluster/type cellulaire.
  src2 <- paste(readLines(file.path("..", "..", "modules", "spatial",
                                    "mod_spatial_niche.R"), warn = FALSE),
                collapse = "\n")
  expect_true(grepl("filename_base = \"spatial_niche_composition\", filter = \"none\",\n                   scroll_x = TRUE",
                    src2, fixed = TRUE),
              info = "the niche composition table must opt into scroll_x = TRUE")
})
