# =============================================================================
# test-sc-velocity-nul-defense.R — W1.1 P0 runtime scVelo NUL/encoding defense
# =============================================================================
# Crash mesuré (STATUS.md §0 #21/#22) : au retour du daemon mirai scVelo
# opt-in, `scvelo_task$result()` peut contenir des octets NUL / mojibake
# (`Vecteurs calculC)s dans l'app (scVelo, opt-in) b\0\024 ...`), tuant le
# process app entier (`Error in execCallbacks(...) : embedded nul in string`).
#
# Contrat (RED avant fix) :
#   1. `modules/sc/mod_sc_velocity.R` définit `.vel_sanitize_text()` (déléguant
#      à `.ts_nul_free_text()` si chargé, sinon filtre octet 0x00 + enc2utf8)
#      et `.vel_sanitize_payload()`.
#   2. Tous les points d'entrée daemon/log sont assainis :
#      - `scvelo_task$result()` (succès et erreur)
#      - `.vel_scvelo_fail()` (`conditionMessage(e)`)
#      - `output$velocity_log <- renderText(...)`
# =============================================================================

test_that("W1.1 P0 : mod_sc_velocity.R câble la défense NUL/UTF-8 sur daemon et log", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "modules/sc/mod_sc_velocity.R"),
                         warn = FALSE), collapse = "\n")
  expect_match(src, "\\.vel_sanitize_text\\s*<-\\s*function")
  expect_match(src, "\\.vel_sanitize_payload\\s*<-\\s*function")
  expect_match(src, "\\.ts_nul_free_text")
  expect_match(src, "\\.vel_sanitize_payload\\(scvelo_task\\$result\\(\\)\\)")
  expect_match(src, "\\.vel_sanitize_text\\(conditionMessage\\(e\\)\\)")
  expect_match(src, "renderText\\(\\{\\s*\\.vel_sanitize_text\\(paste\\(velocity_log_rv\\(\\)")
})

test_that("W1.1 P0 : .vel_sanitize_text et .vel_sanitize_payload assainissent en UTF-8 sans casser les types", {
  src_lines <- readLines(file.path(ts_project_root(),
                                   "modules/sc/mod_sc_velocity.R"),
                         warn = FALSE)
  idx_start <- grep("^\\s*\\.vel_sanitize_text\\s*<-\\s*function", src_lines)[1]
  expect_false(is.na(idx_start))
  idx_end <- grep("^\\s*\\.vel_scvelo_fail\\s*<-\\s*function", src_lines)[1]
  expect_false(is.na(idx_end))

  env <- new.env(parent = globalenv())
  source_project_file("tools/pkg_skeleton/R/io_helpers.R")
  eval(parse(text = src_lines[idx_start:(idx_end - 1L)]), envir = env)

  latin1_str <- iconv("Vecteurs calcul\u00e9s", from = "UTF-8", to = "latin1")
  out_str <- env$.vel_sanitize_text(latin1_str)
  expect_identical(Encoding(out_str), "UTF-8")
  expect_identical(out_str, "Vecteurs calcul\u00e9s")

  payload <- list(
    banner = latin1_str,
    embedding_vectors = matrix(1:4, nrow = 2,
                               dimnames = list(c("c1", "c2"), c("v1", "v2"))),
    provenance = list(note = latin1_str, count = 2L)
  )
  clean <- env$.vel_sanitize_payload(payload)
  expect_identical(clean$banner, "Vecteurs calcul\u00e9s")
  expect_identical(clean$provenance$note, "Vecteurs calcul\u00e9s")
  expect_identical(clean$provenance$count, 2L)
  expect_true(is.matrix(clean$embedding_vectors))
  expect_identical(dim(clean$embedding_vectors), c(2L, 2L))
})
