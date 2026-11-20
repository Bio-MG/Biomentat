# test-bulk-pathways-pk3-notice.R -- W2.3 #3c : reserve methodologique PK-3
# La BH GLOBALE de l'onglet Convergence est calculee sur l'union des genes
# post-cutoff ; la note l'explicite honnetement (texte statique, aucune cle i18n).

source_project_file("modules/bulk/mod_bulk_pathways.R")

test_that("la note de reserve PK-3 existe et dit BH globale + union post-cutoff", {
  expect_true(exists(".bulk_pk3_reserve_note", mode = "function"))
  html <- as.character(.bulk_pk3_reserve_note())
  expect_match(html, "BH", fixed = TRUE)
  expect_match(html, "globale", ignore.case = TRUE)
  expect_match(html, "post-cutoff", fixed = TRUE)
  expect_match(html, "union", fixed = TRUE)
  expect_match(html, "pk3-reserve-note", fixed = TRUE)
})

test_that("convergence_ui rend la note avant le contenu de convergence", {
  src <- paste(readLines(file.path(ts_project_root(),
                                   "modules/bulk/mod_bulk_pathways.R"),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_match(src, "output\\$convergence_ui\\s*<-\\s*renderUI\\(\\s*tagList\\(\\s*\\.bulk_pk3_reserve_note\\(\\)")
})
