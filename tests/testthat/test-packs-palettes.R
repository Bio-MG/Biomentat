# =============================================================================
# test-packs-palettes.R — PK-4 (pilier 5, fiche PROPOSAL_PACKS_ECOSYSTEM §3
# jalon PK-4 ; détail idéation §7) : palettes de packs (hex embarqués) dans
# les resolvers de palettes + sélecteurs Thème existants.
#
# Créé AVANT l'implémentation (règle maison : test d'abord, rouge d'abord).
#
# Contrat : les palettes de packs sont des vecteurs HEX validés au scan (PK-1
# — attribution + license obligatoires pour le type « palette »). Résolution :
# les choix de pack sont des chaînes « <pack_id>::<palette> » ; les resolvers
# (sc_discrete_scale / sc_discrete_colors / bulk_color_scale) consultent le
# registre TS_PACKS_REGISTRY (options() — arbitrage O6 : 1 utilisateur /
# 1 processus, sévérité re-arbitrée 2026-10-18) SEULEMENT pour ces chaînes ;
# les 5 presets legacy ne touchent JAMAIS au registre (zéro changement visuel,
# philosophie theme.R). Repli honnête : registre absent ou choix inconnu =
# warning + comportement « default », JAMAIS une erreur (un mauvais
# selectInput ne doit pas blanchir tous les plots d'un module).
# =============================================================================

source_project_file("tools/pkg_skeleton/R/palettes.R")

# -- fixture ------------------------------------------------------------------

.pack_reg <- function() {
  structure(list(
    structure(list(id = "wes_packs", type = "palette", name = "Wes",
                   version = "1.0", attribution = "Wes Anderson palettes",
                   source = NULL, downloaded = NULL,
                   license = "CC-BY-4.0", files = character(0),
                   palettes = list(
                     GrandBudapest = c("#F1BB7B", "#FD5E53", "#5B7254"),
                     Second = c("#111111", "#EEEEEE")),
                   dir = tempdir(),
                   manifest = file.path(tempdir(), "manifest.json")),
              class = "biomentat_pack"),
    structure(list(id = "chem_packs", type = "libraries", name = "Chem",
                   version = "1.0", attribution = "chem", source = NULL,
                   downloaded = NULL, license = NULL,
                   files = c("libA.gmt"), palettes = NULL,
                   dir = tempdir(),
                   manifest = file.path(tempdir(), "manifest.json")),
              class = "biomentat_pack")
  ))
}

# -- packs_palette_choices / hex / attribution --------------------------------

test_that("packs_palette_choices liste les palettes des packs « palette » uniquement", {
  ch <- packs_palette_choices(.pack_reg())
  expect_setequal(names(ch), c("wes_packs::GrandBudapest", "wes_packs::Second"))
  expect_match(ch[["wes_packs::GrandBudapest"]], "Wes", fixed = TRUE)
  expect_match(ch[["wes_packs::GrandBudapest"]], "GrandBudapest", fixed = TRUE)
  expect_identical(packs_palette_choices(list()), character(0))
  expect_identical(packs_palette_choices(NULL), character(0))
})

test_that("packs_palette_hex rend le vecteur hex verbatim et refuse honnêtement", {
  hex <- packs_palette_hex(.pack_reg(), "wes_packs::GrandBudapest")
  expect_identical(hex, c("#F1BB7B", "#FD5E53", "#5B7254"))
  expect_error(packs_palette_hex(.pack_reg(), "ghost::X"), class = "pathway_error")
  expect_error(packs_palette_hex(.pack_reg(), "wes_packs::Inconnu"),
               class = "pathway_error")
  expect_error(packs_palette_hex(.pack_reg(), "default"), class = "pathway_error")
})

test_that("packs_palette_attribution porte attribution + license (fiche §9)", {
  att <- packs_palette_attribution(.pack_reg(), "wes_packs::GrandBudapest")
  expect_match(att, "Wes Anderson palettes", fixed = TRUE)
  expect_match(att, "CC-BY-4.0", fixed = TRUE)
  expect_error(packs_palette_attribution(.pack_reg(), "ghost::X"),
               class = "pathway_error")
})

# -- intégration resolvers (registre options TS_PACKS_REGISTRY) ---------------

test_that("sc_discrete_colors résout un choix de pack (hex répétés à n, nommés par niveaux)", {
  old <- options(TS_PACKS_REGISTRY = .pack_reg()); on.exit(options(old), add = TRUE)
  out <- sc_discrete_colors(c("A", "B", "C", "D", "E"), "wes_packs::GrandBudapest")
  expect_identical(unname(out),
                   c("#F1BB7B", "#FD5E53", "#5B7254", "#F1BB7B", "#FD5E53"))
  expect_identical(names(out), c("A", "B", "C", "D", "E"))
})

test_that("sc_discrete_scale / bulk_color_scale rendent une échelle manual pour un choix de pack", {
  old <- options(TS_PACKS_REGISTRY = .pack_reg()); on.exit(options(old), add = TRUE)
  s <- sc_discrete_scale("wes_packs::GrandBudapest", aesthetic = "color")
  expect_true(all(c("ScaleDiscrete", "Scale") %in% class(s)))
  b <- bulk_color_scale("wes_packs::GrandBudapest")
  expect_true(all(c("ScaleDiscrete", "Scale") %in% class(b)))
})

test_that("repli honnête : registre absent ou choix inconnu = warning + comportement default (jamais erreur)", {
  options(TS_PACKS_REGISTRY = NULL)
  expect_warning(out <- sc_discrete_colors(c("A", "B"), "wes_packs::GrandBudapest"),
                 "pack")
  # repli = chemin « default » inchangé (Dark 3, zéro changement visuel)
  expect_identical(unname(out),
                   unname(sc_discrete_colors(c("A", "B"), "default")))
  # les presets legacy ne consultent JAMAIS le registre (aucun warning)
  expect_silent(sc_discrete_colors(c("A", "B"), "okabeito"))
  old <- options(TS_PACKS_REGISTRY = .pack_reg())
  on.exit(options(old), add = TRUE)
  expect_silent(sc_discrete_colors(c("A", "B"), "okabeito"))
  # registre posé mais choix de pack inconnu → warning + repli default
  expect_warning(out2 <- sc_discrete_colors(c("A", "B"), "wes_packs::Inconnu"),
                 "pack")
  expect_identical(unname(out2), unname(sc_discrete_colors(c("A", "B"), "default")))
})

# -- câblage (assertions source-level, style freeze tests) --------------------

test_that("boot + les deux sélecteurs référencés portent les palettes de packs", {
  .src <- function(p) paste(readLines(file.path(ts_project_root(), p), warn = FALSE), collapse = "\n")
  app_src <- .src("app.R")
  expect_match(app_src, "TS_PACKS_REGISTRY", fixed = TRUE,
               info = "app.R ne pose pas le registre TS_PACKS_REGISTRY au boot")
  for (m in c("modules/bulk/mod_bulk_filter.R", "modules/sc/mod_sc_viz.R")) {
    src <- .src(m)
    expect_match(src, "packs_palette_choices(", fixed = TRUE,
                 info = paste(m, "— choix de packs absents du sélecteur"))
  }
})
