# =============================================================================
# test-sc-communication-comparative.R — CC-2 (O34/D3) : communication
# comparative inter-conditions (fiche PROPOSAL_CCC_CENTRALITY_COMPARATIVE_
# 2026-11-03.md §CC-2, décisions D2/D4/D5/D6).
# =============================================================================
# GO D3 = renversement explicite du descope V1X:183 (le descope protégeait le
# moteur v1, pas une interdiction de fond). Périmètre v1 : mergeCellChat +
# rankNet(mode="comparison") + computeNetSimilarityPairwise ; liftCellChat et
# identifyCommunicationPatterns = CC-3 (hors périmètre).
#
# FIXTURE : 2 VRAIS objets CellChat construits sur le même patron que
# test-sc-communication-engine.R (240 cellules, 3 groupes A/B/C, nboot=5) —
# la voie comparative exige des OBJETS (les tables extraites ne suffisent pas,
# mergeCellChat/rankNet consomment @netP). Mesures RAM §9.3 : QC/_cc2_ram_log.txt.
# =============================================================================

source_project_file("tools/pkg_skeleton/R/io_helpers.R")
# Même chaîne de sources que test-sc-communication-engine.R:30-38 — le
# moteur + ses helpers (cellchat_input_from_matrix vit dans sc_communication_input.R).
source_project_file("tools/pkg_skeleton/R/provenance.R")
source_project_file("tools/pkg_skeleton/R/sc_velocity.R")
source_project_file("tools/pkg_skeleton/R/sc_communication.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_input.R")
source_project_file("tools/pkg_skeleton/R/sc_communication_engine.R")

if (cellchat_engine_available()) {
  # Fixture partagé des deux objets : même patron génétique (voies communes),
  # échantillons distincts — le comparative exige des IDENTITÉS identiques
  # (liftCellChat différé CC-3 : le mismatch est testé plus bas).
  .cc2_build_object <- function(seed, tag) {
    set.seed(seed)
    cc_env <- new.env(parent = emptyenv())
    utils::data(list = "CellChatDB.human", package = "CellChat", envir = cc_env)
    cc_db <- get("CellChatDB.human", envir = cc_env)
    cc_pw <- c("TGFb", "MIF", "MHC-I", "MHC-II", "GALECTIN", "CD99",
               "ANNEXIN", "VISFATIN", "COMPLEMENT", "CXCL")
    cc_sel <- cc_db$interaction[cc_db$interaction$pathway_name %in% cc_pw,
                                c("ligand", "receptor")]
    genes <- sort(unique(c(cc_sel$ligand, cc_sel$receptor)))
    genes <- genes[nzchar(genes)]
    cells <- 240L
    feats <- c(genes, paste0("FILL", tag, seq_len(60L)))
    grp <- rep(c("A", "B", "C"), length.out = cells)
    mat <- matrix(rpois(length(feats) * cells, 2), nrow = length(feats),
                  dimnames = list(feats, paste0(tag, seq_len(cells))))
    lig <- intersect(cc_sel$ligand, feats)
    rec <- intersect(cc_sel$receptor, feats)
    iA <- which(grp == "A"); iB <- which(grp == "B")
    mat[lig, iA] <- mat[lig, iA] + rpois(length(lig) * length(iA), 25)
    mat[rec, iB] <- mat[rec, iB] + rpois(length(rec) * length(iB), 25)
    input <- cellchat_input_from_matrix(
      data = mat, features = feats, labels = factor(grp), species = "human",
      group_by = "celltype", log_normalize = TRUE
    )
    o <- CellChat::createCellChat(object = input$data, meta = input$meta,
                                  group.by = cellchat_group_by_column(input))
    o@DB <- .cellchat_engine_db("human")$db
    o <- CellChat::subsetData(o)
    o <- CellChat::identifyOverExpressedGenes(o)
    o <- CellChat::identifyOverExpressedInteractions(o)
    o <- CellChat::computeCommunProb(o, nboot = 5L, seed.use = seed)
    o <- CellChat::computeCommunProbPathway(o)
    o <- CellChat::netAnalysis_computeCentrality(o, slot.name = "netP")
    # Le moteur retient l'objet amputé de @data (rétention CC-2) — même état.
    o@data <- o@data.signaling
    o
  }
  .cc2_o1 <- .cc2_build_object(20260914L, "c1")
  .cc2_o2 <- .cc2_build_object(20260915L, "c2")
}

test_that("CC-2 moteur : rankNet + similarité sur 2 objets réels (comparatif)", {
  skip_if_not(cellchat_engine_available(), "CellChat absent (dependance GitHub)")

  res <- run_ccc_comparison(list(ctrl = .cc2_o1, stim = .cc2_o2),
                            labels = c("ctrl", "stim"))
  # Table comparative (flow_table) : présente, bornée par les voies communes,
  # parts relatives sommant à 1 par condition.
  pathways_common <- intersect(.cc2_o1@netP$pathways, .cc2_o2@netP$pathways)
  df <- res$flow_table
  expect_true(is.data.frame(df))
  expect_true(nrow(df) <= length(pathways_common),
              info = paste("flow_table:", nrow(df), "vs voies communes:", length(pathways_common)))
  expect_true(all(c("pathway", "ctrl", "stim") %in% colnames(df)),
              info = paste("colonnes:", paste(colnames(df), collapse = ", ")))
  expect_true("rel_ctrl" %in% colnames(df) && "rel_stim" %in% colnames(df))
  if (nrow(df) > 0L) {
    expect_true(abs(sum(df$rel_ctrl) - 1) < 1e-6,
                info = "les parts relatives somment à 1 (ctrl)")
    expect_true(abs(sum(df$rel_stim) - 1) < 1e-6,
                info = "les parts relatives somment à 1 (stim)")
  }
  # Similarité : DÉFAUT sauté (hang/crash session mesuré sur la version
  # épinglée, même sur voies ≥ 5 — QC/_sim_diag.txt) — NULL + raison tracée ;
  # opt-in TS_CCC_SIMILARITY=1 hors test (l'appel peut tuer la session).
  expect_true(is.null(res$similarity),
              info = "similarité sautée par défaut (hang mesuré)")
  expect_true(any(grepl("computeNetSimilarityPairwise", res$warnings, fixed = TRUE)))
  expect_true(is.null(res$merged) || isS4(res$merged))
})

test_that("CC-2 moteur : refus honnêtes (singleton, netP absent, identités différentes)", {
  skip_if_not(cellchat_engine_available(), "CellChat absent (dependance GitHub)")

  # < 2 objets : l'état insufficient_objects (message i18n propre côté UI).
  err <- tryCatch(run_ccc_comparison(list(ctrl = .cc2_o1)),
                  error = function(e) e)
  expect_s3_class(err, "cellchat_engine_error")
  expect_identical(cellchat_engine_error_state(err), "insufficient_objects")
  err2 <- tryCatch(run_ccc_comparison(NULL), error = function(e) e)
  expect_identical(cellchat_engine_error_state(err2), "insufficient_objects")

  # Objet sans voies netP : invalid_input (S4 stub — prototype vide + slot).
  stub <- methods::new("CellChat")
  methods::slot(stub, "netP") <- list(pathways = character(0))
  err3 <- tryCatch(run_ccc_comparison(list(a = .cc2_o1, b = stub)),
                   error = function(e) e)
  expect_identical(cellchat_engine_error_state(err3), "invalid_input")

  # Identités différentes : liftCellChat différé CC-3 → identity_mismatch.
  o3 <- .cc2_build_object(20260916L, "c3")
  o3@meta$labels <- factor(paste0("X", o3@meta$labels))
  err4 <- tryCatch(run_ccc_comparison(list(a = .cc2_o1, b = o3)),
                   error = function(e) e)
  expect_identical(cellchat_engine_error_state(err4), "identity_mismatch")
})

test_that("CC-2 rétention : run_cellchat attache l'objet amputé de @data (additif)", {
  skip_if_not(cellchat_engine_available(), "CellChat absent (dependance GitHub)")

  # Le patron du test moteur réel : run_cellchat complet sur le fixture A/B/C.
  set.seed(20260914)
  cc_env <- new.env(parent = emptyenv())
  utils::data(list = "CellChatDB.human", package = "CellChat", envir = cc_env)
  cc_db <- get("CellChatDB.human", envir = cc_env)
  cc_pw <- c("TGFb", "MIF", "MHC-I", "MHC-II", "GALECTIN", "CD99",
             "ANNEXIN", "VISFATIN", "COMPLEMENT", "CXCL")
  cc_sel <- cc_db$interaction[cc_db$interaction$pathway_name %in% cc_pw,
                              c("ligand", "receptor")]
  genes <- sort(unique(c(cc_sel$ligand, cc_sel$receptor)))
  genes <- genes[nzchar(genes)]
  cells <- 240L
  feats <- c(genes, paste0("FILL", seq_len(60L)))
  grp <- rep(c("A", "B", "C"), length.out = cells)
  mat <- matrix(rpois(length(feats) * cells, 2), nrow = length(feats),
                dimnames = list(feats, paste0("k", seq_len(cells))))
  lig <- intersect(cc_sel$ligand, feats)
  rec <- intersect(cc_sel$receptor, feats)
  iA <- which(grp == "A"); iB <- which(grp == "B")
  mat[lig, iA] <- mat[lig, iA] + rpois(length(lig) * length(iA), 25)
  mat[rec, iB] <- mat[rec, iB] + rpois(length(rec) * length(iB), 25)
  input <- cellchat_input_from_matrix(
    data = mat, features = feats, labels = factor(grp), species = "human",
    group_by = "celltype", log_normalize = TRUE
  )
  res <- run_cellchat(input, seed = 1L, nboot = 5L)
  # Champ additif (D4) : présent et AMPLI de @data (= @data.signaling, gènes DB).
  expect_true("cellchat_object" %in% names(res))
  if (!is.null(res$cellchat_object)) {
    expect_identical(dim(res$cellchat_object@data),
                     dim(res$cellchat_object@data.signaling))
    expect_true(length(res$cellchat_object@netP$pathways) > 0L)
  }
  # Le contrat canonique n'est PAS réécrit : les champs préexistants intacts.
  expect_true(all(communication_contract_fields() %in% colnames(res$canonical_table)))
})

test_that("CC-2 UI : jumelles drawer + onglet Comparaison + serveur épinglés", {
  src <- paste(readLines(file.path(ts_project_root(), "modules", "sc",
                                   "mod_sc_communication.R"), warn = FALSE),
               collapse = "\n")
  # Les 2 jumelles .comm_drawer_groups portent l'entrée (égalité verrouillée
  # par test-mod-sc-communication-ui.R — ici on vérifie la présence ×2).
  expect_identical(lengths(regmatches(src, gregexpr(
    'list(value = "ccc_tab_comparaison", label = "Comparaison")', src,
    fixed = TRUE))), c(2L))
  # Onglet + widgets + consommateurs du registre.
  expect_true(grepl('nav_panel(i18n$t("Comparaison"), value = "ccc_tab_comparaison"', src, fixed = TRUE))
  expect_true(grepl('output$ccc_comp_selector_ui', src, fixed = TRUE))
  expect_true(grepl("observeEvent(input$btn_ccc_compare", src, fixed = TRUE))
  expect_true(grepl("run_ccc_comparison(objs, labels = pick)", src, fixed = TRUE))
  expect_true(grepl("comm_state$ccc_registry <- reg", src, fixed = TRUE))
  # Purge du registre au clear (session-scoped).
  expect_true(grepl("comm_state$ccc_registry <- NULL", src, fixed = TRUE))
  # Moteur : rétention + comparaison câblées.
  eng_src <- paste(readLines(file.path(ts_project_root(), "tools", "pkg_skeleton", "R",
                                       "sc_communication_engine.R"), warn = FALSE),
                   collapse = "\n")
  expect_true(grepl("result$cellchat_object <- ccc_object", eng_src, fixed = TRUE))
  expect_true(grepl("run_ccc_comparison <- function", eng_src, fixed = TRUE))
})
