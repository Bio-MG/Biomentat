# =============================================================================
# test-mod-sc-pipeline.R
#   tests for modules/sc/mod_sc_pipeline.R
# =============================================================================
# 37e increment de la dette de conventions (2026-09-19, §2cw).
# Cree AVANT la conversion (test d'abord, rouge d'abord) => ce fichier doit
# ECHOUER sur les assertions de classe et sur le verrou C10 tant que le `stop()`
# de la ligne 260 ne porte pas `class = "mod_sc_pipeline_error"`.
#
# 🔑 CLASSE `mod_sc_pipeline_error` — branche (c) de CONVENTIONS.md §7 : le nom
# MECANIQUE `sc_pipeline_error` est DEJA PRIS par `R/sc/sc_pipeline.R:71`, le
# JUMEAU `R/` du meme garde (meme message « Seulement %d cellule(s) apres QC »)
# => on retombe sur la branche (b) : la classe porte le FICHIER, `mod_` inclus.
# Mesure prealable : `mod_sc_pipeline_error` n'apparait dans AUCUN fichier .R du
# depot, et `sc_pipeline_error` n'apparait que dans `R/sc/sc_pipeline.R` + son test.
# ⚠️ La doctrine §7 n'autorise AUCUNE classe partagee `R/` <-> `modules/`
# (mesure) : on ne reutilise donc PAS `sc_pipeline_error`.
#
# ---------------------------------------------------------------------------
# 🔴 VERDICTS — ce que ce lot PROUVE, et ce qu'il NE prouve PAS
# ---------------------------------------------------------------------------
#   (1) ATTEIGNABLE ................ NON PROUVE, et c'est DIT. Le site vit dans
#       le corps d'un `observeEvent` (L220) dont le PROLOGUE exige un objet
#       Seurat : `PercentageFeatureSet()` (L236), `subset()` (L253),
#       `obj@meta.data` (L245). Ce lot ne monte donc PAS le reactif : il execute
#       la GARDE ELLE-MEME, extraite de l'AST — un APPEL `eval()` EST execute
#       (§2ct.2). C'est PLUS FORT qu'un simple verrou source (le jumeau §2ci
#       s'en tenait la), et MOINS FORT qu'une preuve comportementale (§2cr).
#   (2) LA CLASSE S'ECHAPPE ........ NON — le `tryCatch` ouvert L227 se ferme
#       L471 sur un gestionnaire (L469-471) qui NOTIFIE au lieu de RELANCER
#       => il AVALE l'erreur (§2cd.2). => preuve = VERROU SOURCE, rendu
#       FALSIFIABLE par le test du gestionnaire (ci-dessous).
#   (3) CONDITIONNEL ................ OUI, prouve A L'EXECUTION sur la garde
#       extraite : 5 cellules => la condition est TRUE et la garde tire ;
#       50 cellules => la condition est FALSE et la garde ne tire pas.
#   (4) MESSAGE ..................... assere des DEUX bouts (tete ET queue) : la
#       queue est le detecteur de troncature (C16). Le message est a UN SEUL
#       argument (`sprintf`) => aucun `paste0()` a ajouter, C16 sans objet ici.
# =============================================================================

.MSP_FILE <- "modules/sc/mod_sc_pipeline.R"
# P1-6a (2026-11-07) : le CORPS du pipeline vit dans le cœur PUR
# run_sc_manual_pipeline_core (préchargé daemons) — les gardes AST du garde
# QC pointent le CŒUR, plus le module.
.MSCORE_FILE <- "tools/pkg_skeleton/R/sc_pipeline.R"

# Fragments ASCII : le litteral source porte des accents (L262) dont l'encodage
# de lecture n'est pas garanti => on assere des tranches SANS accent.
# Aiguille UNIQUE au garde QC du pipeline manuel — sc_pipeline.R porte AUSSI
# un garde QC pour l'AUTO-pipeline (message « Réduisez les seuils »), le
# fragment « Seulement %d cellule(s) » seul matcherait les deux.
.MSP_MSG_ASCII     <- "Diagnostic par filtre"
# Fragment identifiant le gestionnaire avaleur PAR SON CONTENU (jamais par sa
# position : le fichier porte plusieurs `tryCatch`).
.MSP_HANDLER_ASCII <- 'paste(.tr("Erreur pipeline:")'


# --- Les deux noeuds APPEL extraits de la SOURCE -----------------------------
# ⚠️ On ne recopie RIEN : ces noeuds sont les expressions REELLES du fichier.
#    (P1-6a : le cœur ; les params sont RÉSOLUS — qc_min_gene etc. — plus
#    input$qc_*.)
.msp_guard_call <- function() {
  ts_ast_find(ts_ast_parse(.MSCORE_FILE), function(x) {
    ts_ast_is_call_to(x, "if") &&
      grepl(.MSP_MSG_ASCII, ts_ast_deparse(x), fixed = TRUE)
  })
}

.msp_stop_call <- function() {
  ts_ast_find(ts_ast_parse(.MSCORE_FILE), function(x) {
    ts_ast_is_call_to(x, "stop") &&
      grepl(.MSP_MSG_ASCII, ts_ast_deparse(x), fixed = TRUE)
  })
}


# --- Environnement d'execution de la garde -----------------------------------
# La garde ne lit QUE : `ncol(obj)`, `n_before`, les 3 seuils résolus
# `qc_min_gene/qc_max_gene/qc_mt`, les 4 compteurs `n_ok_*`. On injecte
# exactement cette surface — un maillon de moins et la garde meurt sur
# « object not found » (faux rouge).
# `obj` : seul `ncol()` est lu par la garde => une matrice suffit, et c'est ce
# qui rend le test HERMETIQUE (ni Seurat charge, ni objet reel a construire).
.msp_env <- function(n_cells, n_before = 100L, min_gene = 200L,
                     max_gene = 5000L, mt = 20L) {
  e <- new.env(parent = globalenv())
  e$obj      <- matrix(0, nrow = 1L, ncol = n_cells)
  e$n_before <- n_before
  e$n_ok_min <- as.integer(n_cells)
  e$n_ok_max <- as.integer(n_cells)
  e$n_ok_mt  <- as.integer(n_cells)
  e$n_ok_all <- as.integer(n_cells)
  e$qc_min_gene <- min_gene
  e$qc_max_gene <- max_gene
  e$qc_mt       <- mt
  e
}

# Execute la garde REELLE et rend la condition levee, ou NULL si elle ne tire pas.
.msp_run_guard <- function(envir) {
  tryCatch({
    eval(.msp_guard_call(), envir = envir)
    NULL
  }, error = function(e) e)
}


# ---------------------------------------------------------------------------
# CONDITIONNEL — sous la borne, la condition est TRUE et la garde TIRE
# ---------------------------------------------------------------------------
test_that("le garde 260 est conditionnel : 5 cellules => condition TRUE et il TIRE", {
  e <- .msp_env(n_cells = 5L)

  # La CONDITION du `if`, reelle (2e element du noeud) — evaluee, pas recopiee.
  expect_true(isTRUE(eval(.msp_guard_call()[[2]], envir = e)))

  cond <- .msp_run_guard(e)
  expect_true(inherits(cond, "error"))
  # Tete du message : les compteurs RENDUS (pas le gabarit brut).
  expect_true(grepl("Seulement 5 cellule(s)", conditionMessage(cond), fixed = TRUE))
})


# ---------------------------------------------------------------------------
# CONDITIONNEL — au-dessus de la borne, la condition est FALSE : TEMOIN DE
# TRAVERSEE (§2cj.3). « Notre message est absent » serait vrai AUSSI si la
# garde n'avait pas ete atteinte => on prouve que la condition est FALSE, ce qui
# garantit que l'execution CONTINUE.
# ---------------------------------------------------------------------------
test_that("le garde 260 est conditionnel : 50 cellules => condition FALSE, il ne tire pas", {
  e <- .msp_env(n_cells = 50L)

  expect_false(isTRUE(eval(.msp_guard_call()[[2]], envir = e)))
  expect_null(.msp_run_guard(e))
})


# ---------------------------------------------------------------------------
# CLASSE — executee sur l'ARGUMENT du stop (jamais sur le `stop()` lui-meme :
# `eval()` d'un APPEL l'EXECUTERAIT, §2ct.2)
# ---------------------------------------------------------------------------
test_that("la classe est mod_sc_pipeline_error (execution de l'ARGUMENT du stop)", {
  call <- .msp_stop_call()
  expect_false(is.null(call))

  arg  <- call[[2]]
  cond <- eval(arg, envir = .msp_env(n_cells = 5L))

  expect_true(inherits(cond, "condition"))
  expect_true(inherits(cond, "mod_sc_pipeline_error"))

  msg <- if (inherits(cond, "condition")) conditionMessage(cond) else ""
  # Tete ET QUEUE : la queue est le detecteur de troncature (C16).
  expect_true(grepl("Seulement 5 cellule(s)", msg, fixed = TRUE))
  expect_true(grepl("% Mito (35%).", msg, fixed = TRUE))
})


# ---------------------------------------------------------------------------
# VERROU SOURCE CIBLE : le site porte bien NOTRE classe
# ---------------------------------------------------------------------------
test_that("le site 260 porte errorCondition + class = mod_sc_pipeline_error (cœur P1-6a)", {
  # ⚠️ La classe etant INOBSERVABLE a l'execution (gestionnaire avaleur), la
  # SOURCE est le seul canal ou elle est verifiable — et c'est exactement ce que
  # C10 mesure. P1-6a : le site vit dans le CŒUR (sc_pipeline.R).
  hits <- ts_ast_stop_sites(.MSCORE_FILE, .MSP_MSG_ASCII)
  expect_length(hits, 1L)                                   # exactement UN site
  expect_true(grepl("errorCondition", hits[[1]], fixed = TRUE))
  expect_true(grepl('class = "mod_sc_pipeline_error"', hits[[1]], fixed = TRUE))
})


# ---------------------------------------------------------------------------
# Justification EXECUTABLE du verrou source : le gestionnaire ne RELANCE pas
# ---------------------------------------------------------------------------
test_that("le gestionnaire englobant notifie SANS stop( (falsifiable)", {
  # P1-6a : le gestionnaire d'erreur est .apply_manual_result (module) — il
  # consomme le résultat ok=FALSE du cœur et NOTIFIE, sans jamais relancer.
  lines <- readLines(file.path(ts_project_root(), .MSP_FILE), warn = FALSE)
  i0 <- grep(".apply_manual_result <- function(res)", lines, fixed = TRUE)
  expect_true(length(i0) == 1L, info = "fonction .apply_manual_result introuvable")
  # fin de fonction : première ligne '    }' au niveau 4 après i0
  i1 <- i0 + 1L
  while (i1 < length(lines) && !grepl("^    }\\s*$", lines[i1])) i1 <- i1 + 1L
  d <- paste(lines[i0:i1], collapse = "\n")
  # Il rend la main en notifiant...
  expect_true(grepl("showNotification", d, fixed = TRUE))
  # ...et NE relance PAS. Si un jour il relancait, ce test echouerait et
  # signalerait que le lot devient prouvable a l'execution (§2cd.2).
  expect_false(grepl("stop(", d, fixed = TRUE))
})


# ---------------------------------------------------------------------------
# VERROU SOURCE : mod_sc_pipeline.R ne doit plus contribuer aucun C10
# ---------------------------------------------------------------------------
test_that("mod_sc_pipeline.R ne contribue aucun signalement C10", {
  expect_identical(ts_c10_sites(.MSP_FILE), integer(0))
})

# ---------------------------------------------------------------------------
# P1-7 + P1-10 (fix 2026-11-07, plan post-E2E #7/#10) : point de pose ASR du
# pipeline manuel (le moteur session_checkpoint n'avait aucun publish depuis ce
# module — reload = perte de 30+ min de calcul) + timeline persistante des
# étapes (Progress fermée en fin de course, notifications éphémères). Gardes
# structurelles sur la fenêtre du point de pose (forme ASR-WIRING, contrat
# §3.1 — le garde T8 de test-asr-wiring.R couvre l'appel GARDÉ verbatim).
# ---------------------------------------------------------------------------
test_that("P1-7 : le publish checkpoint est en fin de course, APRÈS l'écriture de sc_obj (P1-6a)", {
  lines <- readLines(file.path(ts_project_root(), .MSP_FILE), warn = FALSE)
  hits <- grep("global_data$sc_obj <- res$obj", lines, fixed = TRUE)
  expect_true(length(hits) >= 1L, info = "site d'écriture sc_obj introuvable (.apply_manual_result)")
  win <- paste(lines[hits[1]:(hits[1] + 16L)], collapse = " ")
  expect_match(win,
               "if (is.function(global_data$session_checkpoint_publish)) global_data$session_checkpoint_publish(\"pipeline_step_complete\")",
               fixed = TRUE,
               info = "point de pose ASR absent de la fenêtre de fin de course")
  # L'écriture PRÉCÈDE le publish (le checkpoint capture l'objet terminal) :
  pub <- grep("session_checkpoint_publish(\"pipeline_step_complete\")", lines, fixed = TRUE)
  expect_true(hits[1] < pub[1], info = "publish AVANT l'écriture de sc_obj — ordre invalide")
  # Le publish est GARDÉ (forme ASR-WIRING, contrat §3.1 — aussi verrouillé
  # par le site T8 de test-asr-wiring.R) :
  expect_match(win, "if (is.function(global_data$session_checkpoint_publish))",
               fixed = TRUE)
})

test_that("P1-10 + P1-6a : timeline persistante alimentée par [TL] (daemon) ET on_progress (repli)", {
  src <- paste(readLines(file.path(ts_project_root(), .MSP_FILE),
                         warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  core <- paste(readLines(file.path(ts_project_root(), .MSCORE_FILE),
                          warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  # UI : la timeline vit SOUS le bouton, dans le panneau (pas une notification) :
  expect_match(src, 'uiOutput(ns("pipeline_timeline"))', fixed = TRUE)
  # Machine d'état : set (courant), end (done/skip préservés, error), skip :
  expect_match(src, ".tl_set <- function(step, status = \"running\")", fixed = TRUE)
  expect_match(src, ".tl_end <- function(ok, failed_step = NULL)", fixed = TRUE)
  expect_match(src, ".tl_skip <- function(step)", fixed = TRUE)
  # Honnêteté : les étapes non exécutées ne deviennent JAMAIS « done »
  expect_match(src, 'ifelse(cur_tl$status == "skip", "skip", "done")', fixed = TRUE,
               info = "les skip ne sont pas écrasés en done à la fin de course")
  # État miroir dans shared_rv (checkpoint/rapport) :
  expect_match(src, 'state_set(shared_rv, "pipeline_timeline", tl)', fixed = TRUE)
  # CÂBLAGE ASYNC (P1-6a) : ExtendedTask + mirai + bind_task_button + tracker
  expect_match(src, "ext_task_cls$new(function(obj, params, log_file, deps = core_deps)",
               fixed = TRUE,
               info = "la tâche ExtendedTask du pipeline manuel est absente")
  expect_match(src, "run_sc_manual_pipeline_core(obj, params, log_file = log_file)",
               fixed = TRUE, info = "le cœur n'est pas exécuté en daemon mirai")
  expect_match(src, 'bslib::bind_task_button(t, "run_pipeline")', fixed = TRUE,
               info = "le bouton n'est pas lié à la tâche (anti double-clic)")
  expect_match(src, "create_reactive_tracker(session, manual_log_file", fixed = TRUE,
               info = "pas de tracker de log — la timeline ne se mettrait pas à jour en daemon")
  expect_match(src, '\\[TL\\] ', fixed = FALSE, info = "parseur de lignes [TL] absent")
  # REPLI SYNCHRONE TRACÉ (jamais muet) :
  expect_match(src, "Worker mirai indisponible", fixed = TRUE)
  expect_match(src, "run_sc_manual_pipeline_core(obj, params, on_progress", fixed = TRUE,
               info = "le repli synchrone n'alimente pas la timeline via on_progress")
  # GARDE RE-ENTRANCE : pas de double lancement
  expect_match(src, 'identical(manual_task$status(), "running")', fixed = TRUE)
  # LE CŒUR émet les marqueurs [TL] (le module ne les fabrique pas) :
  expect_true(length(gregexpr("\\[TL\\] %s\\|%s", core, fixed = TRUE)[[1]]) >= 1L,
              info = "le cœur n'émet pas de lignes [TL] step|status")
  # Le cœur est PUR : aucun appel Shiny dans run_sc_manual_pipeline_core
  # (segment = la fonction seule — vient avant le wrapper Shiny auto-pipeline)
  i0 <- gregexpr("run_sc_manual_pipeline_core <- function", core, fixed = TRUE)[[1]][1]
  i1 <- gregexpr("run_sc_auto_pipeline <- function", core, fixed = TRUE)[[1]][1]
  seg <- substr(core, i0, if (i1 > i0) i1 - 1L else nchar(core))
  expect_false(grepl("showNotification|shiny::Progress|renderUI", seg),
               info = "le cœur contient un appel Shiny — il casserait en daemon")
})

# ---------------------------------------------------------------------------
# P1-6a : falsification EXÉCUTABLE du cœur — fixture seedée 500g × 300c × 3
# clusters, le cœur tourne hors Shiny (comme en daemon : Seurat attaché, pas
# d'i18n) et rend une timeline honnête + 3 clusters retrouvés + 3 réductions.
# ---------------------------------------------------------------------------
test_that("P1-6a : run_sc_manual_pipeline_core exécute le pipeline complet (fixture seedée)", {
  skip_if_not(requireNamespace("Seurat", quietly = TRUE), "Seurat requis")
  skip_if_not(requireNamespace("SeuratObject", quietly = TRUE), "SeuratObject requis")
  # l'app attache Seurat (init différée) et les daemons via everywhere —
  # l'équivalent honnête en test :
  suppressMessages(suppressWarnings(library(Seurat)))
  # mêmes deps que le smoke headless (l'harnais ne charge pas le squelette) :
  for (dep in c("defaults.R", "thresholds.R", "io_helpers.R", "state.R",
                "provenance.R", "validation.R", "memory.R", "sc_bpcells.R",
                "sc_helpers.R", "sc_batch_metrics.R", "sc_qc_doublets.R")) {
    source_project_file(file.path("tools", "pkg_skeleton", "R", dep))
  }
  source_project_file(.MSCORE_FILE)

  set.seed(42)
  ng <- 500L; nc <- 300L
  cl <- sample(rep(1:3, length.out = nc))
  counts <- matrix(rpois(ng * nc, lambda = 2), nrow = ng, ncol = nc)
  counts[sample(seq_len(ng * nc), 2000)] <- 50L
  for (k in 1:3) {
    idx <- which(cl == k)
    genes <- ((k - 1L) * 60L + 1L):(k * 60L)
    counts[genes, idx] <- counts[genes, idx] + rpois(length(genes) * length(idx), 8)
  }
  obj <- SeuratObject::CreateSeuratObject(counts = counts)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  params <- list(qc_min_gene = 50, qc_max_gene = 10000, qc_mt = 40,
                 doublets_enabled = FALSE, bpcells_mode = "never",
                 norm_method = "log", pca_dim = 10, clust_res = 0.4,
                 cluster_algo = "1", reduction_method = "umap", seed = 989L)
  tl <- character(0)
  res <- run_sc_manual_pipeline_core(obj, params,
                                     on_progress = function(stp, sta)
                                       tl[[length(tl) + 1L]] <<- paste(stp, sta, sep = "|"))
  expect_true(isTRUE(res$ok),
              info = paste("le cœur a échoué sur la fixture :",
                           if (inherits(res$error, "error")) conditionMessage(res$error)
                           else as.character(res$error)))
  # skip honnêtes : doublets (désactivé), harmony (non demandée), tsne (déjà
  # présent ? non — tsne tourne sur cette fixture ; skip = doublets + harmony)
  expect_true(all(c("doublets", "harmony") %in% res$skipped_steps),
              info = paste("skipped :", paste(res$skipped_steps, collapse = ",")))
  # timeline honnête : chaque étape exécutée passe running→done, END|ok
  expect_true("END|ok" %in% tl, info = paste("timeline :", paste(tl, collapse = " ")))
  expect_true("qc|done" %in% tl && "cluster|done" %in% tl && "reduction|done" %in% tl)
  expect_false(any(grepl("harmony[|]done|doublets[|]done", tl)),
               info = "une étape sautée est marquée done — honnêteté bafouée")
  # produit : clusters retrouvés + réductions
  expect_true(length(unique(Seurat::Idents(res$obj))) >= 2L)
  expect_true(all(c("pca", "umap") %in% names(res$obj@reductions)))
  expect_s3_class(res$qc_snapshot, "data.frame")
})
