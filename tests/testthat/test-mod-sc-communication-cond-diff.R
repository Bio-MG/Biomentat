# =============================================================================
# test-mod-sc-communication-cond-diff.R — W2.2 §11.8 Gem #1 UI wiring
# =============================================================================

.mscd_root <- if (file.exists("modules/sc/mod_sc_communication.R")) "." else "../.."
.mscd_src <- function() {
  readLines(file.path(.mscd_root, "modules/sc/mod_sc_communication.R"),
            warn = FALSE, encoding = "UTF-8")
}

.mscd_load_views <- function() {
  if (!exists("finalize_communication_result", mode = "function")) {
    source(file.path(.mscd_root, "tools/pkg_skeleton/R/sc_communication.R"), local = FALSE)
  }
  if (!exists("build_communication_condition_diff", mode = "function") ||
      !exists("plot_communication_condition_diff", mode = "function")) {
    source(file.path(.mscd_root, "tools/pkg_skeleton/R/sc_communication_views.R"), local = FALSE)
  }
}

test_that("mod_sc_communication UI wires condition diff plot and table in ccc_tab_comparaison", {
  src <- paste(.mscd_src(), collapse = "\n")

  expect_match(src, 'plotOutput(ns("ccc_cond_diff_plot")', fixed = TRUE)
  expect_match(src, 'DT::dataTableOutput(ns("ccc_cond_diff_table"), height = "auto", fill = FALSE)', fixed = TRUE)

  pos_comp <- regexpr('nav_panel(i18n$t("Comparaison"), value = "ccc_tab_comparaison"', src, fixed = TRUE)[1]
  pos_sim <- regexpr('DT::dataTableOutput(ns("ccc_comp_similarity")', src, fixed = TRUE)[1]
  pos_diff_plot <- regexpr('plotOutput(ns("ccc_cond_diff_plot")', src, fixed = TRUE)[1]
  pos_diff_tbl <- regexpr('DT::dataTableOutput(ns("ccc_cond_diff_table")', src, fixed = TRUE)[1]
  pos_ctx_menu <- regexpr('nav_menu(i18n$t("Contextes multi-modaux")', src, fixed = TRUE)[1]

  expect_gt(pos_diff_plot, pos_comp)
  expect_gt(pos_diff_plot, pos_sim)
  expect_gt(pos_diff_tbl, pos_diff_plot)
  expect_lt(pos_diff_tbl, pos_ctx_menu)
})

test_that("mod_sc_communication server wires build_communication_condition_diff and plot_communication_condition_diff with guards", {
  src <- paste(.mscd_src(), collapse = "\n")

  expect_match(src, 'exists("build_communication_condition_diff", mode = "function")', fixed = TRUE)
  expect_match(src, 'exists("plot_communication_condition_diff", mode = "function")', fixed = TRUE)
  expect_match(src, "build_communication_condition_diff(", fixed = TRUE)
  expect_match(src, "plot_communication_condition_diff(", fixed = TRUE)
  expect_match(src, "output$ccc_cond_diff_plot <- renderPlot", fixed = TRUE)
  expect_match(src, "output$ccc_cond_diff_table <- DT::renderDataTable", fixed = TRUE)
})

test_that("mod_sc_communication server computes condition diff from registry results or LIANA collection and degrades gracefully", {
  .mscd_load_views()

  mod_env <- new.env(parent = globalenv())
  for (helper_f in c("tools/pkg_skeleton/R/defaults.R",
                     "tools/pkg_skeleton/R/error_state.R",
                     "tools/pkg_skeleton/R/validation.R",
                     "tools/pkg_skeleton/R/datatable.R",
                     "tools/pkg_skeleton/R/sc_communication.R",
                     "tools/pkg_skeleton/R/sc_communication_input.R",
                     "tools/pkg_skeleton/R/sc_communication_engine.R",
                     "tools/pkg_skeleton/R/sc_communication_views.R")) {
    sys.source(file.path(.mscd_root, helper_f), envir = mod_env)
  }
  for (subf in c("mod_sc_communication_spatial.R",
                 "mod_sc_communication_trajectory.R",
                 "mod_sc_communication_velocity.R",
                 "mod_sc_communication_perturbation.R",
                 "mod_sc_communication.R")) {
    sys.source(file.path(.mscd_root, "modules", "sc", subf), envir = mod_env)
  }

  .mk_res <- function(scores, src = "cellchat") {
    tbl <- data.frame(
      sender = c("CD4 T", "B"),
      receiver = c("B", "CD8 T"),
      ligand = c("LIG1", "LIG2"),
      receptor = c("REC1", "REC2"),
      interaction = c("LIG1 -> REC1", "LIG2 -> REC2"),
      pathway = c("PW1", "PW2"),
      score = as.numeric(scores),
      p_value = c(0.01, 0.02),
      p_adjusted = c(0.02, 0.04),
      source_method = rep(src, 2L),
      source_file = rep("mock.csv", 2L),
      source_cell_identity_level = rep("celltype", 2L),
      stringsAsFactors = FALSE
    )
    mod_env$finalize_communication_result(canonical_table = tbl, source_method = src)
  }

  res_a <- .mk_res(c(0.2, 0.7), "cellchat")
  res_b <- .mk_res(c(0.8, 0.1), "cellchat")
  res_cpdb <- .mk_res(c(0.5, 0.5), "cellphonedb")

  gd <- shiny::reactiveValues(sc_obj = NULL, language = "fr", i18n = NULL)
  srv <- shiny::reactiveValues()

  shiny::testServer(mod_env$mod_sc_communication_server,
                    args = list(global_data = gd, shared_rv = srv), {
    # 1. Missing results -> NULL (graceful degradation)
    expect_null(ccc_cond_diff_res())

    # 2. Two valid same-source results in ccc_registry -> diff data.frame computed
    comm_state$ccc_registry <- list(
      "Ctrl" = list(label = "Ctrl", object = NULL, result = res_a),
      "Stim" = list(label = "Stim", object = NULL, result = res_b)
    )
    session$setInputs(ccc_comp_pick = c("Ctrl", "Stim"))
    session$flushReact()

    df <- ccc_cond_diff_res()
    expect_s3_class(df, "ccc_condition_diff")
    expect_identical(nrow(df), 2L)
    expect_identical(attr(df, "label_a"), "Ctrl")
    expect_identical(attr(df, "label_b"), "Stim")
    expect_equal(df$delta_score[df$ligand == "LIG1"], 0.6, tolerance = 1e-12)

    # 3. Mismatched source_method -> graceful NULL via tryCatch
    comm_state$ccc_registry <- list(
      "Ctrl" = list(label = "Ctrl", object = NULL, result = res_a),
      "CPDB" = list(label = "CPDB", object = NULL, result = res_cpdb)
    )
    session$setInputs(ccc_comp_pick = c("Ctrl", "CPDB"))
    session$flushReact()
    expect_null(ccc_cond_diff_res())
  })
})

test_that("critical contracts: 0 extra top-level assignments, 8L comm_f_score_min, identical drawer twins, centrality plot before DT", {
  lines <- .mscd_src()
  src <- paste(lines, collapse = "\n")

  expect_identical(
    lengths(regmatches(src, gregexpr("comm_f_score_min", src, fixed = TRUE)))[1],
    8L
  )

  exprs <- parse(text = lines, keep.source = FALSE)
  nm <- vapply(exprs, function(e) {
    if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=")) as.character(e[[2]])[1] else ""
  }, "")
  expect_setequal(
    nm[nzchar(nm)],
    c("mod_sc_communication_ui", "mod_sc_communication_output_ui", "mod_sc_communication_server")
  )

  p1 <- regexpr('plotOutput(ns("comm_centrality_plot"), height = "420px")', src, fixed = TRUE)[1]
  p2 <- regexpr('DT::dataTableOutput(ns("comm_centrality")', src, fixed = TRUE)[1]
  expect_true(p1 > 0L && p2 > p1)

  if (exists("ts_ast_parse", mode = "function") && exists("ts_ast_find_all", mode = "function")) {
    parsed <- ts_ast_parse("modules/sc/mod_sc_communication.R")
    hits <- ts_ast_find_all(parsed, function(x) {
      if (!is.call(x) || !identical(x[[1]], as.name("<-"))) return(FALSE)
      lhs <- tryCatch(ts_ast_deparse(x[[2]]), error = function(e) "")
      identical(lhs, ".comm_drawer_groups")
    })
    expect_identical(length(hits), 2L)
    twins <- lapply(hits, function(h) eval(h[[3]], envir = new.env(parent = baseenv())))
    expect_identical(twins[[1]], twins[[2]])
  }
})
