# =============================================================================
# test-ccc-gems-ui.R — W2.2 Gem #4 (alluvial) & Gem #3 (L-R verification) UI
# =============================================================================

.GEMS_UI_FILE <- "modules/sc/mod_sc_communication.R"

.gems_read_src <- function() {
  repo_root <- if (file.exists(.GEMS_UI_FILE)) "." else "../.."
  paste(readLines(file.path(repo_root, .GEMS_UI_FILE),
                  warn = FALSE, encoding = "UTF-8"),
        collapse = "\n")
}

.gems_source_views <- function() {
  repo_root <- if (file.exists("tools/pkg_skeleton/R/sc_communication_views.R")) "." else "../.."
  if (!exists("finalize_communication_result", mode = "function")) {
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication.R"), local = FALSE)
  }
  if (!exists("plot_communication_alluvial", mode = "function") ||
      !exists("build_communication_lr_expression", mode = "function")) {
    source(file.path(repo_root, "tools/pkg_skeleton/R/sc_communication_views.R"), local = FALSE)
  }
}

test_that("mod_sc_communication UI wires Gem #4 alluvial controls inside ccc_tab_reseau", {
  src <- .gems_read_src()

  expect_match(src, 'selectInput(ns("comm_alluvial_middle_by")', fixed = TRUE)
  expect_match(src, 'c("pathway", "ligand")', fixed = TRUE)
  expect_match(src, 'numericInput(ns("comm_alluvial_top_n")', fixed = TRUE)
  expect_match(src, 'value = 10L, min = 2L, max = 30L', fixed = TRUE)
  expect_match(src, 'plotOutput(ns("comm_alluvial_plot"), height = "420px")', fixed = TRUE)
  expect_match(src, 'downloadButton(ns("dl_comm_alluvial_png")', fixed = TRUE)
  expect_match(src, 'downloadButton(ns("dl_comm_alluvial_pdf")', fixed = TRUE)

  pos_reseau <- regexpr('nav_panel(i18n$t("Réseau"), value = "ccc_tab_reseau"', src, fixed = TRUE)[1]
  pos_circle <- regexpr('plotOutput(ns("comm_circle")', src, fixed = TRUE)[1]
  pos_alluvial <- regexpr('plotOutput(ns("comm_alluvial_plot")', src, fixed = TRUE)[1]
  pos_cent <- regexpr('nav_panel(i18n$t("Centralité"), value = "ccc_tab_centralite"', src, fixed = TRUE)[1]
  expect_gt(pos_circle, pos_reseau)
  expect_gt(pos_alluvial, pos_circle)
  expect_lt(pos_alluvial, pos_cent)
})

test_that("mod_sc_communication UI wires Gem #3 L-R verification controls inside ccc_tab_dotplot", {
  src <- .gems_read_src()

  expect_match(src, 'selectInput(ns("comm_lr_verify_pair")', fixed = TRUE)
  expect_match(src, 'plotOutput(ns("comm_lr_verify_plot"), height = "340px")', fixed = TRUE)

  pos_dotplot_tab <- regexpr('nav_panel(i18n$t("DotPlot"), value = "ccc_tab_dotplot"', src, fixed = TRUE)[1]
  pos_dotplot <- regexpr('plotOutput(ns("comm_dotplot")', src, fixed = TRUE)[1]
  pos_lr_pair <- regexpr('selectInput(ns("comm_lr_verify_pair")', src, fixed = TRUE)[1]
  pos_lr_plot <- regexpr('plotOutput(ns("comm_lr_verify_plot")', src, fixed = TRUE)[1]
  pos_heatmap_tab <- regexpr('nav_panel(i18n$t("Heatmap pathways"), value = "ccc_tab_heatmap"', src, fixed = TRUE)[1]
  expect_gt(pos_dotplot, pos_dotplot_tab)
  expect_gt(pos_lr_pair, pos_dotplot)
  expect_gt(pos_lr_plot, pos_lr_pair)
  expect_lt(pos_lr_plot, pos_heatmap_tab)
})

test_that("strict freeze pins in mod_sc_communication.R are preserved", {
  src <- .gems_read_src()

  n_score_min <- lengths(regmatches(src, gregexpr('comm_f_score_min', src, fixed = TRUE)))[1]
  expect_identical(n_score_min, 8L)

  pos_cent_plot <- regexpr('plotOutput(ns("comm_centrality_plot"), height = "420px")', src, fixed = TRUE)[1]
  pos_cent_tbl <- regexpr('DT::dataTableOutput(ns("comm_centrality")', src, fixed = TRUE)[1]
  expect_gt(pos_cent_plot, 0L)
  expect_lt(pos_cent_plot, pos_cent_tbl)

  parsed <- ts_ast_parse(.GEMS_UI_FILE)
  top_assigns <- vapply(as.list(parsed), function(expr) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-"))) {
      ts_ast_deparse(expr[[2]])
    } else {
      NA_character_
    }
  }, character(1))
  top_assigns <- top_assigns[!is.na(top_assigns)]
  expect_setequal(
    top_assigns,
    c("mod_sc_communication_ui", "mod_sc_communication_output_ui", "mod_sc_communication_server")
  )

  hits <- ts_ast_find_all(parsed, function(x) {
    if (!is.call(x) || !identical(x[[1]], as.name("<-"))) return(FALSE)
    lhs <- tryCatch(ts_ast_deparse(x[[2]]), error = function(e) "")
    identical(lhs, ".comm_drawer_groups")
  })
  expect_identical(length(hits), 2L)
  twins <- lapply(hits, function(h) eval(h[[3]], envir = new.env(parent = baseenv())))
  expect_identical(twins[[1]], twins[[2]])
  vals <- unlist(lapply(twins[[1]], function(g) vapply(g, function(p) p$value, character(1))))
  expect_identical(length(vals), 16L)
})

test_that("mod_sc_communication server renders alluvial and L-R verification plots via W2.2 gems", {
  src <- .gems_read_src()
  expect_match(src, 'output$comm_alluvial_plot <- renderPlot', fixed = TRUE)
  expect_match(src, 'output$dl_comm_alluvial_png <- downloadHandler', fixed = TRUE)
  expect_match(src, 'output$dl_comm_alluvial_pdf <- downloadHandler', fixed = TRUE)
  expect_match(src, 'plot_communication_alluvial(', fixed = TRUE)
  expect_match(src, 'output$comm_lr_verify_plot <- renderPlot', fixed = TRUE)
  expect_match(src, 'build_communication_lr_expression(', fixed = TRUE)
  expect_match(src, 'plot_communication_lr_expression(', fixed = TRUE)
  expect_match(src, 'comm_state$result$identity_column %||% input$comm_identity_column', fixed = TRUE)

  .gems_source_views()

  tbl <- data.frame(
    sender = c("CD4 T", "CD4 T", "B"),
    receiver = c("B", "Mono", "CD4 T"),
    ligand = c("CD40LG", "CCL5", "IL16"),
    receptor = c("CD40", "CCR5", "CD4"),
    interaction = c("CD40LG_CD40", "CCL5_CCR5", "IL16_CD4"),
    pathway = c("CD40", "CCL", "IL16"),
    score = c(0.85, 0.62, 0.41),
    p_value = c(0.001, 0.01, 0.02),
    p_adjusted = c(0.005, 0.02, 0.04),
    source_method = rep("cellchat", 3),
    source_file = rep("mock.csv", 3),
    source_cell_identity_level = rep("cell_type", 3),
    stringsAsFactors = FALSE
  )
  comm_res <- finalize_communication_result(
    canonical_table = tbl,
    source_method = "cellchat",
    identity_column = "cell_type"
  )

  p_alluv_pw <- plot_communication_alluvial(comm_res, tbl, top_n_pathways = 10L, middle_by = "pathway")
  p_alluv_lig <- plot_communication_alluvial(comm_res, tbl, top_n_pathways = 10L, middle_by = "ligand")
  expect_s3_class(p_alluv_pw, "ggplot")
  expect_s3_class(p_alluv_lig, "ggplot")

  expr <- matrix(
    c(2.1, 0.0, 1.8, 0.0,
      0.0, 1.5, 0.0, 2.4),
    nrow = 2,
    dimnames = list(c("CD40LG", "CD40"), c("c1", "c2", "c3", "c4"))
  )
  meta <- data.frame(
    cell_type = c("CD4 T", "B", "CD4 T", "B"),
    row.names = c("c1", "c2", "c3", "c4"),
    stringsAsFactors = FALSE
  )
  mock_obj <- list(expr = expr, meta.data = meta)

  lr_res <- build_communication_lr_expression(
    seurat_obj = mock_obj,
    sender = "CD4 T",
    receiver = "B",
    ligand = "CD40LG",
    receptor = "CD40",
    identity_column = comm_res$identity_column
  )
  p_lr <- plot_communication_lr_expression(lr_res)
  expect_s3_class(p_lr, "ggplot")
})
