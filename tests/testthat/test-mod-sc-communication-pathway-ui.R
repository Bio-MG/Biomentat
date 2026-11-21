# test-mod-sc-communication-pathway-ui.R — W1.3 Lot 2 §11: role + top_n wiring
.pw_src <- function() {
  root <- if (file.exists("modules/sc/mod_sc_communication.R")) "." else "../.."
  paste(readLines(file.path(root, "modules/sc/mod_sc_communication.R"), warn = FALSE), collapse = "\n")
}

test_that("pathway heatmap tab exposes role and top_n controls", {
  src <- .pw_src()
  expect_match(src, 'selectInput(ns("comm_pathway_role")', fixed = TRUE)
  expect_match(src, 'choices = c("both", "outgoing", "incoming")', fixed = TRUE)
  expect_match(src, 'numericInput(ns("comm_pathway_top_n")', fixed = TRUE)
  expect_match(src, 'value = 20L', fixed = TRUE)
})

test_that("role and top_n are passed to the heatmap view (render + exports)", {
  src <- .pw_src()
  n_calls <- lengths(regmatches(src, gregexpr("plot_communication_pathway_heatmap(", src, fixed = TRUE)))
  n_role <- lengths(regmatches(src, gregexpr("role = .pathway_role()", src, fixed = TRUE)))
  n_topn <- lengths(regmatches(src, gregexpr("top_n = .parse_top_n(input$comm_pathway_top_n)", src, fixed = TRUE)))
  expect_gte(n_calls, 3L)
  expect_equal(n_role, n_calls)
  expect_equal(n_topn, n_calls)
})

test_that("module contracts preserved (this worktree base: 7 comm_f_score_min, 8 once sibling lots land)", {
  src <- .pw_src()
  n <- lengths(regmatches(src, gregexpr('"comm_f_score_min"', src, fixed = TRUE)))
  expect_true(n %in% c(7L, 8L))
  a <- regexpr('plotOutput(ns("comm_centrality_plot"), height = "420px")', src, fixed = TRUE)[1]
  b <- regexpr('DT::dataTableOutput(ns("comm_centrality"))', src, fixed = TRUE)[1]
  if (a > 0L && b > 0L) expect_lt(a, b)
  expect_gt(length(parse(text = src, keep.source = FALSE)), 0L)
})