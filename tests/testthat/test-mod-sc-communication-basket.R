# test-mod-sc-communication-basket.R - W2.2 11.8 Gem #5 UI wiring
.mscb_root <- if (file.exists("modules/sc/mod_sc_communication.R")) "." else "../.."
.mscb_src <- function() readLines(file.path(.mscb_root, "modules/sc/mod_sc_communication.R"),
                                  warn = FALSE, encoding = "UTF-8")

test_that("server wires shared_rv$ccc_gene_basket from build_communication_basket_selection", {
  src <- paste(.mscb_src(), collapse = "\n")
  expect_match(src, "shared_rv\\$ccc_gene_basket <-")
  expect_match(src, "build_communication_basket_selection\\(")
  expect_match(src, 'exists\\("build_communication_basket_selection", mode = "function"\\)')
  expect_match(src, "comm_table_filtered_rows_selected")
})

test_that("contracts: score_min count, no new top-level assignments, plot before centrality", {
  lines <- .mscb_src()
  src <- paste(lines, collapse = "\n")
  expect_equal(lengths(regmatches(src, gregexpr('comm_f_score_min', src, fixed = TRUE))), 8L)
  exprs <- parse(text = lines, keep.source = FALSE)
  nm <- vapply(exprs, function(e) {
    if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=")) as.character(e[[2]])[1] else ""
  }, "")
  expect_true(all(nm[nzchar(nm)] %in% c("mod_sc_communication_ui",
    "mod_sc_communication_output_ui", "mod_sc_communication_server")))
  p1 <- regexpr('plotOutput(ns("comm_centrality_plot"), height = "420px")', src, fixed = TRUE)
  p2 <- regexpr('DT::dataTableOutput(ns("comm_centrality")', src, fixed = TRUE)
  expect_true(p1 > 0 && p2 > p1)
})
