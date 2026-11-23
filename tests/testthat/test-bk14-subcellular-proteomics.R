# ==============================================================================
# test-bk14-subcellular-proteomics.R — BK-14 Subcellular Proteomics Engine & UI
# ==============================================================================

.bk14_root <- function(rel) {
  if (file.exists(rel)) return(rel)
  file.path("..", "..", rel)
}

test_that("assign_subcellular_localization classifies proteins by cosine similarity", {
  source(.bk14_root("tools/pkg_skeleton/R/bulk_mnar_de.R"), local = globalenv())

  prof <- rbind(
    M1 = c(0.85, 0.10, 0.05),
    M2 = c(0.80, 0.15, 0.05),
    N1 = c(0.05, 0.85, 0.10),
    N2 = c(0.10, 0.80, 0.10),
    C1 = c(0.05, 0.05, 0.90),
    U_Mito = c(0.82, 0.12, 0.06),
    U_Nuc  = c(0.08, 0.84, 0.08),
    U_Zero = c(0.00, 0.00, 0.00)
  )
  colnames(prof) <- c("Frac1", "Frac2", "Frac3")

  markers <- data.frame(
    protein = c("M1", "M2", "N1", "N2", "C1"),
    organelle = c("Mitochondria", "Mitochondria", "Nucleus", "Nucleus", "Cytosol"),
    stringsAsFactors = FALSE
  )

  res <- assign_subcellular_localization(prof, markers, min_score = 0.7)
  expect_equal(res$n_proteins, 8L)
  expect_equal(res$n_organelles, 3L)
  expect_equal(res$assignments["U_Mito", "assigned_organelle"], "Mitochondria")
  expect_equal(res$assignments["U_Nuc", "assigned_organelle"], "Nucleus")
  expect_equal(res$assignments["U_Zero", "assigned_organelle"], "Unassigned")
})

test_that("test_subcellular_translocation detects organelle shift between conditions", {
  source(.bk14_root("tools/pkg_skeleton/R/bulk_mnar_de.R"), local = globalenv())

  ctrl <- rbind(
    M1 = c(0.9, 0.1, 0.0),
    N1 = c(0.0, 0.9, 0.1),
    C1 = c(0.1, 0.0, 0.9),
    STAT_TF = c(0.1, 0.05, 0.85)
  )
  treat <- rbind(
    M1 = c(0.9, 0.1, 0.0),
    N1 = c(0.0, 0.9, 0.1),
    C1 = c(0.1, 0.0, 0.9),
    STAT_TF = c(0.05, 0.9, 0.05)
  )
  markers <- c(M1 = "Mitochondria", N1 = "Nucleus", C1 = "Cytosol")

  tr <- test_subcellular_translocation(ctrl, treat, markers, min_score = 0.7)
  expect_equal(tr$n_common, 4L)
  expect_equal(tr$n_translocated, 1L)
  expect_equal(tr$table$protein[1L], "STAT_TF")
  expect_true(tr$table$translocated[1L])
  expect_equal(tr$table$organelle_ctrl[1L], "Cytosol")
  expect_equal(tr$table$organelle_treat[1L], "Nucleus")
})

test_that("mod_bulk_multi_ui exposes tab_subcellular and BK-14 controls", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  suppressPackageStartupMessages({
    library(shiny)
    library(bslib)
    library(DT)
  })
  local_i18n <- list(t = function(k) k)
  assign("i18n", local_i18n, envir = globalenv())
  on.exit(if (exists("i18n", envir = globalenv())) rm("i18n", envir = globalenv()), add = TRUE)
  if (!exists(".tr_plain", envir = globalenv())) {
    assign(".tr_plain", function(k) k, envir = globalenv())
  }

  source(.bk14_root("modules/bulk/mod_bulk_multi.R"), local = globalenv())
  ui_html <- as.character(mod_bulk_multi_ui("bm_test"))
  expect_true(grepl("tab_subcellular", ui_html, fixed = TRUE))
  expect_true(grepl("bm_test-subcell_markers_file", ui_html, fixed = TRUE))
  expect_true(grepl("bm_test-run_subcell", ui_html, fixed = TRUE))
})
