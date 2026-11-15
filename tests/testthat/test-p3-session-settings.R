# =============================================================================
# test-p3-session-settings.R — Étape 1 (P3) du proposal
# PROPOSAL_UX_UI_HARMONISATION_WORKFLOW_2026-11-15 : section « ⚙️ Paramètres »,
# reset déménagé, purge « inverse du load ».
# =============================================================================
# WHY SOURCE-LEVEL. The moved blocks live in app.R's Système sidebar — the
# operator-facing placement IS the deliverable (arbitrage « de guide aide ->
# au dessus du bouton guide aide »), so the pins assert ORDER and MOUNTING on
# the app source rather than re-rendering the whole page. The ASR suites
# (test-asr3-ui.R:403) already pin the checkpoint module MOUNT at source
# level the same way; none of the moved input ids change, so no drive/ASR
# pin needs re-basing (verified: no test references the POSITION of
# reset_app/help_btn/ram_*/session_ckpt).
# =============================================================================

.app_src <- function() {
  paste(readLines(file.path("..", "..", "app.R"), warn = FALSE), collapse = "\n")
}

test_that("P3: reset_app sits ABOVE help_btn in the Système sidebar, out of the guide modal", {
  src <- .app_src()
  i_reset <- gregexpr('actionButton("reset_app"', src, fixed = TRUE)[[1]]
  i_help  <- gregexpr('actionButton("help_btn"',  src, fixed = TRUE)[[1]]
  # Exactly ONE sidebar declaration of each (the old modal-footer reset is GONE).
  expect_equal(length(i_reset[i_reset > 0]), 1L,
               info = "reset_app must be declared exactly once (modal footer occurrence removed)")
  expect_true(all(i_help > 0))
  # Positional guarantee: reset comes BEFORE help_btn in the sidebar source.
  expect_true(i_reset[1] < i_help[1],
              info = "reset_app must sit ABOVE the Guide/Aide button (operator arbitrage)")
  # And the guide modal footer no longer carries it.
  expect_false(grepl('footer = tagList\\(\\s*\\n\\s*modalButton\\(.tr_h\\("Fermer"\\)\\),\\s*\\n\\s*actionButton\\("reset_app"', src),
               info = "the guide modal footer must not contain reset_app anymore")
})

test_that("P3: the Paramètres section hosts RAM settings AND the checkpoint toggle", {
  src <- .app_src()
  # Section header exists.
  expect_true(grepl("Param\u00e8tres", src, fixed = TRUE),
              info = "the ⚙️ Paramètres h6 must exist")
  # Positional: the checkpoint module mount moved INTO the Paramètres section —
  # i.e. AFTER the apply_ram_settings button and BEFORE the Objets Chargés h6.
  i_apply <- gregexpr('actionButton("apply_ram_settings"', src, fixed = TRUE)[[1]]
  i_ckpt  <- gregexpr('mod_session_checkpoint_ui("session_ckpt")', src, fixed = TRUE)[[1]]
  i_obj   <- gregexpr('h6(i18n$t("Objets Charg\\u00e9s")', src, fixed = TRUE)[[1]]
  if (i_obj[1] == -1) i_obj <- gregexpr('"Objets Charg', src, fixed = TRUE)[[1]]
  expect_true(i_apply[1] > 0 && i_ckpt[1] > 0 && i_obj[1] > 0,
              info = "all three anchors must exist")
  expect_true(i_ckpt[1] > i_apply[1] && i_ckpt[1] < i_obj[1],
              info = "session_ckpt must live INSIDE the Paramètres section (after RAM apply, before Objets Chargés)")
  # And the OLD Session-section mount is gone: exactly one mount overall.
  expect_equal(length(i_ckpt[i_ckpt > 0]), 1L,
               info = "session_ckpt must be mounted exactly once (moved, not duplicated)")
})

test_that("P3: purge is the INVERSE of load — explicit, settings-preserving, shared body", {
  src <- .app_src()
  # The purge button lives in the Session section, AFTER the load fileInput.
  i_load  <- gregexpr('fileInput("load_session_file"', src, fixed = TRUE)[[1]]
  i_purge <- gregexpr('actionButton("purge_data_btn"', src, fixed = TRUE)[[1]]
  expect_true(i_load[1] > 0 && i_purge[1] > 0 && i_purge[1] > i_load[1],
              info = "purge_data_btn must sit in the Session section, after Charger Session")
  # The clearing body is SHARED (one helper), so reset and purge can never drift.
  expect_true(grepl(".ts_purge_loaded_data <- function", src, fixed = TRUE),
              info = "the purge body must be extracted into .ts_purge_loaded_data()")
  expect_equal(lengths(regmatches(src, gregexpr("\\.ts_purge_loaded_data\\(", src))), 2L,
               info = "exactly 2 CALL sites (confirm_reset + purge_data_confirm); the definition uses '<- function' and does not match this pattern")
  # NO auto-purge on load: the load_session_file observer must not call the purger.
  i_obs <- gregexpr("observeEvent(input$load_session_file", src, fixed = TRUE)[[1]]
  expect_true(i_obs[1] > 0, info = "the load_session_file observer must exist")
  load_obs <- substr(src, i_obs[1], min(i_obs[1] + 3000L, nchar(src)))
  expect_false(grepl(".ts_purge_loaded_data", load_obs, fixed = TRUE),
               info = "loading a session must NEVER auto-purge (multi-load A-then-B preserved)")
  expect_true(grepl("observeEvent(input$purge_data_btn", src, fixed = TRUE) &&
               grepl("actionButton(\"purge_data_confirm\"", src, fixed = TRUE),
              info = "purge must go through an explicit confirmation modal")
})
