# =============================================================================
# test-nav-dropdown-close.R -- W1.2 #3b : close nav_menu Bootstrap dropdowns
# =============================================================================
# Verifies that app.R registers a client-side JS handler in tags$head(...)
# that closes open .nav-item.dropdown menus on Escape keydown and outside click,
# using Bootstrap's Dropdown API when available and cleaning up .show /
# aria-expanded="false".
#
# Pure ASCII (CONVENTIONS C7). Does not source app.R (test-app-sourcing.R).
# =============================================================================

.ndc_app_text <- function() {
  paste(readLines(file.path(ts_project_root(), "app.R"), warn = FALSE),
        collapse = "\n")
}

.ndc_extract_js <- function(src) {
  m <- regexec(
    "tags\\$head\\(\\s*tags\\$script\\(HTML\\(\\s*\"((?:[^\"\\\\]|\\\\.)*)\"\\s*\\)\\)\\s*\\)",
    src,
    perl = TRUE
  )
  parts <- regmatches(src, m)[[1L]]
  if (length(parts) < 2L) return(NULL)
  js <- parts[2L]
  js <- gsub('\\\\"', '"', js)
  js <- gsub("\\\\n", "\n", js)
  js
}

test_that("nav-dropdown-close extractor rejects empty / missing tags$head script", {
  expect_null(.ndc_extract_js("ui_static <- page_navbar(title = 'x')"))
  expect_null(.ndc_extract_js('tags$script(HTML("console.log(1);"))'))
})

test_that("app.R defines nav_menu dropdown close handler in tags$head(tags$script(HTML(...)))", {
  src <- .ndc_app_text()
  js <- .ndc_extract_js(src)
  expect_false(is.null(js),
               info = "app.R must contain tags$head(tags$script(HTML(...))) for nav dropdown close")
  expect_match(js, "document\\.addEventListener\\(\\s*['\"]keydown['\"]", perl = TRUE)
  expect_match(js, "Escape", fixed = TRUE)
  expect_match(js, "document\\.addEventListener\\(\\s*['\"]click['\"]", perl = TRUE)
  expect_match(js, "\\.nav-item\\.dropdown", perl = TRUE)
  expect_match(js, "\\.dropdown-menu\\.show", perl = TRUE)
  expect_match(js, "aria-expanded", fixed = TRUE)
  expect_match(js, "bootstrap\\.Dropdown", perl = TRUE)
})

test_that("nav-dropdown-close JS closes on Escape and outside click, ignores inside click", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "node.exe not available for JS behavioral check")

  src <- .ndc_app_text()
  js <- .ndc_extract_js(src)
  expect_false(is.null(js))

  harness <- paste(
    "const assert = require('assert');",
    "const listeners = {};",
    "function makeEl(opts) {",
    "  const classes = new Set(opts.classes || []);",
    "  const attrs = Object.assign({}, opts.attrs || {});",
    "  return {",
    "    _insideDropdown: !!opts.insideDropdown,",
    "    classList: {",
    "      contains: (c) => classes.has(c),",
    "      add: (c) => classes.add(c),",
    "      remove: (c) => classes.delete(c)",
    "    },",
    "    setAttribute: (k, v) => { attrs[k] = String(v); },",
    "    getAttribute: (k) => attrs[k],",
    "    closest: (sel) => (opts.insideDropdown && sel === '.nav-item.dropdown') ? {} : null",
    "  };",
    "}",
    "let toggleEl, menuEl, hideCalls = 0;",
    "function resetState() {",
    "  hideCalls = 0;",
    "  toggleEl = makeEl({ classes: ['dropdown-toggle', 'show'], attrs: { 'aria-expanded': 'true' }, insideDropdown: true });",
    "  menuEl = makeEl({ classes: ['dropdown-menu', 'show'], insideDropdown: true });",
    "}",
    "global.document = {",
    "  addEventListener: (ev, fn) => { (listeners[ev] = listeners[ev] || []).push(fn); },",
    "  querySelectorAll: (sel) => {",
    "    if (sel.indexOf('aria-expanded') !== -1 || sel.indexOf('dropdown-toggle') !== -1) {",
    "      return toggleEl.getAttribute('aria-expanded') === 'true' || toggleEl.classList.contains('show') ? [toggleEl] : [];",
    "    }",
    "    if (sel.indexOf('.dropdown-menu.show') !== -1) {",
    "      return menuEl.classList.contains('show') ? [menuEl] : [];",
    "    }",
    "    return [];",
    "  }",
    "};",
    "global.window = {",
    "  bootstrap: {",
    "    Dropdown: {",
    "      getInstance: () => ({ hide: () => { hideCalls++; } }),",
    "      getOrCreateInstance: () => ({ hide: () => { hideCalls++; } })",
    "    }",
    "  }",
    "};",
    js,
    "// 1. Non-Escape keydown leaves dropdown open",
    "resetState();",
    "listeners.keydown.forEach((fn) => fn({ key: 'Enter' }));",
    "assert.strictEqual(menuEl.classList.contains('show'), true);",
    "assert.strictEqual(toggleEl.getAttribute('aria-expanded'), 'true');",
    "// 2. Escape keydown closes dropdown via Bootstrap API + class/aria cleanup",
    "listeners.keydown.forEach((fn) => fn({ key: 'Escape' }));",
    "assert.strictEqual(hideCalls, 1);",
    "assert.strictEqual(menuEl.classList.contains('show'), false);",
    "assert.strictEqual(toggleEl.classList.contains('show'), false);",
    "assert.strictEqual(toggleEl.getAttribute('aria-expanded'), 'false');",
    "// 3. Inside click does NOT close dropdown",
    "resetState();",
    "listeners.click.forEach((fn) => fn({ target: makeEl({ insideDropdown: true }) }));",
    "assert.strictEqual(hideCalls, 0);",
    "assert.strictEqual(menuEl.classList.contains('show'), true);",
    "assert.strictEqual(toggleEl.getAttribute('aria-expanded'), 'true');",
    "// 4. Outside click closes dropdown even without window.bootstrap",
    "resetState();",
    "delete global.window.bootstrap;",
    "listeners.click.forEach((fn) => fn({ target: makeEl({ insideDropdown: false }) }));",
    "assert.strictEqual(menuEl.classList.contains('show'), false);",
    "assert.strictEqual(toggleEl.classList.contains('show'), false);",
    "assert.strictEqual(toggleEl.getAttribute('aria-expanded'), 'false');",
    "console.log('OK');",
    sep = "\n"
  )

  tmp <- tempfile(fileext = ".js")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(harness, tmp, useBytes = TRUE)
  out <- system2(node, tmp, stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  expect_true(is.null(status) || status == 0L,
              info = paste(out, collapse = "\n"))
  expect_true(any(grepl("^OK$", out)))
})
