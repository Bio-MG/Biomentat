#!/usr/bin/env Rscript
# =============================================================================
# tools/biomentat_cli.R — Headless Batch CLI Runner (BK-13)
# =============================================================================
# Executes Biomentat pure-R analytical engines outside the Shiny UI:
#   - bulk_de   : Differential expression (deseq2, edger, limma, mnar)
#   - meta_de   : Horizontal multi-study DE meta-analysis (BK-07)
#   - attrition : OMOP / cohort attrition table & Mermaid diagram (BK-11)
#
# Usage:
#   Rscript tools/biomentat_cli.R --task bulk_de --counts counts.csv \
#     --meta meta.csv --condition group --ref ctrl --target case \
#     --engine limma --out de_results.csv
# =============================================================================

.cli_parse_args <- function(args) {
  out <- list()
  i <- 1L
  n <- length(args)
  while (i <= n) {
    a <- args[i]
    if (a %in% c("-h", "--help")) {
      out$help <- TRUE
      i <- i + 1L
    } else if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (grepl("=", key, fixed = TRUE)) {
        parts <- strsplit(key, "=", fixed = TRUE)[[1L]]
        out[[parts[1L]]] <- paste(parts[-1L], collapse = "=")
        i <- i + 1L
      } else if (i < n && !startsWith(args[i + 1L], "--")) {
        out[[key]] <- args[i + 1L]
        i <- i + 2L
      } else {
        out[[key]] <- TRUE
        i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }
  out
}

.cli_read_table <- function(path, rownames_col = NULL) {
  if (is.null(path) || !file.exists(path)) {
    stop(sprintf("Input file not found: %s", as.character(path)), call. = FALSE)
  }
  sep <- if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) "\t" else ","
  if (!is.null(rownames_col)) {
    utils::read.delim(path, sep = sep, row.names = rownames_col,
                      check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    utils::read.delim(path, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

.cli_usage <- function() {
  paste(
    c(
      "Biomentat Headless CLI (BK-13)",
      "Usage: Rscript tools/biomentat_cli.R --task <bulk_de|meta_de|attrition> [options]",
      "",
      "Tasks:",
      "  bulk_de    --counts <csv/tsv> --meta <csv/tsv> --condition <col> --ref <lvl> --target <lvl> [--engine limma|edger|deseq2|mnar] [--out <csv>]",
      "  meta_de    --inputs <f1,f2,...> [--method stouffer|fisher|maxP|minP] [--out <csv>]",
      "  attrition  --steps <csv/tsv> [--unit samples] [--out <csv>] [--mermaid <mmd>]"
    ),
    collapse = "\n"
  )
}

biomentat_cli_main <- function(args = commandArgs(trailingOnly = TRUE),
                               root_dir = ".") {
  opts <- .cli_parse_args(args)
  if (isTRUE(opts$help) || length(args) == 0L || is.null(opts$task)) {
    msg <- .cli_usage()
    cat(msg, "\n")
    return(invisible(list(status = "help", message = msg)))
  }

  task <- tolower(as.character(opts$task)[1L])

  if (identical(task, "bulk_de")) {
    source(file.path(root_dir, "tools/pkg_skeleton/R/batch_correction.R"), local = globalenv())
    source(file.path(root_dir, "tools/pkg_skeleton/R/bulk_mnar_de.R"), local = globalenv())
    source(file.path(root_dir, "tools/pkg_skeleton/R/bulk_helpers.R"), local = globalenv())

    req_keys <- c("counts", "meta", "condition", "ref", "target")
    missing_k <- req_keys[!vapply(req_keys, function(k) !is.null(opts[[k]]) && nzchar(as.character(opts[[k]])), logical(1L))]
    if (length(missing_k) > 0L) {
      stop(sprintf("bulk_de missing required args: %s", paste(missing_k, collapse = ", ")), call. = FALSE)
    }

    counts <- as.matrix(.cli_read_table(opts$counts, rownames_col = 1L))
    meta   <- .cli_read_table(opts$meta, rownames_col = 1L)
    engine <- opts$engine %||% "limma"

    res <- run_de_core(
      counts_matrix = counts,
      metadata = meta,
      engine = engine,
      condition_col = opts$condition,
      group_target = opts$target,
      group_ref = opts$ref
    )
    out_df <- cbind(feature = rownames(res$res), res$res)
    if (!is.null(opts$out)) {
      utils::write.csv(out_df, opts$out, row.names = FALSE)
    }
    return(invisible(list(status = "ok", task = task, engine = engine,
                          n_features = nrow(out_df), table = out_df)))
  }

  if (identical(task, "meta_de")) {
    source(file.path(root_dir, "tools/pkg_skeleton/R/bulk_meta_analysis.R"), local = globalenv())
    if (is.null(opts$inputs) || !nzchar(as.character(opts$inputs))) {
      stop("meta_de requires --inputs <file1,file2,...>", call. = FALSE)
    }
    files <- trimws(strsplit(as.character(opts$inputs), ",", fixed = TRUE)[[1L]])
    studies <- stats::setNames(lapply(files, function(f) .cli_read_table(f, rownames_col = 1L)),
                               basename(files))
    method <- opts$method %||% "stouffer"
    res <- run_horizontal_meta_de(studies, method = method, min_studies = min(2L, length(studies)))
    if (!is.null(opts$out)) {
      utils::write.csv(res$table, opts$out, row.names = FALSE)
    }
    return(invisible(list(status = "ok", task = task, method = method,
                          n_genes = res$n_genes, table = res$table)))
  }

  if (identical(task, "attrition")) {
    source(file.path(root_dir, "tools/pkg_skeleton/R/cohort_attrition.R"), local = globalenv())
    if (is.null(opts$steps) || !nzchar(as.character(opts$steps))) {
      stop("attrition requires --steps <csv/tsv>", call. = FALSE)
    }
    steps_df <- .cli_read_table(opts$steps)
    unit <- opts$unit %||% "samples"
    tbl <- build_attrition_table(steps_df, unit = unit)
    if (!is.null(opts$out)) {
      utils::write.csv(tbl, opts$out, row.names = FALSE)
    }
    mmd <- format_attrition_mermaid(tbl)
    if (!is.null(opts$mermaid)) {
      writeLines(mmd, opts$mermaid, useBytes = TRUE)
    }
    return(invisible(list(status = "ok", task = task, table = tbl, mermaid = mmd)))
  }

  stop(sprintf("Unsupported CLI task: %s", task), call. = FALSE)
}

if (sys.nframe() == 0L) {
  biomentat_cli_main()
}
