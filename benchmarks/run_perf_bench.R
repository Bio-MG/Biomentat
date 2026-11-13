# =============================================================================
# benchmarks/run_perf_bench.R — P1.3b : harnais de benchmark offline (JSON)
# Usage : Rscript benchmarks/run_perf_bench.R [--out=benchmarks/latest_perf.json]
# =============================================================================

.parse_out_arg <- function(args) {
  hit <- grep("^--out=", args, value = TRUE)
  if (length(hit) >= 1L) sub("^--out=", "", hit[1L]) else file.path("benchmarks", "latest_perf.json")
}

.bench_root <- function() {
  if (file.exists(file.path(getwd(), "tools", "pkg_skeleton", "R", "bulk_helpers.R"))) {
    return(getwd())
  }
  f_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f_arg) >= 1L) {
    cand <- dirname(dirname(normalizePath(sub("^--file=", "", f_arg[1L]), winslash = "/", mustWork = FALSE)))
    if (file.exists(file.path(cand, "tools", "pkg_skeleton", "R", "bulk_helpers.R"))) {
      return(cand)
    }
  }
  stop("Executez run_perf_bench.R depuis la racine du depot.", call. = FALSE)
}

run_perf_bench <- function(out_path = file.path("benchmarks", "latest_perf.json")) {
  root <- .bench_root()
  source(file.path(root, "tools", "pkg_skeleton", "R", "io_helpers.R"), local = FALSE)
  source(file.path(root, "tools", "pkg_skeleton", "R", "validation.R"), local = FALSE)
  source(file.path(root, "tools", "pkg_skeleton", "R", "batch_correction.R"), local = FALSE)
  source(file.path(root, "tools", "pkg_skeleton", "R", "memory.R"), local = FALSE)
  source(file.path(root, "tools", "pkg_skeleton", "R", "bulk_helpers.R"), local = FALSE)

  set.seed(20261113)
  n_genes <- 800L
  n_pairs <- 4L
  n_samples <- 2L * n_pairs
  patient <- factor(rep(sprintf("P%d", seq_len(n_pairs)), each = 2L))
  cond <- factor(rep(c("ref", "target"), n_pairs), levels = c("ref", "target"))
  base <- stats::rgamma(n_genes, shape = 2, scale = 50) + 5
  counts <- vapply(seq_len(n_samples), function(i) {
    stats::rpois(n_genes, base * ifelse(i %% 2L == 0L, 1.5, 1.0))
  }, numeric(n_genes))
  rownames(counts) <- sprintf("G%04d", seq_len(n_genes))
  colnames(counts) <- sprintf("S%d", seq_len(n_samples))
  meta <- data.frame(condition = cond, patient = patient, row.names = colnames(counts))

  tasks <- list()

  # Task 1 : limma-voom paired DE
  rss0 <- ts_process_rss_mb()
  t0 <- proc.time()[["elapsed"]]
  res_limma <- run_bulk_de_dispatch(
    engine = "limma", counts_matrix = counts, metadata = meta,
    condition_col = "condition", group_target = "target", group_ref = "ref",
    pair_col = "patient"
  )
  elapsed1 <- as.numeric(proc.time()[["elapsed"]] - t0)
  rss1 <- ts_process_rss_mb()
  tasks[[1L]] <- list(
    name          = "bulk_de_limma_paired",
    status        = "ok",
    elapsed_s     = round(elapsed1, 4),
    rss_before_mb = round(as.numeric(rss0), 2),
    rss_after_mb  = round(as.numeric(rss1), 2),
    n_genes       = nrow(res_limma),
    n_samples     = n_samples
  )

  # Task 2 : dream mixed-effects DE (skipped cleanly if upstream lme4/variancePartition incompat)
  rss2_0 <- ts_process_rss_mb()
  t2 <- proc.time()[["elapsed"]]
  dream_status <- "ok"
  dream_err <- NULL
  res_dream <- tryCatch(
    run_bulk_de_dispatch(
      engine = "dream", counts_matrix = counts, metadata = meta,
      condition_col = "condition", group_target = "target", group_ref = "ref",
      pair_col = "patient"
    ),
    error = function(e) {
      dream_status <<- "skipped"
      dream_err <<- conditionMessage(e)
      NULL
    }
  )
  elapsed2 <- as.numeric(proc.time()[["elapsed"]] - t2)
  rss2_1 <- ts_process_rss_mb()
  tasks[[2L]] <- list(
    name          = "bulk_de_dream_paired",
    status        = dream_status,
    elapsed_s     = round(elapsed2, 4),
    rss_before_mb = round(as.numeric(rss2_0), 2),
    rss_after_mb  = round(as.numeric(rss2_1), 2),
    n_genes       = if (!is.null(res_dream)) nrow(res_dream) else n_genes,
    n_samples     = n_samples,
    note          = dream_err
  )

  git_head <- tryCatch(
    trimws(paste(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE), collapse = "")),
    error = function(e) NA_character_
  )

  report <- list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    git_head  = git_head,
    tasks     = tasks
  )

  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(report, out_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(report)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  out_path <- .parse_out_arg(args)
  run_perf_bench(out_path)
}
