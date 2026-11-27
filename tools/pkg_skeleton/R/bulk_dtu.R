# =============================================================================
# tools/pkg_skeleton/R/bulk_dtu.R — Differential Transcript Usage (DTU) engine
# =============================================================================
# Pure domain logic (zero Shiny symbols). Implements transcript-level isoform
# switching / Differential Transcript Usage (DTU) analysis via
# edgeR::diffSpliceDGE() and edgeR::topSpliceDGE().
#
# Dependencies: edgeR (checked via requireNamespace), stats, utils.
# =============================================================================

#' Run Differential Transcript Usage (DTU) analysis via edgeR::diffSpliceDGE
#'
#' @param counts Numeric/integer matrix (or coercible data.frame) of transcript
#'   counts with rownames = transcript IDs and columns = samples.
#' @param tx2gene data.frame mapping transcripts to genes, containing
#'   `transcript_id` and `gene_id` (also accepts `tx_id` / `gene_name`).
#' @param design_or_groups Either a numeric design matrix (`nrow == ncol(counts)`)
#'   or a factor/character vector of sample group labels (`length == ncol(counts)`,
#'   at least 2 levels).
#' @param coef Integer or character coefficient index/name in the design matrix
#'   to test (default `2L`).
#' @param min_count Minimum count per transcript passed to `edgeR::filterByExpr`.
#' @param min_total_count Minimum total count per transcript passed to
#'   `edgeR::filterByExpr`.
#' @param robust Logical; whether to use robust empirical Bayes squeezing in
#'   `edgeR::glmQLFit`.
#' @param test_method Character; gene-level test method: `"Simes"`, `"F"` (maps
#'   to `"gene"` in `edgeR::topSpliceDGE`), or `"t"` (maps to `"exon"`).
#' @return A list with `gene_table`, `transcript_table`, `n_genes_tested`,
#'   `n_transcripts_tested`, `test_method`, and `provenance`.
run_bulk_dtu <- function(counts,
                         tx2gene,
                         design_or_groups,
                         coef = 2L,
                         min_count = 10L,
                         min_total_count = 15L,
                         robust = TRUE,
                         test_method = c("Simes", "F", "t")) {
  .stop <- function(msg) {
    errorCondition(msg, class = "bulk_dtu_error")
  }

  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop(.stop("Package 'edgeR' is required for DTU analysis but is not installed."))
  }

  test_method <- match.arg(test_method)

  # --- Validate counts ---
  if (is.data.frame(counts)) {
    counts <- as.matrix(counts)
  }
  if (is.null(counts) || !is.matrix(counts) || !is.numeric(counts)) {
    stop(.stop(sprintf(
      "DTU counts must be a numeric/integer matrix or data.frame (got: %s).",
      if (is.null(counts)) "NULL" else paste(class(counts), collapse = "/")
    )))
  }
  if (nrow(counts) < 2L || ncol(counts) < 2L) {
    stop(.stop("DTU counts matrix must have at least 2 transcripts and 2 samples."))
  }
  tx_names <- rownames(counts)
  if (is.null(tx_names) || anyNA(tx_names) || any(!nzchar(tx_names))) {
    stop(.stop("DTU counts matrix must have non-empty rownames corresponding to transcript IDs."))
  }
  if (anyNA(counts) || any(!is.finite(counts)) || any(counts < 0)) {
    stop(.stop("DTU counts matrix must contain finite non-negative values without NAs."))
  }

  # --- Validate tx2gene ---
  if (is.null(tx2gene) || !is.data.frame(tx2gene)) {
    stop(.stop("tx2gene must be a data.frame mapping transcripts to genes."))
  }
  tx_col <- if ("transcript_id" %in% colnames(tx2gene)) {
    "transcript_id"
  } else if ("tx_id" %in% colnames(tx2gene)) {
    "tx_id"
  } else {
    NULL
  }
  gene_col <- if ("gene_id" %in% colnames(tx2gene)) {
    "gene_id"
  } else if ("gene_name" %in% colnames(tx2gene)) {
    "gene_name"
  } else {
    NULL
  }
  if (is.null(tx_col) || is.null(gene_col)) {
    stop(.stop("tx2gene must contain 'transcript_id' (or 'tx_id') and 'gene_id' (or 'gene_name') columns."))
  }

  map_df <- data.frame(
    transcript_id = as.character(tx2gene[[tx_col]]),
    gene_id = as.character(tx2gene[[gene_col]]),
    stringsAsFactors = FALSE
  )
  valid_map <- !is.na(map_df$transcript_id) & nzchar(map_df$transcript_id) &
    !is.na(map_df$gene_id) & nzchar(map_df$gene_id)
  map_df <- map_df[valid_map, , drop = FALSE]
  map_df <- map_df[!duplicated(map_df$transcript_id), , drop = FALSE]

  # --- Validate design_or_groups ---
  if (is.matrix(design_or_groups) && is.numeric(design_or_groups)) {
    design <- design_or_groups
    if (nrow(design) != ncol(counts)) {
      stop(.stop(sprintf(
        "Design matrix row count (%d) does not match sample count (%d).",
        nrow(design), ncol(counts)
      )))
    }
    if (ncol(design) < 2L || anyNA(design)) {
      stop(.stop("Design matrix must have at least 2 columns and no missing values."))
    }
  } else if (is.factor(design_or_groups) || is.character(design_or_groups)) {
    if (length(design_or_groups) != ncol(counts)) {
      stop(.stop(sprintf(
        "Group vector length (%d) does not match sample count (%d).",
        length(design_or_groups), ncol(counts)
      )))
    }
    if (anyNA(design_or_groups) || any(!nzchar(as.character(design_or_groups)))) {
      stop(.stop("Group vector must not contain NA or empty values."))
    }
    grp <- if (is.factor(design_or_groups)) {
      droplevels(design_or_groups)
    } else {
      factor(design_or_groups)
    }
    if (nlevels(grp) < 2L) {
      stop(.stop("Group vector must have at least 2 distinct levels for DTU testing."))
    }
    design <- stats::model.matrix(~ grp)
  } else {
    stop(.stop("design_or_groups must be a numeric design matrix or a factor/character group vector."))
  }

  # --- Intersect transcripts and enforce >= 2 isoforms per gene ---
  common_tx <- intersect(tx_names, map_df$transcript_id)
  if (length(common_tx) < 2L) {
    stop(.stop("Fewer than 2 transcripts overlap between counts rownames and tx2gene."))
  }

  tx_to_gene <- stats::setNames(map_df$gene_id, map_df$transcript_id)
  gene_vec <- unname(tx_to_gene[common_tx])

  gene_counts <- table(gene_vec)
  multi_genes <- names(gene_counts)[gene_counts >= 2L]
  if (length(multi_genes) == 0L) {
    stop(.stop("No multi-isoform genes (>= 2 transcripts per gene) found after intersecting counts with tx2gene."))
  }

  keep_multi <- gene_vec %in% multi_genes
  common_tx <- common_tx[keep_multi]
  gene_vec <- gene_vec[keep_multi]
  counts_sub <- counts[common_tx, , drop = FALSE]

  # --- Build DGEList and filter ---
  dge <- edgeR::DGEList(counts = counts_sub)
  keep_expr <- edgeR::filterByExpr(
    dge,
    design = design,
    min.count = min_count,
    min.total.count = min_total_count
  )
  if (sum(keep_expr) < 2L) {
    stop(.stop("Fewer than 2 transcripts remain after edgeR::filterByExpr filtering."))
  }

  dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
  common_tx <- common_tx[keep_expr]
  gene_vec <- gene_vec[keep_expr]

  # Re-enforce >= 2 transcripts per gene after expression filtering
  gene_counts_post <- table(gene_vec)
  multi_genes_post <- names(gene_counts_post)[gene_counts_post >= 2L]
  if (length(multi_genes_post) == 0L) {
    stop(.stop("No multi-isoform genes (>= 2 transcripts per gene) remain after expression filtering."))
  }
  keep_post <- gene_vec %in% multi_genes_post
  if (!all(keep_post)) {
    dge <- dge[keep_post, , keep.lib.sizes = FALSE]
    common_tx <- common_tx[keep_post]
    gene_vec <- gene_vec[keep_post]
  }

  # --- Fit edgeR QLF + diffSpliceDGE ---
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = isTRUE(robust))
  sp <- edgeR::diffSpliceDGE(
    fit,
    coef = coef,
    geneid = gene_vec,
    exonid = common_tx,
    verbose = FALSE
  )

  top_test <- switch(
    test_method,
    Simes = "Simes",
    F = "gene",
    t = "exon"
  )

  gene_table <- edgeR::topSpliceDGE(sp, test = top_test, number = Inf)
  transcript_table <- edgeR::topSpliceDGE(sp, test = "exon", number = Inf)

  list(
    gene_table = gene_table,
    transcript_table = transcript_table,
    n_genes_tested = length(unique(gene_vec)),
    n_transcripts_tested = length(common_tx),
    test_method = test_method,
    provenance = list(
      engine = "edgeR::diffSpliceDGE",
      edgeR_version = as.character(utils::packageVersion("edgeR"))
    )
  )
}

#' Summarize transcript-level DTU results for a single gene
#'
#' @param dtu_res Result list returned by `run_bulk_dtu()`.
#' @param gene_id Character scalar specifying the gene ID to extract.
#' @return A `data.frame` of transcript-level DTU statistics for `gene_id`,
#'   sorted by `FDR` and `P.Value`.
summarize_dtu_gene <- function(dtu_res, gene_id) {
  .stop <- function(msg) {
    errorCondition(msg, class = "bulk_dtu_error")
  }

  if (!is.list(dtu_res) || is.null(dtu_res$transcript_table) ||
      !is.data.frame(dtu_res$transcript_table)) {
    stop(.stop("dtu_res must be a valid result list from run_bulk_dtu() containing transcript_table."))
  }
  if (!is.character(gene_id) || length(gene_id) != 1L || is.na(gene_id) || !nzchar(gene_id)) {
    stop(.stop("gene_id must be a single non-empty character string."))
  }

  tt <- dtu_res$transcript_table
  gene_col <- if ("GeneID" %in% colnames(tt)) {
    "GeneID"
  } else if ("gene_id" %in% colnames(tt)) {
    "gene_id"
  } else {
    stop(.stop("transcript_table does not contain a 'GeneID' or 'gene_id' column."))
  }

  sub <- tt[tt[[gene_col]] == gene_id, , drop = FALSE]
  if (nrow(sub) > 1L) {
    fdr_vals <- if ("FDR" %in% colnames(sub)) sub$FDR else rep(NA_real_, nrow(sub))
    pval_vals <- if ("P.Value" %in% colnames(sub)) sub$P.Value else rep(NA_real_, nrow(sub))
    ord <- order(fdr_vals, pval_vals, na.last = TRUE)
    sub <- sub[ord, , drop = FALSE]
  }
  rownames(sub) <- NULL
  sub
}

#' Detect transcript or gene ID format from a character vector of feature IDs
#'
#' @param tx_ids Character vector of non-empty feature IDs.
#' @return Character scalar: `"enst"`, `"ensmust"`, `"refseq"`, `"ucsc"`,
#'   `"gene_ensembl"`, `"gene_entrez"`, or `"unknown"`.
.detect_tx_id_kind <- function(tx_ids) {
  base_ids <- sub("\\.[0-9]+$", "", tx_ids)
  n <- length(base_ids)
  if (n == 0L) return("unknown")
  frac <- function(pat) mean(grepl(pat, base_ids, perl = TRUE))
  if (frac("^ENST[0-9]{6,}$") >= 0.2) return("enst")
  if (frac("^ENSMUST[0-9]{6,}$") >= 0.2) return("ensmust")
  if (frac("^[NX][MR]_[0-9]+$") >= 0.2) return("refseq")
  if (frac("^uc[0-9a-z]{6,}$") >= 0.2) return("ucsc")
  if (frac("^ENS(MUS)?G[0-9]{6,}$") >= 0.2) return("gene_ensembl")
  if (frac("^[0-9]+$") >= 0.5) return("gene_entrez")
  "unknown"
}

#' Query an AnnotationDbi OrgDb package for transcript-to-gene mapping
#'
#' @param tx_ids Character vector of transcript IDs (may include `.version`).
#' @param org_pkg Character scalar OrgDb package name (e.g. `"org.Hs.eg.db"`).
#' @param keytype Character scalar keytype (`"ENSEMBLTRANS"` or `"REFSEQ"`).
#' @return `data.frame(transcript_id, gene_id)` preserving original `tx_ids`.
.map_tx_via_orgdb <- function(tx_ids, org_pkg, keytype) {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace(org_pkg, quietly = TRUE)) {
    return(NULL)
  }
  db <- getExportedValue(org_pkg, org_pkg)
  stripped <- sub("\\.[0-9]+$", "", tx_ids)
  sel <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      db, keys = unique(stripped), keytype = keytype,
      columns = c(keytype, "SYMBOL", "ENTREZID")
    )),
    error = function(e) NULL
  )
  if (is.null(sel) || nrow(sel) == 0L) {
    return(data.frame(transcript_id = character(0L), gene_id = character(0L), stringsAsFactors = FALSE))
  }
  sym <- as.character(sel$SYMBOL)
  ent <- as.character(sel$ENTREZID)
  gene_val <- ifelse(!is.na(sym) & nzchar(sym), sym, ent)
  key_col <- as.character(sel[[keytype]])

  # Map stripped keys back to original versioned tx_ids (supports 1:many if duplicate versions)
  split_orig <- split(tx_ids, stripped)
  mapped_list <- lapply(seq_len(nrow(sel)), function(i) {
    orig <- split_orig[[key_col[i]]]
    if (is.null(orig) || is.na(gene_val[i]) || !nzchar(gene_val[i])) return(NULL)
    data.frame(transcript_id = orig, gene_id = rep(gene_val[i], length(orig)), stringsAsFactors = FALSE)
  })
  mapped_list <- Filter(Negate(is.null), mapped_list)
  if (length(mapped_list) == 0L) {
    return(data.frame(transcript_id = character(0L), gene_id = character(0L), stringsAsFactors = FALSE))
  }
  do.call(rbind, mapped_list)
}

#' Build a transcript-to-gene (tx2gene) mapping table in-app
#'
#' @param counts Transcript count matrix/data.frame (uses `rownames(counts)`)
#'   or a character vector of transcript IDs.
#' @param source Character scalar: `"auto"` (detects Ensembl/RefSeq/UCSC offline),
#'   `"org_hs"` (`org.Hs.eg.db`), `"org_mm"` (`org.Mm.eg.db`), `"txdb"`
#'   (installed TxDb/AnnotationDbi package), or `"gtf"` (uploaded GTF/GFF file).
#' @param txdb_pkg Character scalar TxDb package name (default
#'   `"TxDb.Hsapiens.UCSC.hg19.knownGene"`).
#' @param gtf_path Optional path to a GTF/GFF file when `source == "gtf"`.
#' @return Named list with `table` (`data.frame(transcript_id, gene_id)`),
#'   `n_input`, `n_mapped`, `n_unmapped`, `unmapped_ids`, and `source`.
#' @export
build_tx2gene <- function(counts,
                          source = c("auto", "txdb", "gtf", "org_hs", "org_mm"),
                          txdb_pkg = "TxDb.Hsapiens.UCSC.hg19.knownGene",
                          gtf_path = NULL) {
  .stop <- function(msg) errorCondition(msg, class = "bulk_dtu_error")
  source <- match.arg(source)

  tx_ids <- if (is.character(counts)) {
    counts
  } else if (is.matrix(counts) || is.data.frame(counts) || methods::is(counts, "Matrix")) {
    rownames(counts)
  } else {
    NULL
  }
  tx_ids <- unique(trimws(as.character(tx_ids)))
  tx_ids <- tx_ids[!is.na(tx_ids) & nzchar(tx_ids)]
  if (length(tx_ids) == 0L) {
    stop(.stop("counts must provide non-empty transcript IDs (rownames or character vector)."))
  }

  effective_source <- source
  id_kind <- .detect_tx_id_kind(tx_ids)

  if (identical(source, "auto")) {
    if (id_kind %in% c("gene_ensembl", "gene_entrez")) {
      stop(.stop(
        "Gene-level IDs detected (ENSG/Entrez). DTU requires transcript-level counts (>= 2 isoforms per gene, e.g. ENST, ENSMUST, NM_, or uc IDs)."
      ))
    }
    effective_source <- switch(
      id_kind,
      enst    = "org_hs",
      ensmust = "org_mm",
      refseq  = "refseq",
      ucsc    = "txdb",
      "txdb"
    )
  }

  if (effective_source %in% c("org_hs", "org_mm", "refseq")) {
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      stop(.stop("Package 'AnnotationDbi' is required for OrgDb transcript mapping."))
    }
    if (identical(effective_source, "refseq")) {
      tbl <- .map_tx_via_orgdb(tx_ids, "org.Hs.eg.db", "REFSEQ")
      if (is.null(tbl) || nrow(tbl) == 0L) {
        tbl_mm <- .map_tx_via_orgdb(tx_ids, "org.Mm.eg.db", "REFSEQ")
        if (!is.null(tbl_mm) && nrow(tbl_mm) > 0L) tbl <- tbl_mm
      }
      if (is.null(tbl)) {
        stop(.stop("Neither 'org.Hs.eg.db' nor 'org.Mm.eg.db' is installed for RefSeq mapping."))
      }
    } else {
      org_pkg <- if (identical(effective_source, "org_hs")) "org.Hs.eg.db" else "org.Mm.eg.db"
      if (!requireNamespace(org_pkg, quietly = TRUE)) {
        stop(.stop(sprintf("Annotation package '%s' is not installed.", org_pkg)))
      }
      kt <- if (identical(id_kind, "refseq")) "REFSEQ" else "ENSEMBLTRANS"
      tbl <- .map_tx_via_orgdb(tx_ids, org_pkg, kt)
      if ((is.null(tbl) || nrow(tbl) == 0L) && identical(kt, "ENSEMBLTRANS")) {
        tbl_ref <- .map_tx_via_orgdb(tx_ids, org_pkg, "REFSEQ")
        if (!is.null(tbl_ref) && nrow(tbl_ref) > 0L) tbl <- tbl_ref
      }
    }
  } else if (identical(effective_source, "txdb")) {
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      stop(.stop("Package 'AnnotationDbi' is required for source = 'txdb'."))
    }
    if (!is.character(txdb_pkg) || length(txdb_pkg) != 1L || !nzchar(txdb_pkg) ||
        !requireNamespace(txdb_pkg, quietly = TRUE)) {
      stop(.stop(sprintf("Annotation package '%s' is not installed.", as.character(txdb_pkg))))
    }
    txdb <- getExportedValue(txdb_pkg, txdb_pkg)
    sel <- tryCatch(
      suppressMessages(AnnotationDbi::select(
        txdb, keys = tx_ids, keytype = "TXNAME", columns = c("TXNAME", "GENEID")
      )),
      error = function(e) NULL
    )
    if (is.null(sel) || nrow(sel) == 0L) {
      stripped <- sub("\\.[0-9]+$", "", tx_ids)
      sel <- tryCatch(
        suppressMessages(AnnotationDbi::select(
          txdb, keys = unique(stripped), keytype = "TXNAME", columns = c("TXNAME", "GENEID")
        )),
        error = function(e) NULL
      )
      if (!is.null(sel) && nrow(sel) > 0L) {
        idx <- match(sel$TXNAME, stripped)
        sel$TXNAME[!is.na(idx)] <- tx_ids[idx[!is.na(idx)]]
      }
    }
    if (is.null(sel) || nrow(sel) == 0L) {
      tbl <- data.frame(transcript_id = character(0L), gene_id = character(0L), stringsAsFactors = FALSE)
    } else {
      tbl <- data.frame(
        transcript_id = as.character(sel$TXNAME),
        gene_id = as.character(sel$GENEID),
        stringsAsFactors = FALSE
      )
    }
  } else {
    if (is.null(gtf_path) || !is.character(gtf_path) || length(gtf_path) != 1L ||
        !nzchar(gtf_path) || !file.exists(gtf_path)) {
      stop(.stop("gtf_path must point to an existing GTF/GFF file."))
    }
    tbl <- NULL
    if (requireNamespace("rtracklayer", quietly = TRUE)) {
      gr <- tryCatch(rtracklayer::import(gtf_path), error = function(e) NULL)
      if (!is.null(gr)) {
        mc <- as.data.frame(S4Vectors::mcols(gr), stringsAsFactors = FALSE)
        if (all(c("transcript_id", "gene_id") %in% colnames(mc))) {
          tbl <- data.frame(
            transcript_id = as.character(mc$transcript_id),
            gene_id = as.character(mc$gene_id),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    if (is.null(tbl)) {
      lines <- readLines(gtf_path, warn = FALSE)
      lines <- lines[!grepl("^#", lines) & nzchar(trimws(lines))]
      has_both <- grepl('transcript_id\\s+"[^"]+"', lines, perl = TRUE) &
                  grepl('gene_id\\s+"[^"]+"', lines, perl = TRUE)
      if (any(has_both)) {
        sub_l <- lines[has_both]
        tx_v <- sub('^.*transcript_id\\s+"([^"]+)".*$', "\\1", sub_l, perl = TRUE)
        gn_v <- sub('^.*gene_id\\s+"([^"]+)".*$', "\\1", sub_l, perl = TRUE)
        tbl <- data.frame(transcript_id = tx_v, gene_id = gn_v, stringsAsFactors = FALSE)
      } else {
        tbl <- data.frame(transcript_id = character(0L), gene_id = character(0L), stringsAsFactors = FALSE)
      }
    }
    tbl <- tbl[tbl$transcript_id %in% tx_ids, , drop = FALSE]
  }

  keep <- !is.na(tbl$transcript_id) & nzchar(tbl$transcript_id) &
          !is.na(tbl$gene_id) & nzchar(tbl$gene_id)
  tbl <- tbl[keep, , drop = FALSE]
  tbl <- tbl[!duplicated(tbl$transcript_id), , drop = FALSE]
  rownames(tbl) <- NULL

  unmapped <- setdiff(tx_ids, tbl$transcript_id)
  list(
    table = tbl,
    n_input = length(tx_ids),
    n_mapped = nrow(tbl),
    n_unmapped = length(unmapped),
    unmapped_ids = unmapped,
    source = effective_source
  )
}
