#!/usr/bin/env Rscript
# =============================================================================
# tools/docs_hash_manifest.R — garde `docs_hash` (ACT-05) : intégrité de l'arbre
#                             docs/ (gitigné par décision opérateur 2026-10-10)
# =============================================================================
# POURQUOI CE FICHIER EXISTE
# --------------------------
# `docs/` est HORS DE TOUTE PROTECTION GIT (.gitignore:54, décision 2026-10-10 —
# ne jamais compenser par `git add -f`, rappelé en fin de STATUS.md). Un arbre de
# référence entier vit donc sans historique : personne ne peut répondre à « ce
# fichier a-t-il été retouché depuis hier ? ». Le manifeste est cette réponse :
# un sha256 par fichier, écrit une fois, vérifié ensuite. Ce garde est celui qui
# « aurait attrapé les 8 dérives » mesurées par le health check 2026-11-03
# (bannières retournées après coup, pins périmés découverts tard).
#
# ORIGINE — RECONSTRUCTION (comme run_all_gates.R)
# ------------------------------------------------
# L'original (lot 6-7, 2026-10-12) a disparu de l'arbre avec le reste du lot ;
# l'ARTEFACT lui survit (`docs/audits/docs_hash_manifest.txt`, 2026-10-12,
# 34348 o) et fait foi : en-tête, format `<sha256>  <chemin relatif>`, tri par
# chemin (ordre d'octets), auto-exclusion du manifeste lui-même. Le contrat de
# périmètre vient du health check §5 (ACT-05) : `docs/**/*.{md,txt,json}`
# (casse insensible) + `AGENTS.md` racine, INCLUS `docs/archive/` (gelé —
# hashes stables = anti-retouche gratuite), EXCLUS `docs/repomix/` (23 f.,
# régénérables — dérive forcée à chaque commit) et `docs/audits/visual/`
# (sorties d'audit réécrites) en PLUS de l'exclusion historique
# `docs/audits/gates_history/` (artefacts machine).
#
# ⚠️ CONTRAT DE SORTIE AVEC run_all_gates.R. Le meta-gate greffe ce script comme
# garde `docs_hash` (kind="exit") et extrait son `details` par regex
# `(\d+)\s*erreur\(s\)` / `(\d+)\s*avertissement`. L'original portait
# `details: ""` dans gates_history/latest.json — le reproduire oblige à NE PAS
# imprimer ces motifs : les comptes s'écrivent « dérives », « MODIFIÉS »,
# « MANQUANTS », « HORS MANIFESTE ». Exit 0 = VERT, 1 = ROUGE, 2 = usage/infra.
#
# USAGE
#   Rscript tools/docs_hash_manifest.R --build    # (re)génère le manifeste
#   Rscript tools/docs_hash_manifest.R --check    # vérifie (garde, exit 0/1)
#   Rscript tools/docs_hash_manifest.R            # = --check (défaut : le
#                                                 #   meta-gate lance sans arg)
#   Rscript tools/docs_hash_manifest.R --root=DIR # pour les tests
#
# RÈGLE D'USAGE (health check §B4) : `--build` après chaque passe doc délibérée,
# `--check` à chaque clôture de session (et via run_all_gates.R). Un --check
# ROUGE sans modification délibérée = retouche non tracée : arrêter et expliquer.
#
# ⚠️ UN MANIFESTE ABSENT EST UN ÉCHEC, PAS UN Avertissement. `--check` sans
# manifeste rend ROUGE (exit 1) : c'est exactement le trou du 2026-10-12 au
# 2026-11-09 (script perdu, aucun garde ne rougissait) qu'on ne veut plus.
# =============================================================================

# ── Racine du projet : déduite du script, jamais du répertoire courant --------
.ts_docs_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (!length(hit)) return(getwd())
  script <- normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE)
  normalizePath(file.path(dirname(script), ".."), mustWork = FALSE)
}

.ts_docs_args <- function(argv) {
  opt <- list(mode = NULL, root = NULL)
  for (a in argv) {
    if (identical(a, "--build")) { opt$mode <- "build"; next }
    if (identical(a, "--check")) { opt$mode <- "check"; next }
    if (grepl("^--root=", a)) { opt$root <- sub("^--root=", "", a); next }
    stop(sprintf("argument inconnu '%s' (usage : --build | --check | --root=DIR)", a),
         call. = FALSE)
  }
  # Sans mode explicite → --check. Le meta-gate (run_all_gates.R, garde
  # `docs_hash`) lance ce script SANS argument : le défaut DOIT être le mode
  # garde, sinon chaque run du meta-gate rendrait ROUGE exit 2 (mesuré
  # 2026-11-09 avant ce défaut). `--build` reste une action délibérée, jamais
  # un défaut.
  if (is.null(opt$mode)) opt$mode <- "check"
  opt
}

# ── Périmètre (contrat health check §5 / artefact survivant) ------------------
.ts_docs_manifest_path <- function(root) {
  file.path(root, "docs", "audits", "docs_hash_manifest.txt")
}

# Répertoires exclus EN TIER (les 25-93 fichiers qui se réécrivent tout seuls :
# les hacher n'apporterait que des dérives permanentes, donc un garde muet).
# Préfixes RELATIFS À docs/ — `list.files(recursive = TRUE)` retourne des
# chemins relatifs au répertoire scanné, séparés par `/` (normalisés ensuite).
.ts_docs_excluded_prefixes <- c(
  "repomix/",                  # régénérables (O25/G5 — health check §5)
  "audits/gates_history/",     # artefacts machine
  "audits/visual/"             # sorties d'audit réécrites — health check §5
)

.ts_docs_scan <- function(root) {
  docs_dir <- file.path(root, "docs")
  all <- list.files(docs_dir, recursive = TRUE, include.dirs = FALSE,
                    no.. = TRUE, all.files = FALSE)
  rel_docs <- gsub("\\\\", "/", all, perl = TRUE)
  # Casse insensible (contrat d'origine) : on filtre l'extension en minuscule
  # pour être portable (le scan retourne l'extension en l'état).
  keep <- grepl("\\.(md|txt|json)$", tolower(rel_docs), perl = TRUE)
  rel_docs <- rel_docs[keep]
  drop <- vapply(rel_docs, function(p) any(startsWith(p, .ts_docs_excluded_prefixes)),
                 logical(1))
  rel_docs <- rel_docs[!drop]
  # Chemins RELATIFS À LA RACINE du projet (format du manifeste d'origine),
  # + AGENTS.md racine (contrat : haché au même titre que docs/ — retombe si
  # absent, un scan sans AGENTS.md ne doit pas produire de ligne fantôme).
  out <- c(paste0("docs/", rel_docs), "AGENTS.md")
  out <- out[file.exists(file.path(root, out))]
  # Auto-exclusion : le manifeste ne se hache pas lui-même (l'auto-hachage est
  # sans objet et forcerait une dérive à chaque exécution — contrat d'origine).
  out <- out[out != "docs/audits/docs_hash_manifest.txt"]
  # Tri par chemin, ORDRE D'OCTETS (méthode radix = locale C ; le tri par
  # défaut dépend de la locale et rendrait deux runs non comparables).
  sort(out, method = "radix")
}

.ts_docs_hash_file <- function(path) {
  digest::digest(file = path, algo = "sha256")
}

.ts_docs_hashes <- function(root, rel) {
  vapply(rel, function(p) .ts_docs_hash_file(file.path(root, p)), character(1),
         USE.NAMES = FALSE)
}

# ── En-tête du manifeste (mot de l'artefact survivant, étendu au contrat §5) --
.ts_docs_header <- function(n) {
  c(
    "# docs_hash_manifest.txt — manifeste d'intégrité de l'arbre docs/ (gitigné par décision opérateur 2026-10-10).",
    "# Généré par tools/docs_hash_manifest.R — NE PAS ÉDITER À LA MAIN, régénérer :",
    "#   Rscript tools/docs_hash_manifest.R --build",
    sprintf("# Généré le (UTC)  : %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
    "# Algorithme       : sha256 (digest::digest, file=)",
    sprintf("# Fichiers hachés  : %d", n),
    "# Périmètre        : docs/**/*.{md,txt,json} (casse insensible) + AGENTS.md (racine du projet).",
    "# Inclus           : docs/archive/ (gelé — hashes stables = anti-retouche gratuite).",
    "# Auto-exclusion   : ce manifeste S'EXCLUT de son propre balayage (l'auto-hachage",
    "#                    est sans objet et forcerait une dérive à chaque exécution).",
    "# Exclusions aussi : docs/repomix/ (régénérables, dérive à chaque commit — health check §5)",
    "#                    docs/audits/gates_history/ (artefacts machine réécrits à chaque run)",
    "#                    docs/audits/visual/ (sorties d'audit réécrites — health check §5).",
    "# Format           : <sha256>  <chemin relatif à la racine du projet> — tri par chemin (ordre d'octets).",
    "#",
    "# Reconstruction B4 DOC_TIDY (2026-11-09) — l'original (lot 6-7, 2026-10-12) est perdu ;",
    "# l'artefact survivant et le health check §5 (ACT-05) font foi du contrat."
  )
}

.ts_docs_write <- function(root, rel, hashes) {
  lines <- c(.ts_docs_header(length(rel)),
             sprintf("%s  %s", hashes, rel))
  con <- file(.ts_docs_manifest_path(root), open = "wb")
  writeLines(lines, con, sep = "\n", useBytes = TRUE)
  close(con)
}

# ── --check : recompute et compare. Rend les 3 classes de dérive séparées ----
.ts_docs_check <- function(root) {
  man <- .ts_docs_manifest_path(root)
  if (!file.exists(man)) {
    message("ROUGE: manifeste absent : ", man)
    message("        (jamais construit ou perdu — lancer --build APRÈS une passe doc délibérée ;",
            " un manifeste absent doit rougir, jamais passer silencieusement)")
    return(1L)
  }
  ln <- readLines(man, warn = FALSE, encoding = "UTF-8")
  ln <- ln[!startsWith(ln, "#") & nzchar(trimws(ln))]
  parts <- strsplit(ln, "  ", fixed = TRUE)
  ok_fmt <- vapply(parts, function(p) length(p) == 2L, logical(1))
  if (!all(ok_fmt)) {
    message("ROUGE: ", sum(!ok_fmt), " ligne(s) de manifeste mal formée(s) (attendu '<sha256>  <chemin>')")
    return(1L)
  }
  man_hash <- vapply(parts, `[`, character(1), 1)
  man_path <- vapply(parts, `[`, character(1), 2)
  names(man_hash) <- man_path

  rel <- .ts_docs_scan(root)
  hashes <- .ts_docs_hashes(root, rel)
  names(hashes) <- rel

  modified <- character(0)
  for (p in intersect(names(man_hash), names(hashes))) {
    if (!identical(man_hash[[p]], hashes[[p]])) modified <- c(modified, p)
  }
  missing   <- setdiff(names(man_hash), names(hashes))  # au manifeste, absent de l'arbre
  unlisted  <- setdiff(names(hashes), names(man_hash))  # dans le périmètre, pas au manifeste

  nd <- length(modified) + length(missing) + length(unlisted)
  if (nd == 0L) {
    cat(sprintf("docs_hash — vérification : %d fichiers recalculés · 0 dérive · VERT\n",
                length(rel)))
    return(0L)
  }
  cat(sprintf("docs_hash — vérification : %d fichiers recalculés · %d dérive(s) · ROUGE\n",
              length(rel), nd))
  show <- function(label, v) {
    if (length(v)) {
      cat(sprintf("  %s (%d) :\n", label, length(v)))
      for (p in head(v, 20)) cat("    -", p, "\n")
      if (length(v) > 20) cat(sprintf("    … + %d autre(s)\n", length(v) - 20))
    }
  }
  show("MODIFIÉS (hash différent)", modified)
  show("MANQUANTS (au manifeste, absents de l'arbre)", missing)
  show("HORS MANIFESTE (dans le périmètre, non hachés)", unlisted)
  if (length(unlisted)) {
    cat("  → une passe doc délibérée ? lancer --build ; sinon, retouche non tracée à expliquer.\n")
  }
  1L
}

# ── Main --------------------------------------------------------------------
.ts_docs_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opt <- tryCatch(.ts_docs_args(argv), error = function(e) e)
  if (inherits(opt, "condition")) {
    message("ERREUR: ", conditionMessage(opt))
    return(invisible(2L))
  }
  root <- if (is.null(opt$root)) .ts_docs_root() else normalizePath(opt$root, mustWork = FALSE)
  if (!dir.exists(root)) {
    message("ERREUR: racine introuvable : ", root)
    return(invisible(2L))
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    message("ERREUR: le paquet digest est requis (sha256) — il est au renv.lock,",
            " lancer renv::restore() ou install.packages('digest')")
    return(invisible(2L))
  }
  if (!dir.exists(file.path(root, "docs"))) {
    message("ERREUR: pas d'arbre docs/ sous ", root)
    return(invisible(2L))
  }
  if (identical(opt$mode, "build")) {
    rel <- .ts_docs_scan(root)
    hashes <- .ts_docs_hashes(root, rel)
    .ts_docs_write(root, rel, hashes)
    cat(sprintf("docs_hash — manifeste écrit : %d fichiers hachés (%s)\n",
                length(rel), .ts_docs_manifest_path(root)))
    return(invisible(0L))
  }
  invisible(.ts_docs_check(root))
}

if (sys.nframe() == 0L) {
  quit(save = "no", status = .ts_docs_main())
}
