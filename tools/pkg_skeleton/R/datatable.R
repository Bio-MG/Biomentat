# =============================================================================
# R/plotting/datatable.R — PLOT-S3 : wrapper DT harmonise (ts_datatable)
# =============================================================================
# Un seul endroit pour construire les tables de resultats de l'application.
#
# MESURE PREALABLE (2026-09-11) — contrairement a PLOT-S2, IL Y A UN PIEGE :
#   pageLength : 10 (x22), 15 (x15), 8 (x5), 20 (x3), 6 (x1)
#                => AUCUNE valeur par defaut neutre (meme situation que le
#                   base_size de PLOT-S1). C'est pourquoi `page_length` est
#                   OBLIGATOIRE : le helper refuse de choisir a la place de
#                   l'appelant.
#   boutons    : 1 seul site sur 49 declare extensions = "Buttons" + dom =
#                "Bfrtip" (la table pathways). Les 48 autres n'en ont aucun.
#                => `buttons = FALSE` par defaut : activer les boutons est un
#                   OPT-IN explicite, jamais un effet de bord de la migration.
#                [MISE A JOUR jalon DT-EXPORT 2026-09-14] : le defaut est
#                   desormais TS_DT_BUTTONS_DEFAULT = TRUE (contrat §6, option
#                   B) — les tables de resultats sont boutonnees ; les tables
#                   d'APERCU passent buttons = FALSE explicitement.
#
# REGLE : ZERO CHANGEMENT DE COMPORTEMENT. Chaque site migre garde ses valeurs
#         actuelles ecrites explicitement.
#
# Pure R : testable hors de Shiny, erreurs classees, messages FR, call. = FALSE.
# =============================================================================
# Cree le 2026-09-11 (PLOT-S3). Toute modification doit passer SIMULTANEMENT
# par : ce fichier, tests/testthat/test-plot-datatable.R et
# docs/contracts/PLOT_DATATABLE_CONTRACT.md.
# =============================================================================
# PC (2026-11-09) — canal future opt-in du render DT (jalon « app
# non-bloquante ») : le filtrage server-side des tables volumineuses peut
# quitter le thread principal via DT::renderDT(future = TRUE) — DT 0.34
# route alors chaque requete AJAX de filtre sur promises::future_promise(),
# c.-a-d. sur le plan future courant du process. DEUX primitives vivent ICI
# (point de passage unique du jalon PC, fichier exclusif de la lane) :
#   - ts_dt_future()        : resolveur OPT-IN (option/env TS_DT_FUTURE,
#                             default FALSE = comportement stock de DT).
#   - ts_dt_future_install(): pose le defaut du formal `future` de
#     DT::renderDT sur ce resolveur (evaluate A CHAQUE appel de renderDT,
#     donc activation a chaud possible). Idempotent, garde sur la presence
#     du formal (DT >= 0.33) et de future.mirai. Sans cet appel, DT reste
#     100% stock.
# MESURE (2026-11-09, sondes /tmp/pc_recon) : le pont mirai<->future exige
# le paquet CRAN future.mirai (mirai 2.7.2 n'expose aucun backend future ;
# future::plan(mirai::mirai) -> « node stack overflow »). GO operateur du
# 2026-11-09 : future.mirai 1.1.0 ajoute au lock (unique ecart). Le backend
# reel est plan(mirai_cluster) (future.mirai) pose par init_spatial_daemons()
# quand TS_DT_FUTURE est actif — PIDs mesures : reutilisation du pool
# existant, 0 worker nouveau. Sans pool init, le plan de boot F2
# (multisession, global.R:558) sert de repli : 0 worker nouveau non plus
# (workers deja spawnes par le boot). Le formal installe evalue l'OPTION a
# l'appel : la desactivation (options(TS_DT_FUTURE = FALSE)) suffit, sans
# re-install.
# =============================================================================

# --- Constantes ---------------------------------------------------------------

.ts_dt_const <- function(name, fallback) {
  if (exists(name, envir = globalenv(), inherits = TRUE)) {
    get(name, envir = globalenv(), inherits = TRUE)
  } else {
    fallback
  }
}

#' Tailles de page effectivement presentes dans l'application.
#' Il n'y en a AUCUNE par defaut — voir l'en-tete du fichier.
ts_datatable_page_lengths <- function() {
  .ts_dt_const("TS_DT_PAGE_LENGTHS", c(6L, 8L, 10L, 15L, 20L))
}

#' Jeu de boutons d'export propose (utilise uniquement si `filename_base` est
#' fourni ; sans nom de fichier, DT applique ses boutons par defaut).
ts_datatable_buttons <- function(filename_base) {
  list("copy", "print",
       list(extend = "csv",   filename = filename_base),
       list(extend = "excel", filename = filename_base))
}

#' Surface publique figee (verifiee par le test de gel).
ts_datatable_public_api <- function() {
  c("ts_datatable", "ts_datatable_buttons", "ts_datatable_page_lengths",
    "ts_datatable_public_api", "ts_dt_future", "ts_dt_future_install",
    "ts_dt_future_restore")
}

# --- Erreurs classees --------------------------------------------------------

.ts_dt_stop <- function(state, message) {
  stop(errorCondition(
    message,
    state = state,
    class = c("plot_datatable_error", "error", "condition")
  ))
}

# --- Helper principal --------------------------------------------------------

#' Table de resultats harmonisee.
#'
#' Encapsule `DT::datatable()`. **Aucun changement de rendu** par rapport aux
#' appels directs existants, a condition de passer explicitement `page_length`.
#'
#' @param df data.frame (ou coercible) a afficher.
#' @param page_length Nombre de lignes par page. **Obligatoire** : aucune
#'   valeur par defaut n'est neutre (10, 15, 8, 20 et 6 coexistent dans l'app).
#' @param filename_base NULL = aucun bouton d'export nomme. Si fourni, ajoute
#'   les boutons copy / print / csv / excel avec ce nom de fichier.
#' @param filter Emplacement du filtre ("top" / "none" / "bottom").
#' @param rownames Afficher les noms de lignes.
#' @param scroll_x Defilement horizontal. **Defaut = FALSE depuis le 2H (O28,
#'   2026-10-30)** : arbitrage du constat 2H (packet round-7) — la COMBINAISON
#'   `scrollX = TRUE` + `filter = "top"` produisait le bug DataTables du
#'   double header / double filtre (le clonage du conteneur `thead` par le
#'   defilement horizontal duplique les boîtes de sélection du header flottant,
#'   `tools/pkg_skeleton/R/datatable.R:105` à HEAD `88f9a43` ; spéc §2.19).
#'   Les deux correctifs JS/CSS de la spéc étant invérifiables hors navigateur,
#'   l'arbitrage retire la cause du glitch : le filtre par colonne — la
#'   fonctionnalité porteuse d'ergonomie sur les tables de résultats — reste
#'   ACTIF par defaut (`filter = "top"`), le scrollX devient un OPT-IN explicite
#'   (`scroll_x = TRUE`) pour les tables réellement larges.
#' @param buttons Ajouter l'extension Buttons + `dom = "Bfrtip"`.
#'   **Defaut = `TS_DT_BUTTONS_DEFAULT` (TRUE depuis le jalon DT-EXPORT,
#'   2026-09-14 — contrat §6 option B)** : les tables de resultats sont
#'   boutonnees par defaut. Les tables d'APERCU (previews `dom = "t"/"tip"`,
#'   petites tables QC) passent `buttons = FALSE` explicitement.
#' @param dom,extensions Surcharges manuelles (prioritaires sur `buttons`).
#' @param extra_options Liste nommee fusionnee dans `options` APRES tout le
#'   reste (priorite maximale) — canal unique pour les options DT que le
#'   wrapper n'expose pas (`language`, `lengthMenu`, `order`, `autoWidth`...).
#'   Ne jamais passer `options =` via `...` : collision avec l'option interne.
#' @param ... Transmis a `DT::datatable()` (selection, caption, ...).
#' @return Un objet `DT::datatable()` (htmlwidget).
#'
#' @examples
#' \dontrun{
#' ts_datatable(df, page_length = 15)
#' ts_datatable(df, page_length = 10, buttons = TRUE, filename_base = "pathways")
#' }
ts_datatable <- function(df, page_length, filename_base = NULL,
                         filter = "top", rownames = FALSE,
                         scroll_x = FALSE,   # 2H (O28) : était TRUE — voir @param scroll_x
                         buttons = .ts_dt_const("TS_DT_BUTTONS_DEFAULT", TRUE),
                         dom = NULL, extensions = NULL, extra_options = NULL,
                         ...) {
  # --- Gardes ---
  if (missing(page_length) || is.null(page_length) ||
      length(page_length) != 1L || is.na(page_length) ||
      !is.numeric(page_length) || page_length < 1) {
    .ts_dt_stop("invalid_page_length",
                paste0("ts_datatable() : page_length est obligatoire (entier >= 1). ",
                       "Aucune valeur par defaut n'est neutre : l'application ",
                       "utilise aujourd'hui 10, 15, 8, 20 et 6 selon les tables."))
  }
  if (!is.data.frame(df)) {
    if (is.matrix(df) || is.list(df)) {
      df <- as.data.frame(df, check.names = FALSE)
    } else {
      .ts_dt_stop("invalid_data",
                  "ts_datatable() : donnees invalides (data.frame, matrix ou liste attendus).")
    }
  }
  if (!is.logical(buttons) || length(buttons) != 1L || is.na(buttons)) {
    .ts_dt_stop("invalid_buttons",
                "ts_datatable() : buttons doit etre TRUE ou FALSE.")
  }
  if (!is.null(extra_options) &&
      (!is.list(extra_options) || is.null(names(extra_options)) ||
       any(!nzchar(names(extra_options))))) {
    .ts_dt_stop("invalid_extra_options",
                paste0("ts_datatable() : extra_options doit etre une liste ",
                       "nommee (options DT fusionnees dans options)."))
  }

  # --- Options : on ne pose `dom` QUE s'il est demande (zero changement) ---
  # page_length est transmis TEL QUEL (pas de coercion) : les sites historiques
  # ecrivent pageLength = 15 (double) et le JSON produit est identique, mais on
  # ne prend aucun risque de changer la valeur serialisee.
  opts <- list(pageLength = page_length, scrollX = scroll_x)

  ext <- extensions
  if (isTRUE(buttons)) {
    if (is.null(ext)) ext <- "Buttons"
    if (is.null(dom))  dom  <- "Bfrtip"
  }
  if (!is.null(dom)) opts$dom <- dom
  if (!is.null(filename_base)) opts$buttons <- ts_datatable_buttons(filename_base)
  # Priorite maximale : extra_options ecrase les cles posees ci-dessus.
  if (!is.null(extra_options)) opts <- utils::modifyList(opts, extra_options)

  # --- Appel ---
  args <- list(data = df, filter = filter, rownames = rownames,
               options = opts, ...)
  if (!is.null(ext)) args$extensions <- ext

  do.call(DT::datatable, args)
}

# --- PC (2026-11-09) : canal future opt-in du render DT ----------------------

#' Resolveur du canal future des tables DT (jalon PC).
#'
#' Le filtrage server-side DT est OPT-IN : sans activation explicite, le
#' comportement stock de `DT::renderDT` est preserve (filtrage synchrone
#' sur le thread principal — zero changement). Activation : option
#' `TS_DT_FUTURE` (logique ou chaine "true"/"TRUE"/"1") ou variable
#' d'environnement du meme nom.
#'
#' Le backend effectif est le plan future courant du process :
#' `plan(mirai_cluster)` (pool mirai existant, pose par
#' `init_spatial_daemons()` quand ce canal est actif) ou, a defaut de pool
#' init, le plan de boot F2 (`multisession`, global.R:558). Dans les deux
#' cas : 0 worker nouveau (mesure PIDs 2026-11-09).
#'
#' @return Un logique scalaire : TRUE = canal demande.
#' @export
ts_dt_future <- function() {
  v <- getOption("TS_DT_FUTURE", NULL)
  if (!is.null(v)) return(.ts_dt_truthy(v))
  .ts_dt_truthy(Sys.getenv("TS_DT_FUTURE", unset = ""))
}

#' Predicat truthy des reglages opt-in (option/env TS_*).
#'
#' MESURE (2026-11-09) : as.logical("1") vaut NA (les chaines numeriques ne
#' se coercent PAS en logique) — le predicat accepte explicitement les
#' formes logique, "TRUE"/"T"/"1"/"YES" (insensibles a la casse) et refuse
#' tout le reste ("false", "", NA, autres).
.ts_dt_truthy <- function(v) {
  if (is.null(v) || length(v) != 1L || is.na(v)) return(FALSE)
  if (is.logical(v)) return(isTRUE(v))
  if (is.character(v)) {
    toupper(trimws(v)) %in% c("TRUE", "T", "1", "YES")
  } else FALSE
}

#' Installe le defaut `future` du render DT sur le resolveur (jalon PC).
#'
#' Pose `formals(DT::renderDT)$future <- quote(getOption("TS_DT_FUTURE",
#' FALSE))` : l'expression est evaluee dans le namespace DT A CHAQUE appel
#' de `DT::renderDT(...)`, donc (a) le defaut reste FALSE (stock) tant que
#' l'option n'est pas posee, (b) l'activation est possible a chaud (chaque
#' re-enregistrement d'output la lit), (c) `getOption()` est base — aucune
#' dependance a l'environnement d'evaluation.
#'
#' Gardees : DT present, formal `future` existant (DT >= 0.33), paquet
#' `future.mirai` present (pont mesure — voir l'en-tete du fichier). Sans
#' ces gardes, rien n'est modifie et la fonction retourne FALSE : l'app
#' tourne en DT stock.
#'
#' Implantation mesuree (2026-11-09) : `formals(DT::renderDT)$future <- v`
#' (assignment complexe via `::`) echoue avec « object 'DT' not found » —
#' la seule voie fiable est la copie locale du closure + assignInNamespace
#' (restauration symetrique dans les tests). Idempotent : re-executer ne
#' fait rien si le defaut est deja pose. Appelle par
#' `init_spatial_daemons()` (spatial_async.R) apres verification du pool —
#' le seul appel de production ; les tests appellent directement et
#' restaurent le formal precedent.
#'
#' @return invisible(TRUE) si le defaut est (ou devient) pose, sinon
#'   invisible(FALSE) — jamais d'erreur.
ts_dt_future_install <- function() {
  if (!requireNamespace("DT", quietly = TRUE)) return(invisible(FALSE))
  if (!requireNamespace("future.mirai", quietly = TRUE)) return(invisible(FALSE))
  ns <- asNamespace("DT")
  fm <- get("renderDT", envir = ns)
  if (!"future" %in% names(formals(fm))) return(invisible(FALSE))
  target <- quote(getOption("TS_DT_FUTURE", FALSE))
  if (identical(formals(fm)$future, target)) return(invisible(TRUE))
  formals(fm)$future <- target
  assignInNamespace("renderDT", fm, ns = "DT")
  invisible(TRUE)
}

#' Restaure le formal `future` stock de DT::renderDT (compagnon de test).
#'
#' Uniquement pour les tests et les teardowns : pose le formal `future`
#' a `FALSE` (defaut DT stock) via la meme voie que ts_dt_future_install().
#' @return invisible(TRUE) si DT est present, invisible(FALSE) sinon.
ts_dt_future_restore <- function() {
  if (!requireNamespace("DT", quietly = TRUE)) return(invisible(FALSE))
  ns <- asNamespace("DT")
  fm <- get("renderDT", envir = ns)
  if (identical(formals(fm)$future, FALSE)) return(invisible(TRUE))
  formals(fm)$future <- FALSE
  assignInNamespace("renderDT", fm, ns = "DT")
  invisible(TRUE)
}
