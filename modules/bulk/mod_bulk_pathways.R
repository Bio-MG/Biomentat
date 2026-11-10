# mod_bulk_pathways.R — Bulk Child 3 (i18n Phase 3.1)

mod_bulk_pathways_ui <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("enrich_mode"), i18n$t("M\u00e9thode d'enrichissement"),
                choices = stats::setNames(c("ora","gsea","dir","compare"),
                  c(.tr_plain("ORA \u2014 sur g\u00e8nes significatifs (classique)"),
                    .tr_plain("GSEA \u2014 sur tous les g\u00e8nes class\u00e9s (sans seuil)"),
                    .tr_plain("Directionnelle \u2014 goana/kegga (table DE compl\u00e8te, Up/Down par voie)"),
                    .tr_plain("Comparative — compareCluster (inter-contrastes)"))),
                selected = "ora"),
    div(class = "small text-muted mb-2",
        i18n$t("GSEA n'a pas besoin de seuil de significativit\u00e9 arbitraire \u2014 elle classe tous les g\u00e8nes par Log2FC et teste l'enrichissement cumul\u00e9. Plus robuste statistiquement, recommand\u00e9e si peu de g\u00e8nes passent vos seuils.")),

    conditionalPanel(
      condition = "input.enrich_mode == 'ora'", ns = ns,
      selectInput(ns("pathway_source"), i18n$t("Source de g\u00e8nes"),
                  choices = stats::setNames(c("up","down","all_sig","manual"),
                    c(.tr_plain("G\u00e8nes Up (significatifs)"), .tr_plain("G\u00e8nes Down (significatifs)"),
                      .tr_plain("Tous g\u00e8nes significatifs"), .tr_plain("S\u00e9lection manuelle")))),
      conditionalPanel(
        condition = "input.pathway_source == 'manual'", ns = ns,
        selectizeInput(ns("pathway_genes"), i18n$t("G\u00e8nes"), choices = NULL, multiple = TRUE)
      )
    ),

    # PB (GEMS D5) — ORA directionnelle : seuil FDR de la sélection DE (BH
    # interne de goana/kegga sur la p-valeur BRUTE ; le padj du module est le
    # même BH — la cohérence est mesurée au test du moteur).
    conditionalPanel(
      condition = "input.enrich_mode == 'dir'", ns = ns,
      numericInput(ns("dir_fdr"), i18n$t("Seuil FDR (s\u00e9lection DE directionnelle)"),
                   value = 0.05, min = 0.001, max = 0.5, step = 0.01)
    ),

    fluidRow(
      column(6, selectInput(ns("pathway_db"), i18n$t("Base de donn\u00e9es"),
                            choices = c("GO Biological Process" = "GOBP",
                                        "KEGG Pathways"         = "KEGG",
                                        "Reactome"              = "Reactome",
                                        "Packs (GMT libraries)" = "Packs"))),
      column(6, selectInput(ns("pathway_org"), i18n$t("Organisme"),
                            choices = stats::setNames(c("human","mouse"),
                              c(.tr_plain("Humain"), .tr_plain("Souris")))))
    ),
    # PK-2 (pilier 2 local, fiche packs §3/§4.3) : multi-sélection des
    # librairies GMT installées — le sélecteur n'apparaît que sur « Packs » ;
    # les librairies viennent du registre session (scan PK-1 au démarrage).
    uiOutput(ns("packs_libs_ui")),
    # PB v2 (GEMS D5, jalon 2) — compareCluster inter-contrastes : listes
    # = DE results par contraste (shared_rv$contrasts), fond = rownames de la
    # matrice filtrée (M12, un seul fond déclaré).
    conditionalPanel(
      condition = "input.enrich_mode == 'compare'", ns = ns,
      selectInput(ns("comp_source"), i18n$t("Gènes par contraste"),
                  choices = stats::setNames(c("all","up","down"),
                    c(.tr_plain("Tous significatifs"), .tr_plain("Up seulement"), .tr_plain("Down seulement")))),
      numericInput(ns("comp_fdr"), i18n$t("Seuil FDR (sélection DE comparative)"),
                   value = 0.05, min = 0.001, max = 0.5, step = 0.01),
      div(class = "small text-muted",
          i18n$t("Compare les enrichissements des contrastes chargés (au moins 2 contrastes avec des gènes significatifs au seuil requis)."))
    ),
    numericInput(ns("pathway_pval"), i18n$t("P-value cutoff"), value = 0.05, min = 0.001, max = 0.1, step = 0.01),

    # STAT-Q1 — correction pour tests multiples de l'enrichissement. Choix
    # déclarés dans config/defaults.R (TS_PADJ_METHODS) ; "fdr" est
    # volontairement absent (alias de BH dans stats::p.adjust.methods).
    selectInput(ns("pathway_padj_method"), i18n$t("M\u00e9thode de correction (p-adj)"),
                choices  = stats::setNames(TS_PADJ_METHODS, TS_PADJ_METHODS),
                selected = TS_PADJ_METHOD_DEFAULT),

    actionButton(ns("run_pathway"), i18n$t("Lancer Enrichissement"),
                 class = "btn-warning w-100", icon = icon("dna")),

    downloadButton(ns("dl_pathway"), i18n$t("Export CSV"), class = "btn-sm btn-info w-100 mt-2"),

    # STAT-Q4 — meme table, en Excel. Le CSV ci-dessus reste inchange ; le repli
    # CSV du helper n'est utilise que si openxlsx est absent.
    downloadButton(ns("dl_pathway_excel"), i18n$t("Export Excel"), class = "btn-sm btn-success w-100 mt-2"),

    div(class = "small text-muted mt-1", textOutput(ns("pathway_status"))),

    # ── Bulk V2 M2 — scores de voies PAR ÉCHANTILLON (GSVA/ssGSEA/PLAGE/zscore).
    # Complément de l'ORA/GSEA ci-dessus : attribue un score par voie à CHAQUE
    # échantillon (aucun contraste requis). La matrice (voies x échantillons)
    # vit dans bulk_obj$pathways$per_sample + shared_rv$pathway_scores.
    hr(),
    h6(i18n$t("Scores par \u00e9chantillon (GSVA / ssGSEA)"), style = "font-weight:bold;"),
    div(class = "small text-muted mb-2",
        i18n$t("Attribue \u00e0 chaque \u00e9chantillon un score par voie \u2014 PCA et heatmaps par voie m\u00eame sans contraste. Exige la matrice VST (\u00e9tape 1). Les jeux de g\u00e8nes viennent d'une source NATIVE embarqu\u00e9e (MSigDB, PROGENy, DoRothEA \u2014 aucun fichier \u00e0 fournir) ou d'un .gmt fourni. Les jeux dont moins de 20 % des g\u00e8nes sont retrouv\u00e9s sont rejet\u00e9s (d\u00e9calage d'identifiants).")),
    # PLOT-S6b — source NATIVE des jeux de gènes (bulk_gene_sets.R). Le .gmt
    # reste disponible pour les jeux maison : le fileInput n'apparaît que pour
    # cette source (conditionalPanel), il n'est plus un prérequis.
    selectInput(ns("scores_source"), i18n$t("Source de jeux de g\u00e8nes"),
                choices  = bulk_gene_set_choices(.tr_plain),
                selected = "msigdb_hallmark"),
    selectInput(ns("scores_org"), i18n$t("Organisme"),
                choices = stats::setNames(c("human", "mouse"),
                  c(.tr_plain("Humain"), .tr_plain("Souris"))),
                selected = "human"),
    conditionalPanel(
      condition = "input.scores_source == 'file'", ns = ns,
      fileInput(ns("scores_gmt"), i18n$t("Fichier de jeux de g\u00e8nes (.gmt)"),
                accept = c(".gmt", ".txt"), buttonLabel = i18n$t("Parcourir..."))
    ),
    div(class = "small text-muted mb-2", textOutput(ns("scores_source_hint"))),
    selectInput(ns("scores_method"), i18n$t("M\u00e9thode de scoring"),
                choices = stats::setNames(c("ssgsea", "gsva", "plage", "zscore"),
                  c("ssGSEA", "GSVA", "PLAGE", "zscore")),
                selected = "ssgsea"),
    fluidRow(
      column(6, numericInput(ns("scores_min_size"), i18n$t("Taille min voie"),
                            value = TS_BULK_GSVA_MIN_SIZE, min = 1, step = 1)),
      column(6, numericInput(ns("scores_max_size"), i18n$t("Taille max voie"),
                            value = TS_BULK_GSVA_MAX_SIZE, min = 2, step = 1))
    ),
    actionButton(ns("run_scores"), i18n$t("Lancer Scores par \u00e9chantillon"),
                 class = "btn-warning w-100", icon = icon("layer-group")),
    # Les jeux natifs restent exportables en .gmt : reproductibilité hors
    # application, et réutilisation dans un autre outil.
    downloadButton(ns("dl_scores_gmt"), i18n$t("Exporter les jeux de g\u00e8nes (.gmt)"),
                   class = "btn-sm btn-secondary w-100 mt-2"),
    div(class = "small text-muted mt-1", textOutput(ns("scores_status")))
  )
}

mod_bulk_pathways_output_ui <- function(id) {
  ns <- NS(id)
  card(
    full_screen = TRUE, max_height = "900px",
    card_header("Pathway Enrichment"),
    navset_tab(
      id = ns("pathways_tabs"),
      nav_panel(i18n$t("Barplot Top 15"), plotOutput(ns("pathway_barplot"), height = "580px")),
      nav_panel(i18n$t("Dotplot"),        plotOutput(ns("pathway_dotplot"), height = "580px")),
      # STAT-S2 — réseau d'enrichissement (descriptif : similarité de gènes / appartenance)
      nav_panel(i18n$t("Réseau"),         uiOutput(ns("network_ui"))),
      nav_panel(i18n$t("Table"),          DTOutput(ns("pathway_table"))),
      # PK-3 (pilier 3 convergence, fiche PROPOSAL_PACKS_ECOSYSTEM_2026-11-09.md
      # §3 jalon PK-3 ; note PACKS_PK3_CONVERGENCE_DESIGN_2026-11-10, D5 GO
      # opérateur 2026-11-10 = option A) — « convergent evidence » des runs
      # Packs multi-librairies (≥ 2) ; calcul pur côté serveur
      # (packs_split_libraries → packs_converge), aucun moteur re-run.
      nav_panel(i18n$t("Convergence"), value = "tab_convergence",
                uiOutput(ns("convergence_ui"))),
      nav_panel(i18n$t("Courbe GSEA"),    uiOutput(ns("gsea_curve_ui"))),
      # PB (GEMS D5) — ORA directionnelle : table goana/kegga + barcodeplot
      nav_panel(i18n$t("Directionnelle"), value = "tab_dir",
        fluidRow(
          column(6, DTOutput(ns("dir_table"))),
          column(6,
            selectInput(ns("dir_term"), i18n$t("Voie (barcodeplot)"), choices = NULL),
            plotOutput(ns("dir_barcode"), height = "480px"),
            div(class = "small text-muted",
                i18n$t("Rang log2FC complet de la table DE ; la voie s\u00e9lectionn\u00e9e est surlign\u00e9e (limma::barcodeplot)."))
          )
        )),
      # Bulk V2 M2 — scores par échantillon (voies x échantillons)
      # PB v2 (GEMS D5) — comparaison inter-contrastes (compareCluster)
      nav_panel(i18n$t("Comparative"), value = "tab_comp",
        fluidRow(column(12, DTOutput(ns("comp_table"))))),
      nav_panel(i18n$t("Scores par \u00e9chantillon"), value = "tab_scores",
        uiOutput(ns("scores_dropped_ui")),
        fluidRow(
          column(4, selectizeInput(ns("scores_color_by"), i18n$t("Colorer la PCA par"),
                                   choices = NULL,
                                   options = list(placeholder = .tr_placeholder(), allowEmptyOption = TRUE))),
          column(4, numericInput(ns("scores_heat_top_n"), i18n$t("Voies affichées (heatmap)"),
                                 value = 50, min = 5, max = 500, step = 5))
        ),
        plotOutput(ns("scores_pca"), height = "420px"),
        hr(),
        plotOutput(ns("scores_heatmap"), height = "620px"),
        hr(),
        fluidRow(
          column(6, downloadButton(ns("dl_scores_csv"), i18n$t("Export CSV (scores)"),
                                   class = "btn-sm btn-info w-100")),
          column(6, downloadButton(ns("dl_scores_rds"), i18n$t("Export RDS (résultat complet)"),
                                   class = "btn-sm btn-secondary w-100"))
        ),
        DTOutput(ns("scores_table"))
      )
    )
  )
}

mod_bulk_pathways_server <- function(id, global_data, shared_rv) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── DRIVE LIVE CONTROL (docs/DRIVE_LIVE_CONTROL_PLAN.md, grade G2) ─────
    # Two click sites in this module are fired by the file-drop poller through
    # counters it publishes in the session-scoped registry. `updateActionButton()`
    # does NOT click; these second triggers do.
    drive_counter_pathway <- shiny::reactiveVal(0L)
    drive_counter_scores  <- shiny::reactiveVal(0L)
    # READINESS GUARDS (G2). Each observer below opens with its own `req()` —
    # `shared_rv$filtered_counts` for the enrichment, `shared_rv$vst_mat` for
    # the per-sample scores — and `req()` aborts in SILENCE. Without a guard
    # the poller would fire the token, the observer would quietly do nothing,
    # and `result.json` would report `done` for an analysis that never ran.
    # Each guard mirrors its OWN observer's precondition, and
    # `shiny::isolate()` keeps the read out of the poller's dependency set.
    ready_pathway <- function() {
      if (is.null(shiny::isolate(shared_rv$filtered_counts))) {
        "no bulk object loaded (shared_rv$filtered_counts is NULL)"
      } else TRUE
    }
    ready_scores <- function() {
      if (is.null(shiny::isolate(shared_rv$vst_mat))) {
        "Step 1 has not produced a VST matrix (shared_rv$vst_mat is NULL)"
      } else TRUE
    }

    # ── STATE PUBLISHER (jalon « pathway state observability », 2026-09-24) ──
    # 🔴 LE DÉFAUT QUE CE PROBE EXISTE POUR TUER, mesuré en session VIVANTE :
    # `run_pipeline` sur `bulk-pathways-run_pathway` répond `done` AUSSI BIEN
    # quand l'enrichissement a produit 2 voies que quand l'observateur sort par
    # son garde silencieux des `< 10` gènes (l. ~471, `return()` nu). `done` est
    # HONNÊTE sur ce qu'il mesure — « le token a bougé » (spec §2.3) — et MUET
    # sur ce qui compte. Preuve de la cécité : au même instant `result.json`
    # portait `status="done"` et `snapshot$modules` ne contenait AUCUNE entrée
    # `bulk_pathways`, alors que le DOM portait la vérité (`✓ 2 pathways…`).
    #
    # Un `done` de PROTOCOLE n'est pas un `done` de MÉTIER. Ce probe publie le
    # second, par la couture EXISTANTE (`state =`, cf. `bulk_filter`/`bulk_de`).
    #
    # ⚠️ CASSETTE, pas lecture de l'objet vivant. Les deux sorties qui comptent
    # ne se distinguent PAS sur `shared_rv$pathway_results` : le garde `< 10`
    # gènes laisse le slot INCHANGÉ (`return()` nu). Sans enregistreur, « a
    # tourné, trop peu de gènes » et « n'a jamais tourné » rendent le même
    # `NULL`. L'enregistreur est donc le plus petit témoin honnête.
    #
    # ⚠️ `running` N'A PAS DE PRODUCTEUR ICI, et n'en aura pas : le job est
    # SYNCHRONE, donc aucun tick ne tourne pendant qu'il s'exécute. Un drapeau
    # posé dans le corps de l'observateur serait INOBSERVABLE — un `running`
    # sans producteur est un MENSONGE, pas une observabilité. Le producteur est
    # le poller lui-même : `long = TRUE` fait écrire `running` au dispatch. Ce
    # probe se contente donc de LIRE l'état du job, il ne l'invente pas.
    #
    # ⚠️ Un ENVIRONNEMENT et non une liste : écrit depuis les observateurs et lu
    # depuis le probe, il doit avoir une sémantique de RÉFÉRENCE (une liste
    # serait copiée à chaque écriture, donc l'écriture serait perdue).
    drive_last <- local({
      d <- new.env(parent = emptyenv())
      d$action    <- NA_character_
      d$status    <- NA_character_
      d$n_results <- NA_integer_
      d$error     <- ""
      d$elapsed_s <- NA_real_
      d$started   <- NULL
      d$seq       <- NULL
      d
    })

    # Enregistre une issue RÉSOLUE. `started` est posé par l'observateur, pas
    # ici : la durée doit mesurer le TRAVAIL, pas le temps écoulé depuis le
    # dernier appel réussi. Ne JAMAIS y écrire `running` — voir plus haut.
    drive_record <- function(action, status, n_results = NA_integer_, error = "") {
      drive_last$action    <- action
      drive_last$status    <- status
      drive_last$n_results <- n_results
      drive_last$error     <- error
      drive_last$elapsed_s <- if (is.null(drive_last$started)) {
        NA_real_
      } else {
        as.numeric(Sys.time()) - drive_last$started
      }
      # N7 timebox (satellites slice S3, jalon round8, 2026-11-03) : au site de
      # FINALISATION du job (durée MESURÉE ci-dessus), un dépassement du budget
      # déclaré (defaults.R TS_JOB_BUDGET_S_*) est consigné via le canal
      # provenance EXISTANT (ts_job_overrun_entry → new_provenance_entry /
      # provenance_append) — une ligne {job, budget, actual, ratio}, jamais
      # bloquant ; un job dans son budget ne produit AUCUNE ligne.
      elapsed_tb <- drive_last$elapsed_s
      budget_tb <- ts_job_budget(paste0("bulk-pathways-", action))
      if (is.finite(elapsed_tb) && elapsed_tb > as.numeric(budget_tb)) {
        provenance_append(shared_rv, ts_job_overrun_entry(
          job = paste0("bulk-pathways-", action),
          budget_s = as.numeric(budget_tb),
          actual_s = as.numeric(elapsed_tb)
        ))
      }
      invisible(NULL)
    }

    # STATE PROBE (G3) — cf. `state` dans `ts_drive_publish_token()`.
    # Scalaire UNIQUEMENT : jamais la matrice de voies, jamais un tracé, jamais
    # la table biologique, jamais un log. Chaque lecture partagée est
    # isolate()-gardée : ce probe tourne dans le battement réactif du poller
    # (spec §6), donc une lecture nue enrôlerait les données du module dans son
    # jeu de dépendances.
    drive_state <- function() {
      out <- list(
        module    = "bulk_pathways",
        action    = NA_character_,
        seq       = NA_integer_,
        status    = "not_ready",
        elapsed_s = NULL,
        n_results = NULL,
        error     = ""
      )
      # 1) Un job en vol est le SEUL producteur de `running` (cf. plus haut).
      #    N8 (SW-5b, registre round3, J2 2026-10-26) : un job TERMINAL —
      #    `$pending` pose par `ts_drive_job_finish()` (drive_watcher.R:4182) —
      #    n'est PAS en vol, meme s'il reste au registre : le beat qui ecrit le
      #    result terminal projette `ts_drive_snapshot()` AVANT le clear
      #    (drive_watcher.R:3723-3747), donc repondre `running` ici gelait le
      #    mensonge dans le result terminal d'un run_de en erreur. Sur chemin
      #    d'erreur le miroir se reinitialise a l'etat PROPRE du module
      #    (branches 2/3 ci-dessous — le job run_de appartient a bulk_de : ni
      #    son statut ni son action ne sont a publier ici). Un job REELLEMENT
      #    en vol (sans pending) reste le seul producteur de `running`.
      job <- ts_drive_job_state()
      if (!is.null(job) && is.null(job$pending)) {
        btn <- if (is.null(job$button)) "" else as.character(job$button)[1]
        out$action <- if (nzchar(btn)) sub("^bulk-pathways-", "", btn) else NA_character_
        out$seq    <- if (is.null(job$seq)) NA_integer_ else as.integer(job$seq)
        out$status <- "running"
        if (!is.null(job$started)) {
          out$elapsed_s <- round(as.numeric(Sys.time()) - job$started, 1L)
        }
        return(out)
      }
      # 2) Aucun objet de travail ⇒ `not_ready`, que l'objet n'ait jamais existé
      #    ou qu'il ait disparu. Un statut RÉSOLU ne survit pas à sa disparition :
      #    le module ne peut plus ni rejouer ni défendre ce résultat, donc il ne
      #    doit pas le présenter comme `done`. Les DEUX lectures sont isolées.
      ready <- !is.null(shiny::isolate(shared_rv$filtered_counts)) ||
               !is.null(shiny::isolate(shared_rv$vst_mat))
      st <- drive_last$status
      if (!isTRUE(ready) || is.null(st) || length(st) != 1L || is.na(st)) {
        return(out)
      }
      # 3) Issue RÉSOLUE, telle que l'enregistreur l'a MESURÉE. Liste
      #    RECONSTRUITE : `out$champ <- NULL` SUPPRIMERAIT le champ (mesuré), et
      #    le contrat de champs est fermé.
      return(list(
        module    = "bulk_pathways",
        action    = drive_last$action,
        seq       = if (is.null(drive_last$seq)) NA_integer_ else
                    as.integer(drive_last$seq),
        status    = as.character(st),
        elapsed_s = if (is.null(drive_last$elapsed_s) ||
                        is.na(drive_last$elapsed_s)) NULL else
                    round(as.numeric(drive_last$elapsed_s), 1L),
        n_results = if (is.null(drive_last$n_results) ||
                        is.na(drive_last$n_results)) NULL else
                    as.integer(drive_last$n_results),
        # Le seul champ LIBRE du contrat : il traverse le caviardage du badge
        # (un token ne doit pas fuir dans l'UI), et `known` ne lui est JAMAIS
        # appliqué — c'est là que vit la propriété de sécurité.
        error     = if (identical(as.character(st), "error")) {
          ts_drive_badge_sanitize(drive_last$error, 200L)
        } else ""
      ))
    }

    # SONDE DE VOCABULAIRE (Slice 3) — le domaine de `scores_source`, adressé
    # par index (`TS_DRIVE_SESSION_INPUTS`). LES MÊMES valeurs que le select
    # construit plus haut (`bulk_gene_set_choices()`, la source — pas le
    # libellé i18n) : le domaine dépend des paquets LOCALEMENT installés,
    # c'est le cas d'école du vocabulaire publié. LA MÊME closure alimente
    # l'état publié (canal wire) et la résolution d'index de l'applicateur
    # (`vocab =`), donc la rev épinglée par l'agent et les choix résolus au
    # battement d'application ne peuvent pas diverger. Aucune lecture
    # réactive : le catalogue est une fonction pure du disque.
    #
    # `vocab_rev` : compteur MONOTONE posé par empreinte — la sonde est la
    # seule scriptrice, aucun observateur existant n'est touché, l'UI humaine
    # ne change pas.
    drive_vocab_state <- new.env(parent = emptyenv())
    drive_vocab_state$rev <- 0L
    drive_vocab_state$fp  <- NULL
    drive_vocabulary <- function() {
      srcs <- tryCatch(as.character(unname(bulk_gene_set_choices())),
                       error = function(e) character(0))
      fp <- paste(c(length(srcs), srcs), collapse = "\u0001")
      if (!identical(fp, drive_vocab_state$fp)) {
        drive_vocab_state$rev <- drive_vocab_state$rev + 1L
        drive_vocab_state$fp  <- fp
      }
      list(scores_source = srcs, vocab_rev = drive_vocab_state$rev)
    }
    # Publié DANS l'état (le canal wire) : projété par
    # ts_drive_project_vocabulary (choix verbatim D1, keep-set fermé).
    drive_state2 <- function() {
      s <- drive_state()
      s$vocabulary <- drive_vocabulary()
      s
    }

    ts_drive_publish_token(global_data, "bulk-pathways-run_pathway", drive_counter_pathway,
                           ready = ready_pathway, state = drive_state2, long = TRUE,
                           vocab = drive_vocabulary)
    ts_drive_publish_token(global_data, "bulk-pathways-run_scores",  drive_counter_scores,
                           ready = ready_scores, state = drive_state2, long = TRUE,
                           vocab = drive_vocabulary)
    drive_trigger_pathway <- shiny::reactive(list(drive_counter_pathway(), input$run_pathway))
    drive_trigger_scores  <- shiny::reactive(list(drive_counter_scores(),  input$run_scores))

    # PK-2 (pilier 2 local) : sélecteur multi-librairies des packs — rendu
    # seulement sur la base « Packs » ; librairies = registre session
    # (global_data$packs, scan PK-1 au démarrage de chaque session).
    output$packs_libs_ui <- renderUI({
      global_data$language  # re-rendu au changement de langue (motif du module)
      if (!identical(input$pathway_db, "Packs")) return(NULL)
      ch <- packs_library_choices(shiny::isolate(global_data$packs %||% list()))
      if (length(ch) == 0L) {
        return(shiny::tags$div(class = "small text-muted",
                               .tr("Aucun pack de librairies installé — déposez un dossier dans data_sc/packs/ (effectif à la session suivante).")))
      }
      selectizeInput(ns("pathway_packs_libs"), .tr("Librairies de packs"),
                     choices = ch, multiple = TRUE,
                     options = list(placeholder = .tr("Sélectionnez une ou plusieurs librairies")))
    })

    # PLOT-S6b — jeux de gènes réellement chargés (source native ou .gmt).
    # Déclaré ICI, avant tout observe() qui le lit : un observe() évalue son
    # corps immédiatement, donc une définition plus bas serait introuvable.
    scores_sets_rv <- reactiveVal(NULL)

    .tr <- function(key) {
      tr <- global_data$i18n
      if (is.null(tr)) return(key)
      tryCatch(.strip_i18n_html(tr$t(key)), error = function(e) key)
    }

    # ── i18n push on language switch ─────────────────────────────────────
    observeEvent(global_data$language, {
      updateRadioButtons(session, "enrich_mode", label = .tr("M\u00e9thode d'enrichissement"),
        choices = stats::setNames(c("ora","gsea","dir","compare"),
          c(.tr("ORA \u2014 sur g\u00e8nes significatifs (classique)"), .tr("GSEA \u2014 sur tous les g\u00e8nes class\u00e9s (sans seuil)"),
            .tr("Directionnelle \u2014 goana/kegga (table DE compl\u00e8te, Up/Down par voie)"))))
      updateSelectInput(session, "pathway_source", label = .tr("Source de g\u00e8nes"),
        choices = stats::setNames(c("up","down","all_sig","manual"),
          c(.tr("G\u00e8nes Up (significatifs)"), .tr("G\u00e8nes Down (significatifs)"),
            .tr("Tous g\u00e8nes significatifs"), .tr("S\u00e9lection manuelle"))))
      updateSelectizeInput(session, "pathway_genes", label = .tr("G\u00e8nes"))
      updateSelectInput(session, "pathway_db",  label = .tr("Base de donn\u00e9es"))
      updateSelectInput(session, "pathway_org", label = .tr("Organisme"),
        choices = stats::setNames(c("human","mouse"), c(.tr("Humain"), .tr("Souris"))))
      updateNumericInput(session, "pathway_pval", label = .tr("P-value cutoff"))
      updateSelectInput(session, "pathway_padj_method", label = .tr("Méthode de correction (p-adj)"))
      updateActionButton(session, "run_pathway", label = .tr("Lancer Enrichissement"))
      # Bulk V2 M2 — libellés « scores par échantillon »
      # PLOT-S6b — la source de jeux de gènes est un selectInput : ses libellés
      # doivent être repoussés au changement de langue (les choix natifs sont
      # traduits à la construction de l'UI, cf. bulk_gene_set_choices(.tr_plain)).
      updateSelectInput(session, "scores_source", label = .tr("Source de jeux de gènes"),
        choices = bulk_gene_set_choices(.tr))
      updateSelectInput(session, "scores_org", label = .tr("Organisme"),
        choices = stats::setNames(c("human", "mouse"), c(.tr("Humain"), .tr("Souris"))))
      updateSelectInput(session, "scores_method", label = .tr("Méthode de scoring"),
        choices = stats::setNames(c("ssgsea", "gsva", "plage", "zscore"),
                                  c("ssGSEA", "GSVA", "PLAGE", "zscore")))
      updateNumericInput(session, "scores_min_size", label = .tr("Taille min voie"))
      updateNumericInput(session, "scores_max_size", label = .tr("Taille max voie"))
      updateActionButton(session, "run_scores", label = .tr("Lancer Scores par échantillon"))
      updateSelectizeInput(session, "scores_color_by", label = .tr("Colorer la PCA par"))
      updateNumericInput(session, "scores_heat_top_n", label = .tr("Voies affichées (heatmap)"))
      updateActionButton(session, "dl_scores_csv", label = .tr("Export CSV (scores)"))
      updateActionButton(session, "dl_scores_rds", label = .tr("Export RDS (résultat complet)"))
      updateActionButton(session, "dl_scores_gmt", label = .tr("Exporter les jeux de gènes (.gmt)"))
    }, ignoreInit = TRUE)

    # (unchanged) manual gene picker refresh + mirroring
    observeEvent(shared_rv$filtered_counts, {
      req(shared_rv$filtered_counts)
      updateSelectizeInput(session, "pathway_genes",
                           choices = rownames(shared_rv$filtered_counts), server = TRUE)
    })
    observe({
      shared_rv$pathway_db   <- input$pathway_db
      shared_rv$pathway_mode <- input$enrich_mode
      # STAT-Q1 : reprise dans le script R reproductible exporté
      shared_rv$pathway_padj_method <- input$pathway_padj_method %||% TS_PADJ_METHOD_DEFAULT
    })

    .active_de_results <- function() {
      ac <- shared_rv$active_contrast
      if (is.null(ac) || !ac %in% names(shared_rv$contrasts)) return(NULL)
      shared_rv$contrasts[[ac]]
    }

    observe({
      shinyjs::toggleState("run_pathway", condition = !is.null(shared_rv$filtered_counts))
      shinyjs::toggleState("dl_pathway",  condition = !is.null(shared_rv$pathway_results))
      # Bulk V2 M2 — le scoring exige la matrice VST (étape 1), pas les counts.
      shinyjs::toggleState("run_scores", condition = !is.null(shared_rv$vst_mat))
      shinyjs::toggleState("dl_scores_csv", condition = !is.null(shared_rv$pathway_scores))
      shinyjs::toggleState("dl_scores_rds", condition = !is.null(shared_rv$pathway_scores))
      # PLOT-S6b — l'export .gmt ne dépend que des jeux chargés (pas du scoring).
      shinyjs::toggleState("dl_scores_gmt", condition = !is.null(scores_sets_rv()))
    })

    # ── Bulk V2 M2 — scores de voies par échantillon ────────────────────────
    observeEvent(shared_rv$vst_mat, {
      req(global_data$bulk_obj$metadata)
      meta <- global_data$bulk_obj$metadata
      updateSelectizeInput(session, "scores_color_by", choices = colnames(meta), server = FALSE)
    }, ignoreNULL = TRUE)

    output$scores_status <- renderText({
      global_data$language
      sc <- shared_rv$pathway_scores
      if (is.null(sc)) .tr("En attente — lancez d'abord le Filtrage & VST (étape 1).")
      else .t_fmt(.tr("\u2713 {n} voies scorées x {m} échantillons [ {meth} ]"),
                  n = nrow(sc$scores), m = ncol(sc$scores), meth = sc$method)
    })

    # PLOT-S6b — jeux de gènes réellement chargés (natifs ou .gmt) : conservés
    # pour l'export .gmt, qui doit rendre EXACTEMENT ce qui a été scoré.
    # (le reactiveVal lui-même est déclaré en tête de serveur : lu par un
    #  observe() qui s'exécute avant cette section.)

    output$scores_source_hint <- renderText({
      global_data$language
      src <- input$scores_source %||% "msigdb_hallmark"
      if (identical(src, "file")) {
        return(.tr("Format attendu : nom<TAB>description<TAB>gènes..."))
      }
      cat_df <- bulk_gene_set_catalog()
      row <- cat_df[cat_df$source == src, , drop = FALSE]
      if (nrow(row) == 0L) return("")
      if (!isTRUE(row$available)) {
        return(.t_fmt(.tr("\u26a0\ufe0f Source indisponible — installez {pkg} (aucun téléchargement automatique)."),
                      pkg = row$requires))
      }
      .tr(row$description)
    })

    observeEvent(drive_trigger_scores(), {
      # ── DECLARE THE JOB OVER (drive job contract, spec §5) ────────────────
      # Même contrat que l'observateur d'enrichissement ci-dessus : ce jeton est
      # DÉCLARÉ long (le calcul GSVA/ssGSEA est synchrone et long), donc il doit
      # fermer son job sur TOUTE sortie, sinon le contrat B refuserait à jamais
      # le `run_pipeline` suivant et le module deviendrait impilotable — la
      # panne silencieuse la plus coûteuse.
      outcome <- new.env(parent = emptyenv())
      outcome$v <- "refused"
      if (isTRUE(ts_drive_job_busy()) &&
          identical(ts_drive_job_state()$button, "bulk-pathways-run_scores")) {
        # F4 : même contrat que run_pathway — un refus porte sa raison.
        on.exit(ts_drive_job_finish("bulk-pathways-run_scores",
          status = switch(outcome$v, ok = "done", failed = "error", "invalid"),
          error = if (identical(outcome$v, "failed")) {
            "the per-sample scoring raised — see the app notification"
          } else outcome$reason
        ), add = TRUE)
      }

      req(input$run_scores > 0 || shiny::isolate(drive_counter_scores()) > 0)
      req(shared_rv$vst_mat)
      drive_last$seq     <- ts_drive_job_state()$seq
      drive_last$started <- as.numeric(Sys.time())
      # PLOT-S6b — la source décide : jeux NATIFS (aucun fichier à fournir) ou
      # .gmt fourni. Les deux rendent la même forme (nom -> gènes), donc tout
      # l'aval (compute_pathway_scores, exports, vues) est inchangé.
      src <- input$scores_source %||% "msigdb_hallmark"
      if (identical(src, "file")) {
        if (is.null(input$scores_gmt) || is.null(input$scores_gmt$datapath)) {
          showNotification(.tr("\u26a0\ufe0f Fournissez un fichier .gmt (jeux de gènes)."),
                           type = "warning", duration = 5)
          outcome$reason <- "no .gmt file was given"
          drive_record("run_scores", "not_ready")
          return()
        }
        sets <- tryCatch(bulk_parse_gmt(input$scores_gmt$datapath), error = function(e) e)
        if (inherits(sets, "error")) {
          showNotification(paste(.tr("Erreur GMT:"), conditionMessage(sets)),
                           type = "error", duration = 8)
          outcome$reason <- "the .gmt file could not be read"
          drive_record("run_scores", "error", error = conditionMessage(sets))
          return()
        }
      } else {
        sets <- tryCatch(bulk_load_gene_sets(src, input$scores_org %||% "human"),
                         error = function(e) e)
        if (inherits(sets, "error")) {
          showNotification(paste(.tr("Erreur jeux de gènes :"), conditionMessage(sets)),
                           type = "error", duration = 8)
          outcome$reason <- "the gene set source could not be read"
          drive_record("run_scores", "error", error = conditionMessage(sets))
          return()
        }
      }
      scores_sets_rv(sets)
      # `add = TRUE` est PORTEUR : un `on.exit()` nu effacerait la déclaration
      # de job posée en tête d'observateur (cf. le commentaire long ci-dessus).
      p <- shiny::Progress$new(); on.exit(p$close(), add = TRUE)
      p$set(message = .tr("Scores de voies par échantillon..."), value = 0.2)
      tryCatch({
        res <- compute_pathway_scores(
          shared_rv$vst_mat, sets,
          method      = input$scores_method %||% "ssgsea",
          min_size    = input$scores_min_size %||% TS_BULK_GSVA_MIN_SIZE,
          max_size    = input$scores_max_size %||% TS_BULK_GSVA_MAX_SIZE,
          overlap_min = TS_BULK_GSVA_OVERLAP_MIN
        )
        shared_rv$pathway_scores <- res
        # Stockage contractuel : bulk_obj$pathways$per_sample (réaffectation
        # complète pour déclencher la réactivité de global_data).
        bo <- global_data$bulk_obj
        if (!is.list(bo$pathways)) bo$pathways <- list()
        bo$pathways$per_sample <- res$scores
        bo <- bulk_ensure_provenance(bo)
        global_data$bulk_obj <- bo
        if (nrow(res$qc$dropped) > 0L) {
          showNotification(.t_fmt(.tr("\u26a0\ufe0f {n} jeu(x) de gènes rejeté(s) — voir le détail dans l'onglet Scores."),
                                   n = nrow(res$qc$dropped)), type = "warning", duration = 6)
        }
        showNotification(.t_fmt(.tr("\u2713 {n} voies scorées (ssGSEA/GSVA) sur {m} échantillons."),
                                 n = nrow(res$scores), m = ncol(res$scores)), type = "message")
        nav_select(id = "pathways_tabs", selected = "tab_scores", session = session)
        # `empty` reste distinct de `done` : tous les jeux de gènes peuvent être
        # rejetés (recouvrement < seuil), et « 0 voie » n'est pas « N voies ».
        if (nrow(res$scores) == 0L) {
          drive_record("run_scores", "empty", 0L)
        } else {
          drive_record("run_scores", "done", nrow(res$scores))
        }
        outcome$v <- "ok"
      }, error = function(e) {
        showNotification(paste(.tr("Erreur scores:"), conditionMessage(e)),
                         type = "error", duration = 8)
        shared_rv$pathway_scores <- NULL
        drive_record("run_scores", "error", error = conditionMessage(e))
        outcome$v <- "failed"
      })
    })

    output$scores_dropped_ui <- renderUI({
      global_data$language
      sc <- shared_rv$pathway_scores
      if (is.null(sc) || nrow(sc$qc$dropped) == 0L) return(NULL)
      div(class = "alert alert-warning", style = "font-size:0.82em;",
          icon("triangle-exclamation"), " ",
          .t_fmt(.tr("{n} jeu(x) rejeté(s) — recouvrement < {pct} % ou taille hors bornes. Détail :"),
                 n = nrow(sc$qc$dropped), pct = round(100 * TS_BULK_GSVA_OVERLAP_MIN)),
          ts_datatable(sc$qc$dropped, page_length = 5, buttons = FALSE,
                       filter = "none", scroll_x = FALSE))
    })

    output$scores_pca <- renderPlot({
      global_data$language
      req(shared_rv$pathway_scores)
      color_by <- if (nzchar(input$scores_color_by %||% "")) input$scores_color_by else NULL
      tryCatch(
        plot_pathway_scores_pca(shared_rv$pathway_scores,
                                metadata = global_data$bulk_obj$metadata,
                                color_by = color_by, tr = .tr_fn(global_data)),
        error = function(e) {
          ggplot2::ggplot() +
            ggplot2::annotate("text", x = 1, y = 1,
                              label = paste(.tr("Erreur:"), conditionMessage(e)), color = "red") +
            ggplot2::theme_void()
        }
      )
    })

    output$scores_heatmap <- renderPlot({
      global_data$language
      req(shared_rv$pathway_scores)
      h <- plot_pathway_scores_heatmap(shared_rv$pathway_scores,
                                       top_n = input$scores_heat_top_n %||% 50,
                                       tr = .tr_fn(global_data))
      # ComplexHeatmap -> print() explicite ; ggplot de repli s'imprime pareil.
      print(h)
    })

    output$scores_table <- renderDT({
      global_data$language
      req(shared_rv$pathway_scores)
      s <- round(shared_rv$pathway_scores$scores, 4)
      ts_datatable(as.data.frame(s), page_length = 15L,
                   filename_base = "bulk_pathways_scores", filter = "none",
                   rownames = TRUE)
    })

    output$dl_scores_csv <- downloadHandler(
      filename = function() paste0("pathway_scores_", shared_rv$pathway_scores$method,
                                   "_", Sys.Date(), ".csv"),
      content  = function(file) {
        req(shared_rv$pathway_scores)
        write.csv(build_pathway_scores_export(shared_rv$pathway_scores), file, row.names = FALSE)
      }
    )
    output$dl_scores_rds <- downloadHandler(
      filename = function() paste0("pathway_scores_", shared_rv$pathway_scores$method,
                                   "_", Sys.Date(), ".rds"),
      content  = function(file) {
        req(shared_rv$pathway_scores)
        saveRDS(shared_rv$pathway_scores, file)
      }
    )
    # PLOT-S6b — export .gmt des jeux RÉELLEMENT scorés : une source native
    # (MSigDB/PROGENy/DoRothEA) reste ainsi réutilisable hors application, et
    # le .gmt exporté se relit à l'identique (bulk_write_gmt -> bulk_parse_gmt).
    output$dl_scores_gmt <- downloadHandler(
      filename = function() {
        src <- input$scores_source %||% "msigdb_hallmark"
        paste0("gene_sets_", src, "_", Sys.Date(), ".gmt")
      },
      content  = function(file) {
        sets <- scores_sets_rv()
        req(sets)
        bulk_write_gmt(sets, file)
      }
    )

    # ── Enrichment (ORA + GSEA) — strings translated ─────────────────────
    observeEvent(drive_trigger_pathway(), {
      # ── DECLARE THE JOB OVER (drive job contract, spec §5) ────────────────
      # Producer TERMINAL que `long = TRUE` engage. Enregistré AVANT le premier
      # `req()` pour que TOUTE sortie rapporte quelque chose : l'abandon d'un
      # `req()`, le garde GSEA sans DE et le garde `< 10` gènes laissent
      # `outcome$v` à "refused" -> `invalid` (un refus n'est PAS un échec, même
      # raisonnement qu'ailleurs : un refus ne doit pas apprendre à ignorer un
      # badge rouge) ; le chemin nominal pose "ok" -> `done` ; le `catch` de la
      # computation pose "failed" -> `error`.
      #
      # ⚠️ `outcome` est un ENVIRONNEMENT, pas une variable locale. MESURÉ
      # (2026-09-24) : un `job_outcome <- "failed"` écrit DANS un handler
      # `error = function(e) …` se lie à la frame DU HANDLER et ne remonte
      # JAMAIS — la sonde a rendu `refused` sur un `stop()` réel. C'est le motif
      # exact de `mod_bulk_de_run.R` (l. 352), qui annonce donc `invalid` là où
      # la DE a levé. Une cellule d'environnement a une sémantique de RÉFÉRENCE
      # et traverse la frontière de closure.
      outcome <- new.env(parent = emptyenv())
      outcome$v <- "refused"
      # `add = TRUE` est PORTEUR, pas du style : un `on.exit()` nu REMPLACE
      # toute expression en attente, donc il effacerait cette déclaration — le
      # job resterait `running` à jamais et l'agent sonderait une complétion
      # qui ne serait jamais écrite.
      if (isTRUE(ts_drive_job_busy()) &&
          identical(ts_drive_job_state()$button, "bulk-pathways-run_pathway")) {
        # F4 (2026-10-06, mesuré en session réelle) : un refus porte sa raison
        # (`outcome$reason`) sur le wire — un `invalid` à errors[] vide force
        # l'agent à deviner. Les textes restent sanitiser-sure (mots < 8 cars).
        on.exit(ts_drive_job_finish("bulk-pathways-run_pathway",
          status = switch(outcome$v, ok = "done", failed = "error", "invalid"),
          error = if (identical(outcome$v, "failed")) {
            "the pathway enrichment raised — see the app notification"
          } else outcome$reason
        ), add = TRUE)
      }

      req(input$run_pathway > 0 || shiny::isolate(drive_counter_pathway()) > 0)
      req(shared_rv$filtered_counts)
      # `seq` appartient au PROTOCOLE, pas au module : il n'est lisible que tant
      # que le job est en vol. Capturé ici pour que l'issue résolue reste
      # corrélable à la commande qui l'a demandée.
      drive_last$seq     <- ts_drive_job_state()$seq
      drive_last$started <- as.numeric(Sys.time())
      p <- shiny::Progress$new(); on.exit(p$close(), add = TRUE)

      if (input$enrich_mode == "gsea") {
        res_de <- .active_de_results()
        if (is.null(res_de)) {
          showNotification(.tr("\u26a0\ufe0f Lancez d'abord l'\u00e9tape 2 (Analyse Diff\u00e9rentielle)."), type = "warning")
          # Refus : rien n'a tourné. `not_ready` est la réponse honnête — le
          # module ne peut pas produire de résultat faute d'amont.
          outcome$reason <- "no DE result yet: run the DE step first"
          drive_record("run_pathway", "not_ready")
          return()
        }
        p$set(message = .tr("GSEA en cours..."), value = 0.3)
        tryCatch({
          res <- run_gsea_enrichment(res_de, organism = input$pathway_org,
                                     database = input$pathway_db, pval_cutoff = input$pathway_pval,
                                     p_adjust_method = input$pathway_padj_method %||% TS_PADJ_METHOD_DEFAULT)
          if (nrow(res) == 0) {
            showNotification(.tr("\u2139\ufe0f Aucun pathway enrichi trouv\u00e9 (GSEA)."), type = "warning")
            shared_rv$pathway_results <- NULL
            # `empty` et non `done` : l'analyse a TOURNÉ et n'a rien trouvé.
            outcome$reason <- "the GSEA ran and found no pathway"
            drive_record("run_pathway", "empty", 0L)
            return()
          }
          shared_rv$pathway_results <- res
          showNotification(.t_fmt(.tr("\u2705 {n} pathways enrichis (GSEA)"), n = nrow(res)), type = "message")
          shared_rv$active_tab <- "tab_pathway"
          drive_record("run_pathway", "done", nrow(res))
          outcome$v <- "ok"
        }, error = function(e) {
          showNotification(.t_fmt(.tr("\u274c Erreur GSEA: {msg}"), msg = as.character(e$message)[1]),
                           type = "error", duration = 8)
          shared_rv$pathway_results <- NULL
          drive_record("run_pathway", "error", error = as.character(e$message)[1])
          outcome$v <- "failed"
        })
        return()
      }

      # ── PB (GEMS D5) — ORA DIRECTIONNELLE (goana/kegga sur la table DE) ────
      # Même protocole drive que les deux modes existants (not_ready / empty /
      # done / error) : les gardes moteur classés (pathway_error) sortent par
      # la branche error avec leur raison honnête.
      if (input$enrich_mode == "dir") {
        res_de <- .active_de_results()
        if (is.null(res_de)) {
          showNotification(.tr("\u26a0\ufe0f Lancez d'abord l'\u00e9tape 2 (Analyse Diff\u00e9rentielle)."), type = "warning")
          outcome$reason <- "no DE result yet: run the DE step first"
          drive_record("run_pathway", "not_ready")
          return()
        }
        db_dir <- switch(input$pathway_db %||% "GOBP", GOBP = "GO", KEGG = "KEGG", "GO")
        p$set(message = .tr("ORA directionnelle (goana/kegga) en cours..."), value = 0.3)
        tryCatch({
          res <- run_pathway_direction(res_de,
                                       universe_ids = rownames(shiny::isolate(shared_rv$filtered_counts)),
                                       db = db_dir,
                                       organism = input$pathway_org,
                                       fdr_cutoff = input$dir_fdr %||% 0.05)
          tbl <- rbind(res$go, res$kegg)
          if (is.null(tbl) || nrow(tbl) == 0L) {
            showNotification(.tr("\u2139\ufe0f Aucune voie testable (goana/kegga)."), type = "warning")
            shared_rv$pathway_dir_results <- NULL
            outcome$reason <- "the directional ORA ran and found no pathway"
            drive_record("run_pathway", "empty", 0L)
            return()
          }
          shared_rv$pathway_dir_results <- res
          showNotification(.t_fmt(.tr("\u2705 ORA directionnelle : {n} voies (goana/kegga)"), n = nrow(tbl)),
                           type = "message")
          shared_rv$active_tab <- "tab_dir"
          drive_record("run_pathway", "done", nrow(tbl))
          outcome$v <- "ok"
        }, error = function(e) {
          showNotification(.t_fmt(.tr("\u274c Erreur directionnelle: {msg}"), msg = as.character(e$message)[1]),
                           type = "error", duration = 8)
          shared_rv$pathway_dir_results <- NULL
          drive_record("run_pathway", "error", error = as.character(e$message)[1])
          outcome$v <- "failed"
        })
        return()
      }

      # ── PB v2 (GEMS D5, jalon 2) — COMPARATIVE (compareCluster) ────────────
      # Listes = gènes significatifs PAR CONTRASTE (shared_rv$contrasts) ;
      # fond = rownames de la matrice filtrée (M12, un seul fond déclaré) ;
      # même protocole drive que les 3 modes existants.
      if (input$enrich_mode == "compare") {
        fdr_c <- input$comp_fdr %||% 0.05
        direction <- input$comp_source %||% "all"
        lists <- lapply(names(shared_rv$contrasts), function(nm) {
          de <- shared_rv$contrasts[[nm]]
          if (is.null(de) || !is.data.frame(de) ||
              !all(c("gene", "padj") %in% names(de))) return(NULL)
          sel <- !is.na(de$padj) & de$padj < fdr_c
          if (identical(direction, "up"))   sel <- sel & de$log2FoldChange > 0
          if (identical(direction, "down")) sel <- sel & de$log2FoldChange < 0
          unique(as.character(de$gene[sel]))
        })
        names(lists) <- names(shared_rv$contrasts)
        lists <- Filter(function(x) !is.null(x) && length(x) > 0, lists)
        if (length(lists) < 2) {
          showNotification(.tr("⚠️ Comparative : au moins 2 contrastes avec des gènes significatifs au seuil requis — lancez l'étape 2 (Analyse Différentielle) sur plusieurs contrastes."), type = "warning")
          outcome$reason <- "fewer than 2 contrasts with DE genes at threshold"
          drive_record("run_pathway", "not_ready")
          return()
        }
        db_c <- switch(input$pathway_db %||% "GOBP", GOBP = "GO", KEGG = "KEGG", "GO")
        p$set(message = .tr("Comparative (compareCluster) en cours..."), value = 0.3)
        tryCatch({
          res <- run_pathway_compare_cluster(
            gene_lists = lists,
            universe_ids = rownames(shiny::isolate(shared_rv$filtered_counts)),
            db = db_c, organism = input$pathway_org,
            pvalue_cutoff = 0.05, qvalue_cutoff = 0.2)
          tbl <- rbind(res$go, res$kegg)
          if (is.null(tbl) || nrow(tbl) == 0L) {
            showNotification(.tr("ℹ️ Aucune voie comparative sous les seuils."), type = "warning")
            shared_rv$pathway_comp_results <- NULL
            outcome$reason <- "compareCluster ran and found no pathway"
            drive_record("run_pathway", "empty", 0L)
            return()
          }
          shared_rv$pathway_comp_results <- res
          shared_rv$active_tab <- "tab_comp"
          showNotification(.t_fmt(.tr("✅ Comparative (compareCluster) : {n} voies"), n = nrow(tbl)),
                           type = "message")
          drive_record("run_pathway", "done", nrow(tbl))
          outcome$v <- "ok"
        }, error = function(e) {
          showNotification(.t_fmt(.tr("❌ Erreur comparative: {msg}"), msg = as.character(e$message)[1]),
                           type = "error", duration = 8)
          shared_rv$pathway_comp_results <- NULL
          drive_record("run_pathway", "error", error = as.character(e$message)[1])
          outcome$v <- "failed"
        })
        return()
      }

      genes_to_test <- NULL
      if (input$pathway_source %in% c("up","down","all_sig")) {
        res <- .active_de_results()
        if (is.null(res)) {
          showNotification(.tr("\u26a0\ufe0f Lancez d'abord l'\u00e9tape 2 (Analyse Diff\u00e9rentielle)."), type = "warning")
          outcome$reason <- "no DE result yet: run the DE step first"
          drive_record("run_pathway", "not_ready")
          return()
        }
        sig <- res$padj < (shared_rv$padj_thresh %||% 0.05) &
               abs(res$log2FoldChange) > (shared_rv$lfc_thresh %||% 1)
        sig[is.na(sig)] <- FALSE
        genes_to_test <- switch(input$pathway_source,
          up = res$gene[sig & res$log2FoldChange > 0],
          down = res$gene[sig & res$log2FoldChange < 0],
          all_sig = res$gene[sig])
      } else {
        req(input$pathway_genes); genes_to_test <- input$pathway_genes
      }
      genes_to_test <- unique(trimws(genes_to_test)); genes_to_test <- genes_to_test[nchar(genes_to_test) > 0]

      if (length(genes_to_test) < 10) {
        showNotification(.t_fmt(.tr("\u26a0\ufe0f Trop peu de g\u00e8nes ({n}). Minimum 10 requis."), n = length(genes_to_test)),
                         type = "warning", duration = 5)
        # 🔴 LE GARDE SILENCIEUX QUE CE JALON REND VISIBLE. Avant lui, cette
        # sortie et un enrichissement réussi de 2 voies rendaient le MÊME
        # `done` de protocole, et `shared_rv$pathway_results` restait INCHANGÉ
        # (d'où la cassette : l'objet vivant ne peut pas les distinguer).
        # F4 : la raison part aussi sur le wire (mesuré muet en session réelle).
        outcome$reason <- "fewer than 10 genes to test; the min is 10"
        drive_record("run_pathway", "empty", 0L)
        return()
      }

      # Univers ORA = gènes réellement mesurés du jeu chargé (aligne Bulk sur le
      # chemin SC, qui passe déjà rownames(sc_obj)) ; NULL en saisie manuelle
      # sans objet -> comportement inchangé (enrichissement pan-génome).
      fc <- shiny::isolate(shared_rv$filtered_counts)
      universe_genes <- if (is.null(fc)) NULL else rownames(fc)
      # PK-2 (pilier 2 local, fiche packs §3/§4) : la base « Packs » — refus
      # honnêtes AVANT le vol (motif du garde < 10 gènes : une raison porte sa
      # sortie wire, jamais un `done` de protocole ambigu).
      packs_libs <- NULL
      if (identical(input$pathway_db, "Packs")) {
        if (!identical(input$pathway_org, "human")) {
          # D4 (fiche packs §4.4) = décision JOINT GEMS D5 : les librairies
          # Enrichr sont des symboles HUMAINS ; l'orthologie souris attend cet
          # arbitrage — jamais un mapping silencieux.
          showNotification(.tr("Packs : orthologie souris en attente de l'arbitrage D4 (GEMS D5) — librairies symboles humains uniquement."),
                           type = "warning", duration = 8)
          outcome$reason <- "packs: mouse orthology awaits D4 arbitration"
          drive_record("run_pathway", "empty", 0L)
          return()
        }
        packs_libs <- input$pathway_packs_libs
        if (is.null(packs_libs) || length(packs_libs) == 0L) {
          showNotification(.tr("Packs : sélectionnez au moins une librairie."),
                           type = "warning", duration = 6)
          outcome$reason <- "packs: no library selected"
          drive_record("run_pathway", "empty", 0L)
          return()
        }
      }
      p$set(message = .tr("Enrichissement..."), value = 0.3)
      tryCatch({
        res <- if (identical(input$pathway_db, "Packs")) {
          # PK-2 : enricher() en espace SYMBOLES — aucun bitr (les librairies
          # Enrichr consomment les symboles directement, fiche §4.1) ;
          # l'univers reste les gènes mesurés du jeu chargé.
          t2g <- packs_term2gene(shiny::isolate(global_data$packs %||% list()),
                                 packs_libs)
          run_packs_enrichment(genes = genes_to_test,
                               term2gene = t2g$term2gene,
                               term2name = t2g$term2name,
                               universe = universe_genes,
                               pval_cutoff = input$pathway_pval,
                               p_adjust_method = input$pathway_padj_method %||% TS_PADJ_METHOD_DEFAULT)
        } else {
          run_pathway_enrichment(genes = genes_to_test, organism = input$pathway_org,
                                 database = input$pathway_db, pval_cutoff = input$pathway_pval,
                                 p_adjust_method = input$pathway_padj_method %||% TS_PADJ_METHOD_DEFAULT,
                                 universe = universe_genes)
        }
        if (nrow(res) == 0) {
          showNotification(.tr("\u2139\ufe0f Aucun pathway enrichi trouv\u00e9."), type = "warning")
          shared_rv$pathway_results <- NULL
          drive_record("run_pathway", "empty", 0L)
          return()
        }
        shared_rv$pathway_results <- res
        showNotification(.t_fmt(.tr("\u2705 {n} pathways enrichis"), n = nrow(res)), type = "message")
        shared_rv$active_tab <- "tab_pathway"
        drive_record("run_pathway", "done", nrow(res))
        outcome$v <- "ok"
      }, error = function(e) {
        showNotification(.t_fmt(.tr("\u274c Erreur pathway: {msg}"), msg = as.character(e$message)[1]),
                         type = "error", duration = 6)
        shared_rv$pathway_results <- NULL
        drive_record("run_pathway", "error", error = as.character(e$message)[1])
        outcome$v <- "failed"
      })
    })

    output$pathway_status <- renderText({
      global_data$language
      # N7 timebox (satellites slice S3, jalon round8, 2026-11-03) : pendant un
      # job en vol, la surface EXISTANTE drive_state porte déjà elapsed_s —
      # le budget déclaré (defaults.R TS_JOB_BUDGET_S_*) est rendu À CÔTÉ
      # (« en cours — … (typique …) »), zéro spinner nouveau. Le libellé FR
      # passe par l'injecteur i18n du module (clés append-only dans
      # translation.json). Hors vol, la logique existante ci-dessous est
      # INTACTE (additions-only).
      st_tb <- drive_state()
      if (isTRUE(st_tb$status == "running") && !is.null(st_tb$elapsed_s) &&
          !is.na(st_tb$action)) {
        return(ts_job_budget_label(
          elapsed_s = as.numeric(st_tb$elapsed_s),
          budget_s  = as.numeric(ts_job_budget(paste0("bulk-pathways-",
                                                      st_tb$action))),
          tr = .tr
        ))
      }
      if (is.null(shared_rv$pathway_results)) .tr("Aucune analyse en cours")
      else .t_fmt(.tr("\u2713 {n} pathways trouv\u00e9s [ {db} ]"),
                  n = nrow(shared_rv$pathway_results), db = input$pathway_db)
    })

    output$pathway_barplot <- renderPlot({
      global_data$language                     # i18n trigger
      req(shared_rv$pathway_results)
      plot_pathway_barplot(shared_rv$pathway_results, db_label = input$pathway_db, top_n = 15,
                           tr = .tr_fn(global_data),
                           palette = shared_rv$bulk_palette %||% "default")
    })
    output$pathway_dotplot <- renderPlot({
      global_data$language                     # i18n trigger
      req(shared_rv$pathway_results)
      plot_pathway_dotplot(shared_rv$pathway_results, db_label = input$pathway_db, top_n = 20,
                           tr = .tr_fn(global_data),
                           palette = shared_rv$bulk_palette %||% "default")
    })
    output$pathway_table <- renderDT({
      global_data$language                     # i18n trigger
      req(shared_rv$pathway_results)
      build_pathway_dt(shared_rv$pathway_results, tr = .tr_fn(global_data))
    })

    # ── PK-3 (pilier 3 convergence ; note PACKS_PK3_CONVERGENCE_DESIGN_2026-11-
    # 10, D5 GO opérateur 2026-11-10 = option A) — vue « Convergence » ────────
    # Calcul À LA VOLÉE depuis le run Packs en mémoire : packs_split_libraries
    # découpe la table combinée par préfixe pack::fichier:: puis packs_converge
    # recalcule la BH GLOBALE des p brutes portées par le contrat PK-2 et
    # effondre la redondance (enrichplot::pairwise_termsim) — AUCUN moteur
    # re-run (pilier 4 = GEMS PB, jamais double-owner, fiche §4/D1). Le
    # p.adjust PK-2 déjà affiché dans l'onglet Table reste la valeur PAR RUN du
    # contrat (inchangée) ; l'onglet Convergence porte la valeur GLOBALE en
    # colonne dédiée (`p.adjust GLOBAL`) et l'original en audit
    # (`p.adjust PK-2 (run)`). Refus honnêtes (classe pathway_error) affichés
    # tels quels : < 2 librairies, run non-Packs, contrat violé.
    convergence_data <- reactive({
      global_data$language                     # i18n trigger (comme pathway_table)
      res <- shared_rv$pathway_results
      if (is.null(res) || !is.data.frame(res) || nrow(res) == 0L) {
        return(list(ok = FALSE,
                    msg = .tr("La convergence s'applique aux runs Packs multi-librairies (≥ 2 librairies sélectionnées).")))
      }
      tbl <- tryCatch(packs_converge(packs_split_libraries(res)),
                      error = function(e) e)
      if (inherits(tbl, "error")) return(list(ok = FALSE, msg = conditionMessage(tbl)))
      list(ok = TRUE, tbl = tbl)
    })
    output$convergence_ui <- renderUI({
      cd <- convergence_data()
      if (!isTRUE(cd$ok)) {
        return(tags$div(class = "small text-muted", style = "margin:15px;",
                        icon("info-circle"), " ", cd$msg))
      }
      tagList(
        tags$div(class = "small text-muted", style = "margin: 6px 15px;",
                 icon("info-circle"), " ",
                 .tr("BH globale sur l'union des tests (toutes librairies confondues) + effondrement de redondance (similarité de Jaccard enrichplot) — 1 cluster = 1 vote.")),
        DTOutput(ns("convergence_table"))
      )
    })
    output$convergence_table <- renderDT({
      cd <- convergence_data()
      req(isTRUE(cd$ok))
      build_convergence_dt(cd$tbl, tr = .tr_fn(global_data))
    })

    # ── PB (GEMS D5) — résultats ORA directionnelle (goana/kegga) ────────────
    # Tables limma natives (colonne db estampillée par le moteur) empilées en
    # un seul tableau affichage ; ID = rownames (GO:… / hsa…).
    .stack_dir_tables <- function(res) {
      parts <- list()
      if (!is.null(res$go)) {
        t <- res$go
        t$ID <- rownames(t)
        parts[[length(parts) + 1L]] <-
          t[, c("ID", "Term", "Ont", "N", "Up", "Down",
                "P.Up", "FDR.Up", "P.Down", "FDR.Down", "db"), drop = FALSE]
      }
      if (!is.null(res$kegg)) {
        t <- res$kegg
        t$ID <- rownames(t)
        names(t)[names(t) == "Pathway"] <- "Term"
        t$Ont <- "KEGG"
        parts[[length(parts) + 1L]] <-
          t[, c("ID", "Term", "Ont", "N", "Up", "Down",
                "P.Up", "FDR.Up", "P.Down", "FDR.Down", "db"), drop = FALSE]
      }
      if (!length(parts)) return(NULL)
      do.call(rbind, parts)
    }

    output$dir_table <- renderDT({
      global_data$language                     # i18n trigger
      req(shared_rv$pathway_dir_results)
      t <- .stack_dir_tables(shared_rv$pathway_dir_results)
      req(t)
      tr <- .tr_fn(global_data)
      colnames(t) <- c("ID", tr("Voie"), tr("Ontologie"), tr("Nb Gènes"),
                       tr("Up"), tr("Down"), tr("P Up"), tr("FDR Up"),
                       tr("P Down"), tr("FDR Down"), "db")[seq_along(t)]
      ts_datatable(t, page_length = 15L, filename_base = "pathways_dir")
    })

    # PB v2 (GEMS D5, jalon 2) — table comparative (Cluster = contraste)
    output$comp_table <- renderDT({
      global_data$language                     # i18n trigger
      req(shared_rv$pathway_comp_results)
      r <- shared_rv$pathway_comp_results
      t <- rbind(r$go, r$kegg)
      req(t)
      ts_datatable(t, page_length = 15L, filename_base = "pathways_compare")
    })

    # Choix du barcodeplot : un identifiant de voie par entrée (ID conservé
    # séparément du libellé pour le mapping org db).
    observe({
      req(shared_rv$pathway_dir_results)
      t <- .stack_dir_tables(shared_rv$pathway_dir_results)
      if (is.null(t) || !nrow(t)) return()
      labels <- sprintf("%s — %s (%s)", t$ID, t$Term, t$db)
      updateSelectInput(session, "dir_term", choices = stats::setNames(t$ID, labels))
    })

    output$dir_barcode <- renderPlot({
      global_data$language
      res <- shared_rv$pathway_dir_results
      req(res, input$dir_term)
      tryCatch({
        orgpkg <- switch(input$pathway_org %||% "human",
                         human = "org.Hs.eg.db", mouse = "org.Mm.eg.db")
        if (!requireNamespace(orgpkg, quietly = TRUE)) {
          stop(errorCondition(paste0("Package '", orgpkg, "' requis."),
                              class = "pathway_error"))
        }
        res_de <- .active_de_results()
        req(res_de)
        # Rang complet : logFC de TOUTE la table DE, nommé par Entrez
        # (déduplication = même règle que le moteur : min pvalue par Entrez).
        orgdb_obj <- get(orgpkg, envir = asNamespace(orgpkg))
        m <- suppressWarnings(AnnotationDbi::select(
          orgdb_obj, keys = unique(as.character(res_de$gene)),
          keytype = "SYMBOL", column = "ENTREZID"))
        m <- m[!is.na(m$ENTREZID), ]
        p_of <- res_de$pvalue[match(m$SYMBOL, res_de$gene)]
        m <- m[order(m$ENTREZID, p_of), ]
        m <- m[!duplicated(m$ENTREZID), ]
        stats <- stats::setNames(res_de$log2FoldChange[match(m$SYMBOL, res_de$gene)],
                                 m$ENTREZID)
        # Gènes de la voie : GO par GOID, KEGG par PATH (préfixe hsa retiré).
        term_id <- input$dir_term
        keytype <- if (grepl("^GO:", term_id)) "GO" else "PATH"
        if (keytype == "PATH") term_id <- sub("^[a-z]{3}", "", term_id)
        term_map <- suppressWarnings(AnnotationDbi::select(
          orgdb_obj, keys = term_id, keytype = keytype, column = "ENTREZID"))
        term_genes <- unique(term_map$ENTREZID[!is.na(term_map$ENTREZID)])
        plot_pathway_barcode(stats, term_genes,
                             xlab = "log2FC", labels = c("Down", "Up"))
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste(.tr("Erreur:"), conditionMessage(e)),
             col = "red", cex = 0.9)
      })
    })
    output$dl_pathway <- downloadHandler(
      filename = function() paste0("pathways_bulk_", input$pathway_db, "_", Sys.Date(), ".csv"),
      content  = function(file) { req(shared_rv$pathway_results); write.csv(shared_rv$pathway_results, file, row.names = FALSE) }
    )
    # STAT-Q4 — export Excel de la MEME table. Meme source de verite, donc pas de
    # divergence possible entre le CSV et le classeur. L'extension annonce le
    # format reellement ecrit (CSV si openxlsx manque).
    output$dl_pathway_excel <- downloadHandler(
      filename = function() paste0("pathways_bulk_", input$pathway_db, "_", Sys.Date(),
                                   if (requireNamespace("openxlsx", quietly = TRUE)) ".xlsx" else ".csv"),
      content  = function(file) {
        req(shared_rv$pathway_results)
        write_table_excel_or_csv(shared_rv$pathway_results, file)
      }
    )

    # ── GSEA curve (strings translated) ──────────────────────────────────
    output$gsea_curve_ui <- renderUI({
      global_data$language
      gsea_obj <- attr(shared_rv$pathway_results, "gsea_obj")
      if (is.null(gsea_obj)) {
        return(div(class = "alert alert-light", style = "font-size:0.85em;margin:15px;",
                   icon("info-circle"), " ", .tr("Disponible uniquement pour les r\u00e9sultats GSEA \u2014 relancez l'enrichissement en mode GSEA (panneau de gauche).")))
      }
      df <- shared_rv$pathway_results
      choices <- setNames(df$ID, sprintf("%s (NES=%.2f, p.adj=%.1e)", df$Description, df$NES, df$p.adjust))
      tagList(
        fluidRow(
          column(8, selectizeInput(ns("gsea_curve_pathway"), .tr("Pathway (tapez pour rechercher)"),
                                   choices = choices, width = "100%",
                                   options = list(placeholder = .tr("Pathway (tapez pour rechercher)")))),
          column(4, div(style = "margin-top:25px;",
                       downloadButton(ns("dl_gsea_curve_png"), .tr("Export PNG"), class = "btn-sm btn-secondary w-100")))
        ),
        checkboxInput(ns("gsea_curve_pvalue_table"), .tr("Afficher la table p-value sur le graphique"), value = TRUE),
        plotOutput(ns("gsea_curve_plot"), height = "500px")
      )
    })

    .gsea_curve_plot_fn <- function() {
      gsea_obj <- attr(shared_rv$pathway_results, "gsea_obj")
      req(gsea_obj, input$gsea_curve_pathway)
      if (!requireNamespace("enrichplot", quietly = TRUE)) {
        stop(errorCondition(.tr("Package 'enrichplot' requis (BiocManager::install('enrichplot'))."),
                            class = "bulk_pathways_error"))
      }
      enrichplot::gseaplot2(gsea_obj, geneSetID = input$gsea_curve_pathway,
                            title = input$gsea_curve_pathway,
                            pvalue_table = isTRUE(input$gsea_curve_pvalue_table))
    }

    output$gsea_curve_plot <- renderPlot({
      global_data$language  # i18n
      tryCatch(
        .gsea_curve_plot_fn(),
        error = function(e) {
          ggplot() +
            annotate("text", x = 1, y = 1, label = paste(.tr("Erreur:"), conditionMessage(e)), color = "red") +
            theme_void()
        }
      )
    })

    output$dl_gsea_curve_png <- downloadHandler(
      filename = function() paste0("gsea_curve_", input$gsea_curve_pathway, "_", Sys.Date(), ".png"),
      content  = function(file) {
        png(file, width = 9, height = 7, units = "in", res = 300)
        print(.gsea_curve_plot_fn())
        dev.off()
      }
    )

    # ── STAT-S2 : réseau d'enrichissement (emapplot / cnetplot) ─────────────
    output$network_ui <- renderUI({
      global_data$language  # i18n
      if (is.null(attr(shared_rv$pathway_results, "enrich_obj"))) {
        return(div(class = "alert alert-light", style = "font-size:0.85em;margin:15px;",
                   icon("info-circle"), " ",
                   .tr("Disponible uniquement après une analyse de voies — relancez l'enrichissement (panneau de gauche).")))
      }
      tagList(
        fluidRow(
          column(3, radioButtons(ns("network_mode"), .tr("Type de réseau"),
                                 choices = stats::setNames(
                                   c("emap", "cnet"),
                                   c(.tr("Voies ↔ voies (similarité de gènes)"),
                                     .tr("Voies ↔ gènes"))),
                                 inline = TRUE)),
          column(3, numericInput(ns("network_top_n"), .tr("Voies affichées (réseau)"),
                                 value = 30, min = 2, max = 100, step = 1)),
          column(3, div(style = "margin-top:25px;",
                        checkboxInput(ns("network_interactive"),
                                      .tr("Réseau interactif (survol des nœuds)"),
                                      value = FALSE))),
          column(3, div(style = "margin-top:25px;",
                        downloadButton(ns("dl_network_png"), .tr("Export PNG"), class = "btn-sm btn-secondary w-100")))
        ),
        # STAT-S2 V2 : le rendu bascule statique (V1) / plotly interactif —
        # même routeur que pca_interactive (mod_bulk_filter.R) ; défaut =
        # statique, zéro changement de comportement à l'ouverture.
        uiOutput(ns("network_plot_ui"))
      )
    })

    output$network_plot_ui <- renderUI({
      if (isTRUE(input$network_interactive)) {
        plotly::plotlyOutput(ns("network_plot_ly"), height = "600px")
      } else {
        plotOutput(ns("network_plot"), height = "600px")
      }
    })

    .network_plot_fn <- function() {
      req(shared_rv$pathway_results)
      top_n <- input$network_top_n
      if (is.null(top_n) || is.na(top_n)) top_n <- 30
      plot_pathway_network(shared_rv$pathway_results,
                           db_label = input$pathway_db,
                           top_n = top_n, mode = input$network_mode, tr = .tr)
    }

    .network_net_fn <- function() {
      req(shared_rv$pathway_results)
      top_n <- input$network_top_n
      if (is.null(top_n) || is.na(top_n)) top_n <- 30
      build_pathway_network_data(shared_rv$pathway_results,
                                 top_n = top_n, mode = input$network_mode)
    }

    output$network_plot <- renderPlot({
      global_data$language  # i18n
      tryCatch(
        .network_plot_fn(),
        error = function(e) {
          ggplot() +
            annotate("text", x = 1, y = 1, label = paste(.tr("Erreur:"), conditionMessage(e)), color = "red") +
            theme_void()
        }
      )
    })

    output$network_plot_ly <- renderPlotly({
      global_data$language  # i18n
      tryCatch(
        plot_pathway_network_interactive(
          .network_net_fn(),
          title = paste(.tr("Réseau d'enrichissement"), "-", input$pathway_db),
          tr = .tr
        ),
        error = function(e) plotly::plot_ly(type = "scatter", mode = "markers") |>
          plotly::layout(annotations = list(
            text = paste(.tr("Erreur:"), conditionMessage(e)),
            showarrow = FALSE))
      )
    })

    output$dl_network_png <- downloadHandler(
      filename = function() paste0("pathway_network_", input$network_mode, "_", Sys.Date(), ".png"),
      content  = function(file) {
        png(file, width = 10, height = 8, units = "in", res = 300)
        print(.network_plot_fn())
        dev.off()
      }
    )

    # S2c: publish THIS module's drive export route. Same seam shape as the
    # bulk_de one (S2b): a closure over this module's own `shared_rv`, so the
    # route exports the table THIS session produced, never a caller-chosen
    # state. The exporter lives in mod_bulk_pathways_export.R so the offline
    # tests can reach it without a Shiny runtime.
    ts_drive_publish_export(global_data, "bulk_pathways", function() {
      bulk_pathways_export_enrichment_csv(shared_rv, global_data)
    })

  }) # /moduleServer
}
