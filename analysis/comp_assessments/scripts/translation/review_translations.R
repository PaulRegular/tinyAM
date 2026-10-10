# Source-led candidates; baseline recipes are pinned before the workflow refactor.
review_revision <- '29068994feee4b2e078ac911e92e0c4fe1d200af'
review_candidates <- list(
  dfo_cod_2j3kl_2025 = list(
    bh_iid = list(N_settings = list(rec_form = ~ bh(ssb) + iid(year))),
    bh_ar1 = list(N_settings = list(rec_form = ~ bh(ssb) + ar1(year)))
  ),
  afsc_cod_goa_2026 = list(
    iid_recruitment = list(N_settings = list(rec_form = ~ iid(year, sd = 0.44))),
    logistic_q = list(index_settings = list(q_form = ~ logistic(age))),
    catch_age_sd = list(catch_settings = list(sd_form = ~ age + I(age^2)))
  ),
  afsc_pollock_ebs_2024 = list(
    shared_f_rw = list(F_settings = list(process = 'iid', mu_form = ~ factor(age) + rw(year)))
  ),
  afsc_pollock_goa_2024 = list(
    iid_recruitment = list(N_settings = list(rec_form = ~ iid(year, sd = 1.3))),
    adfg_q_rw = list(index_settings = list(q_form = ~ 0 + q_key + environmental_effect + rw(year, by = adfg_effect)))
  ),
  dfo_cod_4t4vn_2019 = list(
    ar1_recruitment = list(N_settings = list(rec_form = ~ ar1(year, sd = 0.5))),
    fixed_sd_m_rw_iid_warm = list(M_settings = list(process = 'off',
      mu_form = ~ 0 + M_formula_group + rw(year, by = M_formula_group, sd = 0.075),
      mu_supplied = ~ M_prior_mean)),
    fixed_sd_m_rw = list(M_settings = list(process = 'off',
      mu_form = ~ 0 + M_formula_group + rw(year, by = M_formula_group, sd = 0.075),
      mu_supplied = ~ M_prior_mean),
      warm_start_settings = list(M_settings = list(process = 'iid', mu_form = NULL)))
  ),
  dfo_herring_4tvn_spring_2024 = list(
    cpue_q_rw = list(index_settings = list(q_form = ~ 0 + mono(age, by = survey) + rw(year, by = cpue_effect))),
    fixed_sd_m_rw = list(M_settings = list(process = 'off',
      mu_form = ~ 0 + M_formula_group + rw(year, by = M_formula_group, sd = 0.075),
      mu_supplied = ~ M_process_center))
  ),
  nefsc_atlantic_mackerel_2018 = list(
    shared_f_rw = list(F_settings = list(process = 'iid', mu_form = ~ factor(age) + rw(year)))
  )
)

# These are source checks, not automatic model-selection rules.
review_inventory <- c(
  afsc_cod_goa_2026 = 'Fixed BH steepness one removes SSB dependence; test fixed-SD IID recruitment and a simpler rising survey curve. Conditional catch residual spread differs between young and middle ages: test quadratic log-SD separately. Only approximate SSB is comparable.',
  afsc_pollock_ebs_2024 = 'Test a shared temporal F mean with IID residuals. This does not reproduce correlated fishery selectivity or the composition likelihood.',
  afsc_pollock_goa_2024 = 'Test source fixed-SD IID recruitment and annual ADF&G q variation. Preserve the logit q link and current age sharing.',
  dfo_cod_2j3kl_2025 = 'Test BH with IID/AR1 residuals at modeled age two and lag two. The source curve is age zero and same-year SSB; survivor recruitment is an approximation.',
  dfo_cod_4t4vn_2019 = 'Test recruitment AR1 with source SD and, separately, age-blocked M mean RWs with fixed source increment SD. Source recruitment-rate times SSB and initial-M priors remain different.',
  dfo_herring_4tvn_spring_2024 = 'Test CPUE-only annual q RW and, separately, age-blocked M mean RWs with fixed source increment SD. Source aggregate CPUE, power catchability and initial-M priors remain different.',
  ices_bluewhiting_northeast_atlantic_2026 = 'Retain RW recruitment and source q sharing. Age-correlated F increments and observation errors require capabilities outside this review.',
  ices_cod_north_sea_2025 = 'Retain RW recruitment and substock q groups. Substock scaling, correlated F increments and M GMRF are not equivalent to extra independent formula effects.',
  ices_cod_northeast_arctic_2026 = 'Retain RW recruitment, independent F RWs and source q groups. Observation correlations and age-specific F innovation SDs remain different.',
  ices_haddock_iceland_2025 = 'Retain flexible survey-age q. Exact native transition and sharing settings and index units remain unresolved; no source-supported curve candidate.',
  ices_haddock_north_sea_2026 = 'Retain RW recruitment and source q/error groups. Density-dependent q powers and age-correlated F increments are not logistic catchability.',
  ices_herring_north_sea_2026 = 'Retain RW recruitment and source q sharing. Larval spawning-component observations remain excluded; extra q effects do not supply that mapping.',
  ices_herring_norwegian_spring_2025 = 'Retain baseline; exact process and q-sharing configuration is unresolved. Do not guess a logistic curve or temporal effect.',
  ices_herring_western_baltic_2026 = 'Retain baseline RW approximation. Hockey-stick recruitment is not BH/Ricker; source q keys are already used.',
  ices_norway_pout_north_sea_2026_benchmark = 'Retain annual approximation as non-converged. Quarterly recruitment and mortality cannot be recovered by adding annual formula effects.',
  ices_plaice_north_sea_2026 = 'Retain RW recruitment and flexible survey-age q. Printed sharing keys are not fully mapped to individual surveys; do not invent their assignments.',
  ices_saithe_north_sea_2026 = 'Retain RW recruitment and source q sharing. Shared terminal F states and correlated innovations are distinct from mean-formula effects.',
  ices_sprat_baltic_2026 = 'Retain RW recruitment and current survey groups. Source age-zero index is already shifted to age one; spawning survival and terminal F sharing remain different.',
  ices_whiting_north_sea_2026 = 'Retain RW recruitment, N and F and source survey-age q. Observation correlations and M GMRF remain distinct unsupported likelihood structures.',
  nefsc_atlantic_mackerel_2018 = 'Test a shared temporal F mean with IID residuals to approximate time-constant fishery selectivity. Logistic q is not fishery selectivity; the aggregate egg index lacks an age allocation.',
  nefsc_summer_flounder_2018 = 'Retain non-converged baseline and flag it. Missing source surveys and annual biological inputs prevent a complete comparison; no undocumented rescue tuning.',
  dfo_cod_3pn4rs_2025 = 'Blocked: annual stock weights and fitted maturity surfaces remain unrecovered, despite framework and supporting biological-source checks.',
  nefsc_haddock_georges_bank_2026 = 'Blocked: current model-grid observations and fitted native files have not been recovered; earlier assessment inputs cannot substitute.'
)

review_decisions <- c(
  'dfo_cod_2j3kl_2025/bh_iid' = 'rejected_numerical',
  'dfo_cod_2j3kl_2025/bh_ar1' = 'rejected_numerical',
  'afsc_cod_goa_2026/iid_recruitment' = 'retained_baseline_tradeoff',
  'afsc_cod_goa_2026/logistic_q' = 'retained_baseline_tradeoff',
  'afsc_cod_goa_2026/catch_age_sd' = 'rejected_numerical',
  'afsc_pollock_ebs_2024/shared_f_rw' = 'retained_baseline_worse_agreement',
  'afsc_pollock_goa_2024/iid_recruitment' = 'retained_baseline_tradeoff',
  'afsc_pollock_goa_2024/adfg_q_rw' = 'retained_baseline_tradeoff',
  'dfo_cod_4t4vn_2019/ar1_recruitment' = 'retained_baseline_worse_agreement',
  'dfo_cod_4t4vn_2019/fixed_sd_m_rw_iid_warm' = 'rejected_numerical',
  'dfo_cod_4t4vn_2019/fixed_sd_m_rw' = 'retained_baseline_tradeoff',
  'dfo_herring_4tvn_spring_2024/cpue_q_rw' = 'retained_baseline_tradeoff',
  'dfo_herring_4tvn_spring_2024/fixed_sd_m_rw' = 'retained_baseline_worse_agreement',
  'nefsc_atlantic_mackerel_2018/shared_f_rw' = 'retained_baseline_worse_agreement'
)
review_decision_notes <- c(
  'dfo_cod_2j3kl_2025/bh_iid' = 'Non-positive-definite Hessian; worse scale agreement; correlated recruitment innovations.',
  'dfo_cod_2j3kl_2025/bh_ar1' = 'Non-positive-definite Hessian; modestly improved trends but worse scale agreement.',
  'afsc_cod_goa_2026/iid_recruitment' = 'Approximate SSB average error improves; terminal error and trend worsen. No matched N/F/recruitment evidence.',
  'afsc_cod_goa_2026/logistic_q' = 'Only a small improvement in approximate SSB; survey residual age-patterns remain and trend weakens. Source length-based selectivity is not this curve.',
  'afsc_cod_goa_2026/catch_age_sd' = 'Hessian fails and abundance scale diverges; do not rescue using undocumented constraints.',
  'afsc_pollock_ebs_2024/shared_f_rw' = 'N/F/recruitment scale and trends worsen; weak process-SD interval and strong fixed-parameter correlations.',
  'afsc_pollock_goa_2024/iid_recruitment' = 'N/recruitment improve but matched SSB scale and trend worsen; q boundary and residual advisories remain.',
  'afsc_pollock_goa_2024/adfg_q_rw' = 'Tiny N/recruitment gains but worse average SSB and trends; adds a parameter-correlation advisory.',
  'dfo_cod_4t4vn_2019/ar1_recruitment' = 'Recruitment error increases from 36% to 87%; N and recruitment trends weaken despite improved F scale.',
  'dfo_cod_4t4vn_2019/fixed_sd_m_rw_iid_warm' = 'Mean RW plus preliminary IID M fails its Hessian check; no final fit attempted. Also test the unchanged baseline preliminary fit.',
  'dfo_cod_4t4vn_2019/fixed_sd_m_rw' = 'With the baseline preliminary fit, all trends and N/SSB scale improve, but F scale and terminal recruitment worsen. Retain baseline and flag this alternative for review.',
  'dfo_herring_4tvn_spring_2024/cpue_q_rw' = 'Trends improve but SSB error increases from 25% to 67%, F from 63% to 192%; 15 active bounds.',
  'dfo_herring_4tvn_spring_2024/fixed_sd_m_rw' = 'Numerical checks pass but abundance scale diverges; AR1 persistence and fixed-parameter correlation advisories. Initial-M priors are not reproduced.',
  'nefsc_atlantic_mackerel_2018/shared_f_rw' = 'F trend improves but N/recruitment/SSB scale and trends worsen; terminal SSB difference doubles.'
)

# Confirmed source features, counted once per assessment. Unknown is not absent.
capability_evidence <- list(
  correlated_F_increments = c('ices_cod_north_sea_2025', 'ices_herring_north_sea_2026',
    'ices_haddock_north_sea_2026', 'ices_bluewhiting_northeast_atlantic_2026',
    'ices_sprat_baltic_2026', 'ices_plaice_north_sea_2026', 'ices_saithe_north_sea_2026',
    'ices_whiting_north_sea_2026', 'ices_herring_western_baltic_2026', 'ices_norway_pout_north_sea_2026_benchmark'),
  correlated_observations = c('ices_cod_north_sea_2025', 'ices_herring_north_sea_2026',
    'ices_cod_northeast_arctic_2026', 'ices_bluewhiting_northeast_atlantic_2026',
    'ices_plaice_north_sea_2026', 'ices_saithe_north_sea_2026',
    'ices_herring_western_baltic_2026', 'ices_whiting_north_sea_2026', 'afsc_pollock_ebs_2024'),
  latent_SD_groups = c('ices_cod_north_sea_2025', 'ices_herring_north_sea_2026',
    'ices_haddock_north_sea_2026', 'ices_cod_northeast_arctic_2026',
    'ices_plaice_north_sea_2026', 'ices_saithe_north_sea_2026',
    'ices_whiting_north_sea_2026', 'ices_norway_pout_north_sea_2026_benchmark'),
  composition_likelihoods = c('dfo_cod_2j3kl_2025', 'dfo_cod_4t4vn_2019',
    'dfo_herring_4tvn_spring_2024', 'afsc_pollock_goa_2024', 'afsc_pollock_ebs_2024',
    'afsc_cod_goa_2026', 'nefsc_summer_flounder_2018', 'nefsc_atlantic_mackerel_2018', 'ices_cod_north_sea_2025'),
  fleet_selectivity = c('afsc_cod_goa_2026', 'afsc_pollock_ebs_2024',
    'afsc_pollock_goa_2024', 'dfo_cod_4t4vn_2019', 'dfo_herring_4tvn_spring_2024',
    'nefsc_atlantic_mackerel_2018', 'nefsc_summer_flounder_2018', 'ices_cod_north_sea_2025'),
  spawning_survival = c('afsc_pollock_ebs_2024', 'afsc_pollock_goa_2024',
    'ices_sprat_baltic_2026', 'ices_herring_western_baltic_2026',
    'ices_haddock_iceland_2025', 'dfo_herring_4tvn_spring_2024'),
  priors = c('dfo_cod_4t4vn_2019', 'dfo_herring_4tvn_spring_2024', 'afsc_pollock_goa_2024'),
  other_recruitment_curves = c('ices_herring_western_baltic_2026', 'dfo_cod_4t4vn_2019')
)

.review_recipe <- function(id, database) {
  path <- paste0('analysis/comp_assessments/scripts/translation/stocks/', id, '.R')
  text <- system2('git', c('show', paste0(review_revision, ':', path)), stdout = TRUE)
  if (!is.null(attr(text, 'status'))) cli::cli_abort('Cannot read the pinned baseline recipe for {id}.')
  e <- new.env(parent = environment())
  eval(parse(text = text), e)
  e$translate_stock(read_assessment(id, database))
}

.review_fit <- function(id, recipe, changes, database) {
  z <- recipe
  if (!is.null(changes$warm_start_settings)) z$warm_start_settings <- changes$warm_start_settings
  changes$warm_start_settings <- NULL
  z$settings <- utils::modifyList(z$settings, changes)
  z$obs$index$adfg_effect <- as.numeric(z$obs$index$survey == 'ADF&G crab/groundfish trawl')
  z$obs$index$cpue_effect <- as.numeric(grepl('CPUE', z$obs$index$survey, ignore.case = TRUE))
  if (id == 'dfo_cod_4t4vn_2019') {
    z$obs$weight$M_formula_group <- factor(ifelse(z$obs$weight$age < 5, '3-4',
      ifelse(z$obs$weight$age < 9, '5-8', '9+')))
  }
  if (id == 'dfo_herring_4tvn_spring_2024') {
    z$obs$weight$M_formula_group <- factor(ifelse(z$obs$weight$age < 7, '2-6', '7+'))
  }
  warnings <- character()
  started <- Sys.time()
  out <- tryCatch(withCallingHandlers({
    args <- c(list(data = z$obs, years = z$years, ages = z$ages, silent = TRUE), z$settings)
    if (!is.null(z$start_par)) args$start_par <- z$start_par
    if (!is.null(z$warm_start_settings)) {
      warm_args <- c(list(data = z$obs, years = z$years, ages = z$ages, silent = TRUE),
                     utils::modifyList(z$settings, z$warm_start_settings))
      warm <- do.call(tinyAM::fit_tam, warm_args)
      if (!isTRUE(warm$is_converged)) cli::cli_abort('The preliminary fit did not converge.')
      args$start_par <- as.list(warm$sdrep, 'Estimate')
    }
    fit <- do.call(tinyAM::fit_tam, args)
    fit <- .assessment_catch_reporting(fit, z$catch_reporting)
    source <- read_assessment(id, database)
    outputs <- if (is.null(z$comparison_outputs)) database$outputs else z$comparison_outputs
    ref <- database_to_tam_ref(
      id, outputs, obs = z$obs, years = z$years, ages = z$ages,
      terminal_year = source$assessment$terminal_year[[1]], age_plus_group = z$age_plus_group,
      comparison_scales = z$comparison_scales, template = fit, assumptions = source$assumptions,
      comparison_aggregates = z$comparison_aggregates, comparison_age_groups = z$comparison_age_groups,
      comparison_definitions = z$comparison_definitions
    )
    list(fit = fit, ref = ref, recipe = z,
         diagnostics = .assessment_diagnostics(id, database,
           if (fit$is_converged) 'converged' else 'not_converged', fit),
         summary = .assessment_comparison_summary(.assessment_percent_differences(ref), id))
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart('muffleWarning')
  }), error = function(e) list(error = conditionMessage(e)))
  out$warnings <- warnings
  out$recipe <- z
  out$elapsed <- as.numeric(difftime(Sys.time(), started, units = 'secs'))
  out
}

review_translations <- function(database = read_database(), assessment_ids = NULL,
                                output_file = file.path(.assessment_root, 'results/translation_review.csv')) {
  if (is.null(assessment_ids)) {
    assessment_ids <- database$assessments$assessment_id[database$assessments$is_current %in% TRUE]
  }
  inputs_unchanged <- identical(system2('git', c('diff', '--name-only', review_revision, '--',
    'analysis/comp_assessments/database'), stdout = TRUE), character())
  cache <- file.path(.assessment_root, 'results/cache/translation_review')
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  attempts <- list()
  for (id in assessment_ids) {
    has_recipe <- file.exists(file.path(.assessment_root, 'scripts/translation/stocks', paste0(id, '.R')))
    recipe <- if (has_recipe) .review_recipe(id, database) else NULL
    options <- c(list(baseline = list()), review_candidates[[id]])
    for (candidate in names(options)) {
      path <- file.path(cache, paste0(id, '_', candidate, '.rds'))
      # Reuse recorded attempts only for unchanged inputs and candidate settings.
      previous <- if (file.exists(path)) readRDS(path) else NULL
      changes <- options[[candidate]]
      warm_settings <- if (is.null(changes$warm_start_settings)) recipe$warm_start_settings else changes$warm_start_settings
      changes$warm_start_settings <- NULL
      same_settings <- !is.null(previous$recipe) && isTRUE(all.equal(
        previous$recipe$settings, utils::modifyList(recipe$settings, changes),
        check.environment = FALSE)) && isTRUE(all.equal(previous$recipe$warm_start_settings,
        warm_settings, check.environment = FALSE))
      original_baseline <- candidate == 'baseline' && inputs_unchanged &&
        identical(previous$diagnostics$database_revision, review_revision)
      if (!has_recipe) {
        out <- list(error = review_inventory[[id]], elapsed = NA_real_)
      } else if (!is.null(previous) && same_settings && (original_baseline ||
          (identical(previous$review_revision, review_revision) &&
          (inputs_unchanged || identical(previous$database_revision, database$commit))))) {
        out <- previous
      } else {
        cat(id, candidate, '\n'); flush.console()
        out <- .review_fit(id, recipe, options[[candidate]], database)
        out$review_revision <- review_revision
        out$database_revision <- database$commit
        saveRDS(out, path)
      }
      s <- data.frame(assessment_id = id, candidate = candidate)
      for (metric in c('ssb', 'recruitment', 'N', 'F')) {
        z <- if (is.null(out$summary)) NULL else out$summary[out$summary$metric == metric, , drop = FALSE]
        fields <- c(status = 'comparison_status', percent = 'mean_absolute_percent_difference',
                    trend = 'trend_correlation', terminal_percent = 'terminal_mean_percent_difference')
        for (field in names(fields)) {
          s[[paste(metric, field, sep = '_')]] <- if (!is.null(z) && nrow(z)) z[[fields[[field]]]][1] else NA
        }
      }
      d <- out$diagnostics
      s$baseline_revision <- review_revision
      s$database_revision <- if (!is.null(out$diagnostics)) out$diagnostics$database_revision else database$commit
      s$is_converged <- if (is.null(d)) FALSE else d$is_converged
      s$optimizer_code <- if (is.null(d)) NA_integer_ else d$optimizer_code
      s$max_abs_gradient <- if (is.null(d)) NA_real_ else d$max_abs_gradient
      s$pdHess <- if (is.null(d)) FALSE else d$positive_definite_hessian
      s$warnings <- paste(unique(out$warnings), collapse = ' | ')
      s$error <- if (is.null(out$error)) '' else out$error
      s$elapsed_seconds <- out$elapsed
      checks <- if (inherits(out$fit, 'tam_fit')) tinyAM::check_tam(out$fit) else NULL
      s$advisories <- if (is.null(checks)) '' else paste(unique(checks$advisories$detail), collapse = ' | ')
      s$review_note <- review_inventory[[id]]
      key <- paste(id, candidate, sep = '/')
      s$decision <- if (!has_recipe) 'blocked' else if (key %in% names(review_decisions))
        review_decisions[[key]] else if (candidate == 'baseline')
          if (isTRUE(s$is_converged)) 'retained_baseline' else 'retained_not_converged' else 'pending_review'
      s$decision_note <- if (key %in% names(review_decision_notes)) review_decision_notes[[key]] else ''
      attempts[[paste(id, candidate)]] <- s
      table <- do.call(rbind, attempts)
      write.csv(table, output_file, row.names = FALSE, na = '')
    }
  }
  invisible(table)
}
