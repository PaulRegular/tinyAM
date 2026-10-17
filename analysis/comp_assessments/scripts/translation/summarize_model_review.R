source("analysis/comp_assessments/R/read_database.R")

## Pinned comparisons and trial records ----
root <- "analysis/comp_assessments"
file <- tempfile(fileext = ".csv")
status <- system2("git", c("show", "86e44bf:analysis/comp_assessments/results/comparison_summary.csv"), stdout = file)
stopifnot(status == 0)
baseline <- read.csv(file)
unlink(file)
ids <- c("afsc_cod_ebs_2024", "afsc_cod_ai_2024", "afsc_pollock_ai_2024",
  "afsc_atka_mackerel_bsai_2024", "afsc_yellowfin_sole_ebs_2024",
  "afsc_northern_rock_sole_ebs_2024", "afsc_alaska_plaice_bsai_2024",
  "afsc_arrowtooth_flounder_bsai_2022", "afsc_pop_bsai_2024",
  "joint_hake_coastal_2025", "dfo_herring_sog_2025", "dfo_pop_bc_2023")
selected <- c(joint_hake_coastal_2025 = "source_M_prior", dfo_pop_bc_2023 = "source_M_prior")
blocked <- c("afsc_cod_ai_2024", "afsc_arrowtooth_flounder_bsai_2022", "dfo_herring_sog_2025")
fixed_M <- c("afsc_cod_ebs_2024", "afsc_atka_mackerel_bsai_2024",
  "afsc_yellowfin_sole_ebs_2024", "afsc_alaska_plaice_bsai_2024")
rows <- list()
attempt_rows <- list()
for (id in ids) {
  directory <- file.path(root, "results/cache/model_review_20261016", id)
  for (trial in list.dirs(directory, recursive = FALSE, full.names = TRUE)) {
    if (basename(trial) == "baseline") next
    specification <- file.path(trial, "specification.rds")
    if (!file.exists(specification)) next
    specification <- readRDS(specification)
    path <- file.path(trial, "result.rds")
    x <- if (file.exists(path)) readRDS(path) else list()
    d <- x$diagnostics
    outcome <- if (length(d$status)) d$status else if (length(x$error)) "fit_error" else "no completed result"
    decision <- if (id %in% names(selected) && basename(trial) == selected[[id]]) "retained" else "not retained"
    attempt_rows[[length(attempt_rows)+1L]] <- data.frame(
      assessment_id=id, candidate=basename(trial), status=outcome,
      optimizer_code=if(length(d$optimizer_code))d$optimizer_code else NA,
      objective=if(length(d$objective))d$objective else NA_real_,
      max_gradient=if(length(d$max_abs_gradient))d$max_abs_gradient else NA_real_,
      elapsed_seconds=if(length(d$elapsed_seconds))d$elapsed_seconds else NA_real_,
      sdreport_created=if(length(d$sdreport_success))d$sdreport_success else NA,
      positive_definite_hessian=if(length(d$positive_definite_hessian))d$positive_definite_hessian else NA,
      decision=decision)
    if (is.null(x$summary)) next
    z <- merge(baseline[baseline$assessment_id==id,], x$summary, by=c("assessment_id","metric"),
               suffixes=c("_before","_after"))
    metrics <- c("ssb","recruitment","abundance","N","F",if(!id %in% fixed_M)"M")
    z <- z[z$metric %in% metrics,]
    rows[[length(rows)+1L]] <- data.frame(assessment_id=id,candidate=basename(trial),
      metric=z$metric,status=outcome,comparison_status=z$comparison_status_after,
      n_before=z$n_before,n_after=z$n_after,
      median_absolute_percent_before=z$median_absolute_percent_difference_before,
      median_absolute_percent_after=z$median_absolute_percent_difference_after,
      correlation_before=z$trend_correlation_before,correlation_after=z$trend_correlation_after,
      decision=decision)
  }
}
comparisons <- do.call(rbind,rows)
attempts <- do.call(rbind,attempt_rows)
write.csv(comparisons,file.path(root,"results/pacific_model_review.csv"),row.names=FALSE,na="")

## Source and model review ----
data_review <- c(
  "VAST is a number index in thousands, not biomass. Active ages and empirical fishery age-length allocation retained; uncovered length mass stays unallocated. Fitted weights/proxy maturity remain explicit.",
  "Final November and preliminary September runs do not match. Original total catch/biomass inputs retained; final likelihood arrays/weights remain unverified.",
  "Original biomass, number proportions and fleet weights retained. Survey weights proxy unavailable population weights; later aggregate-only years are not given artificial ages.",
  "Original kg weights, maturity, number compositions and biomass CVs retained. Actual composition weights are not replaced by raw otolith counts. Whole-million N rounding limits precision.",
  "Joint-sex observed compositions are converted using grams before retaining females. Malformed native tokens use verified observed report echoes, not fitted predictions. Covariates retain native centering.",
  "Joint-sex conversion precedes female restriction; kg weights differ from yellowfin grams. Rounded-derived F stays approximate. Framework confirms the female log-M penalty; male-prior activation is not established.",
  "Joint-sex conversions retain fixed source q calibration and computed growth weights. Only three fishery age years are available; absent compositions remain missing.",
  "ATF-only catch, sex-specific fixed biology and active surveys retained. Numerical length/age likelihood inputs, observed age-length key and survey times remain unverified.",
  "Original grams, ages 40+, two surveys and active biomass likelihoods retained. Fitted maturity is a labelled plug-in; source M density includes its documented bias adjustment.",
  "Native percentage compositions normalize jointly; original weights/fecundity, age-2 index timing and posterior-median definitions retained. No fitted M is used to change observations.",
  "Native survey is aggregate spawning biomass with no observed survey ages. Five-OM ensemble, model-year timing, both-sex SSB and posterior point conventions remain distinct.",
  "Original sex/area biomass conversions and zero-evaluation fixed biology retained. Catch pools only complete three-area years. Historical official N/F/R/SSB and native reproductive-output definition remain unresolved.")
assumptions <- c(
  "Test shared annual F, correlated F and catch SD. Retain the baseline: the shared-level alternative improves F/SSB but weakens N/recruitment trends. Its dashboard remains available as an alternative.",
  "Keep the source-data guard. Do not mix preliminary fitted states with final observations.",
  "Test common annual F, correlated F and the verified M log-penalty. Fixed-q normalization, BH settings and source population weights remain limitations.",
  "Keep baseline: all tested F/precision alternatives weaken important comparisons. Joint composition likelihood, q prior and source selectivity penalties remain different.",
  "Test correlated/common-level F and age-dependent catch SD. Failed numerical checks remain failures; different biomass/composition reference times and source composition likelihood remain material.",
  "Test correlated/common-level F, catch SD, random initial N and verified female M/q penalties. The joint M/q trial still fails (gradient 3.22e10, non-PD Hessian). Random N passes numerical checks but produces implausible abundance/q and is rejected. Source joint composition likelihood remains different.",
  "Test correlated/common-level F and catch SD; preserve the source fixed-M/q meanings. Sex pooling and omitted lengths constrain interpretation.",
  "Keep the source-data guard. Fitted N/F cannot manufacture catch or survey observations.",
  "Test common/correlated F and the source M density. Incomplete/interrupted results are not numerical successes.",
  "Retain one shared log-M effect with native fixed center/SD. M is estimated and its density counted once. A BH alternative gives essentially the same comparisons and a weakly estimated curve, so constant recruitment is retained. Laplace integration differs from native Bayesian reporting; bounds are omitted.",
  "Keep the observation-model guard. A fishery-only fit or fabricated survey composition would not represent the source abundance information.",
  "Retain a moment-matched positive approximation to the native normal M prior. It is not the same density or bound treatment. Only M has a usable historical official comparison.")
decisions <- vapply(ids,function(id) {
  if(id %in% names(selected))paste("Adopt",selected[[id]]) else if(id %in% blocked)"Blocked; no fit invented" else "Retain baseline; alternatives not adopted"
},character(1))
inventory <- data.frame(assessment_id=ids,data_review=data_review,
  assumptions_and_fit_review=assumptions,decision=decisions)
markdown_table <- function(x) {
  clean <- function(z)gsub("[\r\n|]+"," ",as.character(z))
  c(paste0("| ",paste(names(x),collapse=" | ")," |"),
    paste0("| ",paste(rep("---",ncol(x)),collapse=" | ")," |"),
    apply(x,1,function(z)paste0("| ",paste(clean(z),collapse=" | ")," |")))
}
number <- function(x)ifelse(is.finite(x),formatC(x,digits=3,format="g"),"—")
out <- c("# Review of the twelve new Pacific cases","",
  "Baseline: commit 86e44bf. This review checks the canonical inputs, importers, source reviews, translation code and numerical fits. Native prior equations were revisited directly for the M trials. It is not a claim of exhaustive independent verification of every published cell. Canonical observations are not tuned to match fitted outputs.","",
  "Each trial changes a documented assumption or precision treatment while retaining observations, years and ages. Normal package starts run first; one existing source-start construction is allowed after numerical failure. All failed and incomplete trials remain visible. Priors use external source assumptions, not fitted official M values.","",
  "The rejected northern rock sole q-prior trial expresses the index in units relative to the fixed source prior center 1.5, retaining original values in source_obs. Its log-index residuals are mathematically unchanged. This is an explicit unit conversion, not an adjustment using fitted q; the retained recipe and canonical observations do not change.","",
  "Median differences below are median absolute percent differences, with the existing annual/mean-within-age correlation definitions. No agreement threshold is optimized. Retained changes can have trade-offs, which are reported rather than hidden.","",
  "## Case-by-case review","",markdown_table(inventory),"",
  "## Trial diagnostics","",markdown_table(attempts),"",
  "## Retained changes: every comparable core metric","",
  "Hake's M treatment improves all available core scale/trend comparisons. BC POP improves M only; other historical core comparisons remain unavailable or non-equivalent. The EBS cod shared-level F alternative improves F/SSB but weakens N/recruitment trends, so it is not adopted automatically.","")
keep <- comparisons[comparisons$decision=="retained" & comparisons$n_after>0,]
out <- c(out,markdown_table(data.frame(Case=keep$assessment_id,Metric=keep$metric,
  `Median absolute % before`=number(keep$median_absolute_percent_before),
  `After`=number(keep$median_absolute_percent_after),
  `r before`=number(keep$correlation_before),`r after`=number(keep$correlation_after),
  Status=keep$comparison_status,check.names=FALSE)),"",
  "M prior implementation: the existing fixed-SD Gaussian mean effect represents one shared log-M level, with one density contribution across all years/ages. Hake matches the source quadratic log-M penalty; BC POP matches the mean and variance of its natural-scale normal prior with a lognormal approximation. Neither fixes M to an official estimate. Native hard bounds and integration/point-estimate conventions remain different. Mean-only M uncertainty is available for the shared log-M effect; the current full M surface is reported without an automatically propagated interval.","",
  "The full baseline/candidate comparison table is results/pacific_model_review.csv. Fits, exact specifications, initialization histories and rejected trials remain in gitignored results/cache/model_review_20261016/. Selected dashboards are refreshed in the usual per-stock cache. Missing statistics are not passes; summaries of failed fits are exploratory.")
out <- c(out,"","## Validation and boundaries","",
  "Canonical numerical inputs and outputs are unchanged. The northern rock sole source note now distinguishes the verified female M penalty from unresolved male-prior activation. Core package mathematics, APIs and likelihood code are unchanged.","",
  "The normal-start/conditional-retry policy applies separately to each documented model specification. A failed initialization does not justify changing observations. Ten-minute wall limits were used for expensive later trials; missing completed results remain incomplete, not convergence failures inferred from their last gradient.","",
  "The full package suite and all assessment-analysis files were run. Current validation results and existing analysis-fixture failures are recorded in PACIFIC_EXPANSION.md. Targeted tests cover prior scales, one shared density contribution, unchanged source observations, conditional starts, common-definition SSB and projection-safe alternative F design.")
current <- read.csv(file.path(root,"results/comparison_summary.csv"))
existing <- merge(baseline[baseline$assessment_id=="ices_saithe_north_sea_2026",],
  current[current$assessment_id=="ices_saithe_north_sea_2026",],
  by=c("assessment_id","metric"),suffixes=c("_before","_after"))
existing <- existing[existing$metric %in% c("ssb","recruitment","abundance","N","F"),]
out <- c(out,"","## Related existing-case corrections","",
  "North Sea saithe now estimates free initial older-age abundance, with normal package starts. Every comparable core median difference and trend improves; observation conversion and plus-group comparisons are unchanged.","",
  markdown_table(data.frame(Metric=existing$metric,
    `Median absolute % before`=number(existing$median_absolute_percent_difference_before),
    After=number(existing$median_absolute_percent_difference_after),
    `r before`=number(existing$trend_correlation_before),
    `r after`=number(existing$trend_correlation_after),check.names=FALSE)),"",
  "Baltic sprat now applies the source 40% pre-spawning fraction for F and M. The derived official SSB uses official N/F/M and shared biology, with no reconstructed uncertainty. This corrects the comparison definition; its statistics are not presented as a like-for-like improvement over the previous January-1 comparison. Native reported SSB remains separate.")
writeLines(out,file.path(root,"PACIFIC_MODEL_REVIEW.md"),useBytes=TRUE)
