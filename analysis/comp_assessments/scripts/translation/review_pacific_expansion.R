source("analysis/comp_assessments/R/run_assessment.R")
source("analysis/comp_assessments/R/review_flags.R")
source("analysis/comp_assessments/tests/helper_stock.R")

## Inventory and existing results ----
database <- read_database()
current <- database$assessments[database$assessments$is_current,]
new <- c("afsc_cod_ebs_2024","afsc_cod_ai_2024","afsc_pollock_ai_2024",
  "afsc_atka_mackerel_bsai_2024","afsc_yellowfin_sole_ebs_2024","afsc_northern_rock_sole_ebs_2024",
  "afsc_alaska_plaice_bsai_2024","afsc_arrowtooth_flounder_bsai_2022","afsc_pop_bsai_2024",
  "joint_hake_coastal_2025","dfo_herring_sog_2025","dfo_pop_bc_2023")
stopifnot(length(new)==12,all(new %in% current$assessment_id))
summary <- read.csv("analysis/comp_assessments/results/comparison_summary.csv",check.names=FALSE)
diagnostics <- read.csv("analysis/comp_assessments/results/fit_diagnostics.csv",check.names=FALSE)
source_estimated_m <- c("dfo_cod_2j3kl_2025","dfo_cod_4t4vn_2019","ices_cod_north_sea_2025",
  "afsc_cod_goa_2026","dfo_cod_3pn4rs_2025","dfo_herring_4tvn_spring_2024",
  "ices_whiting_north_sea_2026","afsc_pop_bsai_2024","afsc_pollock_ai_2024",
  "afsc_northern_rock_sole_ebs_2024","joint_hake_coastal_2025","dfo_herring_sog_2025",
  "dfo_pop_bc_2023","afsc_cod_ai_2024")
m_mode <- setNames(rep("unresolved",nrow(current)),current$assessment_id)
for(id in names(m_mode)) {
  path <- file.path("analysis/comp_assessments/scripts/translation/stocks",paste0(id,".R"))
  if(!file.exists(path))next
  stock <- .assessment_stock(read_assessment(id,database),FALSE,TRUE)
  if(!is.null(stock$error))next
  spec <- .test_stock(read_assessment(id,database))
  estimated <- isTRUE(ncol(spec$dat$M_modmat)>0) || length(spec$dat$M_terms)>0L ||
    spec$dat$M_settings$process!="off"
  m_mode[id] <- if(estimated)"estimated" else if(id %in% source_estimated_m)"fixed approximation" else "fixed input"
}
flags <- assessment_review_flags(summary,m_mode)
flag_status <- ifelse(flags$metric=="M",paste(flags$comparison_status,flags$M_treatment,sep="; "),flags$comparison_status)
fit_status <- diagnostics$status[match(flags$assessment_id,diagnostics$assessment_id)]
exploratory <- fit_status!="converged"
flag_status[exploratory] <- paste(flag_status[exploratory],fit_status[exploratory],sep="; exploratory ")
limitations <- c(
  dfo_cod_2j3kl_2025="Reduced fleets and M/observation dependence",
  dfo_cod_4t4vn_2019="Mortality priors and observation likelihood differ",
  ices_cod_north_sea_2025="Substocks pooled; source-fitted biology simplified",
  ices_herring_north_sea_2026="Source observation dependence and sharing simplified",
  afsc_pollock_ebs_2024="Composition likelihood and selectivity constraints omitted",
  afsc_cod_goa_2026="Catch-length conversion; source-estimated M fixed",
  ices_haddock_north_sea_2026="Survey precision and observation dependence differ",
  ices_cod_northeast_arctic_2026="Observation covariance and iterative cannibalism M",
  afsc_pollock_goa_2024="Fleet aggregation and source q/selectivity priors",
  dfo_herring_4tvn_spring_2024="CPUE/acoustic conversion and mortality priors",
  ices_bluewhiting_northeast_atlantic_2026="Age pooling and native covariance omitted",
  ices_herring_norwegian_spring_2025="Native process/sharing configuration unresolved",
  ices_sprat_baltic_2026="Simplified observation covariance and biology treatment",
  ices_plaice_north_sea_2026="Fleet aggregation and F sharing/correlation differ",
  ices_saithe_north_sea_2026="Source sharing and observation covariance differ",
  ices_whiting_north_sea_2026="Source GMRF M replaced by supplied inputs",
  ices_herring_western_baltic_2026="Grouped ages and hockey-stick recruitment approximation",
  ices_norway_pout_north_sea_2026_benchmark="Quarterly source model represented annually",
  nefsc_summer_flounder_2018="Aggregation and failed curvature/uncertainty checks",
  nefsc_atlantic_mackerel_2018="Native observation likelihood and selectivity differ",
  ices_haddock_iceland_2025="Native sharing settings not fully recovered",
  afsc_atka_mackerel_bsai_2024="Source selectivity/precision constraints omitted",
  afsc_pollock_ai_2024="Non-convergence and extreme N/q scale trade-off",
  afsc_alaska_plaice_bsai_2024="Sex pooling and selectivity constraints omitted",
  afsc_yellowfin_sole_ebs_2024="Non-convergence; female subset and source timing differ",
  afsc_northern_rock_sole_ebs_2024="Non-convergence; female subset and joint composition likelihood differ",
  joint_hake_coastal_2025="Source M penalty restored; recruitment/selectivity and likelihood differ",
  afsc_cod_ebs_2024="Pooled empirical catch-age key and selectivity simplification",
  dfo_pop_bc_2023="Source-output blocked; historical population outputs not recovered")
case_name <- function(ids) {
  stock <- current$stock_id[match(ids,current$assessment_id)]
  label <- database$stocks$common_name[match(stock,database$stocks$stock_id)]
  paste0(ifelse(ids %in% new,"New: ","Existing: "),label)
}
metric_name <- c(ssb="SSB",recruitment="Recruitment",abundance="Total N",N="N-at-age",F="F-at-age",M="M")
number <- function(v)ifelse(is.finite(v),formatC(v,digits=3,format="g"),"—")
markdown_table <- function(x) {
  if(!nrow(x))return("None.")
  clean <- function(v)gsub("[\r\n|]+"," ",as.character(v))
  c(paste0("| ",paste(names(x),collapse=" | ")," |"),
    paste0("| ",paste(rep("---",ncol(x)),collapse=" | ")," |"),
    apply(x,1,function(v)paste0("| ",paste(clean(v),collapse=" | ")," |")))
}

## Review report ----
out <- c("# Pacific expansion review","",
  "The database has 35 current detailed stock records, including twelve additional Pacific cases. Nine new cases historically received actual optimization attempts. BC POP is now source-output blocked despite a numerical pass. The current status is four numerical passes, three non-converged fits, one interrupted run and four source blockers. Seven active new-case dashboards and the historical BC trial are cached locally. This falls short of twelve runnable translations; blocked guards are not counted as fitted recipes.","",
  "The initial expansion retained existing stock specifications. The subsequent source/model review updates hake and BC POP mortality treatment, North Sea saithe initial abundance and Baltic sprat spawning time; other existing comparisons remain retained results. New diagnostics and comparisons come from the documented default-first attempts. Source additions were committed before their translation/readiness records and pushed separately. No core model mathematics changed.","",
  "Review thresholds are mean absolute percent difference >20, absolute terminal mean percent difference >20, or trend correlation <0.90. Annual aggregate correlation and mean within-age temporal correlation are retained, not squared. Missing statistics are not passes. These are review prompts, not tuning targets or scientific-equivalence tests.","",
  "## Source and model review","",
  "The twelve new cases were reviewed against their cached source inputs, native prior equations and translation code. Per-stock source reviews retain the unit, timing, sex/area, likelihood and completeness findings. Canonical numerical inputs and official outputs were not changed.","",
  "Hake's restored source log-M penalty improves median absolute SSB difference from 69.6% to 44.8% (r 0.619 to 0.751), total abundance from 71.8% to 42.8% (r 0.785 to 0.904), and age-2 recruitment from 77.8% to 55.2% (r 0.959 to 0.969). Hake M differs by 21.5%, with finite shared log-M effect uncertainty. The historical BC POP M-prior trial also had finite uncertainty, but is now retired from comparison because complete matching coastwide population series have not been recovered. Native bounds and integration conventions remain different.","",
  "EBS cod, Atka and Alaska plaice retain their baselines because alternatives fail checks or worsen other core comparisons. Northern rock sole still fails with verified M/q penalties restored; a random-N initializer passes numerical checks but yields implausible states and is rejected. AI pollock and yellowfin sole remain non-converged; later expensive yellowfin/POP trials are incomplete. The original three source-blocked cases remain blocked, and BC POP is now a fourth.","",
  "The related saithe free-initial-abundance change improves all five comparable core metrics; sprat now applies the source 40% pre-spawning fraction. Sprat's revised SSB definition is not a like-for-like agreement improvement. The full one-off trial log and table are retained in [the review checkpoint](https://github.com/PaulRegular/tinyAM/blob/e7b2ac6/analysis/comp_assessments/PACIFIC_MODEL_REVIEW.md); trial fits and initialization histories remain in the local gitignored cache. The experimental runners and duplicate table are retired from the routine workflow.","",
  "## All flagged core comparisons","",
  markdown_table(data.frame(Case=case_name(flags$assessment_id),Metric=metric_name[flags$metric],
    `Mean absolute %`=number(flags$mean_absolute_percent_difference),
    `Terminal %`=number(flags$terminal_mean_percent_difference),r=number(flags$trend_correlation),
    Status=flag_status,
    Limitation=unname(limitations[flags$assessment_id]),
    check.names=FALSE)),"",
  "Constant M series have no defined temporal correlation. A fixed approximation to source-estimated M is labelled explicitly; equality is not evidence that tinyAM estimated M successfully.","",
  "## New-case numerical outcomes","",
  "SD-report creation is separate from usable uncertainty. AI pollock, yellowfin sole and northern rock sole have failed curvature/covariance checks. Hake passes numerical checks with its restored source mortality information. BC POP's historical trial passed too, but its current diagnostic is a source-output blocker, not an optimizer failure. A numerical pass does not establish scientific agreement.","")
z <- diagnostics[match(new,diagnostics$assessment_id),]
out <- c(out,markdown_table(data.frame(Case=case_name(new),Status=z$status,Code=z$optimizer_code,
  Objective=number(z$objective),`Max gradient`=number(z$max_abs_gradient),
  `SD report created`=z$sdreport_success,`PD Hessian`=z$positive_definite_hessian,check.names=FALSE)),"",
  "Initial expansion fits drove M near zero (hake log-M SE 298.1; BC POP 1924.9). Restoring a shared source log-M penalty gives hake M 0.18335 with effect SE 0.09591. The retired BC POP trial used a positive moment approximation to the source normal prior, giving M 0.04053 with effect SE 0.22890. Native hard bounds and integration conventions remain different. Each model specification uses normal starts first and at most one retry after numerical failure; poor agreement alone never triggers an initialization retry. Atka, Alaska plaice and hake used one conditional source-start optimization; EBS cod and the historical BC POP trial passed from defaults. Both allowed attempts failed numerical checks for AI pollock, yellowfin sole and northern rock sole. The interrupted BSAI POP run lacks a completed attempt history; its stage is not guessed.","",
  "## Existing numerical or source blockers","",
  markdown_table(diagnostics[!diagnostics$assessment_id %in% new &
    (!diagnostics$is_converged %in% TRUE | !diagnostics$positive_definite_hessian %in% TRUE),
    c("assessment_id","status","max_abs_gradient","positive_definite_hessian")]),"",
  "## Unavailable or non-equivalent comparisons","",
  "Each entry below lists missing or incompatible core comparisons. Constant-series correlation is also undefined even when scale comparisons are available. New source-blocked cases and the interrupted POP case have no fitted comparison statistics.","")
missing <- lapply(current$assessment_id,function(id) {
  x<-summary[summary$assessment_id==id & summary$metric %in% names(metric_name),]
  absent<-setdiff(names(metric_name),x$metric)
  bad<-x[!x$comparison_status %in% c("matched","approximate") | x$n==0,]
  description<-if(nrow(bad))paste(metric_name[bad$metric],bad$comparison_status,sep=": ") else character()
  if(length(absent))description<-c(description,paste(metric_name[absent],"no fitted statistics",sep=": "))
  if(!length(description))return(NULL)
  data.frame(Case=case_name(id),Unavailable=paste(description,collapse="; "))
})
out<-c(out,markdown_table(do.call(rbind,missing)),"","## Fixed-M input checks","",
  "This checks aligned values only; it does not credit fixed M as successful estimation. Differences can reflect age pooling or biological approximations and require source-definition review. Missing checks remain unavailable.","")
fixed <- summary[summary$metric=="M" & m_mode[summary$assessment_id] %in% "fixed input",]
out<-c(out,markdown_table(data.frame(Case=case_name(fixed$assessment_id),
  `Mean absolute difference / year`=number(fixed$mean_absolute_difference),
  Check=ifelse(!is.finite(fixed$mean_absolute_difference),"Unavailable",
    ifelse(fixed$mean_absolute_difference<1e-10,"No detected difference in aligned values","Review input matching/aggregation")),check.names=FALSE)),"",
  "## Limits and source blockers","",
  "- SoG herring: aggregate spawning-biomass survey with no observed survey ages; density-dependent M and discrete fleets also differ. No artificial age allocation or catch-only substitute.",
  "- AI cod: public September native fit differs from final November production run; final composition/weight inputs and SSB scaling unresolved.",
  "- Arrowtooth: no verified numerical likelihood compositions/age–length key or survey timing; recovered fitted N/F are not observations.",
  "- BC POP: a follow-up recovered regional biomass draws for 3CD and 5DE, but the public 5ABC object ends in 2017. They cannot form the matching coastwide comparison. N/F/recruitment and the native reproductive-output/mature-biomass definition remain unresolved; the case stays blocked. The source review also records the U.S. West Coast petrale alternative.",
  "- BC petrale: an additional partial 2024 summary/output record preserves 1938–2024 SSB posterior summaries, separate published terminal/reference values, three M scalars and 33 assumptions. Original observations, biology and configuration remain missing. It is an applied but non-current summary record, excluded from the detailed-case batch and blocked before optimization. No replacement was fitted and no agreement statistics are available. Mixture intent is documented; component membership of the 2000-draw subset remains unresolved.",
  "- Other new fits omit source composition covariance, ageing error, normalized/sex/fleet selectivity constraints, and in some cases informative M/q priors. Zero composition categories retain source zero values but tinyAM lognormal observations treat them as missing.",
  "- Every new source file is cached with provenance. Completeness remains partial where likelihood arrays, activation, point conventions, biological inputs or joint posterior outputs are unresolved.","",
  "See the per-stock source reviews, compact result CSVs and gitignored results/cache/<assessment_id>/dashboard.html for diagnostics and visual comparisons. Source-blocked records have no current comparison dashboard or invented optimizer; the prior BC trial remains archived as historical work. CAPABILITY_GAPS.md distinguishes supported features, approximations and source-data blockers.")
writeLines(out,"analysis/comp_assessments/PACIFIC_REVIEW.md",useBytes=TRUE)
