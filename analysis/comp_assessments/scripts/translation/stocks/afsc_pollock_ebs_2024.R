## Observations ----
years <- 1964:2024
ages <- 1:15
surveys <- c(
  "NMFS bottom-trawl VAST",
  "NMFS acoustic-trawl",

  # The accepted model treats ATS age 1 as a separate recruitment index.
  # Include that source stream here so database_to_tam_obs() imports it;
  # it is combined back into the main acoustic-trawl survey below for the
  # simpler tinyAM age-specific index representation.
  "NMFS acoustic-trawl age-1 index"
)

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  weight_survey = "",
  index_weight_source = "source",
  maturity_reference_year = 1964,
  maturity_multiplier = 0.5,
  surveys = surveys,
  assumptions = source$assumptions
)

## The accepted assessment removes ATS age 1 from the ordinary age-composition
## likelihood and fits it separately as a recruitment index. For tinyAM, retain
## that information but treat it as age 1 of the same acoustic-trawl survey.
## The source-specific 2024 age-1 observation remains absent because it failed
## the accepted model's uncertainty/exclusion criterion during database import.
obs$index$survey[
  obs$index$survey == "NMFS acoustic-trawl age-1 index"
] <- "NMFS acoustic-trawl"

## Approximate survey selectivity/catchability with survey-by-age q blocks.
## Ages 9+ share q, matching the accepted ATS terminal selectivity treatment
## more closely while keeping the tinyAM representation parsimonious.
q_age_block <- cut_ages(obs$index$age, c(1:9, 15))
levels(q_age_block)[levels(q_age_block) == "9-15"] <- "9+"
obs$index$q_age_block <- as.character(q_age_block)
obs$index$q_key <- interaction(
  obs$index$survey,
  obs$index$q_age_block,
  drop = TRUE
)


## Background and comparisons ----

age_plus_group <- 15

comparison_scales <- c(N = 1e-09, recruitment = 1e-06, ssb = 1e-06, biomass_at_age = 1e-06)

comparison_age_groups <- list(N = list(`10+` = 10:15), biomass_at_age = list(`3+` = 3:15))

comparison_aggregates <- "ssb"

comparison_definitions <- list(ssb = list(status = "approximate", definition = "Female SSB at source spawning time versus tinyAM begin-year SSB, ages 1-15",
    reason = "tinyAM does not reproduce the source spawning-time survival adjustment."))

background <- c("### Eastern Bering Sea pollock: accepted 2024 Model 23.0", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|---|---|---|",
    "| Years | Population model years 1964-2024. | Fit 1964-2024. | The accepted model period is retained. |",
    "| Ages | Ages 1-15, with age 15 as the model plus group. | Ages 1-15, with age 15 as the plus group. | The model age range is retained. Published tables sometimes aggregate ages 10-15 as 10+. |",
    "| Recruitment and N | Age-1 recruitment varies annually. Initial ages 2-15 have a shared mean with regularized age deviations; subsequent cohorts are propagated through F and M. | Recruitment follows tinyAM's temporal process. Older ages follow cohort survival with IID abundance deviations, and initial abundance uses exponential survivorship. | IID N deviations provide additional flexibility relative to the accepted model and materially improve convergence and residual behaviour. The source initial-age parameterization is not reproduced. |",
    "| F | Annual fishing mortality is a scalar level multiplied by time-varying fishery selectivity-at-age. Selectivity changes are themselves regularized through time. | Age-specific F follows independent temporal random walks. | This is a compact approximation to annual fishing intensity combined with changing selectivity rather than a direct reproduction of the source parameterization. |",
    "| M | Fixed age-specific values: 0.9 at age 1, 0.45 at age 2, and 0.3 at ages 3-15. | The same supplied age-specific M values are used with no M process. | The accepted numerical M vector is retained directly. |",
    "| Catch | Total fishery biomass is fitted separately from annual age compositions. Composition likelihoods use annual effective weights and therefore give very little influence to extremely rare age cells. The 2024 fishery composition is unavailable. | Reconstruct catch-at-age in numbers using source catch weights and fit observed age-year cells directly. Log-SD follows a quadratic function of age, allowing greater uncertainty at young and old ages and the lowest uncertainty at intermediate ages. | tinyAM uses an age-specific lognormal observation model rather than separate total-catch and composition likelihoods. The age-dependent SD structure prevents sparse tail observations, especially age 1, from dominating the fit. |",
    "| Bottom-trawl survey | VAST biomass index and age composition are fitted with survey selectivity and a full covariance treatment for the biomass series. | Reconstruct an age-specific mid-year index using source survey weights. Catchability is estimated by survey and age block, with ages 9+ sharing q. | tinyAM uses independent lognormal age-specific observations and does not reproduce the VAST biomass covariance or native composition likelihood. |",
    "| Acoustic-trawl survey | The ATS biomass index is fitted with age composition over ages 2-15. Age 1 is removed from that composition and fitted separately as a recruitment index; the 2024 age-1 observation is excluded by the source uncertainty rule. | Combine the accepted age-1 recruitment-index information with the remaining ATS observations as one age-specific acoustic-trawl survey. Catchability is estimated by age, with ages 9+ sharing q. | This preserves the source age-1 information while simplifying the accepted model's likelihood decomposition into a single tinyAM survey representation. The excluded 2024 age-1 observation remains omitted. |",
    "| Survey q | Bottom-trawl and acoustic surveys use survey-specific catchability together with structured, partly time-varying selectivity. ATS selectivity is constant from age 9 onward in the accepted model. | Estimate survey-by-age q, pooling ages 9+ within each survey. | The age blocks provide a simple approximation to survey selectivity while avoiding weakly identified catchability parameters at sparse old ages. |",
    "| Index error | Survey biomass and age-composition information are represented by separate likelihood components with source-specific variance or covariance structures. | Use one log-scale observation-error parameter per survey for the reconstructed age-specific indices. | This is a simpler independent-error approximation and does not reproduce composition sampling weights or within-survey covariance. |",
    "| Weights and maturity | Annual stock, catch, BTS and ATS weights-at-age are supplied. Maturity is fixed and multiplied by 0.5 for female spawning biomass. | Use stock weights for population biomass, source-specific weights for observation conversion, and source maturity multiplied by 0.5. | The accepted biological inputs and female-SSB convention are retained, but source spawning-time survival is not reproduced. |",
    "| SSB | Female spawning biomass is calculated from ages 1-15 using annual stock weights and the fixed maturity schedule. | Female SSB is calculated over the same modeled ages using the translated stock weights and maturity schedule. | The ages and female convention align, but the different spawning-time survival adjustment makes this an approximate comparison. |",
    "", "The translation deliberately simplifies the accepted likelihood architecture while retaining the same model years, ages, fixed natural mortality and core biological inputs. The final tinyAM specification uses IID abundance deviations, age-specific random-walk F, age-dependent catch observation error, and survey-by-age catchability with ages 9+ pooled.",
    "", "The accepted model treats acoustic age 1 as a separate recruitment index. tinyAM instead includes it as age 1 of the acoustic-trawl survey, while retaining the source exclusion of the 2024 age-1 observation. This keeps the observation structure compact without discarding the recruitment information.",
    "", "Historical fishery CPUE and acoustic-vessel-of-opportunity indices remain omitted because they do not have the age-composition information needed for the current age-specific tinyAM translation.",
    "", "Available outputs are compared on explicit common definitions: N ages 1-9 directly and accepted 10+ against the sum of tinyAM ages 10-15; age-1 recruitment directly; and total abundance across ages 1-15. Female SSB is shown as an approximate comparison because source spawning-time survival differs from tinyAM beginning-year SSB. Reported age-3+ biomass is compared with tinyAM biomass summed over ages 3-15. F-at-age is reconstructed from the pinned fitted parameters and source code; uncertainty is not available. Fixed M matches by construction because the accepted age-specific vector is supplied.")

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "exp"),
    F_settings = list(process = "rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~age + I(age^2), fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey, fill_missing = FALSE),
    silent = silent
  )
}
