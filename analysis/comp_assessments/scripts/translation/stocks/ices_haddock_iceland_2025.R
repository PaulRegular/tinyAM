## Observations ----
years <- 1979:2024
ages <- 1:12
obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  assumptions = source$assumptions
)


## Background and comparisons ----

age_plus_group <- 12

comparison_scales <- c(ssb = 0.001, recruitment = 0.001, N = 0.001, abundance = 0.001, biomass = 0.001,
    biomass_at_age = 0.001)

background <- c("### Icelandic haddock: 2025 benchmark assessment", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|---|---|---|",
    "| Years | Catch and biological inputs cover 1979–2024; population estimates extend to 2025. | Fit 1979–2024. | The 2025 estimate is retained in the reference but lies beyond the detailed input period. |",
    "| Ages | Ages 1–12+, with age 12 as the plus group. | Ages 1–12, with age 12 as the plus group. | The accepted age range is retained. |",
    "| N | SAM estimates annual abundance and age-1 recruitment. | IID abundance process with a free initial age structure. | This is a simpler process than the accepted SAM model. |",
    "| F | Selectivity and fishing mortality vary over time. | AR1 F process around an estimated mean for each age. | Approximates time variation; SAM covariance and parameter sharing are not reproduced. |",
    "| M | Fixed at 0.2 per year for all ages. | Supplied fixed M=0.2. | The accepted assumption is retained. |",
    "| Catch | Commercial catch numbers-at-age and catch weights. | Catch numbers-at-age with one estimated observation-error scale. | tinyAM uses its standard catch likelihood rather than SAM's full likelihood structure. |",
    "| Index | IS-SMB in March and IS-SMH in October, both at age. | Native-scale age indices with series-by-age q and survey-specific observation error. | The published tables do not identify a physical index unit or exact timing fractions. |",
    "| Weights | Stock weights from March survey; catch weights from commercial samples. | Source stock weights for biomass and source catch weights retained in the database. | The input series and units are preserved. |",
    "| Maturity | March survey maturity; pre-1985 vectors use 1985 values. | Source maturity proportions by age and year. | The reported early-year substitution is retained. |",
    "", "For tinyAM, 0.20 and 0.80 are month-midpoint approximations for March and October survey timing. The accepted model's exact timing fractions were not recovered.",
    "", "The source's native SSB applies pre-spawning F and M fractions of 0.4 and 0.3. tinyAM does not apply these fractions. The comparison summary uses a common-definition mature biomass from accepted N and the shared translated weights and maturity; it is not the source's native SSB.",
    "", "The source's reference biomass is defined by fish length (45 cm and larger), which tinyAM cannot match from the available age-based inputs.")

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "iid", init = "free"),
    F_settings = list(process = "ar1", mu_form = ~factor(age)),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + survey, fill_missing = FALSE),
    silent = silent
  )
}
