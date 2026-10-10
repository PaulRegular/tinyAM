## Observations ----
years <- 1947:2025
ages <- 0:8
surveys <- c("HERAS", "IBTS-Q1", "IBTS0", "IBTS-Q3")

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  source$inputs,
  years = years,
  ages = ages,
  index_weight_source = "source",
  surveys = surveys,
  assumptions = source$assumptions
)
obs$weight$M_assumption <- obs$weight$M_assumption + 0.02

obs$catch$sd_block <- cut_ages(obs$catch$age, c(0, 2, 7, 8))
levels(obs$catch$sd_block) <- c("ages_0_1", "ages_2_6", "ages_7_8", "ages_7_8")
heras <- obs$index$survey == "HERAS"
q3 <- obs$index$survey == "IBTS-Q3"
heras_sd <- cut_ages(obs$index$age[heras], c(1:4, 7, 8))
levels(heras_sd) <- c("1", "2", "3", "4_6", "7_8", "7_8")
q3_sd <- cut_ages(obs$index$age[q3], c(0, 2, 5))
levels(q3_sd) <- c("0_1", "2_5")
obs$index$sd_block <- obs$index$survey
obs$index$sd_block[heras] <- paste0("HERAS_", heras_sd)
obs$index$sd_block[q3] <- paste0("IBTS-Q3_", q3_sd)
obs$index$sd_block <- factor(obs$index$sd_block)
heras_q <- cut_ages(obs$index$age[heras], c(1, 3, 8))
levels(heras_q) <- c("1_2", "3_8")
obs$index$q_key <- obs$index$survey
obs$index$q_key[heras] <- paste0("HERAS_", heras_q)
obs$index$q_key[q3] <- paste0("IBTS-Q3_", obs$index$age[q3])
obs$index$q_key <- factor(obs$index$q_key)


model_dat <- tinyAM::prepare_tam(data = obs, years = years, ages = ages, N_settings = list(process = "rw",
    init = "free"), F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:6), M_settings = list(process = "off",
    mu_form = NULL, mu_supplied = ~M_assumption), catch_settings = list(sd_form = ~0 + sd_block, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block, fill_missing = FALSE))
start_par <- tinyAM::make_par(model_dat)
source_N <- .translation_start_surface(source$outputs, "population", "numbers_at_age", years, ages, 1000)
source_F <- .translation_start_surface(source$outputs, "mortality", "fishing_mortality_at_age", years, ages)
start_par$log_r0 <- log(source_N[1, 1])
start_par$log_n0 <- log(source_N[1, -1])
start_par$log_r <- log(source_N[-1, 1])
start_par$log_n <- log(source_N[-1, -1, drop = FALSE])
start_par$log_f <- log(source_F)


## Background and comparisons ----

age_plus_group <- 8

comparison_scales <- c(N = 0.001, recruitment = 0.001, ssb = 0.001, F_bar = 1)

background <- c("### North Sea autumn-spawning herring: accepted 2026 single-fleet assessment", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|------|------|------|",
    "| Years | The population model spans 1947–2026. Catch, weights, maturity and published M inputs end in 2025. | Fit 1947–2025 and compare on those common years. | The accepted 2026 biological inputs are not available in the report tables. |",
    "| Ages | Ages 0–8 winter rings, with age 8 as the plus group. | Ages 0–8, with age 8 as the plus group. | Source age labels are winter-ring classes. |",
    "| N | Recruitment has its own process variance; ages 1–8 share another variance. | Random-walk cohort deviations and free initial older-age abundance. | Variance sharing and initial-state integration differ from FLSAM. |",
    "| F | Ages 0–6 have separate states; ages 7–8 share a state. Innovation variance is shared over ages 0–1, 2–5 and 6–8, with age correlation. | Age-specific random walks through age 8. | tinyAM does not reproduce the source's age-state sharing, innovation groups or age correlation. |",
    "| M | Fixed annual SMS-2023 age surface plus 0.02. | The published 1947–2025 surface plus 0.02, supplied as fixed M. | The numerical input is retained; no M process is fitted. |",
    "| Catch | Direct catch-at-age in thousand fish; closure years 1978–1979 are missing. | Convert to fish and fit positive catch cells with lognormal errors; zero cells are treated as missing. | tinyAM's lognormal catch likelihood does not model exact zero observations. |",
    "| Index | HERAS, IBTS-Q1, IBTS0 and IBTS-Q3 are age-structured. Four larval indices use partial spawning-component likelihoods. | Retain the four age-structured surveys and their q/SD sharing; omit the larval indices. | The larval series are not abundance-at-age observations. Native survey units are unresolved and their scale is absorbed by estimated q; IBTS-Q3's cross-age error correlation is not represented. |",
    "| Weights and maturity | Annual stock and catch weights and annual maturity are tabulated through 2025. | Use the original annual inputs on the fitted period. | This retains the source year-age biology available for 1947–2025. |",
    "| SSB | The accepted source reports SSB using its spawning-time convention. | tinyAM calculates SSB from supplied weights and maturity using its own timing convention. | The source spawning fractions are unavailable, so the definitions may differ. |",
    "", "The assessment's fitted F and N outputs extend to 2026, but 2026 is not included in this tinyAM fit or the common-period comparison. The accepted N and F surfaces are starting values only; tinyAM estimates them freely. Survey timing values are approximate seasonal midpoints because the accepted fleet.txt was unavailable. The retained survey observations stay on their published native scales; estimated catchability absorbs unresolved numerical units, but this does not resolve survey timing.")

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "rw", init = "free"),
    F_settings = list(process = "rw", mu_form = NULL, mean_ages = 2:6),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~0 + sd_block, fill_missing = FALSE),
    index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block, fill_missing = FALSE),
    silent = silent,
    start_par = start_par
  )
}
