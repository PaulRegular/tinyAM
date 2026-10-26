## Observations and biology ----
years <- 2000:2023
ages <- 1:12
inputs <- source$inputs

# Pool fishery conditional age-at-length samples into a fixed empirical key.
# These fishery age samples are biological data, not active source likelihood rows.
key_rows <- inputs[
  inputs$type == "catch" & inputs$measure == "conditional_proportion_at_age", ,
  drop = FALSE
]
lengths <- sort(unique(key_rows$length_bin_lower))
key <- matrix(0, length(lengths), 13, dimnames = list(lengths, 0:12))

for (bin in lengths) {
  z <- key_rows[key_rows$length_bin_lower == bin, , drop = FALSE]
  z$value <- z$value / ave(z$value, z$observation_id, FUN = sum)
  count <- tapply(z$value * z$sample_size, z$age, sum)
  key[as.character(bin), ] <- count[as.character(0:12)] / sum(count)
}

# Convert annual catch-at-length to numbers-at-age using the pooled key.
length_obs <- inputs[
  inputs$type == "catch" & inputs$measure == "proportion_at_length" &
    inputs$year %in% years, , drop = FALSE
]
catch <- lapply(split(length_obs, length_obs$year), function(z) {
  z <- z[order(z$length_bin), , drop = FALSE]
  p <- z$value / sum(z$value)
  covered <- z$length_bin %in% lengths
  weight <- 5.40706e-6 * (z$length_bin + 2.5)^3.19601
  biomass <- inputs$value[
    inputs$type == "catch" & inputs$measure == "total_biomass" &
      inputs$year == z$year[1]
  ]
  stopifnot(length(biomass) == 1L, sum(p * weight) > 0)

  numbers <- biomass * 1000 / sum(p * weight)
  age_p <- as.numeric(p[covered] %*%
    key[as.character(z$length_bin[covered]), , drop = FALSE])
  out <- z[rep(1L, 13L), , drop = FALSE]
  out$age <- 0:12
  out$value <- numbers * age_p
  out$measure <- "numbers_at_age"
  out$basis <- "numbers"
  out$unit <- "fish"
  out$notes <- paste(
    "Translation: pooled supplied-N-weighted empirical fishery age-length key;",
    "uncovered number fraction", signif(1 - sum(p[covered]), 5),
    "left unallocated. Total numbers use original fixed weight-length relation",
    "at 5-cm bin midpoints."
  )
  out
})
catch <- do.call(rbind, catch)

# Retain only active survey marginal-age observations, allocated from VAST total N.
survey <- inputs[
  inputs$type == "index" & inputs$measure == "proportion_at_age" &
    grepl("^Active", inputs$notes) & inputs$year %in% years, ,
  drop = FALSE
]
survey$value <- survey$value / ave(survey$value, survey$observation_id, FUN = sum)
total <- inputs[inputs$type == "index" & inputs$measure == "total_numbers", ]
survey$value <- survey$value * total$value[match(survey$year, total$year)] * 1000
survey$measure <- "numbers_at_age"
survey$basis <- "numbers"
survey$unit <- "fish"
survey$notes <- paste(
  "Translation: active marginal survey number proportions times original VAST",
  "total abundance; no fitted parameter adjusts observations."
)

# Fixed plug-in biology reconstructed from the source's fitted outputs.
outputs <- source$outputs
state_key <- function(z) paste(z$year, z$age)
n <- outputs[outputs$measure == "numbers_at_age" & outputs$year %in% years, ]
w <- outputs[outputs$measure == "weight_at_age" & outputs$year %in% years, ]
w <- w[match(state_key(n), state_key(w)), ]
mat <- plogis(0.132 * ((w$value / 5.40706e-6)^(1 / 3.19601) - 58))
mat[n$age < 1] <- 0

biology <- lapply(split(seq_len(nrow(n)), paste(n$year, pmin(n$age, 12))), function(j) {
  out <- inputs[rep(which(inputs$type == "M")[1], 2L), ]
  out$year <- n$year[j[1]]
  out$year_basis <- "calendar_year"
  out$age <- min(n$age[j])
  out$type <- c("weight", "maturity")
  out$measure <- c("weight_at_age", "maturity_at_age")
  out$value <- c(
    weighted.mean(w$value[j], n$value[j]),
    weighted.mean(mat[j], n$value[j] * w$value[j])
  )
  out$unit <- c("kg/fish", "proportion")
  out$basis <- c("kg_per_fish", "proportion")
  out$source_type <- "derived_source_output"
  out$notes <- paste(
    "Translation only: fitted population weights; maturity proxy from",
    "weight-equivalent length, not official integrated maturity.",
    "Age-12+ biology weighted by official N; no fitted states used",
    "to allocate source observations."
  )
  out
})
biology <- do.call(rbind, biology)

obs <- database_to_tam_obs(
  source$assessment$assessment_id,
  rbind(catch, survey, biology, inputs[inputs$type == "M", ]),
  years = years,
  ages = ages,
  maturity_multiplier = 0.5,
  assumptions = source$assumptions
)
obs$maturity$spawn_time <- 0

## Background and comparisons ----
age_plus_group <- 12
comparison_outputs <- outputs
ssb_age <- n
ssb_age$value <- n$value * w$value * 0.5 * mat
ssb_age$type <- "biomass"
ssb_age$measure <- "mature_biomass_at_age"
ssb_age$unit <- "t"
ssb_age$se <- ssb_age$lwr <- ssb_age$upr <- NA_real_
ssb_age$source_type <- "derived_source_output"
ssb_age$notes <- paste(
  "Common-definition female January biomass using weight-derived maturity proxy;",
  "not official length-integrated female SSB. Thousand fish times kg gives tonnes."
)
comparison_outputs <- rbind(comparison_outputs, ssb_age)

comparison_definitions <- list(
  F = list(
    status = "approximate",
    definition = "F from beginning/midyear official N, pooled at age 12+",
    reason = "Rounded N limits precision; tail F is N-weighted, without covariance."
  ),
  F_bar = list(
    status = "approximate",
    definition = "Population-weighted mean of recovered F",
    reason = "Derived from rounded beginning/midyear N; not native apical F."
  ),
  ssb = list(
    status = "approximate",
    definition = "Common-age female biomass with shared maturity proxy",
    reason = "Original maturity integrates over length; native female SSB and reported SD remain separate."
  ),
  recruitment = list(
    status = "matched",
    definition = "Same-calendar-year age-1 entry from official N",
    reason = "Native recruitment is age 0."
  )
)

background <- c(
  "### Eastern Bering Sea Pacific cod: official 2024 Model 24.1", "", print_sources(source$assessment), "",
  "| Component | Official assessment | tinyAM candidate | Main difference |",
  "|---|---|---|---|",
  "| Coverage | 1977–2024, ages 0–20+, combined sex. | Fit 2000–2023, ages 1–12+. | Translation focuses on active survey age data; 2020 survey composition is unavailable. |",
  "| Catch | Catch tonnes and length compositions; fishery conditional-age observations inactive. | Reconstruct catch-at-age from annual length compositions and a pooled empirical age–length key. | Age allocation, composition weighting and length error are not reproduced by the likelihood. Uncovered length mass remains unallocated. |",
  "| Survey | July 1 VAST numbers index, active marginal survey ages, length-based selectivity, ageing error and scalar q. | Reconstruct age-specific numbers using VAST totals and survey age proportions, sampled at year fraction 0.5. | Native aggregate and composition likelihoods are replaced by age-specific lognormal index observations. |",
  "| Survey q | Estimated scalar catchability with time-varying length-selectivity width. | Rising logistic q-at-age plus a shared annual random walk in log-q (`~ logistic(age) + rw(year)`). | Year RW changes total q level, not age selectivity width; the official model does not estimate an unrestricted annual q process. |",
  "| N and recruitment | Deterministic survival; recruitment steepness fixed 1, log SD 0.6908, and early-regime/bias corrections. | N process off with free initial older-age states; IID age-1 recruitment around an estimated constant mean, with the recruitment log SD now estimated. | Unlike the source, tinyAM estimates rather than fixes recruitment SD; recruitment enters at age 1 instead of age 0, and early-regime/bias treatment differs. |",
  "| F | Annual F with length-based selectivity and temporal selectivity changes. | Age-specific log F has correlated annual random-walk increments (`cor_rw`) across ages. | This flexible process does not reproduce length selectivity or growth-driven changes at age. |",
  "| M | Fixed at 0.3866. | Supplied M with process off. | Equality reflects the input, not independent estimation. |",
  "| Biology / SSB | Estimated growth and length-integrated maturity; female SSB on January 1. | Plug-in fitted stock weights, a length-equivalent maturity proxy and 0.5 female factor. | Age-12+ values use official-N weighting; common-proxy SSB is not native length-integrated female SSB. |",
  "", "The current candidate also estimates the IID recruitment-process SD rather than fixing it at the accepted assessment value of 0.6908. This is additional model flexibility and should not be described as an exact match to the source recruitment regularization.",
  "", "Candidate experiment only: the current correlated-F, year-RW catchability and estimated-recruitment-SD configuration has not been refitted or benchmarked here. Previous dashboard agreement and diagnostics belong to the earlier age-specific-q and independent-F-RW model.",
  "", "The pooled empirical age–length key is fixed before fitting. Survey age-specific observations share a reconstructed total but are treated as independent lognormal errors. Source zero observations remain zero and are omitted from tinyAM's positive lognormal likelihood."
)

## Model ----
fit <- NULL
if (do_fit) {
  fit_stage <- "fit"
  fit_started <- Sys.time()
  fit <- tinyAM::fit_tam(
    data = obs,
    years = years,
    ages = ages,
    N_settings = list(process = "off", init = "free"),
    rec_settings = list(process = "iid", mu_form = ~ 1),
    F_settings = list(process = "cor_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~ 1),
    index_settings = list(
      q_form = ~ logistic(age) + rw(year),
      sd_form = ~ 1
    ),
    silent = silent
  )
}
