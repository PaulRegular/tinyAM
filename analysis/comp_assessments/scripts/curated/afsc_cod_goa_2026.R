## Reconstruct the age-structured observations from official length data ----

# Pool conditional age-at-length samples by source fleet/survey and length bin.
.goa_age_length_key <- function(inputs, years) {
  age <- inputs[
    inputs$type %in% c("catch", "index") &
      inputs$measure == "conditional_proportion_at_age" &
      inputs$year %in% years, , drop = FALSE
  ]
  if (!nrow(age) || anyNA(age$sample_size) || any(age$sample_size <= 0)) {
    stop("GOA age-at-length records need positive supplied sample sizes.", call. = FALSE)
  }
  if (any(abs(aggregate(value ~ observation_id, age, sum)$value - 1) > 0.001)) {
    stop("A GOA conditional age composition does not sum to one.", call. = FALSE)
  }

  dims <- c("type", "fleet", "survey", "length_bin_lower", "length_bin_upper")
  age$weighted <- age$value * age$sample_size
  key <- aggregate(
    weighted ~ type + fleet + survey + length_bin_lower + length_bin_upper + age,
    age, sum
  )
  totals <- aggregate(
    weighted ~ type + fleet + survey + length_bin_lower + length_bin_upper,
    key, sum
  )
  join_key <- function(x) do.call(paste, c(x[dims], sep = "\034"))
  key$age_given_length <- key$weighted /
    totals$weighted[match(join_key(key), join_key(totals))]
  key
}

# Combine each annual length composition with its fleet-specific pooled key.
.goa_age_shares <- function(inputs, key, type, fleet, survey, years, ages) {
  lengths <- inputs[
    inputs$type == type & inputs$measure == "proportion_at_length" &
      inputs$fleet == fleet & inputs$survey == survey & inputs$year %in% years,
    , drop = FALSE
  ]
  keys <- key[key$type == type & key$fleet == fleet & key$survey == survey, , drop = FALSE]
  if (!nrow(lengths) || !nrow(keys)) {
    stop("GOA length and age compositions do not overlap for this stream.", call. = FALSE)
  }

  bins <- aggregate(value ~ year + length_bin, lengths, sum)
  p_key <- keys[c("length_bin_lower", "age", "age_given_length")]
  names(p_key)[names(p_key) == "length_bin_lower"] <- "length_bin"

  shares_by_year <- lapply(sort(unique(bins$year)), function(year) {
    p_length <- bins[bins$year == year, c("length_bin", "value"), drop = FALSE]
    p_length$value <- p_length$value / sum(p_length$value)
    coverage <- sum(p_length$value[p_length$length_bin %in% p_key$length_bin])
    if (!is.finite(coverage) || coverage < 0.95) {
      stop("Pooled GOA age-at-length data cover less than 95% of a length composition.", call. = FALSE)
    }
    joint <- merge(p_length, p_key, by = "length_bin")
    joint$joint <- joint$value * joint$age_given_length
    shares <- aggregate(joint ~ age, joint, sum)
    shares$joint <- shares$joint / sum(shares$joint)
    if (!setequal(shares$age, ages) || any(!is.finite(shares$joint))) {
      stop("GOA age-length-key reconstruction does not cover every model age.", call. = FALSE)
    }
    data.frame(year = year, fleet = fleet, survey = survey,
               age = shares$age, proportion = shares$joint)
  })
  do.call(rbind, shares_by_year)
}

## Observations ----
assessment_id <- source$assessment$assessment_id[[1]]
years <- 2007:2025
ages <- 1:10
survey <- "NMFS bottom-trawl survey"

inputs <- source$inputs[source$inputs$assessment_id == assessment_id, , drop = FALSE]
outputs <- source$outputs[source$outputs$assessment_id == assessment_id, , drop = FALSE]
inputs$fleet[is.na(inputs$fleet)] <- ""
inputs$survey[is.na(inputs$survey)] <- ""

# Fleet-specific age shares; survey shares are derived separately.
age_key <- .goa_age_length_key(inputs, years)
fishery_fleets <- unique(inputs$fleet[
  inputs$type == "catch" & inputs$measure == "proportion_at_length"
])
catch_shares <- do.call(rbind, lapply(fishery_fleets, function(fleet) {
  .goa_age_shares(inputs, age_key, "catch", fleet, "", years, ages)
}))
index_shares <- .goa_age_shares(inputs, age_key, "index", "", survey, years, ages)

## Biological inputs ----
growth <- outputs[
  outputs$type == "biology" &
    outputs$measure %in% c("growth_length_at_age", "growth_rate"),
  , drop = FALSE
]
length_at_age_1 <- growth$value[growth$measure == "growth_length_at_age" & growth$age == 1]
length_at_age_10 <- growth$value[growth$measure == "growth_length_at_age" & growth$age == 10]
growth_rate <- growth$value[growth$measure == "growth_rate"]
if (length(length_at_age_1) != 1L || length(length_at_age_10) != 1L || length(growth_rate) != 1L) {
  stop("GOA growth estimates are incomplete in outputs.csv.", call. = FALSE)
}

growth_factor <- exp(-growth_rate * (max(ages) - min(ages)))
length_infinity <- (length_at_age_10 - length_at_age_1 * growth_factor) / (1 - growth_factor)
length <- length_infinity - (length_infinity - length_at_age_1) *
  exp(-growth_rate * (ages - min(ages)))

assumptions <- source$assumptions[
  source$assumptions$assessment_id == assessment_id, , drop = FALSE
]
parameter <- function(setting) {
  value <- as.numeric(assumptions$value[
    assumptions$component == "biology" & assumptions$setting == setting
  ])
  if (length(value) != 1L || !is.finite(value)) {
    stop("Missing GOA biological parameter: ", setting, call. = FALSE)
  }
  value
}
weight <- parameter("weight_length_coefficient") * length^parameter("weight_length_exponent")
maturity <- plogis(
  -parameter("maturity_logistic_slope") *
    (length - parameter("maturity_length_50_percent_cm"))
)

## Reconstruct catch-at-age ----
catch_totals <- inputs[
  inputs$type == "catch" & inputs$measure == "total_biomass" &
    inputs$year %in% years, , drop = FALSE
]
if (anyDuplicated(paste(catch_totals$year, catch_totals$fleet))) {
  stop("GOA catch biomass has duplicate year-fleet totals.", call. = FALSE)
}

catch_rows <- do.call(rbind, lapply(seq_len(nrow(catch_totals)), function(i) {
  total <- catch_totals[i, , drop = FALSE]
  share <- catch_shares[
    catch_shares$year == total$year & catch_shares$fleet == total$fleet,
    , drop = FALSE
  ]
  if (!setequal(share$age, ages)) stop("GOA catch age shares are incomplete.", call. = FALSE)

  age_shares <- share$proportion[match(ages, share$age)]
  mean_weight <- sum(age_shares * weight)
  tonnes_to_kg <- if (tolower(total$unit) == "t") {
    1000
  } else if (tolower(total$unit) == "kg") {
    1
  } else NA_real_
  if (!is.finite(tonnes_to_kg)) stop("Unsupported GOA catch biomass unit.", call. = FALSE)
  data.frame(
    year = total$year, age = ages,
    numbers = total$value * tonnes_to_kg / mean_weight * age_shares,
    fleet = total$fleet
  )
}))
catch_rows <- aggregate(numbers ~ year + age, catch_rows, sum)
catch_totals_by_year <- aggregate(numbers ~ year, catch_rows, sum)

## Reconstruct bottom-trawl index-at-age ----
index_totals <- inputs[
  inputs$type == "index" & inputs$survey == survey &
    inputs$measure == "total_numbers" & inputs$year %in% years,
  , drop = FALSE
]
if (anyDuplicated(index_totals$year)) {
  stop("GOA bottom-trawl index has duplicate annual totals.", call. = FALSE)
}
index_multiplier <- ifelse(
  grepl("thousand", index_totals$unit, ignore.case = TRUE), 1000,
  ifelse(grepl("fish", index_totals$unit, ignore.case = TRUE), 1, NA_real_)
)
if (any(!is.finite(index_multiplier))) {
  stop("Unsupported GOA survey-index unit.", call. = FALSE)
}
index_rows <- do.call(rbind, lapply(seq_len(nrow(index_totals)), function(i) {
  total <- index_totals[i, , drop = FALSE]
  share <- index_shares[index_shares$year == total$year, , drop = FALSE]
  if (nrow(share) != length(ages) || !setequal(share$age, ages)) {
    stop("GOA bottom-trawl age shares do not match an index year.", call. = FALSE)
  }
  data.frame(
    year = total$year, age = ages,
    numbers = total$value * index_multiplier[[i]] *
      share$proportion[match(ages, share$age)]
  )
}))

## Source M on the tinyAM year–age grid ----
m_outputs <- outputs[
  outputs$type == "mortality" & outputs$measure == "natural_mortality_at_age",
  , drop = FALSE
]
m_values <- m_outputs$value[match(
  paste(rep(years, each = length(ages)), rep(ages, times = length(years))),
  paste(m_outputs$year, m_outputs$age)
)]

## Provenance-preserving reconstructed inputs ----
row <- function(type, measure, basis = "", fleet = "", survey = "",
                year = NA_real_, age = NA_real_, value, unit,
                sampling_time = NA_real_, source_reference, transformation, notes) {
  data.frame(
    assessment_id = assessment_id, type = type, measure = measure,
    basis = basis, fleet = fleet, survey = survey, sex = "", region = "",
    season = "", year = year,
    year_basis = ifelse(is.na(year), "", "calendar_year"),
    age = age, value = value, unit = unit, sampling_time = sampling_time,
    source_type = "reconstructed_source_input", source_reference = source_reference,
    transformation = transformation, notes = notes, stringsAsFactors = FALSE
  )
}
pad <- function(rows) {
  out <- inputs[rep(NA_integer_, nrow(rows)), , drop = FALSE]
  rownames(out) <- NULL
  for (name in names(rows)) out[[name]] <- rows[[name]]
  out
}
references <- function(rows) paste(unique(rows$source_reference), collapse = "; ")

catch_source <- inputs[
  inputs$type == "catch" & inputs$measure %in% c(
    "total_biomass", "proportion_at_length", "conditional_proportion_at_age"
  ), , drop = FALSE
]
index_source <- inputs[
  inputs$type == "index" & inputs$survey == survey &
    inputs$measure %in% c(
      "total_numbers", "log_index_sd", "proportion_at_length",
      "conditional_proportion_at_age"
    ), , drop = FALSE
]
growth_source <- outputs[outputs$measure %in% c("growth_length_at_age", "growth_rate"), ]
biology_reference <- references(assumptions[assumptions$component == "biology", ])

# These are translation-time reconstructions; no tinyAM-fitted quantities
# are used in constructing the observations.
reconstructed <- list(
  row(
    "catch", "total_numbers", year = catch_totals_by_year$year,
    age = NA_real_, value = catch_totals_by_year$numbers, unit = "fish",
    source_reference = references(catch_source),
    transformation = "Biomass converted to fish using reconstructed age proportions and weight-at-age.",
    notes = paste(
      "Fleet catches are aggregated after applying fleet-specific pooled age-length keys;",
      "2025 uses the key pooled from observed ages 2007-2024."
    )
  ),
  row(
    "catch", "proportion_at_age", basis = "proportion_numbers",
    year = catch_rows$year, age = catch_rows$age,
    value = catch_rows$numbers / catch_totals_by_year$numbers[
      match(catch_rows$year, catch_totals_by_year$year)
    ],
    unit = "proportion", source_reference = references(catch_source),
    transformation = paste(
      "P(length) multiplied by sample-size-weighted P(age|length),",
      "then normalized by year and fleet before fleet aggregation."
    ),
    notes = "Native bin labels are matched as categories. Original catch and composition rows remain unchanged."
  ),
  row(
    "index", "total_numbers", survey = survey, year = index_totals$year,
    age = NA_real_, value = index_totals$value, unit = index_totals$unit,
    sampling_time = index_totals$sampling_time,
    source_reference = references(index_source),
    transformation = "Native aggregate index retained before age allocation.",
    notes = "Only the age-sampled NMFS bottom-trawl survey is translated."
  ),
  row(
    "index", "proportion_at_age", basis = "proportion_numbers", survey = survey,
    year = index_rows$year, age = index_rows$age,
    value = index_rows$numbers / (
      index_totals$value[match(index_rows$year, index_totals$year)] *
        index_multiplier[match(index_rows$year, index_totals$year)]
    ),
    unit = "proportion", sampling_time = 0.5,
    source_reference = references(index_source),
    transformation = paste(
      "P(length) multiplied by the sample-size-weighted pooled",
      "P(age|length) and normalized."
    ),
    notes = "Original SS3 composition likelihoods, ageing-error likelihood and covariance are not reproduced."
  ),
  row(
    "index", "log_index_sd", survey = survey,
    year = index_source$year[index_source$measure == "log_index_sd"],
    value = index_source$value[index_source$measure == "log_index_sd"],
    unit = "log scale", source_reference = references(index_source),
    transformation = "Native log-scale index SD retained for every age allocated from the aggregate survey index.",
    notes = "The same aggregate-index SD is applied to each reconstructed age row."
  ),
  row(
    "weight", "weight_at_age", basis = "kg_per_fish", age = ages,
    value = weight, unit = "kg per fish", source_reference = references(growth_source),
    transformation = "Von Bertalanffy curve anchored to reported mean lengths at ages 1 and 10; native weight-length curve applied to mean length.",
    notes = "The reported length-at-age SD is not integrated; age 10 is the plus group."
  ),
  row(
    "maturity", "maturity_at_age", basis = "proportion", age = ages,
    value = maturity, unit = "proportion mature", source_reference = biology_reference,
    transformation = "Native Stock Synthesis length-logistic maturity evaluated at reconstructed mean length-at-age.",
    notes = "The female SSB convention applies a 0.5 multiplier in the tinyAM translation."
  ),
  row(
    "M", "natural_mortality_at_age", basis = "per_year",
    year = rep(years, each = length(ages)),
    age = rep(ages, times = length(years)),
    value = as.numeric(m_values), unit = "per year",
    source_reference = references(outputs[
      outputs$type == "mortality" &
        outputs$measure == "natural_mortality_at_age", , drop = FALSE
    ]),
    transformation = "Official age-invariant M estimates expanded over age and year.",
    notes = "The source estimates M; this approximation fixes its published 2014-2016 block estimates."
  )
)
translation_inputs <- do.call(rbind, lapply(reconstructed, pad))
if (any(!is.finite(translation_inputs$value[translation_inputs$type == "M"]))) {
  stop("GOA official M does not cover the selected fit years and ages.", call. = FALSE)
}

## tinyAM observations ----
obs <- database_to_tam_obs(
  assessment_id, translation_inputs, years = years, ages = ages,
  weight_survey = "", index_weight_source = "stock",
  sampling_times = c("NMFS bottom-trawl survey" = 0.5),
  surveys = survey, maturity_multiplier = 0.5,
  assumptions = source$assumptions
)

## Official comparisons ----
# Keep native source SSB in tonnes. The comparison code handles unit conversion.
comparison_outputs <- outputs[
  outputs$measure %in% c("SSB", "total_biomass", "recruitment") &
    outputs$year %in% years, , drop = FALSE
]
if (!any(comparison_outputs$measure == "SSB")) {
  stop("No GOA female SSB outputs are available for the selected comparison years.", call. = FALSE)
}
biomass_rows <- comparison_outputs$measure %in% c("SSB", "total_biomass")
if (any(tolower(trimws(comparison_outputs$unit[biomass_rows])) != "t") ||
    any(tolower(trimws(comparison_outputs$unit[
      comparison_outputs$measure == "recruitment"
    ])) != "billion fish")) {
  stop("GOA aggregate biomass and recruitment outputs have unexpected units.", call. = FALSE)
}

ssb_note <- paste(
  "Source female SSB is retained in native tonnes;",
  "tinyAM SSB is converted from kg to tonnes for comparison.",
  "The official SSB is age 0+ whereas tinyAM represents ages 1-10+,",
  "so this is an approximate aggregate comparison."
)
ssb_rows <- comparison_outputs$measure == "SSB"
comparison_outputs$notes[ssb_rows] <- ifelse(
  is.na(comparison_outputs$notes[ssb_rows]) |
    !nzchar(trimws(comparison_outputs$notes[ssb_rows])),
  ssb_note,
  paste(comparison_outputs$notes[ssb_rows], ssb_note)
)

## Background and comparisons ----
age_plus_group <- 10

comparison_aggregates <- "ssb"

comparison_definitions <- list(
  ssb = list(
    status = "approximate", definition = "Aggregate female SSB: official ages 0+ versus tinyAM ages 1-10+",
    reason = "The official age-0 contribution cannot be removed from the available aggregate output."
  ),
  biomass = list(
    status = "non_equivalent", definition = "Official total biomass ages 0+ versus tinyAM ages 1-10+",
    reason = "The official age-0 biomass cannot be removed from the reported aggregate total."
  )
)

background <- c(
  "### Gulf of Alaska Pacific cod: official 2026 update, Model 24.0", "", print_sources(source$assessment),
  "", "| Component | Official assessment | tinyAM representation | Reason for difference |", "|---|---|---|---|",
  "| Years | Model 24.0 covers 1977–2025. | Fit 2007–2025. | Fishery conditional age-at-length samples begin in 2007; the shorter window avoids extrapolating the reconstructed age–length relationship across earlier decades. |",
  "| Ages | Population ages 0–10+, with recruitment at age 0; fishery and survey age compositions cover ages 1–10. | Ages 1–10+, with recruitment entering at age 1. | Age 0 is not represented, so official age-0 recruitment and age-0+ population totals are not exactly comparable. |",
  "| Recruitment | Beverton–Holt relationship with steepness fixed at 1 and annual recruitment deviations (source recruitment SD fixed at 0.44). | Default tinyAM log-recruitment random walk; first-year recruitment is a fixed anchor. | Recruitment enters at age 1 instead of age 0, and an RW is not the source's independent annual recruitment-deviation structure. |",
  "| N | Deterministic cohort survival given recruitment, F and M, with population initialization before the translated period. | Deterministic older-age survival (`process = 'off'`) and exponential initial abundance (`init = 'exp'`). | The shorter fit period and initialization cannot replicate SS3's history. Free initial states were less numerically reliable in earlier tests. |",
  "| F — principal approximation | Three fisheries have separate size-based selectivity, temporal selectivity blocks and annual fishing mortality. | The candidate fits correlated age-specific log-F random walks (`process = 'cor_rw'`), allowing persistent annual changes with correlated innovations across ages. | The correlated RW allows F-at-age to evolve without fixing a common logistic age pattern; it differs from the fleet-specific length-selectivity structure. This candidate requires refitting. It remains a single aggregated fishery, without native fleet/length selectivity or its temporal constraints. |",
  "| M | Age-constant M estimated in the source, with a higher block over 2014–2016. | Supply published M estimates as fixed values for each modeled year and age (`process = 'off'`). | tinyAM conditions on source estimates instead of estimating M or propagating its uncertainty. |",
  "| Catch | Three fishery biomass series, with length compositions and conditional age-at-length likelihoods. | Fleet-specific pooled age–length keys reconstruct annual numbers-at-age; aggregate the three fleets and fit catch-at-age with a common lognormal observation SD (`sd_form = ~1`). | Does not reproduce fleet-specific likelihoods, ageing error, length compositions or correlation; fitted catch-at-age is reconstructed rather than directly observed. |",
  "| Index | Bottom-trawl and longline survey indices with distinct catchability/selectivity and composition data. | Reconstruct bottom-trawl index-at-age using pooled age–length data; fit a rising logistic age-selectivity curve plus an annual random walk on log catchability (`~ logistic(age) + rw(year)`). | The longline survey is omitted because age composition is unavailable, and the reconstructed age-index likelihood does not reproduce native aggregate/conditional composition likelihoods. |",
  "| Index error | Aggregate survey observations and length/conditional-age compositions have distinct observation likelihoods and uncertainty. | Repeat the supplied aggregate log-scale survey SD across reconstructed age rows and estimate an additional common multiplicative log-SD (`sd_form = ~1`, `sd_supplied = ~relative_sd`). | Allocated age cells are not statistically independent observations in the source; their within-year dependence and source effective sample sizes are not represented. |",
  "| Weights and maturity | Growth-based weight-at-length and length-logistic maturity; reported female SSB uses one-half of the native spawning output. | Reconstruct mean length-at-age, evaluate the source weight–length and maturity functions at mean length, and multiply maturity by 0.5. | Does not integrate length variability within age; biological surfaces remain approximations. |",
  "| SSB | Reported female spawning biomass covers ages 0+ and is in tonnes. | Calculate female SSB across ages 1–10+ using translated weights and maturity. | The source age-0 component is unavailable separately; the aggregate comparison is approximate even with common biomass units. |",
  "| Total biomass | Reported aggregate biomass includes age 0+. | tinyAM population biomass includes ages 1–10+. | No definition-matched scale difference is reported because the official age-0 contribution cannot be separated. |",
  "", "The pooled age–length keys use native length-bin labels and samples from 2007 onward; they cover at least 95% of each selected annual length composition. Uncovered length mass is omitted and remaining shares normalized. The 2025 catch uses the key based on age-at-length samples through 2024.",
  "", "This is a new candidate combining correlated log-F random walks (`cor_rw`), deterministic cohort survival, and a logistic survey q curve with a log-q year RW. The earlier age-specific-mean/IID F formulation was associated with improved agreement but its fit statistics are not results for this candidate. Correlated F increments permit smooth coherent movements across ages without modelling the original fleet-specific length selectivity.",
  "", "The selected candidate retains deterministic cohort survival and fixed source M; an additional N-deviation process was not needed. Survey q is a logistic function of age with a shared year-RW effect (`~ logistic(age) + rw(year)`). The official age-0 recruitment, fleet-specific catch dynamics, longline survey, conditional-composition likelihoods and native N/F-at-age comparisons are not replicated. The reported SSB agreement should be interpreted as approximate because the modeled age ranges differ."
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
    N_settings = list(process = "off", init = "exp"),
    F_settings = list(process = "cor_rw", mu_form = NULL),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption),
    catch_settings = list(sd_form = ~ 1),
    index_settings = list(
      q_form = ~ logistic(age) + rw(year),
      sd_form = ~ 1,
      sd_supplied = ~ relative_sd
    ),
    silent = silent
  )
}
