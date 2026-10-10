.goa_age_length_key <- function(inputs, years) {
  age <- inputs[
    inputs$type %in% c("catch", "index") &
      inputs$measure == "conditional_proportion_at_age" &
      inputs$year %in% years,
    ,
    drop = FALSE
  ]

  if (!nrow(age) || anyNA(age$sample_size) || any(age$sample_size <= 0)) {
    stop(
      "GOA age-at-length records need positive supplied sample sizes.",
      call. = FALSE
    )
  }

  sums <- aggregate(value ~ observation_id, age, sum)
  if (any(abs(sums$value - 1) > 0.001)) {
    stop(
      "A GOA conditional age composition does not sum to one.",
      call. = FALSE
    )
  }

  age$weighted <- age$value * age$sample_size

  key <- aggregate(
    weighted ~ type + fleet + survey + length_bin_lower +
      length_bin_upper + age,
    age,
    sum
  )

  totals <- aggregate(
    weighted ~ type + fleet + survey + length_bin_lower +
      length_bin_upper,
    key,
    sum
  )

  join_key <- function(data) {
    do.call(
      paste,
      c(
        data[
          c(
            "type", "fleet", "survey",
            "length_bin_lower", "length_bin_upper"
          )
        ],
        sep = "\034"
      )
    )
  }

  key$age_given_length <-
    key$weighted / totals$weighted[match(join_key(key), join_key(totals))]

  key
}


.goa_age_shares <- function(inputs, key, type, fleet, survey, years, ages) {
  lengths <- inputs[
    inputs$type == type &
      inputs$measure == "proportion_at_length" &
      inputs$fleet == fleet &
      inputs$survey == survey &
      inputs$year %in% years,
    ,
    drop = FALSE
  ]

  keys <- key[
    key$type == type &
      key$fleet == fleet &
      key$survey == survey,
    ,
    drop = FALSE
  ]

  if (!nrow(lengths) || !nrow(keys)) {
    stop(
      "GOA length and age compositions do not overlap for this stream.",
      call. = FALSE
    )
  }

  bins <- aggregate(value ~ year + length_bin, lengths, sum)

  output <- lapply(sort(unique(bins$year)), function(year) {
    p_length <- bins[
      bins$year == year,
      c("length_bin", "value"),
      drop = FALSE
    ]
    p_length$value <- p_length$value / sum(p_length$value)

    p_key <- keys[
      ,
      c("length_bin_lower", "age", "age_given_length"),
      drop = FALSE
    ]
    names(p_key)[names(p_key) == "length_bin_lower"] <- "length_bin"

    coverage <- sum(
      p_length$value[p_length$length_bin %in% p_key$length_bin]
    )

    if (!is.finite(coverage) || coverage < 0.95) {
      stop(
        "Pooled GOA age-at-length data cover less than 95% of a length composition.",
        call. = FALSE
      )
    }

    joint <- merge(p_length, p_key, by = "length_bin")
    joint$joint <- joint$value * joint$age_given_length

    shares <- aggregate(joint ~ age, joint, sum)
    shares$joint <- shares$joint / sum(shares$joint)

    if (!setequal(shares$age, ages) || any(!is.finite(shares$joint))) {
      stop(
        "GOA age-length-key reconstruction does not cover every model age.",
        call. = FALSE
      )
    }

    data.frame(
      year = year,
      fleet = fleet,
      survey = survey,
      age = shares$age,
      proportion = shares$joint
    )
  })

  do.call(rbind, output)
}


## Observations ----
assessment_id <- source$assessment$assessment_id[[1]]

years <- 2007:2025
ages <- 1:10
survey <- "NMFS bottom-trawl survey"

inputs <- source$inputs[
  source$inputs$assessment_id == assessment_id,
  ,
  drop = FALSE
]
outputs <- source$outputs[
  source$outputs$assessment_id == assessment_id,
  ,
  drop = FALSE
]

inputs$fleet[is.na(inputs$fleet)] <- ""
inputs$survey[is.na(inputs$survey)] <- ""

## Reconstruct age proportions ----

age_key <- .goa_age_length_key(inputs, years)

catch_shares <- do.call(
  rbind,
  lapply(
    unique(
      inputs$fleet[
        inputs$type == "catch" &
          inputs$measure == "proportion_at_length"
      ]
    ),
    function(fleet) {
      .goa_age_shares(
        inputs, age_key,
        "catch", fleet, "",
        years, ages
      )
    }
  )
)

index_shares <- .goa_age_shares(
  inputs,
  age_key,
  "index",
  "",
  survey,
  years,
  ages
)

## Reconstruct biological inputs ----

growth <- outputs[
  outputs$type == "biology" &
    outputs$measure %in%
    c("growth_length_at_age", "growth_rate"),
  ,
  drop = FALSE
]

length_at_age_1 <- growth$value[
  growth$measure == "growth_length_at_age" &
    growth$age == 1
]

length_at_age_10 <- growth$value[
  growth$measure == "growth_length_at_age" &
    growth$age == 10
]

growth_rate <- growth$value[
  growth$measure == "growth_rate"
]

if (
  length(length_at_age_1) != 1L ||
  length(length_at_age_10) != 1L ||
  length(growth_rate) != 1L
) {
  stop(
    "GOA growth estimates are incomplete in outputs.csv.",
    call. = FALSE
  )
}

growth_factor <-
  exp(-growth_rate * (max(ages) - min(ages)))

length_infinity <-
  (length_at_age_10 - length_at_age_1 * growth_factor) /
  (1 - growth_factor)

length <-
  length_infinity -
  (length_infinity - length_at_age_1) *
  exp(-growth_rate * (ages - min(ages)))

assumptions <- source$assumptions[
  source$assumptions$assessment_id == assessment_id,
  ,
  drop = FALSE
]

parameter <- function(setting) {
  value <- as.numeric(
    assumptions$value[
      assumptions$component == "biology" &
        assumptions$setting == setting
    ]
  )

  if (length(value) != 1L || !is.finite(value)) {
    stop(
      "Missing GOA biological parameter: ",
      setting,
      call. = FALSE
    )
  }

  value
}

weight <-
  parameter("weight_length_coefficient") *
  length^parameter("weight_length_exponent")

maturity <-
  plogis(
    -parameter("maturity_logistic_slope") *
      (
        length -
          parameter("maturity_length_50_percent_cm")
      )
  )

## Reconstruct catch-at-age ----

catch_totals <- inputs[
  inputs$type == "catch" &
    inputs$measure == "total_biomass" &
    inputs$year %in% years,
  ,
  drop = FALSE
]

if (
  anyDuplicated(
    paste(catch_totals$year, catch_totals$fleet)
  )
) {
  stop(
    "GOA catch biomass has duplicate year-fleet totals.",
    call. = FALSE
  )
}

catch_rows <- do.call(
  rbind,
  lapply(seq_len(nrow(catch_totals)), function(i) {
    total <- catch_totals[i, , drop = FALSE]

    share <- catch_shares[
      catch_shares$year == total$year &
        catch_shares$fleet == total$fleet,
      ,
      drop = FALSE
    ]

    if (!setequal(share$age, ages)) {
      stop(
        "GOA catch age shares are incomplete.",
        call. = FALSE
      )
    }

    mean_weight <-
      sum(
        share$proportion[match(ages, share$age)] *
          weight
      )

    tonnes_to_kg <-
      if (tolower(total$unit) == "t") {
        1000
      } else if (tolower(total$unit) == "kg") {
        1
      } else {
        NA_real_
      }

    if (!is.finite(tonnes_to_kg)) {
      stop(
        "Unsupported GOA catch biomass unit.",
        call. = FALSE
      )
    }

    data.frame(
      year = total$year,
      age = ages,
      numbers =
        total$value *
        tonnes_to_kg /
        mean_weight *
        share$proportion[match(ages, share$age)],
      fleet = total$fleet
    )
  })
)

catch_rows <-
  aggregate(numbers ~ year + age, catch_rows, sum)

catch_totals_by_year <-
  aggregate(numbers ~ year, catch_rows, sum)

## Reconstruct bottom-trawl index-at-age ----

index_totals <- inputs[
  inputs$type == "index" &
    inputs$survey == survey &
    inputs$measure == "total_numbers" &
    inputs$year %in% years,
  ,
  drop = FALSE
]

if (anyDuplicated(index_totals$year)) {
  stop(
    "GOA bottom-trawl index has duplicate annual totals.",
    call. = FALSE
  )
}

index_multiplier <- ifelse(
  grepl(
    "thousand",
    index_totals$unit,
    ignore.case = TRUE
  ),
  1000,
  ifelse(
    grepl(
      "fish",
      index_totals$unit,
      ignore.case = TRUE
    ),
    1,
    NA_real_
  )
)

if (any(!is.finite(index_multiplier))) {
  stop(
    "Unsupported GOA survey-index unit.",
    call. = FALSE
  )
}

index_rows <- do.call(
  rbind,
  lapply(seq_len(nrow(index_totals)), function(i) {
    total <- index_totals[i, , drop = FALSE]

    share <- index_shares[
      index_shares$year == total$year,
      ,
      drop = FALSE
    ]

    if (
      nrow(share) != length(ages) ||
      !setequal(share$age, ages)
    ) {
      stop(
        "GOA bottom-trawl age shares do not match an index year.",
        call. = FALSE
      )
    }

    data.frame(
      year = total$year,
      age = ages,
      numbers =
        total$value *
        index_multiplier[[i]] *
        share$proportion[match(ages, share$age)]
    )
  })
)

## Accepted M translated onto tinyAM age-year grid ----

m_outputs <- outputs[
  outputs$type == "mortality" &
    outputs$measure == "natural_mortality_at_age",
  ,
  drop = FALSE
]

m_values <- m_outputs$value[
  match(
    paste(
      rep(years, each = length(ages)),
      rep(ages, times = length(years))
    ),
    paste(
      m_outputs$year,
      m_outputs$age
    )
  )
]

## Build reconstructed translation inputs ----

row <- function(
  type,
  measure,
  basis = "",
  fleet = "",
  survey = "",
  year = NA_real_,
  age = NA_real_,
  value,
  unit,
  sampling_time = NA_real_,
  source_reference,
  transformation,
  notes
) {
  data.frame(
    assessment_id = assessment_id,
    type = type,
    measure = measure,
    basis = basis,
    fleet = fleet,
    survey = survey,
    sex = "",
    region = "",
    season = "",
    year = year,
    year_basis = ifelse(
      is.na(year),
      "",
      "calendar_year"
    ),
    age = age,
    value = value,
    unit = unit,
    sampling_time = sampling_time,
    source_type = "reconstructed_source_input",
    source_reference = source_reference,
    transformation = transformation,
    notes = notes,
    stringsAsFactors = FALSE
  )
}

pad <- function(rows) {
  out <- inputs[
    rep(NA_integer_, nrow(rows)),
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL

  for (name in names(rows)) {
    out[[name]] <- rows[[name]]
  }

  out
}

references <- function(rows) {
  paste(
    unique(rows$source_reference),
    collapse = "; "
  )
}

catch_source <- inputs[
  inputs$type == "catch" &
    inputs$measure %in%
    c(
      "total_biomass",
      "proportion_at_length",
      "conditional_proportion_at_age"
    ),
  ,
  drop = FALSE
]

index_source <- inputs[
  inputs$type == "index" &
    inputs$survey == survey &
    inputs$measure %in%
    c(
      "total_numbers",
      "log_index_sd",
      "proportion_at_length",
      "conditional_proportion_at_age"
    ),
  ,
  drop = FALSE
]

growth_source <- outputs[
  outputs$measure %in%
    c("growth_length_at_age", "growth_rate"),
  ,
  drop = FALSE
]

biology_reference <-
  references(
    assumptions[
      assumptions$component == "biology",
      ,
      drop = FALSE
    ]
  )

reconstructed <- list(
  row(
    "catch",
    "total_numbers",
    year = catch_totals_by_year$year,
    age = NA_real_,
    value = catch_totals_by_year$numbers,
    unit = "fish",
    source_reference = references(catch_source),
    transformation =
      "Biomass converted to fish using reconstructed age proportions and weight-at-age.",
    notes =
      "Fleet catches are aggregated after applying fleet-specific pooled age-length keys; 2025 uses the key pooled from observed ages 2007-2024."
  ),

  row(
    "catch",
    "proportion_at_age",
    basis = "proportion_numbers",
    year = catch_rows$year,
    age = catch_rows$age,
    value =
      catch_rows$numbers /
      catch_totals_by_year$numbers[
        match(
          catch_rows$year,
          catch_totals_by_year$year
        )
      ],
    unit = "proportion",
    source_reference = references(catch_source),
    transformation =
      "P(length) multiplied by sample-size-weighted P(age|length), then normalized by year and fleet before fleet aggregation.",
    notes =
      "Native bin labels are matched as categories. Original catch and composition rows remain unchanged."
  ),

  row(
    "index",
    "total_numbers",
    survey = survey,
    year = index_totals$year,
    age = NA_real_,
    value = index_totals$value,
    unit = index_totals$unit,
    sampling_time = index_totals$sampling_time,
    source_reference = references(index_source),
    transformation =
      "Native aggregate index retained before age allocation.",
    notes =
      "Only the age-sampled NMFS bottom-trawl survey is translated."
  ),

  row(
    "index",
    "proportion_at_age",
    basis = "proportion_numbers",
    survey = survey,
    year = index_rows$year,
    age = index_rows$age,
    value =
      index_rows$numbers /
      (
        index_totals$value[
          match(
            index_rows$year,
            index_totals$year
          )
        ] *
          index_multiplier[
            match(
              index_rows$year,
              index_totals$year
            )
          ]
      ),
    unit = "proportion",
    sampling_time = 0.5,
    source_reference = references(index_source),
    transformation =
      "P(length) multiplied by the sample-size-weighted pooled P(age|length) and normalized.",
    notes =
      "Original SS3 composition likelihoods, ageing-error likelihood and covariance are not reproduced."
  ),

  row(
    "index",
    "log_index_sd",
    survey = survey,
    year =
      index_source$year[
        index_source$measure == "log_index_sd"
      ],
    value =
      index_source$value[
        index_source$measure == "log_index_sd"
      ],
    unit = "log scale",
    source_reference = references(index_source),
    transformation =
      "Native log-scale index SD retained for every age allocated from the aggregate survey index.",
    notes =
      "The same aggregate-index SD is applied to each reconstructed age row."
  ),

  row(
    "weight",
    "weight_at_age",
    basis = "kg_per_fish",
    age = ages,
    value = weight,
    unit = "kg per fish",
    source_reference = references(growth_source),
    transformation =
      "Von Bertalanffy curve anchored to reported mean lengths at ages 1 and 10; native weight-length curve applied to mean length.",
    notes =
      "The reported length-at-age SD is not integrated; age 10 is the plus group."
  ),

  row(
    "maturity",
    "maturity_at_age",
    basis = "proportion",
    age = ages,
    value = maturity,
    unit = "proportion mature",
    source_reference = biology_reference,
    transformation =
      "Native Stock Synthesis length-logistic maturity evaluated at reconstructed mean length-at-age.",
    notes =
      "The female SSB convention applies a 0.5 multiplier in the tinyAM translation."
  ),

  row(
    "M",
    "natural_mortality_at_age",
    basis = "per_year",
    year = rep(
      years,
      each = length(ages)
    ),
    age = rep(
      ages,
      times = length(years)
    ),
    value = as.numeric(m_values),
    unit = "per year",
    source_reference =
      references(
        outputs[
          outputs$type == "mortality" &
            outputs$measure ==
            "natural_mortality_at_age",
          ,
          drop = FALSE
        ]
      ),
    transformation =
      "Accepted age-invariant M estimates expanded over age and year.",
    notes =
      "The source estimates M; this approximation fixes its published 2014-2016 block estimates."
  )
)

translation_inputs <-
  do.call(
    rbind,
    lapply(reconstructed, pad)
  )

if (
  any(
    !is.finite(
      translation_inputs$value[
        translation_inputs$type == "M"
      ]
    )
  )
) {
  stop(
    "GOA accepted M does not cover the selected fit years and ages.",
    call. = FALSE
  )
}

## tinyAM observations ----

obs <- database_to_tam_obs(
  assessment_id,
  translation_inputs,
  years = years,
  ages = ages,
  weight_survey = "",
  index_weight_source = "stock",
  sampling_times =
    c("NMFS bottom-trawl survey" = 0.5),
  surveys = survey,
  maturity_multiplier = 0.5,
  assumptions = source$assumptions
)

## Accepted outputs used for comparison ----
##
## Retain source SSB in its native unit (tonnes). tinyAM produces
## biomass in kg, so comparison_scales converts tinyAM kg -> tonnes.
## Do not convert the source values here.

comparison_outputs <- outputs[
  outputs$measure %in% c("SSB", "total_biomass", "recruitment") &
    outputs$year %in% years,
  ,
  drop = FALSE
]

if (!any(comparison_outputs$measure == "SSB")) {
  stop(
    "No GOA female SSB outputs are available for the selected comparison years.",
    call. = FALSE
  )
}

biomass_rows <- comparison_outputs$measure %in% c("SSB", "total_biomass")
if (any(tolower(trimws(comparison_outputs$unit[biomass_rows])) != "t") ||
    any(tolower(trimws(comparison_outputs$unit[
      comparison_outputs$measure == "recruitment"
    ])) != "billion fish")) {
  stop(
    "GOA aggregate biomass and recruitment outputs have unexpected units.",
    call. = FALSE
  )
}

ssb_note <- paste(
  "Source female SSB is retained in native tonnes;",
  "tinyAM SSB is converted from kg to tonnes for comparison.",
  "The accepted SSB is age 0+ whereas tinyAM represents ages 1-10+,",
  "so this is an approximate aggregate comparison."
)

ssb_rows <- comparison_outputs$measure == "SSB"
comparison_outputs$notes[ssb_rows] <- ifelse(
  is.na(comparison_outputs$notes[ssb_rows]) |
    !nzchar(trimws(comparison_outputs$notes[ssb_rows])),
  ssb_note,
  paste(comparison_outputs$notes[ssb_rows], ssb_note)
)

## Model specification ----


## Background and comparisons ----

age_plus_group <- 10


comparison_scales <- c(ssb = 0.001, biomass = 0.001, recruitment = 1e-09)

comparison_aggregates <- "ssb"

comparison_definitions <- list(ssb = list(status = "approximate", definition = "Aggregate female SSB: accepted ages 0+ versus tinyAM ages 1-10+",
    reason = "The accepted age-0 contribution cannot be removed from the available aggregate output."),
    biomass = list(status = "non_equivalent", definition = "Accepted total biomass ages 0+ versus tinyAM ages 1-10+",
        reason = "The accepted age-0 biomass cannot be removed from the reported aggregate total."))

background <- c("### Gulf of Alaska Pacific cod: accepted 2026 update, Model 24.0", "", print_sources(source$assessment),
    "", "| Component | Accepted assessment | tinyAM representation | Reason for difference |", "|---|---|---|---|",
    "| Years | Model 24.0 covers 1977-2025. | Fit 2007-2025. | Fishery age-at-length samples begin in 2007; the shorter window avoids extrapolating the reconstructed age-length relationship across earlier decades. |",
    "| Ages | The population includes ages 0-10+, with recruitment at age 0; fishery and survey age compositions use ages 1-10. | Model ages 1-10, with age 10 as the plus group and recruitment entering at age 1. | Age 0 is not represented in the translated observations, so recruitment and age-0+ population quantities are not directly comparable. |",
    "| Recruitment | Recruitment enters at age 0 under a Beverton-Holt stock-recruitment model with annual recruitment deviations. | Recruitment enters at age 1 and follows tinyAM's temporal recruitment process. | Recruitment age and process differ, so recruitment is not directly compared. |",
    "| N | Stock Synthesis propagates cohorts through fishing and natural mortality from a population initialized before the assessment period. | Cohort abundance after recruitment follows deterministic survival. Initial abundance-at-age is generated by exponential survivorship from the first recruitment state. | The tinyAM fit begins in 2007 rather than 1977 and cannot reproduce the full SS3 initialization. Freely estimating initial abundance was tested but did not converge reliably, so the exponential initialization is retained. |",
    "| F | Three fisheries have distinct size-based selectivity, temporal selectivity blocks and annual fishing mortality. | Fit an age-specific mean F with independent age-year deviations around those means. | The age-specific mean provides a simple analogue of persistent fishery selectivity and substantially improves the translated fit. Fleet-specific and time-varying selectivity are not reproduced, and the three fisheries are combined. |",
    "| M | Age-constant M is estimated, with a higher 2014-2016 block. | Supply the accepted M estimates as fixed values across modeled ages, including the 2014-2016 block. | This conditions on accepted-assessment M rather than estimating it independently and does not propagate its uncertainty. |",
    "| Catch | Three fishery biomass series are fitted with length compositions and age proportions conditional on length. | Fleet-specific pooled age-length keys reconstruct annual age proportions; reconstructed weight-at-age converts biomass to numbers before fisheries are combined and fitted as aggregate catch-at-age. | tinyAM does not reproduce the native conditional-composition likelihood, ageing error, fleet-specific selectivity or covariance structure. |",
    "| Index | Bottom-trawl and longline surveys have their own catchability and size-dependent selectivity; the source model fits aggregate indices together with composition information. | Use the bottom-trawl survey only. Allocate its aggregate index across age using the pooled age-length key and estimate separate catchability by age. | Age-specific q provides a simple representation of survey selectivity and markedly improves the fit. A monotonic age-q constraint produced a very similar assessment, so the less restrictive age-specific formulation is retained. The longline survey and native conditional-composition likelihood are omitted. |",
    "| Index error | The accepted model fits an aggregate survey index and separate length and conditional age-composition likelihoods. | Apply the supplied aggregate log-scale survey SD to each reconstructed age-specific index observation, with one additional observation-error structure. | Splitting an aggregate index into age-specific observations does not reproduce the source covariance or composition uncertainty, so residual diagnostics should be interpreted accordingly. |",
    "| Weights and maturity | The model uses a combined-sex growth pattern, weight-length relationship and length-logistic maturity; published female SSB is one-half of the native spawning output. | Reconstruct mean length-at-age from reported growth parameters, calculate weight-at-age from the accepted weight-length relationship, and evaluate maturity at mean length. Maturity is multiplied by 0.5 for female SSB. | Length-at-age variation is not integrated, but the accepted biological functions and female-SSB convention are retained. |",
    "| SSB | Published female spawning biomass is age 0+ and reported in tonnes. | Calculate female SSB over ages 1-10+ and convert tinyAM biomass from kg to tonnes for comparison. | Source age-specific SSB contributions are unavailable, so the age-0 component cannot currently be removed. The aggregate comparison is therefore approximate rather than exactly age-matched. |",
    "| Total biomass | The published aggregate covers ages 0+ and is reported in tonnes. | tinyAM reports population biomass over ages 1-10+. | Both totals are shown in tonnes for context, but the age-0 contribution cannot be removed, so no like-for-like percent difference is calculated. |",
    "| Recruitment | Recruitment enters at age 0 and is reported in billions of fish. | tinyAM recruitment enters at age 1 and is shown in billions of fish. | The age definitions differ, so the dashboard can show both trajectories but the comparison is marked non-equivalent. |",
    "", "The pooled age-length keys use native bin labels and source observations from 2007 onward. They cover at least 95% of each selected annual length composition; uncovered length-bin mass is renormalized. The 2025 catch uses the key pooled from age samples through 2024.",
    "", "The translation is intentionally simpler than Model 24.0. In particular, adding persistent age structure to both fishing mortality and survey catchability substantially improves the fit while preserving deterministic cohort survival. An IID abundance process is therefore not required for the current representation. Accepted recruitment is age 0 whereas tinyAM recruitment enters at age 1, and accepted N-at-age and F-at-age surfaces are not yet included in the comparison database.")

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
    F_settings = list(process = "iid", mu_form = ~factor(age)),
    M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption),
    catch_settings = list(sd_form = ~1, fill_missing = FALSE),
    index_settings = list(q_form = ~q_block, sd_form = ~1, sd_supplied = ~relative_sd, fill_missing = FALSE),
    silent = silent
  )
}
