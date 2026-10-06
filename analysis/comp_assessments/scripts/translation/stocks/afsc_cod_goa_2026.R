.goa_age_length_key <- function(inputs, years) {
  age <- inputs[inputs$type %in% c("catch", "index") &
                  inputs$measure == "conditional_proportion_at_age" &
                  inputs$year %in% years, , drop = FALSE]
  if (!nrow(age) || anyNA(age$sample_size) || any(age$sample_size <= 0)) {
    stop("GOA age-at-length records need positive supplied sample sizes.", call. = FALSE)
  }

  sums <- aggregate(value ~ observation_id, age, sum)
  if (any(abs(sums$value - 1) > 0.001)) {
    stop("A GOA conditional age composition does not sum to one.", call. = FALSE)
  }

  age$weighted <- age$value * age$sample_size
  key <- aggregate(weighted ~ type + fleet + survey + length_bin_lower +
                     length_bin_upper + age, age, sum)
  totals <- aggregate(weighted ~ type + fleet + survey + length_bin_lower +
                        length_bin_upper, key, sum)
  join_key <- function(data) do.call(paste, c(data[
    c("type", "fleet", "survey", "length_bin_lower", "length_bin_upper")
  ], sep = "\034"))
  key$age_given_length <- key$weighted / totals$weighted[match(join_key(key), join_key(totals))]
  key
}

.goa_age_shares <- function(inputs, key, type, fleet, survey, years, ages) {
  lengths <- inputs[inputs$type == type & inputs$measure == "proportion_at_length" &
                      inputs$fleet == fleet & inputs$survey == survey &
                      inputs$year %in% years, , drop = FALSE]
  keys <- key[key$type == type & key$fleet == fleet & key$survey == survey, , drop = FALSE]
  if (!nrow(lengths) || !nrow(keys)) {
    stop("GOA length and age compositions do not overlap for this stream.", call. = FALSE)
  }

  bins <- aggregate(value ~ year + length_bin, lengths, sum)
  output <- lapply(sort(unique(bins$year)), function(year) {
    p_length <- bins[bins$year == year, c("length_bin", "value"), drop = FALSE]
    p_length$value <- p_length$value / sum(p_length$value)
    p_key <- keys[, c("length_bin_lower", "age", "age_given_length"), drop = FALSE]
    names(p_key)[names(p_key) == "length_bin_lower"] <- "length_bin"
    coverage <- sum(p_length$value[p_length$length_bin %in% p_key$length_bin])
    if (!is.finite(coverage) || coverage < 0.95) {
      stop("Pooled GOA age-at-length data cover less than 95% of a length composition.",
           call. = FALSE)
    }
    joint <- merge(p_length, p_key, by = "length_bin")
    joint$joint <- joint$value * joint$age_given_length
    shares <- aggregate(joint ~ age, joint, sum)
    shares$joint <- shares$joint / sum(shares$joint)
    if (!setequal(shares$age, ages) || any(!is.finite(shares$joint))) {
      stop("GOA age-length-key reconstruction does not cover every model age.",
           call. = FALSE)
    }
    data.frame(year = year, fleet = fleet, survey = survey, age = shares$age,
               proportion = shares$joint)
  })
  do.call(rbind, output)
}

translate_stock <- function(source) {
  assessment_id <- source$assessment$assessment_id[[1]]
  years <- 2007:2025
  ages <- 1:10
  survey <- "NMFS bottom-trawl survey"
  inputs <- source$inputs[source$inputs$assessment_id == assessment_id, , drop = FALSE]
  outputs <- source$outputs[source$outputs$assessment_id == assessment_id, , drop = FALSE]
  inputs$fleet[is.na(inputs$fleet)] <- ""
  inputs$survey[is.na(inputs$survey)] <- ""

  age_key <- .goa_age_length_key(inputs, years)
  catch_shares <- do.call(rbind, lapply(unique(inputs$fleet[
    inputs$type == "catch" & inputs$measure == "proportion_at_length"
  ]), function(fleet) {
    .goa_age_shares(inputs, age_key, "catch", fleet, "", years, ages)
  }))
  index_shares <- .goa_age_shares(
    inputs, age_key, "index", "", survey, years, ages
  )

  growth <- outputs[outputs$type == "biology" &
                      outputs$measure %in% c("growth_length_at_age", "growth_rate"), , drop = FALSE]
  length_at_age_1 <- growth$value[growth$measure == "growth_length_at_age" & growth$age == 1]
  length_at_age_10 <- growth$value[growth$measure == "growth_length_at_age" & growth$age == 10]
  growth_rate <- growth$value[growth$measure == "growth_rate"]
  if (length(length_at_age_1) != 1L || length(length_at_age_10) != 1L ||
      length(growth_rate) != 1L) {
    stop("GOA growth estimates are incomplete in outputs.csv.", call. = FALSE)
  }
  growth_factor <- exp(-growth_rate * (max(ages) - min(ages)))
  length_infinity <- (length_at_age_10 - length_at_age_1 * growth_factor) /
    (1 - growth_factor)
  length <- length_infinity - (length_infinity - length_at_age_1) *
    exp(-growth_rate * (ages - min(ages)))

  assumptions <- source$assumptions[source$assumptions$assessment_id == assessment_id, , drop = FALSE]
  parameter <- function(setting) {
    value <- as.numeric(assumptions$value[assumptions$component == "biology" &
                                           assumptions$setting == setting])
    if (length(value) != 1L || !is.finite(value)) {
      stop("Missing GOA biological parameter: ", setting, call. = FALSE)
    }
    value
  }
  weight <- parameter("weight_length_coefficient") *
    length^parameter("weight_length_exponent")
  maturity <- plogis(-parameter("maturity_logistic_slope") *
                       (length - parameter("maturity_length_50_percent_cm")))

  catch_totals <- inputs[inputs$type == "catch" &
                           inputs$measure == "total_biomass" &
                           inputs$year %in% years, , drop = FALSE]
  if (anyDuplicated(paste(catch_totals$year, catch_totals$fleet))) {
    stop("GOA catch biomass has duplicate year-fleet totals.", call. = FALSE)
  }
  catch_rows <- do.call(rbind, lapply(seq_len(nrow(catch_totals)), function(i) {
    total <- catch_totals[i, , drop = FALSE]
    share <- catch_shares[catch_shares$year == total$year &
                            catch_shares$fleet == total$fleet, , drop = FALSE]
    if (!setequal(share$age, ages)) stop("GOA catch age shares are incomplete.", call. = FALSE)
    mean_weight <- sum(share$proportion[match(ages, share$age)] * weight)
    tonnes_to_kg <- if (tolower(total$unit) == "t") 1000 else if (tolower(total$unit) == "kg") 1 else NA_real_
    if (!is.finite(tonnes_to_kg)) stop("Unsupported GOA catch biomass unit.", call. = FALSE)
    data.frame(year = total$year, age = ages,
               numbers = total$value * tonnes_to_kg / mean_weight *
                 share$proportion[match(ages, share$age)],
               fleet = total$fleet)
  }))
  catch_rows <- aggregate(numbers ~ year + age, catch_rows, sum)
  catch_totals_by_year <- aggregate(numbers ~ year, catch_rows, sum)

  index_totals <- inputs[inputs$type == "index" & inputs$survey == survey &
                           inputs$measure == "total_numbers" &
                           inputs$year %in% years, , drop = FALSE]
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
    data.frame(year = total$year, age = ages,
               numbers = total$value * index_multiplier[[i]] *
                 share$proportion[match(ages, share$age)])
  }))
  m_outputs <- outputs[outputs$type == "mortality" &
                         outputs$measure == "natural_mortality_at_age", , drop = FALSE]
  m_values <- m_outputs$value[match(
    paste(rep(years, each = length(ages)), rep(ages, times = length(years))),
    paste(m_outputs$year, m_outputs$age)
  )]

  row <- function(type, measure, basis = "", fleet = "", survey = "",
                  year = NA_real_, age = NA_real_, value, unit,
                  sampling_time = NA_real_, source_reference, transformation,
                  notes) {
    data.frame(
      assessment_id = assessment_id, type = type, measure = measure,
      basis = basis, fleet = fleet, survey = survey, sex = "", region = "",
      season = "", year = year,
      year_basis = ifelse(is.na(year), "", "calendar_year"), age = age,
      value = value, unit = unit, sampling_time = sampling_time,
      source_type = "reconstructed_source_input",
      source_reference = source_reference, transformation = transformation,
      notes = notes, stringsAsFactors = FALSE
    )
  }
  pad <- function(rows) {
    out <- inputs[rep(NA_integer_, nrow(rows)), , drop = FALSE]
    rownames(out) <- NULL
    for (name in names(rows)) out[[name]] <- rows[[name]]
    out
  }

  references <- function(rows) paste(unique(rows$source_reference), collapse = "; ")
  catch_source <- inputs[inputs$type == "catch" &
                           inputs$measure %in% c("total_biomass", "proportion_at_length",
                                                 "conditional_proportion_at_age"), , drop = FALSE]
  index_source <- inputs[inputs$type == "index" & inputs$survey == survey &
                           inputs$measure %in% c("total_numbers", "log_index_sd",
                                                 "proportion_at_length",
                                                 "conditional_proportion_at_age"), , drop = FALSE]
  growth_source <- outputs[outputs$measure %in% c("growth_length_at_age", "growth_rate"), , drop = FALSE]
  biology_reference <- references(assumptions[assumptions$component == "biology", , drop = FALSE])

  reconstructed <- list(
    row("catch", "total_numbers", year = catch_totals_by_year$year,
        age = NA_real_, value = catch_totals_by_year$numbers,
        unit = "fish", source_reference = references(catch_source),
        transformation = "Biomass converted to fish using reconstructed age proportions and weight-at-age.",
        notes = "Fleet catches are aggregated after applying fleet-specific pooled age-length keys; 2025 uses the key pooled from observed ages 2007-2024."),
    row("catch", "proportion_at_age", basis = "proportion_numbers",
        year = catch_rows$year, age = catch_rows$age,
        value = catch_rows$numbers /
          catch_totals_by_year$numbers[match(catch_rows$year,
                                             catch_totals_by_year$year)],
        unit = "proportion", source_reference = references(catch_source),
        transformation = "P(length) multiplied by sample-size-weighted P(age|length), then normalized by year and fleet before fleet aggregation.",
        notes = "Native bin labels are matched as categories. Original catch and composition rows remain unchanged."),
    row("index", "total_numbers", survey = survey,
        year = index_totals$year, age = NA_real_,
        value = index_totals$value, unit = index_totals$unit,
        sampling_time = index_totals$sampling_time,
        source_reference = references(index_source),
        transformation = "Native aggregate index retained before age allocation.",
        notes = "Only the age-sampled NMFS bottom-trawl survey is translated."),
    row("index", "proportion_at_age", basis = "proportion_numbers",
        survey = survey, year = index_rows$year, age = index_rows$age,
        value = index_rows$numbers /
          (index_totals$value[match(index_rows$year, index_totals$year)] *
             index_multiplier[match(index_rows$year, index_totals$year)]),
        unit = "proportion", sampling_time = 0.5,
        source_reference = references(index_source),
        transformation = "P(length) multiplied by the sample-size-weighted pooled P(age|length) and normalized.",
        notes = "Original SS3 composition likelihoods, ageing-error likelihood and covariance are not reproduced."),
    row("index", "log_index_sd", survey = survey,
        year = index_source$year[index_source$measure == "log_index_sd"],
        value = index_source$value[index_source$measure == "log_index_sd"],
        unit = "log scale", source_reference = references(index_source),
        transformation = "Native log-scale index SD retained for every age allocated from the aggregate survey index.",
        notes = "The same aggregate-index SD is applied to each reconstructed age row."),
    row("weight", "weight_at_age", basis = "kg_per_fish",
        age = ages, value = weight, unit = "kg per fish",
        source_reference = references(growth_source),
        transformation = "Von Bertalanffy curve anchored to reported mean lengths at ages 1 and 10; native weight-length curve applied to mean length.",
        notes = "The reported length-at-age SD is not integrated; age 10 is the plus group."),
    row("maturity", "maturity_at_age", basis = "proportion",
        age = ages, value = maturity, unit = "proportion mature",
        source_reference = biology_reference,
        transformation = "Native Stock Synthesis length-logistic maturity evaluated at reconstructed mean length-at-age.",
        notes = "The female SSB convention applies a 0.5 multiplier in the tinyAM translation."),
    row("M", "natural_mortality_at_age", basis = "per_year",
        year = rep(years, each = length(ages)),
        age = rep(ages, times = length(years)),
        value = as.numeric(m_values),
        unit = "per year", source_reference = references(outputs[
          outputs$type == "mortality" & outputs$measure == "natural_mortality_at_age", , drop = FALSE
        ]),
        transformation = "Accepted age-invariant M estimates expanded over age and year.",
        notes = "The source estimates M; this approximation fixes its published 2014-2016 block estimates." )
  )
  translation_inputs <- do.call(rbind, lapply(reconstructed, pad))
  if (any(!is.finite(translation_inputs$value[translation_inputs$type == "M"]))) {
    stop("GOA accepted M does not cover the selected fit years and ages.", call. = FALSE)
  }

  obs <- database_to_tam_obs(
    assessment_id, translation_inputs, years = years, ages = ages,
    weight_survey = "", index_weight_source = "stock",
    sampling_times = c("NMFS bottom-trawl survey" = 0.5),
    surveys = survey, maturity_multiplier = 0.5,
    assumptions = source$assumptions
  )

  comparison_outputs <- outputs[outputs$type == "biomass" &
                                  outputs$measure == "SSB" &
                                  outputs$year %in% years, , drop = FALSE]
  if (any(comparison_outputs$unit != "t")) {
    stop("GOA SSB conversion expects source tonnes.", call. = FALSE)
  }
  comparison_outputs$value <- comparison_outputs$value * 1000
  comparison_outputs$unit <- "kg"
  comparison_outputs$notes <- paste(
    comparison_outputs$notes,
    "Converted tonnes to kg for comparison with tinyAM; both use female SSB convention."
  )

  list(
    years = years,
    ages = ages,
    age_plus_group = 10,
    obs = obs,
    comparison_outputs = comparison_outputs,
    comparison_scales = c(ssb = 1e-6),
    comparison_aggregates = "ssb",
    settings = list(
      N_settings = list(process = "off", init = "exp"),
      F_settings = list(process = "iid", mu_form = NULL),
      M_settings = list(process = "off", mu_form = NULL,
                        mu_supplied = ~ M_assumption),
      catch_settings = list(sd_form = ~ 1, fill_missing = FALSE),
      index_settings = list(q_form = ~ 1,
                            sd_form = ~ 1,
                            sd_supplied = ~ relative_sd,
                            fill_missing = FALSE)
    ),
    background = c(
      "### Gulf of Alaska Pacific cod: accepted 2026 update, Model 24.0",
      "",
      "| Component | Accepted assessment | tinyAM representation | Reason for difference |",
      "|---|------|------|------|",
      "| Years | Model 24.0 covers 1977-2025. | Fit 2007-2025. | Fishery age-at-length samples begin in 2007; the short window avoids extrapolating the age-length relationship across earlier decades. |",
      "| Ages | Ages 1-10, with age 10 as the plus group. | Ages 1-10, with age 10 as the plus group. | The assessment's age range is retained. |",
      "| N | Annual age-1 recruitment and an estimated population process. | Exponential initialization and no N process. | This is a simpler abundance model and does not reproduce SS3's initial-equilibrium or recruitment structure. |",
      "| F | Three fisheries with fleet selectivity and time-varying parameters. | Aggregate catch observations and independent age-year F process states. | Fleet-specific selectivity and the source F parameterization are not reproduced. IID is retained as a simple illustrative fit; an RW alternative also converged, so this does not identify a preferred process. |",
      "| M | Age-constant M is estimated, with a higher 2014-2016 block. | The published 0.50 and 0.84 block estimates are supplied as fixed M. | This conditions on the accepted estimates and ignores their uncertainty. |",
      "| Catch | Fishery biomass, length proportions, and age proportions conditional on length. | Fleet-specific sample-size-weighted age-length keys reconstruct age proportions; weight-at-age converts biomass to numbers before fisheries are combined. | The source observations are retained; tinyAM fits the reconstructed aggregate catch-at-age with one error scale. |",
      "| Index | Bottom-trawl and longline aggregate indices with length and conditional age compositions. | Bottom-trawl index only, allocated by a pooled age-length key; supplied aggregate log SD is copied to each age. | Longline has no age-composition series. The source likelihood covariance and conditional-composition likelihoods are not represented; month 7 is mapped to midyear (0.5) because within-month timing is unavailable. |",
      "| Weights and maturity | One-sex growth, weight-length and length-logistic maturity parameters. | Reported lengths at ages 1 and 10 anchor a von Bertalanffy curve; the native weight and maturity functions are evaluated at mean length-at-age. | Growth variation is ignored; female SSB is represented by multiplying maturity by 0.5. |",
      "",
      "The pooled age-length keys use native bin labels and only source observations from 2007 onward. They cover at least 95% of each selected annual length composition; missing-bin mass is renormalized. The 2025 catch uses the key pooled from age samples through 2024.",
      "",
      "Only female spawning biomass is compared. Accepted total biomass includes age 0, and the accepted recruitment is age 0 while tinyAM recruitment is age 1, so those metrics are not compared. Accepted N-at-age and F-at-age surfaces were not recoverable from the compact fitted report."
    )
  )
}
