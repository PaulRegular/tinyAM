# Stock-specific conversion: retain the joint-sex denominator before selecting females.
.canary_female_numbers <- function(inputs, type) {
  stopifnot(type %in% c("catch", "index"))
  p <- inputs[inputs$type == type & inputs$measure == "proportion_at_age", ]
  totals <- inputs[inputs$type == type & inputs$measure == "total_biomass", ]
  w <- inputs[inputs$type == if (type == "catch") "catch_weight" else "weight", ]
  out <- lapply(split(p, p$observation_id), function(sample) {
    year <- unique(sample$year)
    stream <- unique(if (type == "catch") sample$fleet else sample$survey)
    if (length(year) != 1L || length(stream) != 1L ||
      !setequal(sample$sex, c("female", "male")) ||
      any(sample$basis != "proportion_numbers") || anyDuplicated(paste(sample$sex, sample$age)) ||
      any(!is.finite(sample$value) | sample$value < 0) ||
      abs(sum(sample$value) - 1) > 1e-5) {
      cli::cli_abort("Canary conversion needs one valid joint-sex composition per sample.")
    }
    weights <- if (type == "catch") w[!is.na(w$fleet) & w$fleet == stream, ] else
      w[!is.na(w$survey) & w$survey == stream, ]
    key <- function(x) paste(x$sex, x$age)
    if (anyDuplicated(key(weights))) cli::cli_abort("Canary source weights have duplicate sex-age cells.")
    weight <- weights$value[match(key(sample), key(weights))]
    if (any(!is.finite(weight) | weight <= 0)) {
      cli::cli_abort("Canary conversion needs positive weights for every sex-age cell.")
    }
    total <- totals[totals$year == year, ]
    if (type == "index") total <- total[total$survey == stream, ]
    expected <- if (type == "catch") 2L else 1L
    if (nrow(total) != expected || any(total$unit != "t") ||
      (type == "catch" && !setequal(total$fleet, c("TRAWL_FISHERY", "OTHER_FISHERY"))) ||
      any(!is.finite(total$value) | total$value < 0)) {
      cli::cli_abort("Canary conversion needs matching original catch or survey biomass totals.")
    }
    female <- sample[sample$sex == "female", ]
    female$value <- 1000 * sum(total$value) * female$value / sum(sample$value * weight)
    female$measure <- "numbers_at_age"
    female$basis <- "numbers"
    female$unit <- "fish"
    if (type == "catch") female$fleet <- "POOLED_CATCH"
    female$source_type <- "derived_tinyAM_observation"
    female$transformation <- paste(
      "Biomass(t)*1000*p(female,age)/sum_over_both_sexes_and_all_source_ages(p*weight_kg).",
      if (type == "catch") "Pooled trawl+other biomass uses the sampled trawl joint-sex profile as an explicit allocation approximation." else
        "Use the matching survey biomass, joint-sex composition and source fleet weights.")
    female$notes <- paste(female$notes, female$transformation,
      "Only sampled years are converted; no composition, sample size or population state is invented.")
    female
  })
  do.call(rbind, out)
}
