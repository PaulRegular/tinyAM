fit_stage <- "source_readiness"
background <- c(
  "### Southern Gulf of St. Lawrence yellowtail flounder: numeric inputs are unavailable", "",
  print_sources(source$assessment), "",
  "The official summary identifies age-based assessment inputs, but numeric catch, survey and biological age schedules have not been recovered from the detailed assessment source.",
  "The available age-5+ biomass output is not a spawning-stock-biomass series. Without the source arrays and biological schedules, the observations and comparison definitions cannot be reconstructed.",
  "No fit or source-based starting values are attempted."
)
cli::cli_abort(c(
  "Yellowtail flounder is blocked before fitting because the numeric age inputs and biology are unavailable.",
  "x" = "Summary descriptions and age-5+ biomass do not provide the required observation arrays or SSB.",
  "i" = "Reopen this recipe if the detailed report or source data become accessible."
))
