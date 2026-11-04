fit_stage <- "source_readiness"
background <- c(
  "### Oregon Black Rockfish: source inputs do not support an age-structured fit", "",
  print_sources(source$assessment), "",
  "The Oregon shore fleet has positive removals but no age or length samples. The other catch-age information is conditional on length, and no compatible catch weights at age are available.",
  "The assessment does not provide age-specific survey observations. Its maturity-times-fecundity surface is a reproductive-output measure, not a maturity proportion.",
  "SS3-generated biological surfaces are not substituted for missing observations. No fit is attempted."
)
cli::cli_abort(c(
  "Oregon Black Rockfish is blocked before fitting because required age-specific catches and survey observations are unavailable.",
  "x" = "The shore catch is unsampled by age and length; age-specific survey data and compatible catch weights are also missing.",
  "i" = "See the 2023 Oregon assessment model inputs and the translation notes for the source limitations."
))
