fit_stage <- "source_readiness"
background <- c(
  "### Greenland halibut: the 2026 advice event has no population fit", "",
  print_sources(source$assessment), "",
  "The current event is a survey-based management procedure, not an applied population assessment with fitted age-specific states.",
  "The 2026 advice event is a survey-based management procedure, not a fitted population model. Separate 2023 SCAA and SSM operating models were developed for MSE and were not used for annual TAC advice; their inputs and outputs are not combined with this event.",
  "No observations or fit are constructed from these mismatched sources."
)
cli::cli_abort(c(
  "The selected 2026 event is a management procedure, not an age-structured population fit.",
  "x" = "The separately imported 2023 SCAA and SSM MSE runs are distinct assessment events, not substitutes for current annual advice.",
  "i" = "No current-event age-structured states are available to translate."
))
