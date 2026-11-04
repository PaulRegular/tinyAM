fit_stage <- "source_readiness"
background <- c(
  "### Greenland halibut in NAFO Divisions 2+3KLMNO: no population fit is available", "",
  print_sources(source$assessment), "",
  "The current event is a survey-based management procedure, not an applied population assessment with fitted age-specific states.",
  "A separate non-applied 2023 MSE state-space model does not provide the stock weights, maturity and numeric survey timing needed for a defensible tinyAM fit. Its inputs are not combined with the current management-procedure event.",
  "No observations or fit are constructed from these mismatched sources."
)
cli::cli_abort(c(
  "Greenland halibut is blocked before fitting because the selected event is a management procedure, not a population fit.",
  "x" = "The separate non-applied MSE model lacks required biological schedules and numeric survey timing.",
  "i" = "See the selected assessment and source review before considering another event."
))
