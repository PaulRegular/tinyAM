fit_stage <- "source_readiness"
background <- c(
  "### British Columbia yellowmouth rockfish: age observations are unavailable", "",
  print_sources(source$assessment), "",
  "The database contains aggregate indices and maturity information, but the native numeric age-composition arrays have not been recovered.",
  "Aggregate series and biological schedules do not determine catch-at-age or age-specific abundance indices. No age structure is inferred from fitted model outputs.",
  "No observations or fit are constructed while the source age arrays are unavailable."
)
cli::cli_abort(c(
  "Yellowmouth rockfish is blocked before fitting because the source age-composition arrays are unavailable.",
  "x" = "Aggregate indices and maturity do not supply catch-at-age or survey abundance-at-age.",
  "i" = "Reopen the recipe when numeric source age compositions are available."
))
