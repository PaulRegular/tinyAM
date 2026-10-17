## Comparison readiness ----
background <- c(
  "### British Columbia Pacific ocean perch: official 2023 assessment", "",
  print_sources(source$assessment), "",
  "The original three-area inputs and fixed biology remain in the database. Historical official N, F, recruitment and SSB series have not been recovered from the public production files or report tables.",
  "Available outputs are constant sex-specific M and a January 2024 biomass summary. These do not support the intended historical population comparison; the native spawning-output definition also remains unresolved.",
  "This case is blocked from fitting and comparison until matching numerical historical outputs are available. Earlier exploratory fits remain archived locally and in Git history. See source_reviews/dfo_pop_bc.md."
)
cli::cli_abort(c(
  "BC Pacific ocean perch comparison is blocked: historical official population outputs have not been recovered.",
  "i" = "Constant M and a projection-boundary biomass summary are insufficient. See source_reviews/dfo_pop_bc.md."
))
