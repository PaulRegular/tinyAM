## Comparison readiness ----
background <- c(
  "### British Columbia Pacific ocean perch: official 2023 assessment", "",
  print_sources(source$assessment), "",
  "The original three-area inputs and fixed biology remain in the database. Historical official N, F, recruitment and a complete matching coastwide SSB series have not been recovered.",
  "Available outputs include constant sex-specific M, January 2024 biomass summaries and historical biomass draws for 3CD and 5DE. The public 5ABC posterior object ends in 2017. These do not support the intended coastwide comparison; the native spawning-output definition also remains unresolved.",
  "This case is blocked from fitting and comparison until matching numerical historical outputs are available. Earlier exploratory fits remain archived locally and in Git history. See source_reviews/dfo_pop_bc.md."
)
cli::cli_abort(c(
  "BC Pacific ocean perch comparison is blocked: matching coastwide historical official population outputs have not been recovered.",
  "i" = "Partial regional biomass, constant M and boundary summaries are insufficient. See source_reviews/dfo_pop_bc.md."
))
