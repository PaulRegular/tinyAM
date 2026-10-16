## Source readiness ----
background <- c(
  "### BSAI arrowtooth flounder: official 2022 Model 18.9","",print_sources(source$assessment),"",
  "Original sex-specific biology and historical N, biomass, recruitment and recoverable cohort F are retained. The 2024 document is an off-cycle projection update.",
  "Numerical fishery length compositions and their original age–length key, survey age arrays and timing fractions are not verified. Fitted N/F are not used to manufacture observations."
)
cli::cli_abort(c(
  "Arrowtooth translation is blocked: numerical likelihood compositions, an observed age–length key and survey timing are missing.",
  "i"="Recovered fitted N/F are assessment outputs, not catch/survey data. See source_reviews/afsc_arrowtooth_flounder_bsai.md."
))
