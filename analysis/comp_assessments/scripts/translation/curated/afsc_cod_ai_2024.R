## Source readiness ----
background <- c(
  "### Aleutian Islands Pacific cod: official 2024 Model 24.1","",print_sources(source$assessment),"",
  "The public native fit predates and differs from the final November production model. Final length/conditional-age arrays and biological weight surfaces are not verified.",
  "Final published catch, survey biomass, M and aggregate outputs are preserved. Preliminary fitted N/F/weights are not used to fill missing final inputs or states. Reported spawning-biomass scaling also needs confirmation."
)
cli::cli_abort(c(
  "AI cod translation is blocked pending matching final production inputs and biological weights.",
  "i"="Public September M24_1 output differs from the final November assessment. See source_reviews/afsc_cod_ai.md."
))
