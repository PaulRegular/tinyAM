fit_stage <- "source_readiness"
background <- c(
  "### Chilipepper rockfish: source inputs do not support an age-structured fit", "",
  print_sources(source$assessment), "",
  "The official likelihood uses catch biomass, aggregate indices, conditional age-at-length and length-composition data. It does not supply observed catch-at-age or an age-specific abundance index.",
  "Model-estimated age surfaces are not used as observations, and length-based biological surfaces are not converted to age without verified source mappings."
)
cli::cli_abort(c(
  "Chilipepper is blocked before fitting because its observed likelihood inputs are not age-specific.",
  "x" = "Model-estimated catch numbers and population surfaces cannot replace observed age data.",
  "i" = "See the 2025 SAIO source review and Step 5 translation-readiness note."
))
