## Source readiness ----
fit_stage <- "source_readiness"
background <- c(
  "### Irish Sea herring: annual-reference review required", "",
  print_sources(source$assessment), "",
  "Catch, ISAS abundance, stock weights and maturity are available for winter-ring ages 1–8+, but their annual origin must be reconciled before assigning survey or spawning fractions.",
  "The source rows are indexed as model_year and the associated assumption describes April–March. The HAWG SSB footnote instead describes January-to-September survival. These may refer to different catch and population conventions; they are not silently treated as one calendar year.",
  "The aggregate ISSS index has no age composition. No observations, timing values, initial parameters or optimizer result are invented while the annual-reference question remains open."
)
cli::cli_abort(c(
  "Irish Sea herring needs verification of its catch-year and population-state reference dates.",
  "i" = "See the annual-reference follow-up in SAIO's discovery/source_reviews/ices_irish_sea_herring_2026.md."
))
