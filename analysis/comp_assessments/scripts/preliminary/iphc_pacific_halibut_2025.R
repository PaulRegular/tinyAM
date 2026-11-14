fit_stage <- "source_readiness"
background <- c(
  "### IPHC Pacific halibut: public inputs do not support an age-structured fit", "",
  print_sources(source$assessment), "",
  "The public assessment tables provide aggregate indices and population outputs, but not the exact catch-age likelihood vectors, age-specific survey observations, or run-matched numeric weight and maturity surfaces.",
  "The native model input files are permission-restricted. No age allocation is inferred from aggregate indices or fitted outputs."
)
cli::cli_abort(c(
  "Pacific halibut is blocked before fitting because public age-specific observations and biology are unavailable.",
  "x" = "Aggregate index series cannot be allocated to age without source data.",
  "i" = "See the 2025 SAIO source review and Step 5 translation-readiness note."
))
