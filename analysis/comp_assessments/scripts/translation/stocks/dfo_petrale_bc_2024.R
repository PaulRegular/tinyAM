## Source readiness ----
background <- c(
  "### British Columbia petrale sole: official 2024 assessment", "",
  print_sources(source$assessment), "",
  "The official model is a two-sex, coastwide Stock Synthesis assessment combining three fixed-M models. Its catches end in 2023; the beginning-2024 state is reported separately.",
  "The database contains public female SSB posterior summaries for 1938-2024, published terminal/reference values, and the three reported M scalars. The public subset summaries are distinct from exact published full-ensemble quantiles.",
  "Original catch/survey age compositions, biological weights and maturity, survey timing, model ages and the final configuration remain unavailable from the searched public sources. Fitted SSB cannot replace these inputs.",
  "This is a partial summary/output record, not a current detailed assessment. Translation is blocked before preparation or fitting. No start parameters, catch-only substitute or borrowed U.S. inputs are used. See source_reviews/dfo_petrale_bc.md."
)
cli::cli_abort(c(
  "BC petrale translation is blocked: original catch/survey age compositions and biological inputs have not been recovered.",
  "i" = "Weight, maturity, timing and the model age grid must be verified before preparing a fit. Available SSB output and component M scalars are insufficient.",
  "i" = "See source_reviews/dfo_petrale_bc.md for public-source coverage and the missing-input inventory."
))
