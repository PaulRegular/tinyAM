fit_stage <- "source_readiness"
background <- c(
  "### West Greenland offshore cod: run-matched model inputs are unavailable", "",
  print_sources(source$assessment), "",
  "The named 2026 WOSC SAM fitted object could not be retrieved from the source host. Numeric observation arrays and survey timing therefore remain unverified.",
  "Earlier advice is not substituted for the selected 2026 assessment run.",
  "No observations, source-based starts or fit are created without the run-matched object."
)
cli::cli_abort(c(
  "West Greenland offshore cod is blocked before fitting because the named 2026 SAM object and its numeric inputs are unavailable.",
  "x" = "Earlier advice cannot replace the selected assessment run.",
  "i" = "Reopen when the named object and survey timing can be verified."
))
