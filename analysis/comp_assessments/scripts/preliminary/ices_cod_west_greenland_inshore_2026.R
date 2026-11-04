fit_stage <- "source_readiness"
background <- c(
  "### West Greenland inshore cod: run-matched model inputs are unavailable", "",
  print_sources(source$assessment), "",
  "The named 2026 N-WISC SAM fitted object could not be retrieved from the source host. Numeric observation arrays and survey timing therefore remain unverified.",
  "The southern West Greenland stock and earlier advice are separate assessment inputs; they are not substituted for the selected inshore run.",
  "No observations, source-based starts or fit are created without the run-matched object."
)
cli::cli_abort(c(
  "West Greenland inshore cod is blocked before fitting because the named 2026 SAM object and its numeric inputs are unavailable.",
  "x" = "The distinct southern stock or earlier advice cannot replace the selected assessment run.",
  "i" = "Reopen when the named object and survey timing can be verified."
))
