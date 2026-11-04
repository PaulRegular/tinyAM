fit_stage <- "source_readiness"
background <- c(
  "### West Coast Dover sole: source inputs do not support an age-structured fit", "",
  print_sources(source$assessment), "",
  "The assessment's removals include retained catch and discards, while the available catch-age samples describe retained fish. Discard age composition is not available to reconstruct total removals at age.",
  "The recorded mean-body-weight observations refer to discard samples, not compatible retained catch weights at age. The SS3 weight-at-age surface is model-generated and is not used as an observation input.",
  "The input includes a source-derived marginal WCGBTS age composition, but using it would still be an approximate survey translation and would not resolve the catch reconstruction. No partial fit is attempted."
)
cli::cli_abort(c(
  "Dover sole is blocked before fitting because total removals cannot be reconstructed at age from the verified source inputs.",
  "x" = "The source age samples do not represent discards, and compatible catch weights at age are unavailable.",
  "i" = "See the 2021 assessment inputs and the translation notes for the specific source gap."
))
