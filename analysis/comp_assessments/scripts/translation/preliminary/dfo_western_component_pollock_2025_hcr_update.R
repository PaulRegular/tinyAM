fit_stage <- "source_readiness"
background <- c(
  "### Western component pollock: the current event is an HCR update", "",
  print_sources(source$assessment), "",
  "The selected 2025 event updates harvest-control-rule information rather than providing a fitted population assessment suitable for translation.",
  "The detailed 2022 J_M4B model arrays and biological inputs have not been curated into the pinned database. The current survey index is not treated as model-estimated SSB.",
  "Inputs from the earlier event are not mixed into the 2025 update. No fit is attempted."
)
cli::cli_abort(c(
  "Western component pollock is blocked before fitting because the selected event is an HCR update and the detailed model inputs are not available in the pinned data.",
  "x" = "The 2025 survey index cannot stand in for model SSB.",
  "i" = "A translation requires a run-matched detailed assessment and its biological inputs."
))
