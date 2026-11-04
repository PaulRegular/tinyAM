fit_stage <- "source_readiness"
background <- c(
  "### West Coast yelloweye rockfish: source inputs do not support a complete fit", "",
  print_sources(source$assessment), "",
  "The assessment has California and Oregon/Washington components. Positive California trawl removals lack age samples, and a numbers-based recreational catch stream has an unresolved multiplier.",
  "Conditional age-at-length samples do not fill those missing catch ages. Compatible catch weights at age are not available; the weight surfaces in SS3 reports are model-generated.",
  "The source reports spawning output in eggs rather than SSB. These definition differences would need separate treatment even if complete catch inputs were recovered. No partial fit is attempted."
)
cli::cli_abort(c(
  "Yelloweye rockfish is blocked before fitting because source catch-at-age is incomplete and catch units are unresolved.",
  "x" = "The active California trawl removals have no age composition; one numbers stream has no verified multiplier.",
  "i" = "See the 2025 source assessment and the translation notes for the documented limits."
))
