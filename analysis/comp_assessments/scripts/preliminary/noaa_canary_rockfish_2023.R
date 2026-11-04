fit_stage <- "source_readiness"
background <- c(
  "### Pacific coast Canary Rockfish: source inputs do not support an age-structured fit", "",
  print_sources(source$assessment), "",
  "The native fishery age samples are conditional age-at-length records. The recorded length-bin keys do not provide a verified mapping sufficient to reconstruct marginal age observations and their sampling weights for tinyAM.",
  "An active recreational catch fleet has removals but no age samples, and compatible catch weights at age are unavailable. SS3-generated weight and maturity surfaces are not treated as observations.",
  "Native spawning output is in eggs, not SSB. No source-based conversion or partial fit is attempted."
)
cli::cli_abort(c(
  "Canary Rockfish is blocked before fitting because catch-at-age and catch-weight inputs cannot be reconstructed from the verified records.",
  "x" = "Conditional ages and unsampled catch cannot be replaced with fitted SS3 states or biology.",
  "i" = "See the 2023 assessment model inputs and the translation notes for the unresolved data mapping."
))
