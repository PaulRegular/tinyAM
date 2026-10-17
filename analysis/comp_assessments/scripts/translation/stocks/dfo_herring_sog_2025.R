## Source readiness ----
background <- c(
  "### Strait of Georgia Pacific herring: official 2025 DDM ensemble","",
  print_sources(source$assessment),"",
  "The native combined survey measures total spawning biomass. Survey age samples are absent, and tinyAM currently requires an abundance observation at each age.",
  "Original fishery samples, weights and survey biomass are preserved. No source-fitted N is used to invent survey ages. This case is blocked before model preparation or fitting.",
  "The detailed DDM ensemble ends in 2023. The separately reported interim 2025 DIM is not substituted for it."
)
cli::cli_abort(c(
  "SoG herring translation is blocked: the survey is aggregate spawning biomass with no observed survey ages.",
  "i"="tinyAM's current age-specific index likelihood cannot represent this observation. See source_reviews/dfo_herring_sog.md."
))
