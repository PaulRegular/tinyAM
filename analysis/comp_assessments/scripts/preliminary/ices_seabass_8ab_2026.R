fit_stage <- "source_readiness"
background <- c(
  "### Seabass in ICES divisions 8.a-b: detailed model inputs are unavailable", "",
  print_sources(source$assessment), "",
  "Aggregate official outputs are available, but the run-matched age or length observations, biological schedules and survey timing have not been recovered.",
  "Those quantities are needed to construct age-specific catch and survey observations and define the spawning biomass comparison. No inputs are inferred from the aggregate outputs.",
  "No fit is attempted."
)
cli::cli_abort(c(
  "Seabass in divisions 8.a-b is blocked before fitting because run-matched observations, biology and survey timing are unavailable.",
  "x" = "Aggregate assessment outputs do not supply the age-structured inputs needed by tinyAM.",
  "i" = "Reopen if the detailed run-matched source files are recovered."
))
