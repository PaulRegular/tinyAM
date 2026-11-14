fit_stage <- "source_readiness"
background <- c(
  "### Gulf of Alaska rex sole: source inputs do not support an age-structured fit", "",
  print_sources(source$assessment), "",
  "The final-run numerical age-composition inputs are unavailable; the report-linked composition appendix is marked draft. The report does not provide an age-specific survey index, and compatible run-matched age biology is incomplete.",
  "Reconstructing marginal ages from conditional age-at-length samples would require unsupported mappings. No observations or fit are constructed."
)
cli::cli_abort(c(
  "Rex sole is blocked before fitting because final-run age observations and biology are incomplete.",
  "x" = "The draft composition appendix and conditional age-at-length data do not verify final-run age-specific observations.",
  "i" = "See the 2026 SAIO source review and Step 5 translation-readiness note."
))
