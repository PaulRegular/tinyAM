## Source readiness ----
background <- c(
  "### Georges Bank cod (US assessment): official 2024 WHAM run 4", "",
  print_sources(source$assessment), "",
  "The database preserves the current combined US/Canadian catch totals, fixed M and reported 2014–2023 SSB, recruitment and fully selected F. These summaries do not provide the original age-composition observations or population state surfaces.",
  "Original fishery and three survey age inputs, current biological weight and maturity matrices, survey timing and fitted N/full F remain unrecovered. A fit cannot be prepared or compared from the available records.",
  "This source gap is documented in source_reviews/nefsc_cod_georges_bank.md. No older assessment inputs, fitted states or assumed biology are used to fill it."
)
cli::cli_abort(c(
  "Georges Bank cod translation is blocked: original catch/survey age inputs and biological matrices are unavailable.",
  "i" = "The aggregate catch/status table and fixed M do not supply the WHAM likelihood data required by tinyAM.",
  "i" = "See source_reviews/nefsc_cod_georges_bank.md for the recovered inputs and missing-source inventory."
))
