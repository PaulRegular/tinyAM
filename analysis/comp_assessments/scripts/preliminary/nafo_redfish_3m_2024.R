fit_stage <- "source_readiness"
background <- c(
  "### NAFO Division 3M beaked redfish: age-specific inputs are not available", "",
  print_sources(source$assessment), "",
  "The latest recoverable full age-based assessment is the 2024 XSA for beaked redfish (*Sebastes mentella* and *S. fasciatus*), using ages 4–19+. The detailed assessment report SCR 24/032 is not accessible in the catalogue, and the available summary does not provide the annual age-specific commercial catch, survey index, stock weight or maturity matrices.",
  "Published total catch and short-term scenario outputs are preserved in SAIO, but they do not replace the missing age-specific observations or biology. Historical annual SSB, N-at-age and F-at-age are not numerically available; values shown only in figures are not digitized.",
  "No age-based tinyAM data, fit or Official–tinyAM comparison can be constructed without inventing inputs. Revisit if the detailed assessment tables or model files become accessible."
)
cli::cli_abort(c(
  "NAFO 3M beaked redfish is blocked before fitting because annual age-specific inputs are unavailable.",
  "x" = "Published total catch and scenario projections do not supply catch-at-age, survey indices, weights or maturity.",
  "i" = "No values are inferred from plots or fitted outputs; see the SAIO source review."
))
