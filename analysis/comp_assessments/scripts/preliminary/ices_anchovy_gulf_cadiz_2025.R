## Source readiness ----
background <- c(
  "### Gulf of Cadiz anchovy: official 2025 SS3 assessment", "",
  print_sources(source$assessment), "",
  "The source has quarterly fishing fleets, spawning on 1 April and age-0 settlement on 1 July. Reported recruitment is not January recruitment.",
  "Current catches and survey tables, fixed M/maturity and annual SSB/recruitment/age-3 F are retained. Exact empirical stock weights, complete quarterly compositions and their adjusted sample sizes remain missing.",
  "The older public MSE reference differs from the current official estimates and is not used to fill those inputs. No 2025 fishery ages or catch weights are manufactured."
)
cli::cli_abort(c(
  "Southern anchovy translation is blocked: current empirical stock weights and native quarterly compositions are missing.",
  "i" = "Raw 2024 catch weights and older MSE states do not supply the missing production inputs. See source_reviews/ices_anchovy_gulf_cadiz.md."
))
