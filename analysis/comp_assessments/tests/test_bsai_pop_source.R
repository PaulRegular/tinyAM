source("analysis/comp_assessments/R/read_database.R")
x <- read_assessment("afsc_pop_bsai_2024")
i <- x$inputs
o <- x$outputs
stopifnot(x$stocks$ocean == "Northeast Pacific", x$stocks$charbonneau_id == "AFSC_alutus_BSAI",
  x$assessment$terminal_year == 2024, x$assessment$is_current,
  sum(i$type == "catch" & i$measure == "total_biomass") == 65L,
  setequal(unique(i$survey[!is.na(i$survey)]), c("AI_survey", "EBS_survey")),
  all(is.na(i$sampling_time)),
  sum(i$type == "weight") == 44L,
  i$value[i$type == "weight" & i$age == 3] == 94,
  all(i$unit[i$type == "weight"] == "g/fish"),
  !any(i$type == "M"), !any(i$type == "maturity"),
  sum(i$measure == "maturity_sample_size") == 44L,
  o$value[o$measure == "SSB" & o$year == 2024] == 359278,
  o$value[o$measure == "total_biomass" & o$year == 2024] == 864800,
  all(o$value[o$measure == "natural_mortality_at_age"] == 0.0505762),
  all(o$se[o$measure == "SSB"] > 0),
  all(is.na(o$se[o$measure == "fishing_mortality_at_age"])))
n <- o[o$measure == "numbers_at_age", ]
stopifnot(identical(sort(unique(n$age)), 3:40),
          all(n$age_group[n$age == 40] == "40+"))
cat("BSAI POP source values, units, groups and active-input distinctions passed.\n")
