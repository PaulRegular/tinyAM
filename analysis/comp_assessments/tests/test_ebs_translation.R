source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

database <- read_database()
source_data <- read_assessment("afsc_pollock_ebs_2024", database)
translated <- .test_stock(source_data)
obs <- translated$obs
documentation <- print_sources(source_data$assessment)
stopifnot(length(documentation) >= 2L,
          all(documentation %in% translated$background),
          any(grepl(source_data$assessment$assessment_url[[1L]],
                    translated$background, fixed = TRUE)))
settings <- translated$settings

tinyAM::check_obs(obs)
dat <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = translated$years, ages = translated$ages),
  settings
))

survey_names <- sort(unique(as.character(obs$index$survey)))
age1_ats <- obs$index[
  obs$index$survey == "NMFS acoustic-trawl" & obs$index$age == 1,
  , drop = FALSE
]
ats_old_age1 <- age1_ats[age1_ats$year < 2024, , drop = FALSE]
high_age <- obs$index[obs$index$age >= 9, , drop = FALSE]
high_age_keys <- unique(high_age[c("survey", "age", "q_key")])
keys_by_age <- aggregate(
  q_key ~ survey + age, high_age_keys,
  function(x) length(unique(as.character(x)))
)
keys_by_survey <- aggregate(
  q_key ~ survey, high_age_keys,
  function(x) length(unique(as.character(x)))
)
m_values <- unique(obs$weight[c("age", "M_assumption")])
m_values <- m_values[order(m_values$age), ]

stopifnot(
  identical(translated$years, 1964:2024),
  identical(translated$ages, 1:15),
  translated$age_plus_group == 15,
  identical(dat$years, 1964:2024),
  identical(dat$ages, 1:15),
  identical(settings$N_settings$process, "iid"),
  identical(settings$N_settings$init, "exp"),
  identical(settings$F_settings$process, "rw"),
  identical(settings$M_settings$process, "off"),
  identical(settings$M_settings$mu_form, NULL),
  grepl("M_assumption",
        paste(deparse(settings$M_settings$mu_supplied), collapse = "")),
  identical(settings$catch_settings$fill_missing, FALSE),
  identical(settings$index_settings$fill_missing, FALSE),
  grepl("q_key", paste(deparse(settings$index_settings$q_form), collapse = "")),
  grepl("survey", paste(deparse(settings$index_settings$sd_form), collapse = "")),
  grepl("age", paste(deparse(settings$catch_settings$sd_form), collapse = ""), fixed = TRUE),
  grepl("I(age^2)",
        paste(deparse(settings$catch_settings$sd_form), collapse = ""), fixed = TRUE),
  identical(survey_names,
           sort(c("NMFS bottom-trawl VAST", "NMFS acoustic-trawl"))),
  !"NMFS acoustic-trawl age-1 index" %in% survey_names,
  nrow(ats_old_age1) > 0L,
  !any(age1_ats$year == 2024),
  nrow(high_age_keys) > 0L,
  all(high_age$q_age_block == "9+"),
  all(keys_by_age$q_key == 1L),
  all(keys_by_survey$q_key == 1L),
  length(unique(as.character(high_age_keys$q_key))) == length(survey_names),
  all(m_values$age == 1:15),
  isTRUE(all.equal(m_values$M_assumption, c(0.9, 0.45, rep(0.3, 13)))),
  identical(translated$comparison_age_groups$N[["10+"]], 10:15),
  identical(translated$comparison_age_groups$biomass_at_age[["3+"]], 3:15),
  identical(translated$comparison_definitions$ssb$status, "approximate")
)

cat("EBS pollock translation tests passed.\n")
