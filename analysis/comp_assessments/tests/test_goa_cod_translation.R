source("analysis/comp_assessments/tests/helper_stock.R")
root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

database <- read_database()
source_data <- read_assessment("afsc_cod_goa_2026", database)
inputs_before <- source_data$inputs
translated <- .test_stock(source_data)
obs <- translated$obs
settings <- translated$settings
stock <- .assessment_stock(source_data, do_fit = FALSE, silent = TRUE)$env

annual_years <- sort(unique(stats::na.omit(stock$age_key$year)))
pooled_key <- stock$age_key[is.na(stock$age_key$year), , drop = FALSE]
pooled_catch <- stock$.goa_age_shares(
  stock$inputs, pooled_key, "catch", "Trawl fishery", "",
  stock$years, stock$ages
)
pooled_index <- stock$.goa_age_shares(
  stock$inputs, pooled_key, "index", "", stock$survey,
  stock$years, stock$ages
)
share_key <- function(x) paste(x$year, x$fleet, x$age)
trawl_shares <- stock$catch_shares[
  stock$catch_shares$fleet == "Trawl fishery", , drop = FALSE
]
pooled_catch_match <- match(share_key(trawl_shares), share_key(pooled_catch))

tinyAM::check_obs(obs)
stopifnot(
  identical(source_data$inputs, inputs_before),
  identical(translated$years, 2007:2025),
  identical(translated$ages, 1:10),
  translated$age_plus_group == 10,
  identical(annual_years, 2007:2024),
  nrow(pooled_key) > 0L,
  !anyNA(pooled_catch_match),
  any(abs(
    trawl_shares$proportion -
      pooled_catch$proportion[pooled_catch_match]
  ) > 1e-6),
  isTRUE(all.equal(
    trawl_shares$proportion[trawl_shares$year == 2025],
    pooled_catch$proportion[pooled_catch$year == 2025]
  )),
  isTRUE(all.equal(stock$index_shares, pooled_index)),
  nrow(obs$catch) == 19L * 10L,
  nrow(obs$index) == 10L * 10L,
  all(is.finite(obs$catch$obs)),
  all(is.finite(obs$index$obs)),
  identical(settings$N_settings$process, "off"),
  identical(settings$N_settings$init, "exp"),
  identical(settings$F_settings$process, "iid"),
  grepl("age", paste(deparse(settings$F_settings$mu_form), collapse = "")),
  identical(settings$M_settings$process, "off"),
  identical(settings$M_settings$mu_form, NULL),
  grepl(
    "M_assumption",
    paste(deparse(settings$M_settings$mu_supplied), collapse = "")
  ),
  grepl(
    "q_block",
    paste(deparse(settings$index_settings$q_form), collapse = "")
  ),
  length(unique(obs$index$q_block)) > 1L,
  all(is.finite(obs$pop_covariates$M_assumption)),
  all(obs$pop_covariates$M_assumption[obs$pop_covariates$year %in% 2014:2016] == 0.84),
  all(obs$pop_covariates$M_assumption[!obs$pop_covariates$year %in% 2014:2016] == 0.50),
  setequal(
    unique(translated$comparison_outputs$measure),
    c("SSB", "total_biomass", "recruitment")
  ),
  all(translated$comparison_outputs$unit[
    translated$comparison_outputs$measure %in% c("SSB", "total_biomass")
  ] == "t"),
  all(translated$comparison_outputs$unit[
    translated$comparison_outputs$measure == "recruitment"
  ] == "billion fish"),
  is.null(translated$comparison_scales),
  identical(sort(unique(translated$comparison_outputs$year)), 2007:2025),
  identical(translated$comparison_definitions$ssb$status, "approximate")
)

index_totals <- source_data$inputs[
  source_data$inputs$type == "index" &
    source_data$inputs$survey == "NMFS bottom-trawl survey" &
    source_data$inputs$measure == "total_numbers" &
    source_data$inputs$year %in% translated$years, ,
  drop = FALSE
]
index_sum <- aggregate(obs ~ year, obs$index, sum)
expected_index <- index_totals$value
stopifnot(isTRUE(all.equal(
  index_sum$obs[match(index_totals$year, index_sum$year)],
  expected_index,
  tolerance = 1e-8
)))

dat <- do.call(tinyAM::prepare_tam, c(
  list(data = obs, years = translated$years, ages = translated$ages),
  settings
))
stopifnot(
  identical(dat$years, translated$years),
  identical(dat$ages, translated$ages)
)

reference <- database_to_tam_ref(
  "afsc_cod_goa_2026", database$outputs,
  obs = obs,
  years = translated$years, ages = translated$ages,
  age_plus_group = translated$age_plus_group,
  comparison_scales = translated$comparison_scales
)
m_reference <- reference$pop$M[reference$pop$M$year %in% translated$years & reference$pop$M$age %in% translated$ages, , drop = FALSE]
m_expected <- expand.grid(year = translated$years, age = translated$ages)
stopifnot(
  nrow(m_reference) == 19L * 10L,
  all(is.finite(m_reference$est)),
  setequal(
    paste(m_reference$year, m_reference$age),
    paste(m_expected$year, m_expected$age)
  ),
  all(m_reference$est ==
    obs$pop_covariates$M_assumption[
      match(
        paste(m_reference$year, m_reference$age),
        paste(obs$weight$year, obs$weight$age)
      )
    ])
)

m_fit <- list(
  dat = dat,
  pop = list(M = data.frame(
    year = obs$weight$year,
    age = obs$weight$age,
    est = obs$pop_covariates$M_assumption
  ))
)
m_fit$rep <- m_fit$obs_pred <- m_fit$random_par <- list()
m_fit$fixed_par <- data.frame()
m_reference <- database_to_tam_ref(
  "afsc_cod_goa_2026", database$outputs,
  obs = obs,
  years = translated$years, ages = translated$ages,
  age_plus_group = translated$age_plus_group,
  comparison_scales = translated$comparison_scales, template = m_fit
)
m_comparison <- .assessment_percent_differences(m_reference)
m_comparison <- m_comparison[m_comparison$metric == "M", ]
stopifnot(
  nrow(m_comparison) == 19L * 10L,
  all(m_comparison$comparison_status == "matched"),
  all(m_comparison$source == m_comparison$tinyAM)
)
cat("GOA Pacific cod translation tests passed.\n")
