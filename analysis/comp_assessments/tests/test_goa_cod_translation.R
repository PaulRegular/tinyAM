root <- file.path("analysis", "comp_assessments")
pkgload::load_all(".", quiet = TRUE)
source(file.path(root, "R", "run_assessment.R"))

database <- read_database()
source_data <- read_assessment("afsc_cod_goa_2026", database)
inputs_before <- source_data$inputs
stock <- new.env(parent = globalenv())
sys.source(file.path(root, "scripts", "translation", "stocks",
                     "afsc_cod_goa_2026.R"), stock)
translated <- stock$translate_stock(source_data)
obs <- translated$obs
settings <- translated$settings

tinyAM::check_obs(obs)
stopifnot(
  identical(source_data$inputs, inputs_before),
  identical(translated$years, 2007:2025),
  identical(translated$ages, 1:10),
  translated$age_plus_group == 10,
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
  grepl("M_assumption",
        paste(deparse(settings$M_settings$mu_supplied), collapse = "")),
  grepl("q_block",
        paste(deparse(settings$index_settings$q_form), collapse = "")),
  length(unique(obs$index$q_block)) > 1L,
  all(is.finite(obs$weight$M_assumption)),
  all(obs$weight$M_assumption[obs$weight$year %in% 2014:2016] == 0.84),
  all(obs$weight$M_assumption[!obs$weight$year %in% 2014:2016] == 0.50),
  setequal(unique(translated$comparison_outputs$measure),
           c("SSB", "total_biomass", "recruitment")),
  all(translated$comparison_outputs$unit[
    translated$comparison_outputs$measure %in% c("SSB", "total_biomass")
  ] == "t"),
  all(translated$comparison_outputs$unit[
    translated$comparison_outputs$measure == "recruitment"
  ] == "billion fish"),
  identical(unname(translated$comparison_scales[["ssb"]]), 1e-3),
  identical(unname(translated$comparison_scales[["biomass"]]), 1e-3),
  identical(unname(translated$comparison_scales[["recruitment"]]), 1e-9),
  identical(sort(unique(translated$comparison_outputs$year)), 2007:2025),
  identical(translated$comparison_definitions$ssb$status, "approximate")
)

index_totals <- source_data$inputs[
  source_data$inputs$type == "index" &
    source_data$inputs$survey == "NMFS bottom-trawl survey" &
    source_data$inputs$measure == "total_numbers" &
    source_data$inputs$year %in% translated$years, , drop = FALSE
]
index_sum <- aggregate(obs ~ year, obs$index, sum)
expected_index <- index_totals$value * 1000
stopifnot(isTRUE(all.equal(
  index_sum$obs[match(index_totals$year, index_sum$year)],
  expected_index, tolerance = 1e-8
)))

dat <- do.call(tinyAM::make_dat, c(
  list(obs = obs, years = translated$years, ages = translated$ages),
  settings
))
stopifnot(identical(dat$years, translated$years),
          identical(dat$ages, translated$ages))

reference <- database_to_tam_ref(
  "afsc_cod_goa_2026", database$outputs, obs = obs,
  years = translated$years, ages = translated$ages,
  age_plus_group = translated$age_plus_group,
  comparison_scales = translated$comparison_scales
)
m_reference <- reference$pop$M[reference$pop$M$year %in% translated$years & reference$pop$M$age %in% translated$ages, , drop = FALSE]
m_expected <- expand.grid(year = translated$years, age = translated$ages)
stopifnot(
  nrow(m_reference) == 19L * 10L,
  all(is.finite(m_reference$est)),
  setequal(paste(m_reference$year, m_reference$age),
           paste(m_expected$year, m_expected$age)),
  all(m_reference$est ==
        obs$weight$M_assumption[
          match(paste(m_reference$year, m_reference$age),
                paste(obs$weight$year, obs$weight$age))
        ])
)

m_fit <- list(
  dat = dat,
  pop = list(M = data.frame(
    year = obs$weight$year,
    age = obs$weight$age,
    est = obs$weight$M_assumption
  ))
)
m_comparison <- .assessment_percent_differences(
  m_fit, reference, scales = c(M = 1)
)
m_comparison <- m_comparison[m_comparison$metric == "M", ]
stopifnot(
  nrow(m_comparison) == 19L * 10L,
  all(m_comparison$comparison_status == "matched"),
  all(m_comparison$source == m_comparison$tinyAM)
)
cat("GOA Pacific cod translation tests passed.\n")
