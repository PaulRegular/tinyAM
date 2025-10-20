
## Exploratory analysis of retro patterns of a model with cohort deviations or
## M deviations

library(RTMB)
library(tinyAM)
library(plotly)

## Terminal fit ----

## More age blocks appear to work better than fewer
cod_obs <- tinyAM::cod_obs
cod_obs$catch$F_a_block <- cut_ages(cod_obs$catch$age, c(0:14))
cod_obs$weight$collapse <- ifelse(cod_obs$weight$year %in% 1991:1994, 1, 0)

## Couple F and M processes across the first few years to stabilize
## early estimates and reduce confounding with initial abundance
year_breaks <- c(1983, 1986:2024)

sam_style1 <- fit_tam(
  cod_obs,
  years = 1983:2024,
  ages = 2:14,
  N_settings = list(
    process = "iid",
    init_N0 = FALSE
  ),
  F_settings = list(
    process = "approx_rw",
    mu_form = NULL,
    mean_ages = 5:14,
    year_breaks = year_breaks
  ),
  M_settings = list(
    process = "off",
    assumption = ~I(0.3),
    mean_ages = 5:14,
    year_breaks = year_breaks
  ),
  obs_settings = list(
    q_form = ~ q_block,
    sd_form = ~ sd_obs_block,
    fill_missing = TRUE
  ),
  proj_settings = list(
    n_proj = 5,
    n_mean = 3,
    F_mult = 1
  )
)
sam_style1$sdrep
sam_style1$opt$objective

sam_style2 <- update(
  sam_style1,
  F_settings = list(
    process = "ar1",
    mu_form = ~F_a_block,
    mean_ages = 5:14,
    year_breaks = year_breaks
  )
)
sam_style2$sdrep

ncam_style1 <- update(
  sam_style1,
  N_settings = list(
    process = "off",
    init_N0 = TRUE
  ),
  F_settings = list(
    process = "approx_rw",
    mu_form = NULL,
    mean_ages = 5:14,
    year_breaks = year_breaks
  ),
  M_settings = list(
    process = "ar1",
    mu_form = ~ 0 + collapse,
    assumption = ~I(0.3),
    age_breaks = c(2:8, 14),
    mean_ages = 5:14,
    year_breaks = year_breaks
  )
)
ncam_style1$sdrep
ncam_style1$opt$objective

ncam_style2 <- update(
  ncam_style1,
  F_settings = list(
    process = "approx_rw",
    mu_form = NULL,
    mean_ages = 5:14,
    year_breaks = year_breaks
  ),
  M_settings = list(
    process = "approx_rw",
    mu_form = NULL,
    assumption = ~I(0.3),
    age_breaks = c(2, 14),
    mean_ages = 5:14,
    year_breaks = year_breaks
  )
)
ncam_style2$sdrep


vis_tam(
  list(
    sam_style1 = sam_style1,
    sam_style2 = sam_style2,
    ncam_style1 = ncam_style1,
    ncam_style2 = ncam_style2
  ),
  output_file = "analysis/comp_retro/001_models.html"
)



