
## TODO: think about supplying sd_index values

library(RTMB)
library(tinyAM)
library(plotly)

## Terminal fit ----

cod_obs <- tinyAM::cod_obs
cod_obs$weight$collapse <- ifelse(cod_obs$weight$year %in% 1991:1994, 1, 0)

## Questions:
## How should the random processes be modeled?
## - Should N deviations be iid or ar1?
## - Should M deviations be iid or ar1?
## - Should F deviations be approx rw or ar1?

## Conventions:
## Naming by deviation_process_deviation_process, e.g.:
## - N_iid_F_rw ~ SAM style model
## - M_ar1_F_rw ~ NCAM style model

N_iid_F_rw <- fit_tam(
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
    mean_ages = 5:14
  ),
  M_settings = list(
    process = "off",
    assumption = ~I(0.3),
    mean_ages = 5:14
  ),
  obs_settings = list(
    q_form = ~ q_block,
    sd_catch_form = ~ 1,
    sd_index_form = ~ 1,
    fill_missing = TRUE
  ),
  proj_settings = list(
    n_proj = 5,
    n_mean = 3,
    F_mult = 1
  )
)
N_iid_F_rw
N_iid_F_rw$opt$objective

N_iid_F_ar1 <- update(
  N_iid_F_rw,
  F_settings = list(
    process = "ar1",
    mu_form = ~F_a_block + F_y_block,
    mean_ages = 5:14
  )
)
N_iid_F_ar1

N_ar1_F_rw <- update(
  N_iid_F_rw,
  N_settings = list(
    process = "ar1",
    init_N0 = FALSE
  )
)
N_ar1_F_rw

M_ar1_F_rw <- update(
  N_iid_F_rw,
  N_settings = list(
    process = "off",
    init_N0 = TRUE
  ),
  F_settings = list(
    process = "approx_rw",
    mu_form = NULL,
    mean_ages = 5:14
  ),
  M_settings = list(
    process = "ar1",
    mu_form = NULL,
    assumption = ~I(0.3),
    age_breaks = c(3, 5, 7, 9, 11, 14), # Coupling ages to avoid convergence issues
    mean_ages = 5:14
  )
)
M_ar1_F_rw

M_ar1_F_ar1 <- update(
  M_ar1_F_rw,
  F_settings = list(
    process = "ar1",
    mu_form = ~F_a_block + F_y_block,
    mean_ages = 5:14
  )
)
M_ar1_F_ar1

M_iid_F_rw <- update(
  M_ar1_F_rw,
  M_settings = list(
    process = "iid",
    mu_form = NULL,
    assumption = ~I(0.3),
    age_breaks = c(3, 5, 7, 9, 11, 14),
    mean_ages = 5:14
  )
)
M_iid_F_rw

model_names <- c(
  "N_iid_F_rw",
  "N_iid_F_ar1",
  "N_ar1_F_rw",
  "M_ar1_F_rw",
  "M_ar1_F_ar1",
  "M_iid_F_rw"
)
models <- mget(model_names)

saveRDS(models, file = "analysis/comp_retro/outputs/001_models.rds")

vis_tam(
  models,
  output_file = "analysis/comp_retro/outputs/001_models.html"
)



