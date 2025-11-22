
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
    mu_supplied = ~I(0.3),
    mean_ages = 5:14
  ),
  catch_settings = list(
    sd_form = ~ 1,
    fill_missing = TRUE
  ),
  index_settings = list(
    sd_form = ~ 1,
    q_form = ~ q_block,
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
    mu_supplied = ~I(0.3),
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
    mu_supplied = ~I(0.3),
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

## To age 20
models_20plus <- lapply(seq_along(models), function(i) {
  m <- models[[i]]
  m$call$M_settings$age_breaks = c(3, 5, 7, 9, 11, 13, 15, 17, 20)
  update(m, ages = 2:20)
})
names(models_20plus) <- names(models)

saveRDS(models_20plus, file = "analysis/comp_retro/outputs/001_models_20plus.rds")

vis_tam(
  models_20plus,
  output_file = "analysis/comp_retro/outputs/001_models_20plus.html"
)

