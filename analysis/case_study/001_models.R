

pkgload::load_all(".", quiet = TRUE)

dir.create("analysis/case_study/outputs", recursive = TRUE, showWarnings = FALSE)

F_dev <- fit_tam(
  data,
  years = 1983:2025,
  ages = 2:14,
  N_settings = list(
    process = "off",
    init = "exp"
  ),
  F_settings = list(
    process = "ar1",
    mu_form = ~ logistic(age) + F_y_block,
    mean_ages = 5:14
  ),
  M_settings = list(
    process = "off",
    mu_supplied = ~ I(0.3),
    mean_ages = 5:14
  ),
  catch_settings = list(
    sd_form = ~ 1
  ),
  index_settings = list(
    sd_form = ~ 1,
    q_form = ~ mono(q_block)
  )
)
F_dev

F_N_dev <- update(
  F_dev,
  N_settings = list(
    process = "ar1"
  )
)

F_M_dev <- update(
  F_dev,
  M_settings = list(
    process = "ar1",
    age_breaks = c(2, 5, 8, 14)
  )
)


AIC(F_dev, F_N_dev, F_M_dev)

cAIC(F_dev, F_N_dev, F_M_dev)


fits <- list(F_dev, F_N_dev, F_M_dev)
hind <- lapply(fits, fit_hindcast, folds = 15)
lapply(hind, `[[`, "hindcast_rmse")
hind[[1]]$mohns_rho |> subset(metric == "F_bar")
hind[[2]]$mohns_rho |> subset(metric == "F_bar")
hind[[3]]$mohns_rho |> subset(metric == "F_bar")





















pkgload::load_all(".", quiet = TRUE)

dir.create("analysis/case_study/outputs", recursive = TRUE, showWarnings = FALSE)

data <- tinyAM::cod_obs
data$pop_covariates <- data$pop_covariates |>
  tidyr::replace_na(list(log_capelin_per_cod = 0)) # i.e., capelin effect will not be applied if absent


## Build up the full model first that attempts to capture key changes in the fishery,
## changes in survey catchability, and variation in natural mortality. Then compare
## alternative choices.

full <- fit_tam(
  data,
  years = 1983:2025,
  ages = 2:14,
  N_settings = list(
    process = "off",
    init = "exp"
  ),
  F_settings = list(
    process = "ar1",
    mu_form = ~ logistic(age) + F_y_block,
    mean_ages = 5:14
  ),
  M_settings = list(
    process = "ar1",
    mu_supplied = ~ I(0.3),
    mu_form = ~ 0 + M_a_block:log_capelin_per_cod,
    age_breaks = c(2, 5, 8, 14),
    mean_ages = 5:14
  ),
  index_settings = list(
    sd_form = ~ 1,
    q_form = ~ mono(q_block) + smith_sound_effect
  )
)
full


# no_F_blocks <- update(
#   full,
#   F_settings = list(
#     mu_form = ~ 1
#   )
# )
# no_F_blocks
#
# no_smith_sound_effect <- update(
#   full,
#   index_settings = list(
#     q_form = ~ mono(q_block)
#   )
# )

no_capelin_effect <- update(
  full,
  M_settings = list(
    process = "ar1",
    mu_form = NULL
  )
)

# constant_M <- update(
#   full,
#   M_settings = list(
#     process = "off",
#     mu_form = NULL
#   )
# )

future::plan(future::multisession, workers = 10)
full_hind <- fit_hindcast(full, folds = 20)
no_cap_hind <- fit_hindcast(no_capelin_effect, folds = 20)
full_hind$hindcast_rmse
no_cap_hind$hindcast_rmse


cAIC(full, no_F_blocks, no_smith_sound_effect, no_capelin_effect, constant_M)

vis_tam(full, no_F_blocks, no_smith_sound_effect, no_capelin_effect, constant_M)









ar1_M <- update(
  base,
  M_settings = list(
    process = "ar1",
    age_breaks = c(2, 5, 8, 14),
    mean_ages = 5:14
  )
)

base_ss_q <- update(
  base,
  index_settings = list(
    q_form = ~mono(q_block) + smith_sound_effect
  )
)
base_ss_q

ar1_M_ss_q <- update(
  ar1_M,
  index_settings = list(
    q_form = ~mono(q_block) + smith_sound_effect
  )
)
ar1_M_ss_q

ar1_cap_M_ss_q <- update(
  ar1_M_ss_q,
  M_settings = list(
    mu_form = ~ 0 + log_capelin_per_cod
  )
)
ar1_cap_M_ss_q

base_cap_M_ss_q <- update(
  base_ss_q,
  M_settings = list(
    mu_form = ~ 0 + log_capelin_per_cod
  )
)
base_cap_M_ss_q




ar1_cap_M <- update(
  ar1_M,
  M_settings = list(
    mu_form = ~ 0 + log_capelin_per_cod
  )
)
ar1_cap_M

base_cap_M <- update(
  base,
  M_settings = list(
    mu_form = ~ 0 + log_capelin_per_cod
  )
)
base_cap_M




cAIC(base, ar1_M, base_ss_q, ar1_M_ss_q, base_cap_M_ss_q, ar1_cap_M_ss_q, base_cap_M, ar1_cap_M, details = TRUE)




model_names <- c(
  "Base" = "base",
  "AR1 M" = "ar1_M",
  "RV q ~ Smith Sound B" = "ss_q",
  "M ~ Capelin" = "cap_M"
)
models <- mget(model_names)
names(models) <- names(model_names)

saveRDS(models, file = "analysis/case_study/outputs/001_models.rds")

vis_tam(
  models,
  output_file = "analysis/case_study/outputs/001_models.html",
  open_file = FALSE
)

future::plan(future::multisession, workers = 10)
hind <- lapply(models, fit_hindcast, folds = 20)

lapply(hind, `[[`, "hindcast_rmse")

## Smith Sound improves cAIC but does little to improve forecast skill (makes sense).
## The capelin covariate improve forecast still, though cAIC is slightly poorer.


