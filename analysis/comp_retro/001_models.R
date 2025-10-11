
## Exploratory analysis of retro patterns of a model with cohort deviations or
## M deviations

library(RTMB)
library(tinyAM)
library(plotly)

## Terminal fit ----

sam_style <- fit_tam(
  tinyAM::cod_obs,
  years = 1983:2024,
  ages = 2:14,
  N_settings = list(process = "iid", init_N0 = FALSE),
  F_settings = list(process = "approx_rw", mu_form = NULL, mean_ages = 5:14),
  M_settings = list(process = "off", assumption = ~M_assumption, mean_ages = 5:14),
  obs_settings = list(q_form = ~ q_block, sd_form = ~ sd_obs_block)
)

cohort_dev <- update(sam_style,
                     N_settings = list(process = "iid", init_N0 = FALSE),
                     F_settings = list(process = "approx_rw", mu_form = ~F_a_block, mean_ages = 5:14),
                     M_settings = list(process = "off", assumption = ~I(0.3), mean_ages = 5:14))

ncam_style <- update(sam_style,
                     N_settings = list(process = "off", init_N0 = TRUE),
                     F_settings = list(process = "approx_rw", mu_form = NULL, mean_ages = 5:14),
                     M_settings = list(process = "ar1",
                                       mu_form = NULL,
                                       assumption = ~M_assumption,
                                       age_breaks = c(2:9, 14),
                                       mean_ages = 5:14))

m_dev <- update(ncam_style,
                F_settings = list(process = "approx_rw", mu_form = ~F_a_block),
                M_settings = list(process = "ar1",
                                  mu_form = NULL,
                                  assumption = ~I(0.3),
                                  age_breaks = seq(2, 14, by = 2),
                                  mean_ages = 5:14))

sam_style$opt$objective # 887.8872
cohort_dev$opt$objective # 987.7389
ncam_style$opt$objective # 811.0894
m_dev$opt$objective # 819.4535

