library(tinyAM)

## Recruitment assumptions ----
base <- fit_tam(cod_obs, years = 1983:2024, ages = 2:14,
  N_settings = list(process = "off"), rec_settings = list(process = "rw", mu_form = NULL),
  F_settings = list(process = "rw", mu_form = NULL),
  M_settings = list(process = "off", mu_supplied = ~ I(.3)), silent = TRUE)
start <- as.list(base$sdrep, "Estimate")
start$rec_beta <- c("(Intercept)" = mean(log(base$rep$recruitment)))
independent <- update(base, N_settings = list(process = "off"), rec_settings = list(process = "iid", mu_form = ~ 1),
                      start_par = start)
beverton_holt <- update(base, N_settings = list(process = "off"), rec_settings = list(process = "iid", mu_form = ~ bh(ssb)),
                       start_par = start)
ricker_ar1 <- update(base, N_settings = list(process = "off"), rec_settings = list(process = "ar1", mu_form = ~ ricker(ssb)),
                     start_par = start)

## An illustrative covariate; replace with measured annual values ----
data <- cod_obs
data$pop_covariates$temperature <- NA_real_
i <- data$maturity$age == 2
data$pop_covariates$temperature[i] <- sin((data$pop_covariates$year[i] - 1983) / 5)
start$rec_beta <- c(start$rec_beta, temperature = 0)
covariate <- update(base, data = data,
  N_settings = list(process = "off"), rec_settings = list(process = "iid", mu_form = ~ temperature), start_par = start)

## Inspect support, not just fitted curves ----
fits <- list(RW = base, IID = independent, BH = beverton_holt,
             Ricker = ricker_ar1, Covariate = covariate)
lapply(fits, check_tam)
tidy_recruitment(beverton_holt)$pairs
if (interactive()) vis_tam(model_list = fits)
