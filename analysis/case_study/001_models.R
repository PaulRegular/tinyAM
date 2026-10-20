
pkgload::load_all(".", quiet = TRUE)

cod_obs <- tinyAM::cod_obs

dir.create("analysis/case_study/outputs", recursive = TRUE, showWarnings = FALSE)

## Include Smith Sound data to account for testing its effect on RV survey catchability
ss <- NCAM::inputs$SS.bms[, c("Year", "Biomass")]
names(ss) <- c("year", "smith_sound_biomass")
ss <- merge(data.frame(year = seq(min(ss$year), max(ss$year))), ss, by = "year", all = TRUE)
interp <- approx(ss$year, ss$smith_sound_biomass)
ind <- which(is.na(ss$smith_sound_biomass))
ss$smith_sound_biomass[ind] <- interp$y[ind] # fill two missing values using linear interpolation
ss$smith_sound_biomass <- ss$smith_sound_biomass / 1000
cod_obs$index <- merge(cod_obs$index, ss, by = "year", all = TRUE)
cod_obs$index$smith_sound_biomass[is.na(cod_obs$index$smith_sound_biomass)] <- 0 # assumed to be zero outside the Smith Sound survey years

## Include data on capelin for testing the effect of prey availability on survival
cap_per_cod <- NCAM::inputs$covariates[, c("year", "cap_per_cod")]
cod_obs$pop_covariates$age_group <- cut_ages(cod_obs$pop_covariates$age, c(0, 2, 5, 8, 20)) # c(0, 2, 4, 6, 8, 20))
cod_obs$pop_covariates <- merge(cod_obs$pop_covariates, cap_per_cod, by = "year", all.x = TRUE)
cod_obs$pop_covariates$log_cap_per_cod <- ifelse(
  is.na(cod_obs$pop_covariates$cap_per_cod),
  0,
  log(cod_obs$pop_covariates$cap_per_cod)
)


base <- fit_tam(
  cod_obs,
  years = 1983:2025,
  ages = 2:14,
  N_settings = list(
    process = "off"
  ),
  F_settings = list(
    process = "ar1",
    mu_form = ~1,
    mean_ages = 5:14
  ),
  M_settings = list(
    process = "ar1",
    mu_supplied = ~ I(0.3),
    age_breaks = c(2, 5, 8, 14), # c(2, 4, 6, 8, 14),
    mean_ages = 5:14
  ),
  index_settings = list(
    sd_form = ~ 1,
    q_form = ~ mono(q_block)
  )
)
base

ss_q <- update(
  base,
  index_settings = list(
    q_form = ~mono(q_block) + smith_sound_biomass
  )
)
ss_q

cap_M <- update(
  ss_q,
  M_settings = list(
    mu_form = ~ 0 + age_group:log_cap_per_cod
  )
)
cap_M

cAIC(base, ss_q, cap_M, details = TRUE)


model_names <- c(
  "Base" = "base",
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
hind <- lapply(models, fit_hindcast)

lapply(hind, `[[`, "hindcast_rmse")

## Smith Sound improves cAIC but does little to improve forecast skill (makes sense).
## The capelin covariate improve forecast still, though cAIC is slightly poorer.


