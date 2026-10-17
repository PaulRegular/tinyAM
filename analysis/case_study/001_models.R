
pkgload::load_all(".", quiet = TRUE)

cod_obs <- tinyAM::cod_obs

dir.create("analysis/case_study/outputs", recursive = TRUE, showWarnings = FALSE)

## Accounting for Smith Sound effect via covariate
ss <- NCAM::inputs$SS.bms[, c("Year", "Biomass")]
names(ss) <- c("year", "smith_sound_biomass")
ss <- merge(data.frame(year = seq(min(ss$year), max(ss$year))), ss, by = "year", all = TRUE)
interp <- approx(ss$year, ss$smith_sound_biomass)
ind <- which(is.na(ss$smith_sound_biomass))
ss$smith_sound_biomass[ind] <- interp$y[ind] # fill two missing values using linear interpolation
ss$smith_sound_biomass <- ss$smith_sound_biomass / 1000
cod_obs$index <- merge(cod_obs$index, ss, by = "year", all = TRUE)
cod_obs$index$smith_sound_biomass[is.na(cod_obs$index$smith_sound_biomass)] <- 0 # assumed to be zero outside the Smith Sound survey years

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
    age_breaks = c(2, 4, 6, 8, 14),
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

model_names <- c(
  "Base" = "base",
  "Smith Sound q" = "ss_q"
)
models <- mget(model_names)
names(models) <- names(model_names)

saveRDS(models, file = "analysis/case_study/outputs/001_models.rds")

vis_tam(
  models,
  output_file = "analysis/case_study/outputs/001_models.html",
  open_file = FALSE
)
