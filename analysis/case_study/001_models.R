
pkgload::load_all(".", quiet = TRUE)

dir.create("analysis/case_study/outputs", recursive = TRUE, showWarnings = FALSE)

scaa <- fit_tam(
  cod_obs,
  years = 1983:2025,
  ages = 2:14,
  N_settings = list(
    process = "off",
    init = "exp"
  ),
  F_settings = list(
    process = "off",
    mu_form = ~ factor(year) + factor(age),
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
scaa

sam <- update(
  scaa,
  N_settings = list(
    process = "iid"
  ),
  F_settings = list(
    process = "cor_rw",
    mu_form = NULL,
    mean_ages = 5:14
  )
)
sam

ncam <- update(
  scaa,
  F_settings = list(
    process = "ar1",
    mu_form = ~ 1,
    mean_ages = 5:14
  ),
  M_settings = list(
    process = "ar1",
    age_breaks = c(2, 5, 8, 14)
  )
)
ncam


future::plan(future::multisession, workers = 10)
hind <- fit_hindcast(
  scaa,
  sam,
  ncam,
  folds = 15,
  proj_settings = list(
    n_proj = 1,
    n_mean = 1,
    catch = "observed"
  ),
  start_from_fit = TRUE
)
hind

vis_tam(scaa, sam, ncam)

