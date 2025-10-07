
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

# fit <- sam_style
# fit <- ncam_style
# fit <- cohort_dev
fit <- m_dev

dat <- fit$dat
rep <- fit$rep
obs_pred <- fit$obs_pred
pop <- fit$pop

obs_pred$catch |>
  mutate(age = as.character(age)) |>
  plot_ly(x = ~year, color = ~age, legendgroup = ~age, colors = viridis::viridis(100)) |>
  add_markers(y = ~obs, showlegend = FALSE) |>
  add_lines(y = ~pred) |>
  layout(yaxis = list(type = "log"))

obs_pred$index |>
  mutate(age = as.character(age)) |>
  plot_ly(x = ~year, color = ~age, legendgroup = ~age, colors = viridis::viridis(100)) |>
  add_markers(y = ~obs, showlegend = FALSE) |>
  add_lines(y = ~pred) |>
  layout(yaxis = list(type = "log"))

obs_pred$catch |>
  mutate(age = as.character(age)) |>
  plot_ly(x = ~year, color = ~age, legendgroup = ~age, colors = viridis::viridis(100)) |>
  add_markers(y = ~std_res)

obs_pred$index |>
  mutate(age = as.character(age)) |>
  plot_ly(x = ~year, color = ~age, legendgroup = ~age, colors = viridis::viridis(100)) |>
  add_markers(y = ~std_res)

plot_ly(x = dat$years, y = dat$ages, z = t(rep$mu_F)) |>
  add_heatmap() |>
  layout(xaxis = list(title = "Year"),
         yaxis = list(title = "Age"),
         title = "log_mu_F")

plot_ly(x = dat$years, y = dat$ages, z = t(rep$F)) |>
  add_heatmap() |>
  layout(xaxis = list(title = "Year"),
         yaxis = list(title = "Age"),
         title = "F")

plot_ly(x = dat$years, y = dat$ages, z = t(rep$M)) |>
  add_heatmap() |>
  layout(xaxis = list(title = "Year"),
         yaxis = list(title = "Age"),
         title = "M")

pop$ssb |>
  plot_ly(x = ~year, showlegend = FALSE) |>
  add_ribbons(ymin = ~lwr, ymax = ~upr, color = I("steelblue"),
              line = list(width = 0)) |>
  add_lines(y = ~est, color = I("steelblue"))



## Retro ----

future::plan(future::multisession, workers = 6)

ncam_style_retros <- fit_retro(ncam_style, folds = 20, hindcast = TRUE)
sam_style_retros <- fit_retro(sam_style, folds = 20, hindcast = TRUE)
cohort_dev_retros <- fit_retro(cohort_dev, folds = 20, hindcast = TRUE)
m_dev_retros <- fit_retro(m_dev, folds = 20, hindcast = TRUE)

m_dev1 <- m_dev2 <- m_dev3 <- m_dev4 <- m_dev

# m_dev1$call$M_settings$age_breaks <- seq(2, 14, by = 1)
# m_dev1 <- update(m_dev1) # did not converge
m_dev3$call$M_settings$age_breaks <- seq(2, 14, by = 3)
m_dev3 <- update(m_dev3)
m_dev4$call$M_settings$age_breaks <- seq(2, 14, by = 4)
m_dev4 <- update(m_dev4)

m_dev3_retros <- fit_retro(m_dev3, folds = 20, hindcast = TRUE)
m_dev4_retros <- fit_retro(m_dev4, folds = 20, hindcast = TRUE)

m_dev_f_shift <- m_dev
m_dev_f_shift$call$F_settings$mu_form <- ~F_y_block + ~F_a_block
m_dev_f_shift <- update(m_dev_f_shift)
m_dev_f_shift_retros <- fit_retro(m_dev_f_shift, folds = 20, hindcast = TRUE)

m_dev_f_null <- m_dev
m_dev_f_null$call$F_settings$mu_form <- NULL
m_dev_f_null <- update(m_dev_f_null)
m_dev_f_null_retros <- fit_retro(m_dev_f_null, folds = 20, hindcast = TRUE, grad_tol = 0.01)

m_dev_f_ar1 <- m_dev
m_dev_f_ar1$call$F_settings$process <- "ar1"
m_dev_f_ar1 <- update(m_dev_f_ar1)
m_dev_f_ar1_retros <- fit_retro(m_dev_f_ar1, folds = 20, hindcast = TRUE, grad_tol = 0.01)

ncam_style_retros$hindcast_rmse
sam_style_retros$hindcast_rmse
m_dev_retros$hindcast_rmse
cohort_dev_retros$hindcast_rmse

ncam_style_retros$mohns_rho |> subset(is.na(age))
sam_style_retros$mohns_rho |> subset(is.na(age))
cohort_dev_retros$mohns_rho |> subset(is.na(age))
m_dev_retros$mohns_rho |> subset(is.na(age))

m_dev_retros$hindcast_rmse
m_dev3_retros$hindcast_rmse
m_dev4_retros$hindcast_rmse

m_dev_retros$mohns_rho |> subset(is.na(age))
m_dev3_retros$mohns_rho |> subset(is.na(age))
m_dev4_retros$mohns_rho |> subset(is.na(age))

m_dev_retros$hindcast_rmse
m_dev_f_shift_retros$hindcast_rmse
m_dev_f_null_retros$hindcast_rmse
m_dev_f_ar1_retros$hindcast_rmse

m_dev_retros$mohns_rho |> subset(is.na(age))
m_dev_f_shift_retros$mohns_rho |> subset(is.na(age))
m_dev_f_null_retros$mohns_rho |> subset(is.na(age))
m_dev_f_ar1_retros$mohns_rho |> subset(is.na(age))

## initial m_dev structure appears best so far


retros <- m_dev_retros
# retros <- ncam_style_retros

retros$pop$ssb |>
  mutate(fold = as.character(fold)) |>
  plot_ly(x = ~year, y = ~est, color = ~fold,
          colors = viridis::viridis(100)) |>
  add_lines()

retros$pop$abundance |>
  mutate(fold = as.character(fold)) |>
  plot_ly(x = ~year, y = ~est, color = ~fold,
          colors = viridis::viridis(100)) |>
  add_lines()

retros$pop$F_bar |>
  mutate(fold = as.character(fold)) |>
  plot_ly(x = ~year, y = ~est, color = ~fold,
          colors = viridis::viridis(100)) |>
  add_lines()

retros$pop$M_bar |>
  mutate(fold = as.character(fold)) |>
  plot_ly(x = ~year, y = ~est, color = ~fold,
          colors = viridis::viridis(100)) |>
  add_lines()

retros$pop$N |>
  mutate(fold = as.character(fold)) |>
  plot_ly(x = ~year, y = ~est, color = ~fold, frame = ~age,
          colors = viridis::viridis(100)) |>
  add_lines()


## M deviations coupled for every second age appears to have the lowest retro
## and best forecast skill


## Simulations ----

future::plan(multisession, workers = 6)
sims <- sim_tam(fit, n = 100, obs_only = TRUE)

sims$total_catch |>
  group_by(sim) |>
  plot_ly(x = ~year, y = ~est) |>
  add_lines(line = list(width = 0.5), alpha = 0.1) |>
  plotly::toWebGL() |>
  layout(yaxis = list(type = "log"))

sims$index |>
  group_by(sim, year) |>
  summarise(obs = sum(obs, na.rm = TRUE)) |>
  plot_ly(x = ~year, y = ~obs) |>
  add_lines(line = list(width = 0.5), alpha = 0.1) |>
  plotly::toWebGL() |>
  layout(yaxis = list(type = "log"))

sims$ssb |>
  group_by(sim) |>
  plot_ly(x = ~year, y = ~est) |>
  add_lines(line = list(width = 0.5), alpha = 0.1)


## TODO
## - Make a sim_test_tam function



