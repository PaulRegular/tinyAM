
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
  F_settings = list(process = "approx_rw", mu_form = NULL),
  M_settings = list(process = "off", assumption = ~M_assumption),
  obs_settings = list(q_form = ~ q_block, sd_form = ~ sd_obs_block)
)

ncam_style <- update(sam_style,
                     N_settings = list(process = "off", init_N0 = TRUE),
                     F_settings = list(process = "approx_rw", mu_form = NULL),
                     M_settings = list(process = "ar1",
                                       mu_form = NULL,
                                       assumption = ~M_assumption,
                                       age_breaks = c(2:9, 14)))

sam_style$opt$objective # 887.8872
ncam_style$opt$objective # 811.0894

# fit <- sam_style
fit <- ncam_style

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
retros <- fit_retro(fit, folds = 20)
retros$mohns_rho |> subset(is.na(age))

plot_ly(retros$pop$ssb, x = ~year, y = ~est, color = ~fold,
        colors = viridis::viridis(100)) |>
  add_lines()

plot_ly(retros$pop$abundance, x = ~year, y = ~est, color = ~fold,
        colors = viridis::viridis(100)) |>
  add_lines()

plot_ly(retros$pop$N, x = ~year, y = ~est, color = ~fold, frame = ~age,
        colors = viridis::viridis(100)) |>
  add_lines()


## Simulations ----

plot(fit$dat$years, fit$rep$ssb, xlab = "Year", ylab = "SSB", type = "l")
for (i in 1:10) {
  sims <- sim_tam(fit)
  lines(fit$dat$years, sims$ssb, lwd = 0.5, col = "grey")
}

