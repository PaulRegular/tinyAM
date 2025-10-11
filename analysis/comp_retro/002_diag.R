
source("analysis/comp_retro/001_models.R")

# fit <- sam_style
# fit <- ncam_style
# fit <- cohort_dev
fit <- sam_style

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


