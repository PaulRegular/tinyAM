
source("analysis/comp_retro/001_models.R")

fit <- ncam_style

fit <- update(fit, proj_settings = list(n_proj = 10, n_mean = 5, F_mult = 1))

future::plan(multisession, workers = 2)
sims <- sim_tam(fit, n = 100, par_uncertainty = "joint")

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

sims$abundance |>
  group_by(sim) |>
  plot_ly(x = ~year, y = ~est) |>
  add_lines(line = list(width = 0.5), alpha = 0.1)

sims$ssb |>
  group_by(sim) |>
  plot_ly(x = ~year, y = ~est) |>
  add_lines(line = list(width = 0.5), alpha = 0.5)


## TODO
## - Make a sim_test_tam function
## - Check effect of lower and higher obs errors
## - Check effect of projecting declining pop vs increasing pop




