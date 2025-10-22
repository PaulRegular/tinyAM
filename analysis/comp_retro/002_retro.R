
library(tinyAM)

# source("analysis/comp_retro/001_models.R")
models <- readRDS("analysis/comp_retro/001_models.rds")

future::plan(future::multisession, workers = 10)

retros <- lapply(names(models), function(nm) {
  try(fit_retro(models[[nm]], folds = 20, hindcast = TRUE, grad_tol = 1e-2))
})

saveRDS(retros, file = "analysis/comp_retro/002_retros.rds")











retros <- vector("list", length(models))
names(retros) <-
for (nm in names(models)) {
  retros[[nm]] <- fit_retro(models[[nm]])
}


ncam_style_retros <- fit_retro(ncam_style, folds = 20, hindcast = TRUE)
vis_tam(ncam_style_retros$fits, output_file = "analysis/comp_retro/003_ncam_style_retros.html")

sam_style_retros <- fit_retro(sam_style, folds = 20, hindcast = TRUE)
vis_tam(sam_style_retros$fits, output_file = "analysis/comp_retro/003_sam_style_retros.html")

cohort_dev_retros <- fit_retro(cohort_dev, folds = 20, hindcast = TRUE)
m_dev_retros <- fit_retro(m_dev, folds = 20, hindcast = TRUE)

ncam_style1 <- ncam_style
ncam_style1$call$M_settings$age_breaks <- seq(2, 14, by = 1)
ncam_style1 <- update(ncam_style1)
ncam_style1_retros <- fit_retro(ncam_style1, folds = 20, hindcast = TRUE)


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
# retros <- sam_style_retros
# retros <- ncam_style1_retros

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

## Also, there are some signs of positive retro with SAM style and negative retro
## with NCAM style



