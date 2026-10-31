fleet_test_data <- function() {
  d <- cod_obs
  d$catch <- rbind(transform(d$catch, fleet = "Trawl", obs = obs * .6),
    transform(d$catch, fleet = "Gillnet", obs = obs * .4))
  d
}

test_that("one labeled fleet preserves the single-fishery likelihood and derivatives", {
  obs <- cod_obs; obs$catch$fleet <- "Fishery"
  for (process in c("off", "iid", "rw", "cor_rw", "ar1")) {
    settings <- list(process=process,mu_form=if(process %in% c("rw","cor_rw")) NULL else ~factor(age))
    a <- make_test_dat(years=2000:2005,ages=2:6,F_settings=settings)
    b <- make_test_dat(data=obs,years=2000:2005,ages=2:6,F_settings=settings)
    p <- make_par(a); q <- make_par(b)
    if (process != "off") {
      set.seed(321); p$log_f[] <- rnorm(length(p$log_f),-.5,.2); q$log_f[] <- p$log_f
    }
    x <- RTMB::MakeADFun(function(p)nll_fun(p,a),p,silent=TRUE)
    y <- RTMB::MakeADFun(function(p)nll_fun(p,b),q,silent=TRUE)
    expect_equal(x$fn(),y$fn(),tolerance=1e-11)
    expect_equal(x$gr(),y$gr(),tolerance=1e-11,ignore_attr=TRUE)
    expect_equal(x$report()$log_pred,y$report()$log_pred,tolerance=1e-11)
  }
})

test_that("fleet inputs have complete explanatory grids and explicit routing", {
  obs <- fleet_test_data()
  expect_true(check_tam_data(obs))
  missing <- obs; missing$catch <- missing$catch[-1L, ]
  expect_error(check_tam_data(missing), "coverage is incomplete")
  duplicate <- obs; duplicate$catch <- rbind(duplicate$catch, duplicate$catch[1L, ])
  expect_error(check_tam_data(duplicate), "duplicate")
  obs$catch$gear <- factor(ifelse(obs$catch$age < 5, "small", "large"))
  obs$catch$obs[obs$catch$year == 2002] <- NA
  d <- make_test_dat(data = obs, years = 2000:2005, ages = 2:6,
    F_settings = list(process = "iid", mu_form = ~ fleet + gear))
  expect_equal(nrow(d$F_modmat), 60L)
  expect_error(make_test_dat(data = obs, F_settings = list(process = "iid", mu_form = ~obs)), "response obs")
  expect_error(make_test_dat(F_settings = list(process = "iid", mu_form = ~M_a_block)), "data\\$catch")
})

test_that("all fleet processes preserve Baranov catches and independent densities", {
  for (process in c("off", "iid", "rw", "cor_rw", "ar1")) {
    d <- make_test_dat(data = fleet_test_data(), years = 2000:2005, ages = 2:6,
      F_settings = list(process = process, mu_form = if (process %in% c("rw", "cor_rw")) NULL else ~fleet))
    p <- make_par(d)
    if (process != "off") p$log_f[] <- matrix(seq_along(p$log_f) / 30, nrow(p$log_f))
    obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
    expect_true(is.finite(obj$fn()))
    expect_true(all(is.finite(obj$gr())))
    r <- obj$report()
    expect_equal(r$F, apply(r$F_fleet, c(1, 2), sum))
    expected <- r$N * r$F_fleet[, , 1L] / r$Z * -expm1(-r$Z)
    rows <- d$obs_map$type == "catch" & d$obs_map$fleet == d$fleets[1L]
    expect_equal(unname(exp(r$log_pred[rows])), as.vector(expected))
    if (process != "off") {
      mean <- .F_matrix(0, d)
      expected <- sum(vapply(d$fleets, function(f) {
        columns <- d$F_states$fleet == f
        x <- p$log_f[, columns, drop = FALSE]
        switch(process, iid = sum(dnorm(x, log = TRUE)), rw = dprocess_rw(x),
          cor_rw = dprocess_cor_rw(x), ar1 = dprocess_ar1(x, phi = c(.5, .5)))
      }, numeric(1)))
      expect_equal(as.numeric(.fleet_F_process(p$log_f, mean, p, d)), expected)
      expect_equal(dim(nll_fun(p, d, simulate = TRUE)$log_f), dim(p$log_f))
    } else expect_null(p$log_f)
  }
})

test_that("fleet blocks and correlation groups describe only applicable states", {
  obs <- fleet_test_data()
  obs$catch$obs[obs$catch$fleet == "Trawl" & obs$catch$age < 3] <- 0
  d <- make_test_dat(data = obs, years = 2000:2005, ages = 2:6,
    F_settings = list(process = "ar1", mu_form = ~fleet,
      age_breaks = list(Trawl = c(3, 7), Gillnet = NULL),
      cor_form = list(age = ~fleet, year = ~1)))
  p <- make_par(d)
  expect_equal(ncol(p$log_f), 6L)
  expect_identical(names(p$logit_phi_f), c("age / Gillnet", "year"))
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  expect_true(all(is.finite(obj$gr())))
  expect_equal(obj$report()$F_fleet[, 1L, "Trawl"], rep(0, 6), ignore_attr = TRUE)
  expect_error(make_test_dat(data = obs, F_settings = list(age_breaks = list(Trawl = NULL))), "every catch fleet")
  expect_error(make_test_dat(data = obs, F_settings = list(cor_form = ~age)), "constant within each fleet")
  expect_error(make_test_dat(data = obs, F_settings = list(cor_form = ~F_y_block)), "vary within fleet")
  expect_error(make_test_dat(data = obs, ages = 2:6,
    F_settings = list(process = "off", mu_form = ~fleet + logistic(age, by = fleet), age_breaks = c(2,4,7))),
    "F age block.*mu_form")
})

test_that("typed fleet weights calculate distinct yields and projection targets", {
  obs <- fleet_test_data()
  obs$weight <- rbind(transform(obs$weight, type = "stock", fleet = NA_character_),
    transform(obs$weight, type = "catch", fleet = "Trawl", obs = 2 * obs),
    transform(obs$weight, type = "catch", fleet = "Gillnet", obs = 3 * obs))
  expect_true(check_tam_data(obs))
  d <- make_test_dat(data = obs, years = 2000:2005, ages = 2:6,
    proj_settings = list(n_proj = 2, n_mean = 2, catch = 1))
  p <- make_par(d)
  obj <- RTMB::MakeADFun(function(p) nll_fun(p, d), p, silent = TRUE)
  expect_true(all(is.finite(obj$gr())))
  r <- obj$report()
  expect_equal(r$total_yield_pred[d$is_proj], rep(1, 2), ignore_attr = TRUE)
  expect_equal(r$W_catch_fleet[, , "Trawl"], 2 * r$W, ignore_attr = TRUE)
  expect_equal(r$W_catch_fleet[, , "Gillnet"], 3 * r$W, ignore_attr = TRUE)
  expect_equal(r$total_yield_pred, rowSums(r$yield_pred_by_fleet), ignore_attr = TRUE)
  bad <- obs; bad$weight$fleet[bad$weight$type == "catch"][1L] <- NA
  expect_error(check_tam_data(bad), "mix shared")
})

test_that("fleet-specific correlation groups agree with independent dense covariances", {
  obs <- fleet_test_data()
  obs$catch$sd_group <- factor(ifelse(obs$catch$age < 4, "young", "old"))
  for (process in c("ar1", "cor_rw")) {
    d <- make_test_dat(data = obs, years = 2000:2003, ages = 2:6,
      F_settings = list(process = process, mu_form = if (process == "ar1") ~fleet else NULL,
        age_breaks = c(2, 4, 7), cor_form = ~fleet, sd_form = ~fleet + sd_group))
    p <- make_par(d)
    p$sd_beta_f[] <- c(log(.3), log(2), log(.6))
    set.seed(45); p$log_f[] <- rnorm(length(p$log_f), sd = .2)
    if (process == "ar1") p$logit_phi_f[] <- qlogis(c(.2, .4, .6, .7)) else p$atanh_rho_f[] <- atanh(c(-.3,.5))
    sd <- .process_sd(p, d, "F")
    expected <- 0
    for (f in seq_along(d$fleets)) {
      columns <- which(d$F_states$fleet == d$fleets[f])
      x <- p$log_f[, columns, drop = FALSE]
      if (process == "cor_rw") {
        x <- x[-1L, ] - x[-nrow(x), ]
        phi <- c(tanh(p$atanh_rho_f[f]), 0); normalization <- 1
      } else {
        phi <- .fleet_phi(p, d, f); normalization <- prod(1 - phi^2)
      }
      covariance <- kronecker(outer(sd[columns], sd[columns]) * phi[1]^abs(outer(1:2,1:2,"-")),
        phi[2]^abs(outer(seq_len(nrow(x)),seq_len(nrow(x)),"-"))) / normalization
      expected <- expected - .5 * (length(x)*log(2*pi) + as.numeric(determinant(covariance)$modulus) + sum(x*solve(covariance, c(x))))
    }
    expect_equal(as.numeric(.fleet_F_process(p$log_f, .F_matrix(0,d),p,d)), expected, tolerance=1e-10)
    objective <- RTMB::MakeADFun(function(p) nll_fun(p,d),p,silent=TRUE)
    expect_true(all(is.finite(objective$gr())))
  }
})

test_that("plus-group F covariates are explicit and partial catch totals retain known fleets", {
  obs <- fleet_test_data()
  obs$catch$predictor <- obs$catch$age
  expect_error(make_test_dat(data=obs,ages=2:6,F_settings=list(process="iid",mu_form=~predictor)),
    "predictor conflicts across biological ages")
  obs$catch$obs[obs$catch$fleet == "Trawl" & obs$catch$year == 2002] <- NA
  d <- make_test_dat(data=obs,years=2000:2003,ages=2:6)
  objective <- RTMB::MakeADFun(function(p) nll_fun(p,d),make_par(d),silent=TRUE)
  r <- objective$report()
  expected <- sum(d$obs$catch$obs[d$obs$catch$year == 2002],na.rm=TRUE)
  expect_equal(unname(r$total_catch[which(d$years == 2002)]),expected)
})

test_that("fleet process plots retain distinct state coordinates", {
  f <- expand.grid(year=2000:2003,age_block=c("2-3","4-6"),fleet=c("A","B"),stringsAsFactors=FALSE)
  f$est <- sin(seq_len(nrow(f))); f$model <- "fit"
  out <- structure(list(deviations=list(F=f), recruitment=list(),
    process_residuals=list(F=f,status=data.frame(model="fit",status="ok"))),label="model")
  before <- serialize(out,NULL)
  d <- .process_nfm(out$deviations,out,"fit")
  d <- d[d$component == "F", ]
  expect_equal(length(unique(d$coordinate)),4L)
  expect_false(anyDuplicated(d[c("year","coordinate")]) > 0)
  expect_s3_class(plot_process_deviations(out),"plotly")
  expect_s3_class(plot_process_residuals(out),"plotly")
  expect_identical(serialize(out,NULL),before)
})

test_that("fleet hindcast targets use actual removals and training weights stay separate", {
  g <- expand.grid(year=2000:2002,age=2:4)
  c <- rbind(transform(g,fleet="A",obs=age,unit="fish"),
    transform(g,fleet="B",obs=2*age,unit="fish"))
  w <- rbind(transform(g,type="stock",fleet=NA_character_,obs=1,unit="kg"),
    transform(g,type="catch",fleet="A",obs=year-1998,unit="kg"),
    transform(g,type="catch",fleet="B",obs=year-1995,unit="kg"))
  obs <- list(catch=c,weight=w,maturity=transform(g,obs=.5),
    index=transform(g,obs=100,survey="RV",samp_time=.5))
  d <- prepare_tam(obs,years=2000:2001,ages=2:3,
    index_settings=list(q_form=~1),
    proj_settings=list(n_proj=1,n_mean=1,F_mult=1))
  f <- list(dat=d,refit_args=list(data=obs))
  target <- .hindcast_controls(f,2002,.retro_projection(NULL))
  expect_equal(target$target,162)
  expect_equal(d$W_catch_fleet[d$is_proj,d$F_columns$fleet == "A"],rep(3,2),ignore_attr=TRUE)
  expect_equal(d$W_catch_fleet[d$is_proj,d$F_columns$fleet == "B"],rep(6,2),ignore_attr=TRUE)
  expect_equal(d$W_catch_plus_fleet$A[d$is_proj,],rep(3,2),ignore_attr=TRUE)
  expect_equal(d$W_catch_plus_fleet$B[d$is_proj,],rep(6,2),ignore_attr=TRUE)
  f$refit_args$data$catch$obs[which(c$year == 2002 & c$fleet == "B")[1L]] <- NA
  expect_match(.hindcast_controls(f,2002,.retro_projection(NULL))$error,"incomplete")
})
