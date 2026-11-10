curve_term <- function(type) {
  .make_q_term(list(type = type, variable = "age", by = NULL),
    data.frame(age = 1:12, obs = 1, is_proj = FALSE))
}

test_that("double logistic matches the pinned WHAM raw product", {
  term <- curve_term("double_logistic")
  # Independent translation of WHAM get_selAA(), model 3, revision
  # ffcbfd093161e8210aa20716933975191b661807, src/selectivity.hpp (GPL-3).
  wham <- function(ages, a50up, kup, a50down, kdown)
    1/(1+exp(-(ages-a50up)/kup)) * 1/(1+exp((ages-a50down)/kdown))
  for (mid in list(c(3,8),c(8,3))) {
    p <- .curve_parameters(term)
    p[[term$a50_up_parameter]][] <- mid[1]
    p[[term$a50_down_parameter]][] <- mid[2]
    p[[term$scale_up_parameter]][] <- log(.8)
    p[[term$scale_down_parameter]][] <- log(1.2)
    expect_equal(exp(.curve_log_selectivity(1:12,p,term,1)),wham(1:12,mid[1],.8,mid[2],1.2),ignore_attr=TRUE)
    # WHAM counters are 1:na; biological ages 0:(na-1) shift both midpoints.
    p[[term$a50_up_parameter]][] <- mid[1]-1
    p[[term$a50_down_parameter]][] <- mid[2]-1
    expect_equal(exp(.curve_log_selectivity(0:11,p,term,1)),wham(1:12,mid[1],.8,mid[2],1.2),ignore_attr=TRUE)
    expect_true(max(wham(1:12,mid[1],.8,mid[2],1.2)) < 1)
  }
  p <- .curve_parameters(term)
  theta <- unlist(p,use.names=FALSE)
  f <- RTMB::MakeTape(function(x) {
    pars <- stats::setNames(lapply(seq_along(x),function(i)x[i]),names(p))
    .curve_log_selectivity(1:12,pars,term,1)
  },theta)
  up <- (1:12-theta[1])/exp(theta[2]); down <- (1:12-theta[3])/exp(theta[4])
  independent <- cbind(-(1-plogis(up))/exp(theta[2]),-up*(1-plogis(up)),
    plogis(down)/exp(theta[4]),down*plogis(down))
  expect_equal(unname(f$jacobian(theta)),independent,tolerance=1e-12)
})

test_that("normal selectivity matches a peak-normalized Gaussian", {
  term <- curve_term("normal")
  par <- .curve_parameters(term)
  peak <- 6.5
  width <- 2.25
  par[[term$peak_parameter]][] <- peak
  par[[term$scale_parameter]][] <- log(width)
  age <- c(1:12, peak)
  log_selectivity <- .curve_log_selectivity(age, par, term, 1)
  expect_equal(unname(exp(log_selectivity)),
    stats::dnorm(age, peak, width) / stats::dnorm(peak, peak, width))
  expect_equal(unname(tail(log_selectivity, 1)), 0)
  expect_true(is.finite(.curve_log_selectivity(1000, par, term, 1)))
  theta <- unlist(par, use.names = FALSE)
  tape <- RTMB::MakeTape(function(x) {
    parameters <- stats::setNames(list(x[1], x[2]), names(par))
    .curve_log_selectivity(age, parameters, term, 1)
  }, theta)
  distance <- age - peak
  expected <- cbind(distance / width^2, distance^2 / width^2)
  expect_equal(unname(tape$jacobian(theta)), expected, tolerance = 1e-12)
  expect_error(normal(age), "formula marker")
})

test_that("selectivity curves compile independently for F and q and reject invalid roles", {
  for (form in list(~logistic(age), ~normal(age), ~double_logistic(age))) {
    d <- make_test_dat(years=2000:2005,ages=2:8,
      F_settings=list(process="off",mu_form=form),index_settings=list(q_form=form))
    p <- make_par(d)
    expect_equal(anyDuplicated(names(p)),0L)
    expect_null(d$F_terms[[1]]$parameter)
    expect_null(d$q_terms[[1]]$parameter)
    o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
    expect_true(is.finite(o$fn()))
    expect_true(all(is.finite(o$gr())))
    expect_true(all(is.finite(o$he())))
  }
  expect_error(make_test_dat(F_settings=list(process="rw",mu_form=~double_logistic(age))),"cancel")
  expect_error(make_test_dat(M_settings=list(mu_form=~double_logistic(age))),"F only")
  expect_error(make_test_dat(index_settings=list(sd_form=~double_logistic(age))),"not catchability")
  expect_error(.make_q_term(list(type = "normal", variable = "age", by = NULL),
    data.frame(age = 1:2, obs = 1)), "needs at least 3 observed")
})

test_that("tiny selectivities retain finite likelihoods and residuals without floors", {
  d <- make_test_dat(years=2000:2005,ages=2:8,
    F_settings=list(process="off",mu_form=~double_logistic(age)),
    index_settings=list(q_form=~double_logistic(age),q_link="logit"))
  p <- make_par(d)
  for (term in c(d$F_terms,d$q_terms)) {
    p[[term$a50_up_parameter]][] <- 4.3
    p[[term$a50_down_parameter]][] <- 6.3
    p[[term$scale_up_parameter]][] <- p[[term$scale_down_parameter]][] <- log(.001)
  }
  o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
  expect_true(is.finite(o$fn()))
  expect_true(all(is.finite(o$gr())))
  expect_true(all(is.finite(o$he())))
  r <- o$report()
  expect_true(any(r$F == 0)) # numerical underflow is not a structural zero
  expect_true(all(is.finite(r$log_pred[d$is_observed])))
  expect_true(all(is.finite(r$S)))
  fit <- structure(list(dat=d,rep=r,obj=o),class="tam_fit")
  tab <- tidy_obs_pred(fit)
  expect_true(all(is.finite(tab$catch$std_res[which(d$obs$catch$obs > 0)])))
  expect_true(all(is.finite(tab$index$std_res[which(d$obs$index$obs > 0)])))
  log_res <- residuals(fit,type="log")
  expect_true(all(is.finite(log_res$residual[which(log_res$obs > 0)])))
})

test_that("scoped fleet curves retain states, structural zeros and projection multipliers", {
  obs <- cod_obs
  obs$catch <- rbind(transform(obs$catch,fleet="Trawl",obs=.6*obs),
    transform(obs$catch,fleet="Gillnet",obs=.4*obs))
  obs$catch$obs[obs$catch$fleet == "Gillnet" & obs$catch$age < 4] <- 0
  obs$index <- rbind(transform(obs$index,survey="RV"),transform(obs$index,survey="Sentinel",obs=.8*obs))
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:10,
    F_settings=list(process="ar1",mu_form=~0+for_fleet("Trawl",~logistic(age))+
      for_fleet("Gillnet",~double_logistic(age)),
      age_breaks=list(Trawl=NULL,Gillnet=c(4:10,11))),
    index_settings=list(q_form=~0+for_survey("RV",~double_logistic(age))+
      for_survey("Sentinel",~double_logistic(age))),
    proj_settings=list(n_proj=2,F_mult=c(1,.5)))
  p <- make_par(d)
  expect_equal(lengths(.formula_random_parameters(d)),integer())
  expect_true(all(p$log_f == .f_state_mean(matrix(drop(d$F_modmat %*% p$log_mu_f)+
    .mean_effects(p,d,"F")$contribution,length(d$years)),d)[!d$is_proj,]))
  o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
  r <- o$report()
  expect_true(all(is.finite(o$gr())))
  expect_true(all(r$F_fleet[,as.character(2:3),"Gillnet"] == 0))
  expect_equal(r$F[7,],r$F[6,],ignore_attr=TRUE)
  expect_equal(r$F_fleet[8,,],r$F_fleet[6,,]*.5,ignore_attr=TRUE)
  expect_equal(r$F_full_fleet,apply(r$F_fleet,c(1,3),max),ignore_attr=TRUE)
  expect_true(all(vapply(d$F_terms,function(term) all(startsWith(unlist(term[c("a50_parameter","a50_up_parameter")]),"F_")),logical(1))))
  set.seed(902)
  draw <- nll_fun(p,d,simulate=TRUE)
  expect_true(all(is.finite(draw$log_f)))
  expect_true(all(is.finite(draw$log_obs[d$is_observed])))
})
