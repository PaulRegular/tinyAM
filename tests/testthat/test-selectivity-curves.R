curve_term <- function(type, b = NULL) {
  .make_q_term(list(type = type, variable = "age", by = NULL, b = b),
    data.frame(age = 1:12, obs = 1, is_proj = FALSE))
}

test_that("warped normal has the specified half-selection, slope, peak and tail", {
  for (b in c(0, .3, 2, Inf)) {
    term <- curve_term("warped_normal", b)
    p <- .curve_parameters(term)
    p[[term$a50_parameter]][] <- 3
    p[[term$scale_parameter]][] <- log(2)
    f <- function(x) exp(.curve_log_selectivity(x, p, term, 1))
    expect_equal(f(3), .5, ignore_attr = TRUE)
    h <- 1e-5
    expect_equal(as.numeric((f(3+h)-f(3-h))/(2*h)), log(2)/2, tolerance = 1e-8)
    if (b > 0) {
      peak <- if (is.infinite(b)) 5 else 3 + 2*(1+b)*log1p(1/b)
      expect_equal(f(peak), 1, ignore_attr = TRUE, tolerance = 1e-12)
    }
    if (b == 0) expect_equal(f(1e4), 1, ignore_attr = TRUE)
    if (is.finite(b)) expect_equal(f(1e4), 2^(-b^2), ignore_attr = TRUE)
    if (is.infinite(b)) expect_equal(f(1:12), exp(-log(2)*((1:12-5)/2)^2), ignore_attr = TRUE)
  }
})

test_that("warped endpoints have finite exact first and second AD derivatives", {
  tape <- RTMB::MakeTape(function(p) .warped_g(p[1],p[2]), c(2,0))
  expect_equal(as.numeric(tape(c(2,0))), -1)
  expect_equal(unname(tape$jacobian(c(2,0))), matrix(c(-1,2),1))
  expect_equal(unname(tape$jacfun()$jacobian(c(2,0))), matrix(c(0,2,2,-8/3),2))
  for (shape in c(0, 1e-10, .01, .5, 1)) {
    expect_true(all(is.finite(tape$jacobian(c(2,shape)))))
    expect_true(all(is.finite(tape$jacfun()$jacobian(c(2,shape)))))
  }
  # One tape crosses both branches; it is not frozen at its starting endpoint.
  for (shape in c(.001,.2,1)) {
    expected <- 1+expm1(-2*shape)/shape
    expect_equal(as.numeric(tape(c(2,shape))),expected,tolerance=1e-12)
  }
  for (u in c(-3,0,2)) for (shape in c(0,.2,.8,1)) {
    e <- exp(-shape*u)
    first <- c(-e,if (shape == 0) u^2/2 else (1-e*(1+shape*u))/shape^2)
    second <- matrix(c(shape*e,u*e,u*e,
      if (shape == 0) -u^3/3 else (e*(2+2*shape*u+shape^2*u^2)-2)/shape^3),2)
    expect_equal(unname(tape$jacobian(c(u,shape))),matrix(first,1),tolerance=1e-10)
    expect_equal(unname(tape$jacfun()$jacobian(c(u,shape))),second,tolerance=1e-9)
  }
})

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

test_that("new curves compile independently for F and q and reject invalid roles", {
  for (form in list(~warped_normal(age),~warped_normal(age,b=0),
    ~warped_normal(age,b=Inf),~double_logistic(age))) {
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
  for (b in list(-1,NA_real_,c(0,1),"dome")) {
    form <- substitute(~warped_normal(age,b=VALUE),list(VALUE=b))
    form <- stats::as.formula(form)
    expect_error(make_test_dat(index_settings=list(q_form=form)),"nonnegative")
  }
  expect_error(make_test_dat(F_settings=list(process="rw",mu_form=~warped_normal(age))),"cancel")
  expect_error(make_test_dat(M_settings=list(mu_form=~warped_normal(age))),"F only")
  expect_error(make_test_dat(index_settings=list(sd_form=~double_logistic(age))),"not catchability")
})

test_that("tiny selectivities retain finite likelihoods and residuals without floors", {
  d <- make_test_dat(years=2000:2005,ages=2:8,
    F_settings=list(process="off",mu_form=~warped_normal(age,b=Inf)),
    index_settings=list(q_form=~warped_normal(age,b=Inf),q_link="logit"))
  p <- make_par(d)
  for (term in c(d$F_terms,d$q_terms)) {
    p[[term$a50_parameter]][] <- 4.3
    p[[term$scale_parameter]][] <- log(.01)
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
      for_fleet("Gillnet",~warped_normal(age)),
      age_breaks=list(Trawl=NULL,Gillnet=c(4:10,11))),
    index_settings=list(q_form=~0+for_survey("RV",~double_logistic(age))+
      for_survey("Sentinel",~warped_normal(age,b=Inf))),
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
  expect_true(all(vapply(d$F_terms,function(term) startsWith(term$a50_parameter,"F_"),logical(1))))
  set.seed(902)
  draw <- nll_fun(p,d,simulate=TRUE)
  expect_true(all(is.finite(draw$log_f)))
  expect_true(all(is.finite(draw$log_obs[d$is_observed])))
})

test_that("Gaussian uncertainty draws respect the warped shape bounds", {
  d <- make_test_dat(years=2000:2005,ages=2:8,
    F_settings=list(process="off",mu_form=~warped_normal(age)))
  p <- make_par(d); name <- d$F_terms[[1]]$shape_parameter
  o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
  fit <- structure(list(dat=d,obj=o,rep=o$report()),class="tam_fit")
  for (shape in c(-.2,1.2)) {
    par <- p; par[[name]][] <- shape
    observed <- .sim_obs(fit,par_fun=function(fit)par,redraw_random=FALSE)
    par[[name]][] <- pmin(1,pmax(0,shape))
    expected <- RTMB::MakeADFun(function(p)nll_fun(p,d),par,silent=TRUE)$report()
    expect_equal(observed$F$est,as.vector(expected$F))
  }
})

test_that("warped derived quantities use joint covariance and expose infinite endpoints", {
  for (b in list(NULL,0,Inf)) {
    d <- make_test_dat(years=2000:2005,ages=2:8,index_settings=list(q_form=~warped_normal(age,b=b)))
    p <- make_par(d); term <- d$q_terms[[1L]]
    p[[term$a50_parameter]][] <- 3
    p[[term$scale_parameter]][] <- log(2)
    if (is.null(b)) p[[term$shape_parameter]][] <- .5
    o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
    covariance <- diag(.01,length(o$par))
    sdr <- list(pdHess=TRUE,par.fixed=o$par,cov.fixed=covariance)
    fit <- structure(list(dat=d,obj=o,rep=o$report(),parameter_values=o$par,sdrep=sdr),class="tam_fit")
    # The helper reads stored modes when given this minimal uncertainty fixture.
    fit$sdrep <- NULL
    tab <- .tidy_selectivity_curves(fit,.95)
    if (is.null(b)) {
      expect_equal(tab$est[tab$par == "q_curve_b"],1)
      expect_equal(tab$est[tab$par == "q_curve_peak"],3+4*log(2))
      expect_equal(tab$est[tab$par == "q_curve_tail"],.5)
    } else {
      infinite <- if (b == 0) "q_curve_peak" else "q_curve_b"
      expect_equal(tab$est[tab$par == infinite],Inf)
      expect_true(is.na(tab$se[tab$par == infinite]))
    }
    if (is.null(b)) {
      # Independently propagate a non-diagonal covariance through an analytic
      # Jacobian; the production helper uses RTMB derivatives.
      active <- which(names(o$par) %in% c(term$a50_parameter,term$scale_parameter,term$shape_parameter))
      covariance[active[1],active[2]] <- covariance[active[2],active[1]] <- .002
      actual <- RTMB::sdreport(o,par.fixed=o$par,hessian.fixed=solve(covariance),getReportCovariance=FALSE)
      fit$sdrep <- actual
      result <- .tidy_selectivity_curves(fit,.95)
      jac <- c(1,4*log(2),8*(1-log(2)))
      expected <- sqrt(drop(jac %*% covariance[active,active] %*% jac))
      expect_equal(result$se[result$par == "q_curve_peak"],expected,tolerance=1e-9)
    }
  }
})
