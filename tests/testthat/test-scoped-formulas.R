.scope_inputs <- function() {
  d <- cod_obs
  d$catch <- rbind(transform(d$catch,fleet="Trawl",obs=.6*obs),transform(d$catch,fleet="Gillnet",obs=.4*obs))
  d$index <- rbind(transform(d$index,survey="RV"),transform(d$index,survey="Sentinel",obs=.8*obs))
  d
}

test_that("scoped fixed effects use only their fleet/survey and preserve inputs", {
  obs <- .scope_inputs(); before <- serialize(obs,NULL)
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:6,
    F_settings=list(process="off",mu_form=~0+
      for_fleet("Trawl",~1+I(year-2000))+for_fleet("Gillnet",~1)),
    index_settings=list(q_form=~0+for_survey("RV",~1+age)+for_survey("Sentinel",~1+I(age^2))))
  p <- make_par(d)
  p$log_mu_f[] <- c(log(.2),.1,log(.3))
  p$log_q[] <- c(log(.1),.05,log(.2),-.01)
  o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
  r <- o$report()
  expected <- rep(.2*exp(.1*(0:5)),5)
  expect_equal(as.vector(r$F_fleet[, ,"Trawl"]),expected)
  expect_equal(as.vector(r$F_fleet[, ,"Gillnet"]),rep(.3,30))
  index <- d$obs$index
  q <- ifelse(index$survey == "RV",.1*exp(.05*index$age),.2*exp(-.01*index$age^2))
  expect_equal(unname(exp(r$log_q_obs)),q)
  expect_true(all(is.finite(o$gr())))
  expect_identical(serialize(obs,NULL),before)
  expect_true(all(grepl("fleet[",names(p$log_mu_f),fixed=TRUE)))
})

test_that("scoped curves agree with the existing categorical-by curves", {
  obs <- .scope_inputs()
  a <- make_test_dat(data=obs,years=2000:2005,ages=2:6,
    F_settings=list(process="off",mu_form=~fleet+logistic(age,by=fleet)),
    index_settings=list(q_form=~survey+logistic(age,by=survey)))
  b <- make_test_dat(data=obs,years=2000:2005,ages=2:6,
    F_settings=list(process="off",mu_form=~0+for_fleet("Trawl",~logistic(age))+for_fleet("Gillnet",~logistic(age))),
    index_settings=list(q_form=~0+for_survey("RV",~logistic(age))+for_survey("Sentinel",~logistic(age))))
  p <- make_par(a); q <- make_par(b)
  x <- RTMB::MakeADFun(function(p)nll_fun(p,a),p,silent=TRUE)
  y <- RTMB::MakeADFun(function(p)nll_fun(p,b),q,silent=TRUE)
  expect_equal(x$fn(),y$fn())
  expect_equal(x$report()$log_pred,y$report()$log_pred)
  expect_equal(x$report()$F,y$report()$F)
  expect_true(all(is.finite(y$gr())))
  expect_true(all(startsWith(vapply(b$q_terms,`[[`,character(1),"slope_parameter"),"log_q_slope_")))
})

test_that("scopes combine different curve types and retain sparse survey support", {
  obs <- .scope_inputs()
  obs$index <- obs$index[obs$index$survey == "RV" | obs$index$age == 5,]
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:6,
    F_settings=list(process="ar1",mu_form=~0+for_fleet("Trawl",~logistic(age))+for_fleet("Gillnet",~mono(age))),
    index_settings=list(q_form=~0+for_survey("RV",~logistic(age))+for_survey("Sentinel",~1)))
  p <- make_par(d)
  o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
  expect_true(all(is.finite(o$gr())))
  expect_length(d$q_terms,1L)
  expect_identical(unique(d$obs$index$age[d$obs$index$survey == "Sentinel"]),5)
  expect_true(all(grepl("fleet[Gillnet]",names(p$dmu_f),fixed=TRUE)))
  expect_true(all(grepl("fleet[Trawl]",names(p)[startsWith(names(p),"F_a50_")],fixed=TRUE)))
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:6,
    F_settings=list(process="off",mu_form=~0+fleet:factor(year)+
      for_fleet("Trawl",~0+logistic(age))+for_fleet("Gillnet",~0+mono(age))))
  expect_equal(ncol(d$F_modmat),12L)
  expect_true(is.finite(nll_fun(make_par(d),d)))
})

test_that("scoped q starts retain compatible coefficients but reject changed contrasts", {
  obs <- .scope_inputs()
  obs$index$group <- factor(ifelse(obs$index$age < 4,"young","old"),levels=c("young","old"))
  form <- ~0+for_survey("RV",~group)+for_survey("Sentinel",~group)
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:6,index_settings=list(q_form=form))
  p <- make_par(d); p$log_q[] <- seq_along(p$log_q)/10
  fake <- structure(list(dat=d,parameter_values=unlist(p),obj=list(env=list(parList=function(par)p))),class="tam_fit")
  peel <- make_test_dat(data=obs,years=2000:2004,ages=2:6,index_settings=list(q_form=form))
  expect_equal(.start_parameters(fake,peel)$log_q,p$log_q)
  contrasts(obs$index$group) <- contr.treatment(c("young","old"),base=2)
  changed <- make_test_dat(data=obs,years=2000:2004,ages=2:6,index_settings=list(q_form=form))
  expect_message(start <- .start_parameters(fake,changed),"Catchability mean starts.*basis changed")
  expect_null(start$log_q)
})

test_that("scoped formula environments and shared global terms survive compilation", {
  obs <- .scope_inputs()
  formula <- local({
    shift <- function(x) x-2
    trawl <- ~shift(age)
    selected <- "Trawl"
    ~0+for_fleet(selected,trawl)+for_fleet("Gillnet",~1)
  })
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:6,F_settings=list(process="off",mu_form=formula))
  expect_equal(unname(d$F_modmat[,"fleet[Trawl]:shift(age)"]),ifelse(d$F_data$fleet == "Trawl",d$F_data$age-2,0))
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:6,
    F_settings=list(process="off",mu_form=~0+age+for_fleet("Trawl",~1)+for_fleet("Gillnet",~1)))
  expect_equal(ncol(d$F_modmat),3L)
  expect_identical(colnames(d$F_modmat)[1],"age")
})

test_that("scope safeguards reject unsupported names, roles and redundancy", {
  obs <- .scope_inputs()
  for (form in list(~0+for_fleet("Missing",~1),~0+for_survey("RV",~1),~0+for_fleet("Trawl",~obs))) {
    expect_error(make_test_dat(data=obs,F_settings=list(process="off",mu_form=form)),"Unknown|only supported|response obs")
  }
  expect_error(make_test_dat(data=obs,F_settings=list(process="off",mu_form=~1+for_fleet("Trawl",~1)+for_fleet("Gillnet",~1))),"redundant")
  expect_error(make_test_dat(data=obs,F_settings=list(process="off",mu_form=~fleet+logistic(age,by=fleet)+for_fleet("Trawl",~0+logistic(age)))),"Overlapping structural")
  expect_error(make_test_dat(data=obs,F_settings=list(process="off",mu_form=~0+for_fleet("Trawl",~factor(age)+logistic(age))+for_fleet("Gillnet",~1))),"duplicates")
  expect_error(make_test_dat(data=obs,F_settings=list(process="off",mu_form=~0+for_fleet("Trawl",~for_fleet("Gillnet",~1)))),"nest")
  expect_error(make_test_dat(data=obs,N_settings=list(mu_form=~for_fleet("Trawl",~age))),"ordinary fixed")
  expect_error(make_test_dat(data=obs,catch_settings=list(sd_form=~for_fleet("Trawl",~age))),"only supported")
  expect_error(make_test_dat(data=obs,index_settings=list(q_form=~for_fleet("Trawl",~age))),"only supported")
})

test_that("scoped Gaussian effects simulate once per coordinate and keep independent densities", {
  obs <- .scope_inputs()
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:6,
    F_settings=list(process="off",mu_form=~0+for_fleet("Trawl",~1+iid(age,sd=.2))+for_fleet("Gillnet",~1)))
  p <- make_par(d); term <- d$F_terms[[1L]]
  expect_true(startsWith(term$parameter,"eta_mu_F_"))
  set.seed(77); effects <- .mean_effects(p,d,"F",simulate=TRUE)
  values <- split(effects$contribution[d$F_data$fleet == "Trawl"],d$F_data$age[d$F_data$fleet == "Trawl"])
  expect_true(all(vapply(values,function(x)length(unique(x)) == 1L,logical(1))))
  expect_equal(effects$contribution[d$F_data$fleet == "Gillnet"],rep(0,30))
  expected <- -sum(dnorm(effects$parameters[[term$parameter]],sd=.2,log=TRUE))
  expect_equal(effects$nll,expected)
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:6,
    F_settings=list(process="iid",mu_form=~0+for_fleet("Trawl",~1+rw(year,sd=.2))+for_fleet("Gillnet",~1+ar1(year,sd=.2,phi=.4))))
  expect_length(d$F_terms,2L)
  o <- RTMB::MakeADFun(function(p)nll_fun(p,d),make_par(d),silent=TRUE)
  expect_true(all(is.finite(o$gr())))
})

test_that("scoped q_power terms keep their abundance-time semantics", {
  d <- make_test_dat(data=.scope_inputs(),years=2000:2005,ages=2:6,
    index_settings=list(q_form=~0+for_survey("RV",~1+q_power(N))+for_survey("Sentinel",~1+q_power(N))))
  p <- make_par(d); p$log_q_power[] <- log(c(.8,1.2))
  o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
  r <- o$report(); rows <- d$obs_map$type == "index"; map <- d$obs_map[rows,]
  at_time <- log(r$N[cbind(match(map$year,d$years),match(map$age,d$ages))])-
    r$Z[cbind(match(map$year,d$years),match(map$age,d$ages))]*map$samp_time
  expected <- ifelse(map$survey == "RV",.8,1.2)*at_time
  expect_equal(unname(r$log_pred[rows]),expected)
  expect_true(all(is.finite(o$gr())))
  expect_error(make_test_dat(data=.scope_inputs(),index_settings=list(q_form=~q_power(N)+for_survey("RV",~0+q_power(N)))),"Overlapping q_power")
  q_only <- ~1+q_power(N)
  obs <- .scope_inputs(); obs$catch$N <- 1
  expect_error(make_test_dat(data=obs,F_settings=list(process="off",
    mu_form=~0+for_fleet("Trawl",q_only)+for_fleet("Gillnet",~1))),"only supported in index_settings")
})

test_that("unused scopes and future tail covariates do not influence historical preparation", {
  obs <- .scope_inputs()
  obs$catch$x <- ifelse(obs$catch$fleet == "Trawl",1,NA_real_)
  obs$catch$x[obs$catch$year > 2005 & obs$catch$fleet == "Trawl"] <- obs$catch$age[obs$catch$year > 2005 & obs$catch$fleet == "Trawl"]
  d <- make_test_dat(data=obs,years=2000:2005,ages=2:6,
    F_settings=list(process="off",mu_form=~0+for_fleet("Trawl",~0+x)+for_fleet("Gillnet",~1)))
  expect_true(all(is.finite(d$F_modmat)))
  expect_equal(d$F_modmat[d$F_data$fleet == "Trawl",1],rep(1,30),ignore_attr=TRUE)
})
