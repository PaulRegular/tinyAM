smooth_data <- function() transform(expand.grid(year=2000:2007,age=1:12,
  group=factor(c("B","A"),levels=c("B","A"))),obs=1,is_proj=FALSE)

test_that("univariate smooth bases and penalties match mgcv reference matrices", {
  d <- smooth_data()
  for (bs in c("tp","cr")) {
    form <- if (bs == "tp") ~s(age,k=5,bs="tp") else ~s(age,k=5,bs="cr")
    actual <- .parse_q_formula(form,d)
    sm <- mgcv::smoothCon(mgcv::s(age,k=5,bs=bs),d,absorb.cons=TRUE,modCon=3)[[1]]
    reference <- mgcv::smooth2random(sm,names(d),type=2)
    expect_equal(unname(actual$q_modmat[,-1,drop=FALSE]),unname(reference$Xf))
    expect_equal(unname(actual$q_terms[[1]]$design),unname(reference$rand[[1]]),ignore_attr=TRUE)
    p <- .q_term_parameters(actual$q_terms)
    term <- actual$q_terms[[1]]
    p[[term$sd_parameter]][] <- log(.4)
    p[[term$parameter]][] <- seq_len(ncol(term$design))/10
    effect <- .q_effects(p,list(obs=list(index=d),q_terms=actual$q_terms))
    expect_equal(effect$nll,-sum(dnorm(p[[term$parameter]],sd=.4,log=TRUE)))
    expect_equal(effect$contribution,drop(reference$rand[[1]] %*% p[[term$parameter]]))
  }
})

test_that("factor by smooths retain groups, contrasts and independent penalties", {
  d <- smooth_data(); before <- d
  actual <- .parse_q_formula(~group+s(age,k=5,bs="cr",by=group),d)
  bases <- mgcv::smoothCon(mgcv::s(age,k=5,bs="cr",by=group),d,absorb.cons=TRUE,modCon=3)
  expect_length(actual$q_terms,2L)
  expect_equal(vapply(actual$q_terms,function(x)x$groups[[1]]$group,character(1)),c("B","A"),ignore_attr=TRUE)
  for (i in seq_along(bases)) {
    reference <- mgcv::smooth2random(bases[[i]],names(d),type=2)
    expect_equal(unname(actual$q_terms[[i]]$design),unname(reference$rand[[1]]),ignore_attr=TRUE)
    other <- as.character(d$group) != bases[[i]]$by.level
    expect_true(all(actual$q_terms[[i]]$design[other,] == 0))
  }
  expect_identical(d,before)
})

test_that("tensor conversion includes every penalty and predicts from retained bases", {
  d <- smooth_data()
  actual <- .parse_q_formula(~t2(age,year,k=c(4,4),bs=c("cr","tp")),d)
  sm <- mgcv::smoothCon(mgcv::t2(age,year,k=c(4,4),bs=c("cr","tp")),
    d,absorb.cons=TRUE,modCon=3)[[1]]
  reference <- mgcv::smooth2random(sm,names(d),type=2)
  expect_length(actual$q_terms,length(reference$rand))
  expect_equal(unname(actual$q_modmat[,-1,drop=FALSE]),unname(reference$Xf))
  for (j in seq_along(reference$rand)) {
    expect_equal(unname(actual$q_terms[[j]]$design),unname(reference$rand[[j]]),ignore_attr=TRUE)
  }
  new <- d[1:12,]; new$age <- seq(1,12,length.out=12)+.1; new$year <- 2003.5
  prediction <- .smooth_prediction(sm,reference,new)
  # Independently reconstruct coefficients on the original mgcv basis.
  beta <- seq_len(ncol(reference$Xf))/10
  u <- lapply(reference$rand,function(x)seq_len(ncol(x))/7)
  index <- reference$pen.ind
  index[reference$rind] <- index[index > 0]
  reordered <- numeric(length(index))
  for (j in seq_along(u)) reordered[index == j] <- u[[j]]
  coefficients <- numeric(length(index))
  coefficients[reference$pen.ind == 0] <- beta
  coefficients[reference$pen.ind != 0] <- reordered[reference$rind]
  coefficients <- reference$trans.D*coefficients
  if (!is.null(reference$trans.U)) coefficients <- drop(reference$trans.U %*% coefficients)
  expected <- drop(mgcv::PredictMat(sm,new) %*% coefficients)
  value <- drop(prediction$fixed %*% beta)
  for (j in seq_along(u)) value <- value + drop(prediction$random[[j]] %*% u[[j]])
  expect_equal(value,expected,tolerance=1e-12)
})

test_that("smooth marginal Gaussian likelihood agrees with independent covariance", {
  d <- smooth_data()[1:96,]
  for (form in list(~s(age,k=5),~t2(age,year,k=c(4,4)))) {
    design <- .parse_q_formula(form,d)
    terms <- design$q_terms
    par <- c(list(beta=seq_len(ncol(design$q_modmat))/10),.q_term_parameters(terms))
    for (j in seq_along(terms)) par[[terms[[j]]$sd_parameter]][] <- log(.2+j/10)
    y <- sin(d$age/3)+cos((d$year-2000)/4)
    obj <- RTMB::MakeADFun(function(p) {
      effect <- .q_effects(p,list(obs=list(index=d),q_terms=terms))
      effect$nll-sum(RTMB::dnorm(y,drop(design$q_modmat %*% p$beta)+effect$contribution,.3,log=TRUE))
    },par,random=vapply(terms,`[[`,character(1),"parameter"),silent=TRUE)
    dense <- function(theta) {
      p <- obj$env$parList(x=theta)
      V <- diag(.3^2,nrow(d))
      for (term in terms) V <- V + exp(2*p[[term$sd_parameter]])*tcrossprod(term$design)
      residual <- y-drop(design$q_modmat %*% p$beta)
      .5*(length(y)*log(2*pi)+as.numeric(determinant(V,logarithm=TRUE)$modulus)+
        drop(crossprod(residual,solve(V,residual))))
    }
    expect_equal(as.numeric(obj$fn()),dense(obj$par),tolerance=1e-9)
    numerical <- vapply(seq_along(obj$par),function(i) {
      up <- down <- obj$par; up[i] <- up[i]+1e-5; down[i] <- down[i]-1e-5
      (dense(up)-dense(down))/(2e-5)
    },numeric(1))
    expect_equal(as.numeric(obj$gr()),numerical,tolerance=1e-6)
  }
})

test_that("temporal smooth forecasts hold training-year values without future bases", {
  d <- smooth_data()
  future <- transform(d[d$year == 2007,],year=2008,is_proj=TRUE)
  full <- rbind(d,future)
  for (form in list(~s(year,k=5),~t2(age,year,k=c(4,4)))) {
    historical <- .parse_q_formula(form,d)
    actual <- .parse_q_formula(form,full)
    for (j in seq_along(actual$q_terms)) {
      expect_false(is.null(.smooth_basis_transform(historical$q_terms[[j]],actual$q_terms[[j]])))
      expect_equal(actual$q_terms[[j]]$design[!full$is_proj,,drop=FALSE],historical$q_terms[[j]]$design)
      expect_equal(unname(actual$q_terms[[j]]$design[full$is_proj,,drop=FALSE]),
        unname(actual$q_terms[[j]]$design[full$year == 2007,,drop=FALSE]))
    }
    expect_equal(unname(actual$q_modmat[full$is_proj,,drop=FALSE]),
      unname(actual$q_modmat[full$year == 2007,,drop=FALSE]))
  }
})

test_that("smooth options, roles, missing predictors and redundant designs fail clearly", {
  d <- smooth_data()
  for (form in list(~s(age,bs="re"),~s(age,fx=TRUE),~s(age,m=3),~s(age,id=1),
    ~s(age,sp=.1),~s(age,k=2),~s(log(age)),~s(age,year),~te(age,year),~ti(age,year),
    ~t2(age,group),~s(age,by=year))) {
    expect_error(.parse_q_formula(form,d),"support|Supported|Unsupported|basis|bare|numeric|categorical|univariate|age and year|unnamed")
  }
  expect_error(.parse_q_formula(~age+s(age,k=5),d),"null-space")
  d$age[1] <- NA
  expect_error(.parse_q_formula(~s(age,k=5),d),"complete finite")
  expect_error(make_test_dat(N_settings=list(mu_form=~s(age,k=5))),"ordinary fixed")
  expect_error(make_test_dat(M_settings=list(mu_form=~s(age,k=5))),"F and q only")
  expect_error(make_test_dat(F_settings=list(process="rw",mu_form=~s(age,k=5))),"intrinsic RW")
  form <- local({basis_size <- 4L; ~s(age,k=basis_size,bs="cr")})
  expect_identical(.formula_input_variables(form),"age")
  expect_length(.parse_q_formula(~s(age,k=4,by=NULL),smooth_data())$q_terms,1L)
})

test_that("smooth simulations draw one coherent function per coefficient realization", {
  d <- smooth_data()
  design <- .parse_q_formula(~group+s(age,k=5,by=group),d)
  p <- .q_term_parameters(design$q_terms)
  for (term in design$q_terms) p[[term$sd_parameter]][] <- log(.3)
  set.seed(6117)
  effect <- .q_effects(p,list(obs=list(index=d),q_terms=design$q_terms),simulate=TRUE)
  expect_length(effect$parameters,2L)
  values <- split(effect$contribution,interaction(d$age,d$group,drop=TRUE))
  expect_true(all(vapply(values,function(x)length(unique(x)) == 1L,logical(1))))
  p[names(effect$parameters)] <- effect$parameters
  repeated <- .q_effects(p,list(obs=list(index=d),q_terms=design$q_terms))
  expect_equal(effect$contribution,repeated$contribution)
  expect_equal(effect$nll,repeated$nll)
})

test_that("smooth warm starts retain compatible age bases and reject changed bases", {
  settings <- list(years=2000:2007,ages=2:10,
    F_settings=list(process="off",mu_form=~factor(year)+s(age,k=5)),
    index_settings=list(q_form=~s(age,k=4)))
  d <- do.call(make_test_dat,settings)
  p <- make_par(d)
  for (term in c(d$F_terms,d$q_terms)) {
    p[[term$parameter]][] <- seq_along(p[[term$parameter]])/10
    p[[term$sd_parameter]][] <- log(.3)
  }
  fake <- structure(list(dat=d,parameter_values=unlist(p),
    obj=list(env=list(parList=function(par)p))),class="tam_fit")
  before <- fake$parameter_values
  settings$years <- 2000:2006
  peel <- do.call(make_test_dat,settings)
  start <- .start_parameters(fake,peel)
  merged <- .merge_start_par(make_par(peel),start)
  for (component in c("F", "q")) {
    source_data <- if (component == "F") d$F_data else d$obs$index
    target_data <- if (component == "F") peel$F_data else peel$obs$index
    key <- function(x) paste(x$year, x$age)
    rows <- match(key(target_data), key(source_data))
    field <- if (component == "F") "F_terms" else "q_terms"
    coefficient <- if (component == "F") "log_mu_f" else "log_q"
    matrix <- paste0(component, "_modmat")
    old_effect <- if (component == "F") .mean_effects(p,d,"F") else .q_effects(p,d)
    new_effect <- if (component == "F") .mean_effects(merged,peel,"F") else .q_effects(merged,peel)
    old_function <- drop(d[[matrix]] %*% p[[coefficient]]) + old_effect$contribution
    new_function <- drop(peel[[matrix]] %*% merged[[coefficient]]) + new_effect$contribution
    expect_equal(new_function, old_function[rows], ignore_attr = TRUE)
    expect_equal(new_effect$nll, old_effect$nll)
    for (term in peel[[field]]) {
      expect_false(is.null(start[[term$parameter]]))
      expect_equal(sum(start[[term$parameter]]^2),sum(p[[term$parameter]]^2))
      expect_equal(start[[term$sd_parameter]],p[[term$sd_parameter]])
    }
  }
  settings$F_settings$mu_form <- ~factor(year)+s(age,k=6)
  changed <- do.call(make_test_dat,settings)
  start <- suppressMessages(.start_parameters(fake,changed))
  expect_null(start[[changed$F_terms[[1]]$parameter]])
  expect_null(start[[changed$F_terms[[1]]$sd_parameter]])
  changed_par <- .merge_start_par(make_par(changed),start)
  q_rows <- match(paste(changed$obs$index$year,changed$obs$index$age),
    paste(d$obs$index$year,d$obs$index$age))
  old_q <- .q_effects(p,d)
  new_q <- .q_effects(changed_par,changed)
  expect_equal(drop(changed$q_modmat %*% changed_par$log_q)+new_q$contribution,
    (drop(d$q_modmat %*% p$log_q)+old_q$contribution)[q_rows],ignore_attr=TRUE)
  expect_equal(new_q$nll,old_q$nll)
  for (term in changed$q_terms) {
    expect_equal(start[[term$sd_parameter]],p[[term$sd_parameter]])
  }
  expect_identical(fake$parameter_values,before)
  active <- settings
  active$F_settings <- list(process="iid",mu_form=~factor(year)+s(age,k=5))
  active <- do.call(make_test_dat,active)
  merged <- .merge_start_par(make_par(active),.start_parameters(fake,active))
  initialized <- .initialize_F_states(merged,active)
  expected <- matrix(drop(active$F_modmat %*% initialized$log_mu_f)+
    .mean_effects(initialized,active,"F")$contribution,length(active$years))
  expect_equal(unname(initialized$log_f),unname(.f_state_mean(expected,active)))
})

test_that("smooth warm starts preserve functions and penalties under explicit rotations", {
  # Deliberately change coordinates instead of relying on a platform's choice
  # of eigenvector signs. Multiple tensor penalties must transfer together.
  for (tensor in c(FALSE, TRUE)) {
    d <- make_test_dat(years=2000:2007, ages=2:10,
      F_settings=list(process="off",mu_form=~factor(year)+s(age,k=5)),
      index_settings=list(q_form=if (tensor) ~t2(age,year,k=c(4,4)) else ~s(age,k=4)))
    old <- d
    p <- make_par(d)
    rotation <- function(n) {
      out <- diag(n)
      if (n == 1L) return(-out)
      out[1:2,1:2] <- matrix(c(cos(.6),sin(.6),-sin(.6),cos(.6)),2)
      out
    }
    for (component in c("F", "q")) {
      field <- if (component == "F") "F_terms" else "q_terms"
      coefficient <- if (component == "F") "log_mu_f" else "log_q"
      pieces <- d[[field]]
      signature <- pieces[[1]]$basis_signature
      fixed <- rotation(ncol(signature$design$fixed))
      random <- lapply(signature$design$random, function(x) rotation(ncol(x)))
      signature$design$fixed <- signature$design$fixed %*% fixed
      signature$design$random <- Map(function(x,r) x %*% r,signature$design$random,random)
      columns <- colnames(pieces[[1]]$smooth_fixed)
      old[[paste0(component,"_modmat")]][,columns] <-
        d[[paste0(component,"_modmat")]][,columns,drop=FALSE] %*% fixed
      p[[coefficient]][] <- seq_along(p[[coefficient]])/7
      for (j in seq_along(pieces)) {
        term <- pieces[[j]]
        old[[field]][[j]]$basis_signature <- signature
        old[[field]][[j]]$smooth_fixed <- term$smooth_fixed %*% fixed
        colnames(old[[field]][[j]]$smooth_fixed) <- columns
        old[[field]][[j]]$design <- term$design %*% random[[term$penalty]]
        p[[term$parameter]][] <- seq_along(p[[term$parameter]])/10
        p[[term$sd_parameter]][] <- log(.3+j/10)
      }
    }
    fake <- structure(list(dat=old,parameter_values=unlist(p),
      obj=list(env=list(parList=function(par)p))),class="tam_fit")
    original <- fake
    start <- .start_parameters(fake,d)
    merged <- .merge_start_par(make_par(d),start)
    for (component in c("F", "q")) {
      coefficient <- if (component == "F") "log_mu_f" else "log_q"
      field <- if (component == "F") "F_terms" else "q_terms"
      old_effect <- if (component == "F") .mean_effects(p,old,"F") else .q_effects(p,old)
      new_effect <- if (component == "F") .mean_effects(merged,d,"F") else .q_effects(merged,d)
      expect_equal(drop(old[[paste0(component,"_modmat")]] %*% p[[coefficient]])+
        old_effect$contribution, drop(d[[paste0(component,"_modmat")]] %*%
        merged[[coefficient]])+new_effect$contribution,ignore_attr=TRUE,tolerance=1e-10)
      expect_equal(old_effect$nll,new_effect$nll,tolerance=1e-10)
      for (term in d[[field]]) {
        expect_false(is.null(start[[term$parameter]]))
        expect_equal(sum(p[[term$parameter]]^2),sum(start[[term$parameter]]^2),tolerance=1e-10)
        expect_equal(p[[term$sd_parameter]],start[[term$sd_parameter]])
      }
    }
    expect_identical(fake,original)
    # A non-orthogonal penalty change invalidates the entire tensor, including
    # its null space and every smoothing SD, even if other blocks still match.
    bad <- old
    for (j in seq_along(bad$q_terms)) {
      bad$q_terms[[j]]$basis_signature$design$random[[1]] <-
        2*bad$q_terms[[j]]$basis_signature$design$random[[1]]
    }
    fake$dat <- bad
    rejected <- suppressMessages(.start_parameters(fake,d))
    for (term in d$q_terms) {
      expect_null(rejected[[term$parameter]])
      expect_null(rejected[[term$sd_parameter]])
    }
    expect_false(any(colnames(d$q_terms[[1]]$smooth_fixed) %in% names(rejected$log_q)))
    changed <- d$q_terms[[1]]
    changed$basis_signature$grid$age <- changed$basis_signature$grid$age+.01
    expect_null(.smooth_basis_transform(d$q_terms[[1]],changed))
    changed <- d$q_terms[[1]]
    changed$scope <- list(type="survey",name="different",id="survey:different")
    expect_null(.smooth_basis_transform(d$q_terms[[1]],changed))
    fake$dat <- old
    fake$dat$q_terms <- fake$dat$q_terms[-1]
    rejected <- suppressMessages(.start_parameters(fake,d))
    for (term in d$q_terms) expect_null(rejected[[term$parameter]])
  }
})

test_that("scoped F and q smooths report complete functions and signed null-space coefficients", {
  obs <- cod_obs
  obs$catch <- rbind(transform(obs$catch,fleet="Trawl",obs=.6*obs),
    transform(obs$catch,fleet="Gillnet",obs=.4*obs))
  obs$index <- rbind(transform(obs$index,survey="RV"),transform(obs$index,survey="Sentinel",obs=.8*obs))
  d <- make_test_dat(data=obs,years=2000:2007,ages=2:10,N_settings=list(process="off"),
    rec_settings=list(process="off",mu_form=~1),
    F_settings=list(process="off",mu_form=~0+for_fleet("Trawl",~s(age,k=5))+
      for_fleet("Gillnet",~t2(age,year,k=c(4,4)))),
    index_settings=list(q_form=~0+for_survey("RV",~s(age,k=4,bs="cr"))+
      for_survey("Sentinel",~logistic(age))))
  p <- make_par(d)
  for (term in c(d$F_terms,d$q_terms)) if (term$type == "smooth") {
    p[[term$parameter]][] <- seq_along(p[[term$parameter]])/20
    name <- if (startsWith(term$id,"F_")) "log_mu_f" else "log_q"
    p[[name]][colnames(term$smooth_fixed)] <- if (name == "log_mu_f") .3 else -.4
  }
  o <- RTMB::MakeADFun(function(p)nll_fun(p,d),p,silent=TRUE)
  expect_true(all(is.finite(o$gr())))
  r <- o$report()
  expect_equal(log(as.vector(r$F_fleet)),rowSums(r$eta_mu_F_smooth),ignore_attr=TRUE)
  expect_length(.formula_random_parameters(d),5L) # 1 + 3 F penalties and 1 q penalty
  expect_null(p$log_f)
  fit <- structure(list(dat=d,obj=o,rep=r,parameter_values=o$par),class="tam_fit")
  functions <- .tidy_formula_effects(fit,.95)$smooths
  expect_length(functions,3L)
  for (tab in functions) {
    expect_equal(length(unique(if (unique(tab$component) == "F") tab$fleet else tab$survey)),1L)
    expect_true(all(is.finite(tab$est)))
    expect_true(all(is.na(tab$se)))
    expect_true(is.character(.formula_effect_label(tab)))
  }
  fixed <- tidy_par(fit)$fixed
  expect_true(all(fixed$est[fixed$par == "mu_F_smooth"] == .3))
  expect_true(all(fixed$est[fixed$par == "q_smooth_fixed"] == -.4))
  expect_true(all(fixed$se_scale[fixed$par %in% c("mu_F_smooth","q_smooth_fixed")] == "reported"))
  expect_false(any(grepl("eta_.*smooth",names(tidy_rep(fit)))))
  # A prescribed fixed-parameter covariance isolates the RTMB delta calculation.
  covariance <- diag(.01,length(o$par))
  fit$sdrep <- RTMB::sdreport(o,par.fixed=o$par,hessian.fixed=solve(covariance),getReportCovariance=FALSE)
  full <- .tidy_smooth_functions(fit,.95,"F")
  expect_true(all(vapply(full,function(x)all(is.finite(x$se)),logical(1))))
  term <- d$F_terms[[1]]
  X <- cbind(term$smooth_fixed,term$design)
  expected <- sqrt(.01*rowSums(X^2))
  tab <- full[[1]]
  rows <- term$groups[[1]]$rows
  expect_equal(tab$se,expected[rows],ignore_attr=TRUE,tolerance=1e-10)
  stacked <- tidy_tam(fit,process_residuals=FALSE)
  expect_length(stacked$formula_effects$smooths,3L)
})

test_that("smooth and residual process cautions stay separate from convergence", {
  d <- make_test_dat(years=2000:2007,ages=2:10,
    F_settings=list(process="ar1",mu_form=~s(age,k=5)))
  expect_true("smooth_F_residuals" %in% .smooth_design_advisories(d)$issue)
  term <- d$F_terms[[1]]
  parameters <- setNames(c(0,0),c(term$sd_parameter,d$process_sd$F$parameter))
  fit <- list(dat=d,sdrep=list(pdHess=TRUE,par.fixed=parameters,
    cov.fixed=matrix(c(1,.95,.95,1),2)))
  expect_equal(.smooth_fit_advisories(fit)$issue,"smooth_F_separation")
  fit$sdrep$pdHess <- FALSE
  expect_equal(nrow(.smooth_fit_advisories(fit)),0L)
  expect_error(make_test_dat(F_settings=list(process="iid",age_breaks=c(2,5,15),
    mu_form=~s(age,k=5))),"conflicting")
})
