sam_fixture <- function() {
  years <- 2000:2002
  ages <- 1:3
  grid <- expand.grid(year = years, age = ages)
  aux <- rbind(cbind(grid, fleet = 1),
               cbind(expand.grid(year = years, age = 1:2), fleet = 2),
               cbind(expand.grid(year = 2001:2002, age = 1), fleet = 3))
  aux <- as.matrix(aux[c("year", "fleet", "age")])
  bio <- function(x) matrix(rep(x, each = 3), 3, dimnames = list(year = years, age = ages))
  data <- list(years = years, aux = aux,
    logobs = log(c(10, 0, NA, 20, 21, NA, 30, 31, NA, 5, NA, 7, 0, 8, 9, 14, 15)),
    fleetTypes = c(0, 2, 2), sampleTimes = c(0, .125, .5), weight = rep(NA_real_, 17),
    stockMeanWeight = bio(c(1, 2, 4)), catchMeanWeight = bio(c(2, 2, 2)),
    propMat = bio(c(0, .5, 1)), natMor = bio(c(.2, .3, .4)), propF = bio(c(0, 0, 0)), propM = bio(c(0, 0, 0)))
  attr(data, "fleetNames") <- c("Catch", "Survey A", "Survey B")
  conf <- list(minAge = 1, maxAge = 3, maxAgePlusGroup = c(1, 0, 0),
    keyLogFsta = rbind(0:2, c(-1, -1, -1), c(-1, -1, -1)), corFlag = 0,
    keyVarF = rbind(c(0, 0, 0), c(-1, -1, -1), c(-1, -1, -1)),
    keyLogFpar = rbind(c(-1, -1, -1), c(0, 1, -1), c(0, -1, -1)),
    keyVarLogN = c(0, 1, 1), keyVarObs = rbind(c(0, 1, 1), c(2, 3, -1), c(4, -1, -1)),
    obsCorStruct = factor(rep("ID", 3), levels = c("ID", "AR", "US")),
    obsLikelihoodFlag = factor(rep("LN", 3), levels = c("LN", "ALN")),
    keyCorObs = matrix(-1, 3, 2), fbarRange = c(2, 3), stockRecruitmentModelCode = 0,
    noScaledYears = 0, keyScaledYears = numeric(), keyParScaledYA = matrix(numeric(), 0, 0),
    initState = 0, stockWeightModel = 0, catchWeightModel = 0, matureModel = 0, mortalityModel = 0,
    fixVarToWeight = matrix(numeric(), 0, 3), keyQpow = matrix(-1, 3, 3))
  pl <- list(logN = log(matrix(1:9, 3)), logF = log(matrix(rep(c(.1, .2, .3), 3), 3)), logFpar = log(c(.2, .6)))
  plsd <- lapply(pl, function(x) { x[] <- .1; x })
  structure(list(data = data, conf = conf, pl = pl, plsd = plsd,
    opt = list(objective = 10, convergence = 0),
    rep = list(predObs = rep(log(12), 17), obsCov = list(diag(c(.1, .3, .3)^2), diag(c(.2, .4)^2), matrix(.5^2, 1, 1))),
    sdrep = list(value = setNames(log(1:9), rep(c("logssb", "logR", "logfbar"), each = 3)), sd = rep(.1, 9))), class = "sam")
}
