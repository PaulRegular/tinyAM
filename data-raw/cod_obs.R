
library(tinyAM) # for cut_years and cut_ages functions

## Re-work small sub-set of inputs for the Northern Cod Assessment Model into TAM structure
inputs <- NCAM::inputs

catch <- as.data.frame.table(inputs$catch.num, responseName = "obs", stringsAsFactors = FALSE)
catch$age <- as.numeric(catch$age)
catch$year <- as.numeric(catch$year)
catch$unit <- "thousand fish"

index <- inputs$index |>
  subset(select = c("year", "age", "index"))
names(index)[names(index) == "index"] <- "obs"
index$survey <- "Fall 2J3KL"
index$samp_time <- 0.8
index$q_block <- cut_ages(index$age, c(min(index$age):6, max(index$age)))
index$obs <- index$obs * 12048556 / 1000 # use number of trawl units to bump up to total abundance
index$unit <- "thousand fish"

## Add Smith Sound biomass
ss <- NCAM::inputs$SS.bms[, c("Year", "Biomass")]
names(ss) <- c("year", "smith_sound_effect")
ss <- merge(data.frame(year = seq(min(ss$year), max(ss$year))), ss, by = "year", all = TRUE)
interp <- approx(ss$year, ss$smith_sound_effect)
ind <- which(is.na(ss$smith_sound_effect))
ss$smith_sound_effect[ind] <- interp$y[ind] # fill two missing values using linear interpolation
ss$smith_sound_effect <- ss$smith_sound_effect / 1000
index <- merge(index, ss, by = "year", all = TRUE)
index$smith_sound_effect[is.na(index$smith_sound_effect)] <- 0 # assumed to be zero outside the Smith Sound survey years

weight <- inputs$pwts[, c("year", "age", "weight")]
names(weight)[names(weight) == "weight"] <- "obs"
weight$unit <- "kg"

maturity <- inputs$mats[
  inputs$mats$year <= max(catch$year),
  c("year", "age", "mat")
]
names(maturity)[names(maturity) == "mat"] <- "obs"


# Label pre-moratorium, limited fishing gear + rec, closed, stewardship + rec / post-moratorium blocks
pop_covariates <- catch[c("year", "age")]
pop_covariates$F_y_block <- cut_years(catch$year, c(min(catch$year), 1992, 1998, 2002, 2006, 2023, max(catch$year)))
pop_covariates$F_a_block <- cut_ages(catch$age, seq(min(catch$age), max(catch$age), 2))


## Allometric M: sensu Cadigan, Rideout & Regular (2022), NAFO SCR 22/039
## Ocean-fish weight exponent from Lorenzen (1996).
allometric_b <- -0.305
reference_M <- 0.30
reference_ages <- 5:14

i_ref <- weight$age %in% reference_ages &
  is.finite(weight$obs) &
  weight$obs > 0

M0 <- reference_M /
  mean(weight$obs[i_ref]^allometric_b)

pop_covariates$allometric_M_assumption <-
  (M0 * weight$obs^allometric_b)[match(
    paste(pop_covariates$year, pop_covariates$age), paste(weight$year, weight$age)
  )]
# plot(pop_covariates$age, pop_covariates$allometric_M_assumption, log = "y", xlab = "age", ylab = "Allometric M")

## Log capelin per cod index (generated using a random walk)
cap_per_cod <- NCAM::inputs$covariates[, c("year", "cap_per_cod")]
pop_covariates$M_a_block <- cut_ages(pop_covariates$age, c(0, 2, 5, 8, 20))
pop_covariates <- merge(pop_covariates, cap_per_cod, by = "year", all.x = TRUE)
pop_covariates$log_capelin_per_cod <- log(pop_covariates$cap_per_cod)
pop_covariates$cap_per_cod <- NULL

## Impose consistent year range
pop_covariates <- pop_covariates[pop_covariates$year >= min(catch$year) &
                                   pop_covariates$year <= max(catch$year), ]
maturity <- maturity[maturity$year >= min(catch$year) &
                       maturity$year <= max(catch$year), ]
weight <- weight[weight$year >= min(catch$year) &
                   weight$year <= max(catch$year), ]
index <- index[index$year >= min(catch$year) &
                 index$year <= max(catch$year), ]

cod_obs <- as_tam_data(
  list(
    catch = catch,
    index = index,
    weight = weight,
    maturity = maturity,
    pop_covariates = pop_covariates
  )
)

usethis::use_data(cod_obs, overwrite = TRUE)
