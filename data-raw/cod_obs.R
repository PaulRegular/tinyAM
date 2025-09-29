
library(TAM) # for cut_years and cut_ages functions

## Re-work small sub-set of inputs for the Northern Cod Assessment Model into TAM structure
inputs <- NCAM::inputs

catch <- as.data.frame.table(inputs$catch.num, responseName = "obs", stringsAsFactors = FALSE)
catch$age <- as.numeric(catch$age)
catch$year <- as.numeric(catch$year)
catch$sd_obs_block <- "catch"

# Label pre-moratorium, limited fishing gear + rec, closed, stewardship + rec / post-moratorium blocks
catch$F_y_block <- cut_years(catch$year, c(min(catch$year), 1992, 1998, 2002, 2006, max(catch$year)))
catch$F_a_block <- cut_ages(catch$age, seq(min(catch$age), max(catch$age), 2))

index <- inputs$index |>
  subset(select = c("year", "age", "index"))
names(index)[names(index) == "index"] <- "obs"
index$samp_time <- 0.8
index$q_block <- cut(index$age, c(min(index$age):6, max(index$age) + 1), right = FALSE)
index$sd_obs_block <- "index"

weight <- inputs$pwts[, c("year", "age", "weight")]
names(weight)[names(weight) == "weight"] <- "obs"

maturity <- inputs$mats[inputs$mats$year <= max(catch$year),
                        c("year", "age", "mat")]
names(maturity)[names(maturity) == "mat"] <- "obs"

nm <- as.data.frame.table(inputs$nm, stringsAsFactors = FALSE)
names(nm) <- c("age", "year", "M_assumption")
nm$age <- as.numeric(nm$age)
nm$year <- as.numeric(nm$year)
weight <- merge(weight, nm, by = c("year", "age"))

cod_obs <- list(catch = catch,
                          index = index,
                          weight = weight,
                          maturity = maturity)

usethis::use_data(cod_obs, overwrite = TRUE)
