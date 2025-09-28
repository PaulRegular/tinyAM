
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

## Set-up projection rows
n_proj <- 3
max_year <- max(catch$year)
proj_years <- seq(max_year + 1, max_year + n_proj)
mean_years <- seq(max_year - n_proj + 1, max_year)
mean_catch <- aggregate(obs ~ age, catch, mean, subset = year %in% mean_years) |>
  cbind(catch[catch$year == max_year, c("sd_obs_block", "F_y_block", "F_a_block")])
proj_catch <- replicate(n_proj, mean_catch, simplify = FALSE)
for (i in seq_along(proj_years)) {
  proj_catch[[i]] <- cbind(data.frame(year = proj_years[i], proj_catch[[i]]))
}
proj_catch <- do.call(rbind, proj_catch)
proj_catch$is_proj <- TRUE
catch$is_proj <- FALSE
catch <- rbind(catch, proj_catch)



index <- inputs$index |>
  subset(select = c("year", "age", "index"))
names(index)[names(index) == "index"] <- "obs"
index$samp_time <- 0.8
index$q_block <- cut(index$age, c(min(index$age):6, max(index$age) + 1), right = FALSE)
index$sd_obs_block <- "index"

weight <- inputs$pwts[, c("year", "age", "weight")]
names(weight)[names(weight) == "weight"] <- "obs"

maturity <- inputs$mats[, c("year", "age", "mat")]
names(maturity)[names(maturity) == "mat"] <- "obs"

nm <- as.data.frame.table(inputs$nm, stringsAsFactors = FALSE)
names(nm) <- c("age", "year", "M_assumption")
nm$age <- as.numeric(nm$age)
nm$year <- as.numeric(nm$year)
weight <- merge(weight, nm, by = c("year", "age"))

northern_cod_data <- list(catch = catch,
                          index = index,
                          weight = weight,
                          maturity = maturity)

usethis::use_data(northern_cod_data, overwrite = TRUE)
