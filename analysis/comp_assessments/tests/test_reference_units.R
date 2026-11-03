pkgload::load_all(".",quiet=TRUE)
source("analysis/comp_assessments/R/database_to_tam_ref.R")
units <- c("thousand fish (inferred; Table 4.6.14 labels units NA)",
  "thousand fish; inferred from Table 4.6.14","million fish (inferred from source)",
  "fish per tow","native survey index units (unit unresolved)",
  "fish per tow (inferred from source)","unresolved (inferred from source)")
data <- data.frame(year=2000,age=seq_along(units),est=seq_along(units),
  se=.1,lwr=.5,upr=10,unit=units,notes="Original source note")
before <- data
clean <- .reference_count_units(data)
stopifnot(identical(data,before),identical(clean$unit[1:3],c("thousand fish","thousand fish","million fish")),
  identical(clean$unit[4:7],units[4:7]),
  identical(clean[setdiff(names(clean),c("unit","notes"))],data[setdiff(names(data),c("unit","notes"))]),
  identical(.reference_count_units(clean),clean))
for(i in 1:3)stopifnot(grepl(units[i],clean$notes[i],fixed=TRUE))
stopifnot(identical(tinyAM:::.unit_label(clean[1:2,],"Abundance"),"Abundance (thousand fish)"))
error <- tryCatch(tinyAM:::.unit_label(clean[1:3,],"Abundance"),error=function(e)e)
stopifnot(inherits(error,"error"),grepl("incompatible units",conditionMessage(error),fixed=TRUE))
outputs <- transform(data[1:2,],assessment_id="inferred_fixture",type="population",measure="numbers_at_age",
  age_group=NA_character_,value=est,source_type="source_table",source_reference="Table 4.6.14")
reference <- database_to_tam_ref("inferred_fixture",outputs,years=2000:2001,ages=1:2)
stopifnot(identical(outputs$unit,units[1:2]),all(reference$pop$N$unit[reference$pop$N$year == 2000] == "thousand fish"))
patterns <- tinyAM:::.tidy_patterns(reference$pop,reference=TRUE,target_units=list(numbers="thousand fish"))$population
stopifnot(sum(is.finite(patterns$est)) == 2L)
widget <- tinyAM::plot_matrix(patterns,metric=~metric)
spec <- widget$jsHooks$render[[1L]]$data
stopifnot(identical(spec$metrics[[1L]]$unit,"thousand fish"),
  any(grepl("Table 4.6.14",unlist(spec$metrics[[1L]]$cells[[1L]]$text),fixed=TRUE)))
cat("Reference count units, provenance, unchanged values and matrix labels passed.\n")
