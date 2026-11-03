source("analysis/comp_assessments/R/database_to_tam_ref.R")
row <- function(type,measure,age,value,fleet=NA_character_)data.frame(
  assessment_id="fleet_fixture",type=type,measure=measure,year=2000,age=age,
  age_group=NA_character_,value=value,se=NA_real_,lwr=NA_real_,upr=NA_real_,
  unit="per year",fleet=fleet,source_type="fixture",source_reference="fixture",
  notes="Unchanged source rates")
outputs <- rbind(row("population","numbers_at_age",1:3,c(100,60,40)),
  row("mortality","fishing_mortality_at_age",1:3,c(.1,.2,.4),"A"),
  row("mortality","fishing_mortality_at_age",1:3,c(.1,.3,.1),"B"))
reference <- database_to_tam_ref("fleet_fixture",outputs,years=2000:2001,
  ages=1:2,age_plus_group=2)
stopifnot(is.null(reference$pop$F),nrow(reference$pop$F_fleet) == 8L)
d <- reference$pop$F_fleet
stopifnot(isTRUE(all.equal(d$est[d$year == 2000 & d$age == 2 & d$fleet == "A"],.28)),
  isTRUE(all.equal(d$est[d$year == 2000 & d$age == 2 & d$fleet == "B"],.22)),
  all(is.na(d$est[d$year == 2001])))
array <- array(NA_real_,dim=c(2,2,2),dimnames=list(year=2000:2001,age=1:2,fleet=c("A","B")))
filled <- .fill_tam_report_array(array,d,2000:2001,1:2)
stopifnot(isTRUE(all.equal(filled["2000","2","A"],.28)),
  isTRUE(all.equal(filled["2000","2","B"],.22)),all(is.na(filled["2001",,])))
shuffled <- .fill_tam_report_array(array,d[rev(seq_len(nrow(d))),],2000:2001,1:2)
stopifnot(identical(filled,shuffled))
years <- 2000:2001; ages <- 1:2
template <- list(dat=list(years=years,ages=ages,
  has_fleet=TRUE,fleets=c("A","C")),pop=list(),rep=list(),obs_pred=list(),
  fixed_par=data.frame(),random_par=list())
error <- tryCatch(database_to_tam_ref("fleet_fixture",outputs,years=years,ages=ages,
  template=template),error=function(e)e)
stopifnot(inherits(error,"error"),grepl("Source F fleet names",conditionMessage(error),fixed=TRUE))
cat("Fleet reference identities, missing source cells, plus ages and 3-D report arrays passed.\n")
