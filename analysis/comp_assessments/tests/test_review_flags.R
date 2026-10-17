source("analysis/comp_assessments/R/review_flags.R")
x <- data.frame(assessment_id=letters[1:9],metric=c("N","N","ssb","recruitment","F","M","M","N","M"),
  comparison_status=c(rep("matched",7),"unavailable","matched"),n=c(rep(5,7),0,5),
  mean_absolute_percent_difference=c(20,20.01,0,0,NA,99,99,99,NA),
  terminal_mean_percent_difference=c(-20,0,-20.01,0,NA,0,0,99,NA),
  trend_correlation=c(.90,.95,.95,.899,NA,NA,NA,.1,NA))
mode <- setNames(c(rep("fixed input",5),"fixed input","estimated","estimated","estimated"),letters[1:9])
z <- assessment_review_flags(x,mode)
stopifnot(identical(z$assessment_id,c("b","c","d","g")))
mode["f"] <- "fixed approximation"
stopifnot("f" %in% assessment_review_flags(x,mode)$assessment_id)
cat("Review thresholds exclude equality, missing/unavailable statistics and fixed-M input checks; source-estimated M approximations remain explicit.\n")
