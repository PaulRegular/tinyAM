root <- 'analysis/comp_assessments'
cache <- file.path(root,'source_cache/afsc_pollock_ebs_2024')
ctl <- trimws(readLines(file.path(cache,'control.dat'),warn=FALSE))
read_control <- function(label,n=1) {
 i <- which(ctl==paste0('#',label)); stopifnot(length(i)==1)
 as.numeric(ctl[i+seq_len(n)])
}
m <- read_control('natmort_in',15)
stopifnot(identical(m,c(.9,.45,rep(.3,13))),read_control('phase_natmort')==-6,
          read_control('switch_pred_mort')==0)
report <- paste(readLines(file.path(cache,'2024_SAFE_text.txt'),warn=FALSE),collapse='\n')
stopifnot(grepl('constant natural mortality rates at age',report,fixed=TRUE))
p <- file.path(root,'database/inputs.csv'); x <- read.csv(p,stringsAsFactors=FALSE,check.names=FALSE)
id <- 'afsc_pollock_ebs_2024'; x <- x[!(x$assessment_id==id & x$type=='M'),]
for(age in 1:15) {
 r <- as.list(setNames(rep('',ncol(x)),names(x)))
 r$assessment_id <- id; r$type <- 'M'; r$measure <- 'natural_mortality_at_age'; r$basis <- 'per_year'
 r$year <- 1964; r$year_basis <- 'calendar_year'; r$age <- age; r$value <- m[age]; r$unit <- 'per year'
 r$source_type <- 'native_model'; r$source_reference <- 'https://github.com/noaa-afsc/EBS_pollock; accepted 2024 control.dat, natmort_in; SAFE section 5.3.1'
 r$notes <- 'Fixed constant age vector valid throughout 1964–2024, stored once; age 15 is 15+. phase_natmort=-6 and switch_pred_mort=0. Accepted report explicitly confirms age-constant-through-time M; predation extensions are not active.'
 x <- rbind(x,as.data.frame(r,stringsAsFactors=FALSE))
}
write.csv(x,p,row.names=FALSE,na='')
cat('\n## Fixed natural mortality resolved\n\nSAFE section 5.3.1 explicitly states constant natural mortality rates at age for M23. Native control.dat supplies 0.9 at age 1, 0.45 at age 2 and 0.3 at ages 3–15, fixes the mortality-scaling parameter at phase -6, and sets switch_pred_mort=0. The implementation copies this vector into annual M when optional alternative mortality switches are inactive. The 15-value supplied vector is now canonical input, stored once for the full modeled period, with age 15 the plus group. It is not a fitted output or an annual predation estimate. scripts/026_import_ebs_pollock_mortality.R checks the source controls and report wording.\n',file=file.path(root,'source_reviews/afsc_pollock_ebs.md'),append=TRUE)
message('Fixed 15-age EBS pollock mortality input recovered from native controls and report.')
