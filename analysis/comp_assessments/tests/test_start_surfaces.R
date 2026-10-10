source('analysis/comp_assessments/R/run_assessment.R')
rows <- data.frame(type = 'population', measure = 'numbers_at_age',
                   year = c(2001, 2000, 2001, 2000), age = c(2, 1, 1, 2),
                   value = c(4, 1, 3, 2))
original <- rows
surface <- .translation_start_surface(rows, 'population', 'numbers_at_age',
                                       2000:2001, 1:2, 1000)
stopifnot(identical(rows, original), identical(dimnames(surface), list(c('2000','2001'), c('1','2'))),
          identical(unname(surface), matrix(c(1000,3000,2000,4000), 2)))
fails <- function(x) inherits(try(.translation_start_surface(x, 'population',
  'numbers_at_age', 2000:2001, 1:2), silent = TRUE), 'try-error')
stopifnot(fails(rows[-1,]), fails(rbind(rows,rows[1,])))
rows$value[1] <- 0
stopifnot(fails(rows))
rows$value[1] <- Inf
stopifnot(fails(rows))
cat('Source starting-surface conversion tests passed.\n')
