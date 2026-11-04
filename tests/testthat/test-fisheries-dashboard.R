test_that("Fisheries pages preserve declared fleets and put the aggregate last", {
  expect_equal(.fishery_pages(data.frame(year=2000))$title,"Aggregate fishery")
  expect_equal(.fishery_pages(data.frame(fleet=character()))$kind,"aggregate")
  expect_equal(.fishery_pages(data.frame(fleet=NA_character_))$kind,"aggregate")
  one <- .fishery_pages(data.frame(fleet=c("A",NA)))
  expect_equal(one$title,c("A","Aggregate fishery"))
  catch <- data.frame(fleet=factor(c("A",NA,"B","A"),levels=c("B","A","Unused")))
  before <- catch
  pages <- .fishery_pages(catch)
  expect_equal(pages$title,c("B","A","Aggregate fishery"))
  expect_equal(pages$id,c("fishery-fleet-1","fishery-fleet-2","fishery-aggregate"))
  expect_identical(catch,before)
  expect_equal(.fishery_pages(data.frame(fleet=c('A "trawl"','B / C')))$title,
    c('A "trawl"','B / C',"Aggregate fishery"))
})

test_that("fleet yield uses reported ages and preserves unavailable references", {
  g <- expand.grid(year=2000:2002,age=1:2,model=c("fit","Official"),fleet=c("B","A"))
  g$est <- ifelse(g$model == "fit",g$age,NA_real_)
  g$est[g$fleet == "A"] <- 100*g$est[g$fleet == "A"]
  g$est[g$year == 2001 & g$fleet == "B" & g$model == "fit"] <- 0
  g$unit <- "tonnes"; g$is_proj <- g$year == 2002
  g$fold <- 2001; g$horizon <- pmax(0,g$year-g$fold)
  p <- g; p$est[p$model == "fit"] <- 3*p$est[p$model == "fit"]
  pop <- list(yield_by_fleet=g,yield_pred_by_fleet=p); before <- pop
  d <- .fleet_yield_data(pop,"B")
  expect_equal(d$obs[d$model == "fit" & d$year == 2000],3)
  expect_equal(d$pred[d$model == "fit" & d$year == 2000],9)
  expect_equal(d$obs[d$model == "fit" & d$year == 2001],0)
  expect_true(all(is.na(d$obs[d$is_proj | d$model == "Official"])))
  expect_true(all(is.na(d$pred[d$model == "Official"])))
  expect_equal(d$horizon[d$year == 2002],c(1,1))
  expect_true(all(d$unit == "tonnes"))
  expect_null(.fleet_yield_data(pop,"Unknown"))
  expect_true(all(is.na(.fleet_yield_data(list(yield_pred_by_fleet=p),"B")$obs)))
  partial <- pop
  partial$yield_pred_by_fleet$est[which(partial$yield_pred_by_fleet$model == "fit" &
    partial$yield_pred_by_fleet$fleet == "B" & partial$yield_pred_by_fleet$year == 2000)[1L]] <- NA_real_
  d <- .fleet_yield_data(partial,"B")
  expect_true(is.na(d$pred[d$model == "fit" & d$year == 2000]))
  expect_identical(pop,before)
})

test_that("Fisheries renders full fleet diagnostics and an unlabeled aggregate", {
  # Reporting-only fixture: split known outputs, with missing reference intervals.
  # This tests presentation without fitting an additional scientific model.
  ref <- default_fit[c("dat","pop","obs_pred","fixed_par","random_par","formula_effects")]
  class(ref) <- c("tam_ref","list")
  ref$dat$has_fleet <- TRUE; ref$dat$fleets <- c('B "trawl"',"A / gillnet")
  fractions <- c(.6,.4)
  ref$obs_pred$catch <- do.call(rbind,lapply(seq_along(fractions),function(i) {
    d <- default_fit$obs_pred$catch
    d$fleet <- ref$dat$fleets[i]; d$obs <- d$obs*fractions[i]; d$pred <- d$pred*fractions[i]
    d$std_res <- NA_real_; d
  }))
  ref$dat$obs$catch <- ref$obs_pred$catch
  ref$pop$F_fleet <- do.call(rbind,lapply(seq_along(fractions),function(i) {
    d <- default_fit$pop$F; d$fleet <- ref$dat$fleets[i]; d$est <- d$est*fractions[i]
    d$lwr <- d$upr <- NA_real_; d
  }))
  # Available means retain their dotted overlay; other references can lack them.
  ref$pop$mu_F_fleet <- ref$pop$F_fleet
  ref$pop$yield_by_fleet <- ref$obs_pred$catch[c("year","age","fleet","is_proj")]
  ref$pop$yield_by_fleet$est <- ref$obs_pred$catch$obs
  ref$pop$yield_by_fleet$unit <- "tonnes"
  ref$pop$yield_pred_by_fleet <- ref$pop$yield_by_fleet
  ref$pop$yield_pred_by_fleet$est <- ref$obs_pred$catch$pred
  ref$fixed_par$est <- NA_real_
  parameters <- default_fit$parameter_values
  reference_before <- ref
  path <- tempfile(fileext=".html"); on.exit(unlink(path))
  vis_tam(Unlabeled=default_fit,Fleet_reference=ref,output_file=path,
    open_file=FALSE,render_args=list(quiet=TRUE))
  expect_true(file.exists(path))
  html <- paste(readLines(path,warn=FALSE),collapse="\n")
  expect_match(html,'id="fishery-fleet-1"[^>]*data-navmenu="Fisheries"')
  expect_match(html,'id="fishery-fleet-2"[^>]*data-navmenu="Fisheries"')
  expect_match(html,'id="fishery-aggregate"[^>]*data-navmenu="Fisheries"')
  first <- regexpr('id="fishery-fleet-1"',html,fixed=TRUE)[1]
  second <- regexpr('id="fishery-fleet-2"',html,fixed=TRUE)[1]
  aggregate <- regexpr('id="fishery-aggregate"',html,fixed=TRUE)[1]
  expect_true(first < second && second < aggregate)
  fleet_html <- substr(html,first,aggregate-1L)
  expect_match(fleet_html,"Numbers-at-age",fixed=TRUE)
  expect_match(fleet_html,"F-at-age",fixed=TRUE)
  expect_match(fleet_html,"Residuals (age, year)",fixed=TRUE)
  expect_match(fleet_html,"(cohort)",fixed=TRUE)
  expect_match(fleet_html,"(expected)",fixed=TRUE)
  expect_match(fleet_html,'"dash":"dot"',fixed=TRUE)
  pages <- gregexpr('<div id="[^"]+" class="section level1"',html)[[1L]]
  aggregate_html <- substr(html,aggregate,min(pages[pages > aggregate])-1L)
  expect_match(aggregate_html,"Total yield",fixed=TRUE)
  expect_false(grepl("Residuals (age, year)",aggregate_html,fixed=TRUE))
  expect_identical(default_fit$parameter_values,parameters)
  expect_identical(ref,reference_before)
})
