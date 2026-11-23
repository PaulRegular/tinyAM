# Pinned SimSurvey fishery_engine.R (GPL-3). This study simulates one fishery
# and year at a time and draws monitoring independently from raw removals.
delta_fast_fishery_engine <- function(engine) {
  source <- list(
    setup = quote(access <- lapply(controls, .fishery_access, cells = cells, year = year)),
    available = quote(eligible <- remaining > 0 & selected[[f]] > 0 & death_remaining[a_index] > 0),
    cells = quote(eligible_cells <- access[[f]] & seq_len(nrow(cells)) %in% c_index[eligible]),
    rows = quote(rows <- which(c_index == cell)),
    debit = quote(death_remaining[a] <- death_remaining[a] - caught_at_age),
    monitoring = quote(result <- .monitor_fisheries(truth, controls, sim$years, sim,
      stats::setNames(kg, lengths))))
  replacement <- list(
    setup = quote({
      access <- lapply(controls, .fishery_access, cells = cells, year = year)
      select_age <- controls[[1L]]$selection$selected[, 1L]
      stopifnot(length(sim$years) == 1L, length(controls) == 1L, !biomass_target[1L],
        identical(as.numeric(selected[[1L]]), as.numeric(select_age[a_index])))
      cell_rows <- vector("list", nrow(cells))
      rows_by_cell <- split(seq_len(nrow(stock)), c_index)
      cell_rows[as.integer(names(rows_by_cell))] <- rows_by_cell
      age_cell <- c_index + (a_index - 1L) * nrow(cells)
      grouped <- rowsum(remaining, age_cell, reorder = FALSE)
      stock_age <- matrix(0, nrow(cells), length(sim$ages))
      stock_age[as.integer(rownames(grouped))] <- grouped[, 1L]
    }),
    available = quote(age_open <- death_remaining > 0 & select_age > 0),
    cells = quote(eligible_cells <- access[[f]] & rowSums(stock_age[, age_open, drop = FALSE]) > 0),
    rows = quote(rows <- as.integer(cell_rows[[cell]])),
    debit = quote({
      death_remaining[a] <- death_remaining[a] - caught_at_age
      stock_age[cell, a] <- stock_age[cell, a] - caught_at_age
    }),
    monitoring = quote(result <- truth))
  seen <- setNames(integer(length(source) + 2L), c(names(source), "density", "residual_deaths"))
  density <- list(quote(density <- numeric(nrow(cells))),
    quote(totals <- rowsum(remaining * selected[[f]] * (death_remaining[a_index] > 0) / cell_area,
      c_index, reorder = FALSE)), quote(density[as.integer(rownames(totals))] <- totals[, 1]))
  locate <- function(entries, code) which(vapply(entries, identical, logical(1), y = code))
  revise <- function(code) {
    hit <- locate(source, code)
    if (length(hit)) {
      seen[hit] <<- seen[hit] + 1L
      return(replacement[[hit]])
    }
    if (!is.call(code)) return(code)
    entries <- as.list(code)
    if (identical(entries[[1L]], quote(`{`))) {
      at <- locate(entries, density[[2L]])
      if (length(at)) {
        stopifnot(length(at) == 1L, identical(entries[[at - 1L]], density[[1L]]),
          identical(entries[[at + 1L]], density[[3L]]))
        seen[["density"]] <<- seen[["density"]] + 1L
        entries[[at - 1L]] <- quote(density <- drop(stock_age %*% (select_age * age_open / cell_area)))
        entries <- entries[-(at:(at + 1L))]
      }
      start <- locate(entries, quote(fishery_deaths <- starting - remaining))
      if (length(start)) {
        end <- locate(entries, quote(ledgers[[length(ledgers) + 1L]] <- ledger))
        stopifnot(length(start) == 1L, length(end) == 1L, end - start == 11L)
        seen[["residual_deaths"]] <<- seen[["residual_deaths"]] + 1L
        # The main study uses the upstream cohort and age-specific budget.
        # This terminal draw only allocates the leftover deaths by cell and length.
        entries <- entries[-(start:end)]
      }
    }
    for (i in seq_along(entries)[-1L]) {
      if (is.call(entries[[i]])) entries[i] <- list(revise(entries[[i]]))
    }
    as.call(entries)
  }
  body(engine) <- revise(body(engine))
  if (any(seen != 1L)) stop("The pinned SimSurvey fishery engine changed. Audit the study adapter.")
  engine
}
