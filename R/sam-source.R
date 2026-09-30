# Ordinary values from a SAM fit; no parsing or serialized optimizer calls.
.sam_source <- function(sam_fit) {
  if (!is.list(sam_fit) || !all(c("data", "conf", "pl", "opt") %in% names(sam_fit)) ||
      !is.list(sam_fit$opt) || length(sam_fit$opt$objective) != 1L ||
      !is.finite(sam_fit$opt$objective)) cli::cli_abort("Supply a fitted SAM object with data, conf, pl and a finite optimizer objective.")
  dat <- sam_fit$data
  conf <- sam_fit$conf
  bounds <- c(conf$minAge, conf$maxAge)
  if (length(bounds) != 2L || anyNA(bounds) || any(bounds != trunc(bounds)) ||
      bounds[1] < 0 || bounds[2] <= bounds[1]) cli::cli_abort("Invalid SAM modeled age range.")
  ages <- seq.int(bounds[1], bounds[2])
  years <- dat$years
  if (length(years) < 2L || anyNA(years) || any(years != trunc(years)) || any(diff(years) != 1)) cli::cli_abort("SAM years must be consecutive.")
  aux <- dat$aux
  if (!is.matrix(aux) || !all(c("year", "fleet", "age") %in% colnames(aux)) ||
      nrow(aux) != length(dat$logobs)) cli::cli_abort("SAM observation rows and logobs do not align.")
  raw <- data.frame(year = aux[, "year"], fleet_id = aux[, "fleet"], age = aux[, "age"], obs = exp(dat$logobs))
  if (anyDuplicated(raw[c("year", "fleet_id", "age")])) cli::cli_abort("Duplicate SAM year/fleet/age observation rows.")
  nf <- length(dat$fleetTypes)
  if (!nf || anyNA(raw[c("year", "fleet_id", "age")]) ||
      any(!raw$fleet_id %in% seq_len(nf)) || any(!raw$year %in% years)) cli::cli_abort("Invalid SAM observation identifiers.")
  fleet_names <- attr(dat, "fleetNames")
  if (length(fleet_names) != nf) fleet_names <- paste("Fleet", seq_len(nf))
  if (length(dat$sampleTimes) != nf) cli::cli_abort("SAM sampling times do not match fleets.")
  if (any(!raw$age %in% ages)) cli::cli_abort("SAM observation ages fall outside the modeled age range.")
  if (any(!seq_len(nf) %in% raw$fleet_id)) cli::cli_abort("Every SAM fleet must contain observation rows.")
  fleet_ages <- lapply(seq_len(nf), function(i) range(raw$age[raw$fleet_id == i]))
  fleets <- data.frame(fleet_id = seq_len(nf), fleet_name = fleet_names, fleet_type = dat$fleetTypes,
    min_age = vapply(fleet_ages, min, numeric(1)), max_age = vapply(fleet_ages, max, numeric(1)),
    samp_time = dat$sampleTimes)
  mats <- lapply(seq_len(nf), function(i) {
    d <- raw[raw$fleet_id == i, ]
    ys <- seq.int(min(d$year), max(d$year))
    aa <- seq.int(fleets$min_age[i], fleets$max_age[i])
    m <- matrix(NA_real_, length(ys), length(aa), dimnames = list(year = ys, age = aa))
    m[cbind(match(d$year, ys), match(d$age, aa))] <- d$obs
    # Effort was already applied by SAM and cannot be recovered from logobs.
    attr(m, "effort") <- rep(NA_real_, length(ys))
    m
  })
  names(mats) <- fleet_names
  data <- list(catch = mats[dat$fleetTypes == 0], surveys = mats[dat$fleetTypes == 2])
  fields <- c(sw = "stockMeanWeight", cw = "catchMeanWeight", mo = "propMat",
              nm = "natMor", pf = "propF", pm = "propM")
  for (nm in names(fields)) {
    m <- dat[[fields[[nm]]]]
    if (length(dim(m)) == 3L && dim(m)[3] == 1L) m <- m[, , 1]
    if (is.matrix(m)) {
      if (is.null(rownames(m))) rownames(m) <- utils::head(years, nrow(m))
      if (is.null(colnames(m))) colnames(m) <- utils::head(ages, ncol(m))
    }
    data[[nm]] <- m
  }
  list(data = data, fleets = fleets, conf = conf, years = years, raw = raw,
       observation_weights = dat$weight)
}
