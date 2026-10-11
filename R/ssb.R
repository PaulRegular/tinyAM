# The same biological definition supplies reports and recruitment parents.
.ssb_at_age <- function(N, W, P, Z, spawn_time = 0) {
  mature_biomass <- N * W * P
  if (all(spawn_time == 0)) mature_biomass else mature_biomass * exp(-spawn_time * Z)
}

.ssb_contributions <- function(N, W, P, Z, dat, N_plus = NULL, year = NULL) {
  "[<-" <- RTMB::ADoverload("[<-")
  time <- if (is.null(year)) dat$spawn_time else dat$spawn_time[year, ]
  out <- .ssb_at_age(N, W, P, Z, time)
  if (!is.null(N_plus)) {
    A <- length(dat$ages)
    plus_time <- if (is.null(year)) dat$spawn_time_plus else dat$spawn_time_plus[year, ]
    terminal_time <- if (is.null(year)) time[, A] else time[A]
    # Uniform timing can use effective W/P exactly as before.
    if (!all(plus_time == terminal_time)) {
      plus_W <- if (is.null(year)) dat$W_plus_input else dat$W_plus_input[year, ]
      plus_P <- if (is.null(year)) dat$P_plus_input else dat$P_plus_input[year, ]
      plus_Z <- if (is.null(year)) Z[, A] else Z[A]
      plus_ssb <- .ssb_at_age(N_plus, plus_W, plus_P, plus_Z, plus_time)
      if (is.null(year)) out[, A] <- rowSums(plus_ssb) else out[A] <- sum(plus_ssb)
    }
  }
  out
}

.ssb_timing_label <- function(dat) {
  time <- c(dat$spawn_time, dat$spawn_time_plus)
  if (all(time == 0)) return("Start of year")
  if (length(unique(time)) == 1L) return(paste("Spawning at year fraction", time[1L]))
  surface <- if (is.null(dat$spawn_time_plus)) dat$spawn_time else
    cbind(dat$spawn_time[, -ncol(dat$spawn_time), drop = FALSE], dat$spawn_time_plus)
  if (all(apply(surface, 1L, function(x) length(unique(x)) == 1L))) {
    "Year-specific spawning times"
  } else "Age-specific spawning times (may vary by year)"
}

.ssb_timing_table <- function(dat) {
  surface <- if (is.null(dat$spawn_time_plus)) dat$spawn_time else
    cbind(dat$spawn_time[, -ncol(dat$spawn_time), drop = FALSE], dat$spawn_time_plus)
  data.frame(year = rep(dat$years, ncol(surface)),
    age = rep(as.numeric(colnames(surface)), each = nrow(surface)),
    spawn_time = as.vector(surface), is_proj = rep(dat$is_proj, ncol(surface)))
}
