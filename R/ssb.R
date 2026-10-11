# The same biological definition supplies reports and recruitment parents.
.ssb_at_age <- function(N, W, P, Z, spawn_time = 0) {
  mature_biomass <- N * W * P
  if (spawn_time == 0) mature_biomass else mature_biomass * exp(-spawn_time * Z)
}
