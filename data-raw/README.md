# Simulated cod teaching data

`sim_cod_obs.R` generates one population and simulates three fisheries and two
surveys on it. Run the script from the tinyAM repository root with tinyAM
available in your R session. Use a clean
SimSurvey checkout at
[`42e5def`](https://github.com/PaulRegular/SimSurvey/commit/42e5def3819a5eec443deaabd5d268c7f126438d)
or install that revision into your R library before running it.

```r
# With a local checkout of the pinned SimSurvey revision
pkgload::load_all("../SimSurvey")
Sys.setenv(TINYAM_SIMSURVEY_REVISION = system2("git",
  c("-C", "../SimSurvey", "rev-parse", "HEAD"), stdout = TRUE))
source("data-raw/sim_cod_obs.R")
```

With SimSurvey installed by `remotes::install_github()` at that revision,
run the script directly. It stages the prepared inputs, truth and diagnostics at
`data-raw/cache/sim_cod/f-full-age/prepared-data.rds`. Set
`TINYAM_SIM_COD_OUTPUT` to a new file for another run. Existing staged files
are kept. Generation runs separately from package checks and documentation.
The seed-438 run took about 43 minutes on the development workstation.

Review the staged run before updating the packaged data. Check catch and survey
coverage, target shortfalls, realized M across years and ages, and the age-20
survivor losses. After reviewing the run, load its two objects and save them
as the packaged datasets.

```r
staged <- readRDS("data-raw/cache/sim_cod/f-full-age/prepared-data.rds")
sim_cod_obs <- staged$sim_cod_obs
sim_cod_truth <- staged$sim_cod_truth
tinyAM::check_tam_data(sim_cod_obs)

save(sim_cod_obs, file = "data/sim_cod_obs.rda", compress = "xz", version = 2)
save(sim_cod_truth, file = "data/sim_cod_truth.rda", compress = "xz", version = 2)
```

The [multi-fleet guide](../vignettes/articles/multiple-fisheries-surveys.Rmd)
fits ages 1–14 with an age-14 plus group. Its fixed-M and AR1-M fits use the
same reference of 0.3 per year. The website renders the dashboard using the
same fits as the guide. Previous simulations and fitted objects remain in
their original caches and worktrees. [Earlier exploration](https://github.com/PaulRegular/tinyAM/tree/a98322f8/analysis/simulated_cod)
