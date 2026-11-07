# Simulated cod teaching data

`sim_cod_obs.R` generates one population and simulates three fisheries and two
surveys on it. Run the script from the tinyAM repository root with tinyAM
available in your R session. Use a clean
SimSurvey checkout at
[`b6c4914`](https://github.com/PaulRegular/SimSurvey/commit/b6c491473e887ca0c913df1f3d155267515ee99f)
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
`data-raw/cache/sim_cod/seed17-simsurvey-estimators/prepared-data.rds`. Set
`TINYAM_SIM_COD_OUTPUT` to a new file for another run. Existing staged files
are kept. Generation runs separately from package checks and documentation.

SimSurvey's `run_catch()` estimates catch at age from fleet monitoring.
`sim_surveys() |> run_strat()` estimates indices for both survey programs,
and `strat_error(scope = "survey_frame")` supplies the matching truth.

Seed 17 generates the population. A separate seed, 1701, generates annual
log-scale random walks for the three fleet targets. Starting in year 11, the
trawl target is 20% of its underlying trajectory to represent bycatch after
directed fishing stops. Gillnet and longline continue on their own trajectories.
The reference M centre is 0.28 per year. A first run with a centre of 0.30
gave mean realized M of 0.312 before the cut and 0.337 after it.

Review the staged run before updating the packaged data. Check catch and survey
coverage, target shortfalls, realized M across years and ages, and the age-20
survivor losses. The diagnostics show realized M over ages 5–10 before and after
the reduction. After reviewing the run, load its two objects and save them as
the packaged datasets.

```r
staged <- readRDS("data-raw/cache/sim_cod/seed17-simsurvey-estimators/prepared-data.rds")
sim_cod_obs <- staged$sim_cod_obs
sim_cod_truth <- staged$sim_cod_truth
tinyAM::check_tam_data(sim_cod_obs)

save(sim_cod_obs, file = "data/sim_cod_obs.rda", compress = "xz", version = 2)
save(sim_cod_truth, file = "data/sim_cod_truth.rda", compress = "xz", version = 2)
```

The [multi-fleet guide](../vignettes/articles/multiple-fisheries-surveys.Rmd)
compares a fleet-specific fit with one based on aggregate catch. Both fit ages
1–14 with an age-14 plus group and fix M at 0.3 per year. The website renders
the dashboard using the same fits as the guide. Previous simulations and
fitted objects remain in their original caches and worktrees. [Earlier exploration](https://github.com/PaulRegular/tinyAM/tree/a98322f8/analysis/simulated_cod)
