# Simulated cod teaching data

`sim_cod_obs.R` is the entry point. Its visible sequence is one population,
distribution and rounding pipeline, followed by fishery and survey branches.
`sim_cod_helpers.R` contains the access, mortality-allocation, monitoring-expansion
and tidy-table calculations. SimSurvey is used only for regeneration.

Use a clean SimSurvey checkout at
`282e3c14a2fbb5d312a2fde79ee3e8d8511b4513`:

```r
source("data-raw/sim_cod_obs.R")
result <- build_sim_cod_obs("../SimSurvey")
# Inspect result, warnings and checks before replacing packaged inputs:
result <- build_sim_cod_obs("../SimSurvey", write_data = TRUE)
```

Run `Rscript data-raw/validate_sim_cod.R` for deterministic monitoring, allocation
and survey-capacity checks. The optional representation validator uses the pinned
SimSurvey checkout to compare sparse and dense fishing inputs.

The selected defaults retain the paper's recruitment and demographic centres,
with a 20×20 grid, one spatial field shared across ages/years, Z innovation SD
0.1, and length-at-age log-SD 0.03. RV/Sentinel set densities are 2/4 per 1,000
km². Cells are sampled without replacement; infeasible allocations fail before
fishing starts. Narrower length variation improves the approximation of
length-integrated selectivity by a curve evaluated at mean length.

Fleet allocations calibrate the geometric reference M baseline to 0.3, not
realized annual M. Inputs supply this known baseline (`M0`) and known population
mean length; annual latent states remain separate in `sim_cod_truth`. Source
papers, SimSurvey revision, settings, hashes and warnings are retained there.

The [multi-fleet guide](../vignettes/articles/multiple-fisheries-surveys.Rmd)
owns the final model specifications and exports its exact fits for the dashboard.
Both models have N process off and AR1 fleet F; they differ only in fixed versus
AR1 M. All 20 years enter the primary SSB comparison:

| Model | Mean SSB RE | Mean absolute RE | Terminal RE | Mean abundance RE | Mean recruitment RE |
|---|---:|---:|---:|---:|---:|
| Reference M | −6.6% | 8.8% | −13.7% | −12.0% | −13.7% |
| AR1 M | −3.7% | 6.7% | −13.1% | −21.8% | −28.9% |

Both passed numerical checks without changed optimizer controls or fitted-value
starts. This case was selected after exploration to meet an absolute mean annual
SSB error below 10%; these are single-realization errors, not Monte Carlo bias.
Abundance/recruitment recovery is weaker. The generator also couples fleet F and
remaining M through an existing death budget, unlike the fitted independent
processes. Agreement in SSB does not establish general accuracy or coverage.

Generation remains a substantial manual job. Run it separately from fitting and
checks. Configuration-specific caches under `data-raw/cache/sim_cod/` preserve
completed stages, RNG states and warnings. Original 80×80 and intermediate 40×40
data, fits and failed attempts remain in local ignored caches and earlier
worktrees; refreshing this example does not overwrite those files. Exploration
code and results are retained in branch history rather than active entry points.
