# Does RV sampling correlation decline with age distance?

We compared three models for the positive RV observations in the
SimSurvey mortality and sampling study in the `observation-error-study` worktree.
All used the same calibration means, marginal log-SDs and held-out
observations. One assumed independence, one used a positive age AR1, and
one assigned the same positive correlation to every age pair within a
survey draw. The last is an exchangeable correlation model. The two
correlated models each estimated one parameter per scenario.

The benchmark used all 1,600 saved RV draws at standard effort, with
22,348 positive observations at modeled ages 1–14+. We split the 20 draws
for each fixed population and year into two halves, calibrated on each
half and scored the other. Means and SDs came from the calibration draws
for each population, year and age. The correlation parameter pooled
calibration information across the four histories and two years in each
spatial and mortality scenario. We scored the joint density of the
observed age vector and divided by the number of included ages. Zeros
and unavailable values remained outside this positive-log benchmark.

## Held-out score differences

Negative differences favor the model listed first. Each estimate averages
the two sampling years and the three main mortality regimes within each
history. The 95% paired t intervals use four histories, so they describe
variation in this pilot with limited precision. The separate temporal
mortality contrast appears in the generated results.

| Spatial field | AR1 minus independent | Exchangeable minus independent | Exchangeable minus AR1 |
|---|---:|---:|---:|
| Shared across ages | −0.321 [−0.806, +0.164] | −0.451 [−0.981, +0.079] | −0.130 [−0.179, −0.081] |
| Default | −0.135 [−0.319, +0.048] | −0.175 [−0.283, −0.067] | −0.039 [−0.125, +0.046] |
| Independent across ages | +0.001 [−0.004, +0.005] | 0.000 [0.000, 0.000] | −0.001 [−0.005, +0.004] |

When ages shared a spatial field, estimates at distant ages tended to move
together. Exchangeable correlation fit this pattern better than AR1 in
each mortality regime. The size of its gain over independence varied
among histories. Under default spatial fields, exchangeable correlation
improved on independence in all four mortality scenarios, including the
temporal sensitivity. Its average gain over AR1 was smaller and the
interval included zero. With independent age fields, every exchangeable
fit selected zero correlation.

Carry both covariance shapes into the controlled recovery experiment for
RV. This benchmark measures prediction of joint sampling errors at known
stock states. It does not establish whether an assessment can estimate
either correlation separately from population-process covariance or improve
stock estimates. No tinyAM model was refitted and no package option was added.

## Reproduce

From the package root, point the runner at the saved Z-study directory
containing `manifest.rds` and `checkpoints/`.

```sh
Rscript analysis/formula_expansion/benchmark_rv_covariance.R path/to/study-grid40
```

The runner checks both covariance densities against dense calculations
with unequal SDs, missing ages and singleton draws. It writes the
per-regime, per-history and per-draw results and source fingerprints to
the ignored `results/observation_error_z_covariance/study-grid40/rv_covariance_benchmark/`
directory.
