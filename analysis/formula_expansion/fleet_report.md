# Fleet and catch-weight validation

`fleet_validation.R` retains two reproducible trials. Run after loading the
package from the repository root:

```r
source('analysis/formula_expansion/fleet_validation.R')
recovery <- fleet_recovery()
vis_tam(attr(recovery, 'example_fit'),
  output_file = 'analysis/formula_expansion/results/fleet_dashboard.html',
  open_file = FALSE)
trial <- fleet_assessment_trial()
```

Attempts, warnings, failures, fits and dashboards remain under ignored
`analysis/formula_expansion/results/`. No retained stock recipe or aggregate
assessment export is changed.

## Equations and safeguards

Each fleet has an independent F process conditional on the shared parameters
and formula means. Age correlation applies within a fleet; it never connects
the last age of one fleet to the first age of another. Mortality uses the sum
of fleet F, and fleet catches use their own F divided by common total mortality.
SD formulas control sharing; process type and correlations are shared.

Catch weights are common across fleets. They change catch biomass and
catch-conditioned projections, not stock biomass or SSB. Hidden plus ages
retain their separate catch weights before abundance-weighted aggregation.
An absent catch-weight column uses stock weights.

Tests compare correlated-RW densities with independent dense multivariate
Normal calculations at positive, negative and zero correlation, including
unequal SDs. A 4,000-row simulation checks increment variance, age correlation,
cross-fleet independence and the unpenalized first row. Other tests cover
Baranov catches, plus groups, biological averaging, units, named warm starts,
uncertainty, evaluation, fleet-specific scores and information-criterion
observation matching. Changed factor bases or fleet definitions do not inherit
fixed coefficients whose coding has changed.

## Two-fleet recovery

Seed 332, 20 simulated 40-year assessments, seven ages, two logistic F curves,
deterministic cohort survival, RW recruitment and lognormal catch/index errors.
F residual variation is off. No rescue fitting is used. All 20 passed numerical
checks, including a positive-definite Hessian.

| Parameter | Fleet A truth / median | Fleet B truth / median |
|---|---:|---:|
| F plateau | 0.25 / 0.2507 | 0.40 / 0.4033 |
| Midpoint age | 3 / 2.9893 | 4 / 4.0537 |
| Slope | 1 / 1.0099 | 0.8 / 0.7966 |

These are point-recovery results in an informed generating model, not general
coverage claims or proof that complex fleet structures will be estimable.

## Source-supported summer flounder trial

The [2018 ASAP source review](../comp_assessments/source_reviews/nefsc_summer_flounder.md)
records four fisheries: commercial landings/discards and recreational
landings/discards. The local trial uses their numbers-at-age once, excluding
the separately published aggregate. Fleet and aggregate tables differ in some
cells, so comparison with the retained aggregate recipe is not a controlled
identical-observation experiment.

The trial retains 1982-2016, ages 0-7+, IID N, exponential initialization,
published fixed M, the two Bigelow surveys and their q structure, and static
common weights. It adds four AR1 F surfaces with fleet-age mean levels and
fleet-specific process SDs. It does not reproduce ASAP's composition likelihood,
selectivity structures or fleet-specific biology. Common stock weights remain
the catch-weight fallback. No extra model features or rescue adjustments are used.

The optimizer returned code 0, objective 2065.485 and maximum gradient 0.00155,
but the Hessian was not positive definite and fixed-effect uncertainty was not
finite (59 fixed, 1,392 random parameters). AR1 persistence exceeded 0.98.
This is **not a successful assessment fit**. Its scale/trajectory disagreements
are large; its dashboard is exploratory. No claim of improvement or recipe
replacement is warranted. Source coverage, biology and likelihood differences
remain separate from numerical support for this parameter-rich approximation.

## Default parity and remaining limits

Unlabeled single-fishery models retain their parameter layouts. At saved
default starts, objectives, gradients and reports are exactly unchanged. At
non-default parameter values, objectives remain unchanged, maximum gradient
differences are below 4e-15, and reports agree within 1e-10.

All comparative-assessment tests were attempted. The reduced-age Western Baltic
herring test remains numerically unstable: a direct reduced-age fit already
failed the Hessian check with preceding code; current runs can also fail its
recipe-based reduced-age fit. The recipe is unchanged, and its convergence
assertions were not weakened. This is a fitting limitation, not a reason to
claim a successful translation.

Local HTML payload checks found valid animation trace references in both the
49-plot simulated example and the 76-plot Official/tinyAM trial comparison.
Browser-level interaction checking was unavailable because the browser tool
blocked local-file navigation; no bypass was attempted. The separately reported
single-model blocked-logistic F blank-panel problem could not be reproduced
with the supplied formula and packaged data. At the user's request it remains deferred until a reproducible example is available.

No fleet-specific process types/correlations or catch weights, TAC allocation,
double-logistic selectivity or new composition likelihoods are introduced.
Predictive-interval coverage, process-integrated forecast scoring and uncertainty
in differences between forecast scores remain future work.
