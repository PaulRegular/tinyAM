# Observation-family studies

The [SimSurvey report](delta_validation_report.Rmd) compares six observation
families on paired 30-year SimSurvey histories. The two fisheries target 90%
or 10% of total mortality with flat selectivity. Both see the same stock and
start-of-year survey with catchability fixed to one. Each share is sampled
at half, standard, and double monitoring effort.

The primary fits use age–year AR(1) fishing mortality and a fixed mean
natural mortality calculated from the realized removal ledger. The
standard-effort sensitivity also estimates age–year AR(1) deviations in M
for the low-fishing scenario. The initial high-fishing attempts remain
available for comparison of numerical diagnostics.
Two figures summarize how the sampling design affects the frequency of zeros
and when improved observation predictions translate into better estimates of
older fish. The report includes paired results for each family and setting.

## Run and render

Run from the package root with tinyAM's development dependencies and Pandoc.
The SimSurvey source must be at [29a1b03](https://github.com/PaulRegular/SimSurvey/commit/29a1b03350c274f88dac108fb95f9446eac79233).

```powershell
$study = './analysis/family_expansion'
$output = "$study/results/default-z-fishing-share"
& "$study/run_delta_study.ps1" -Rscript $rscript -SimSurveySource $simsurvey -Phase spatial -Results $output
& "$study/run_delta_study.ps1" -Rscript $rscript -SimSurveySource $simsurvey -Phase variable -Results $output
& $rscript "$study/simulate_delta.R" report $output
```

The spatial phase audits the first two histories before fitting and saves
each generated history and fit attempt. Use the same command and sources
to resume. Before the two SimSurvey phases finish, rendering shows the
available progress. The final HTML, two summary figures and fit checkpoints
are in the selected output directory. The optional `exact` phase tests fits
to data generated from a specified observation family. Its existing checkpoints
are retained separately from the evidence summarized in the report.
The Windows launcher uses
one below-normal, single-thread worker and records CPU and RAM. It pauses
after two consecutive readings below 2 GiB available memory.

The former Poisson-link comparison and selective-survey pilot remain in
their original result directories. The original study sources are in Git
at `b161dbcd`.

[Sampled-state quantile residual validation](sampled_quantile_validation.md)
is a separate study of state-space residuals.
