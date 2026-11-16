# Observation-family studies

The [PLD recovery report](pld_validation_report.Rmd) compares ordinary lognormal,
Poisson-link delta (PLD) lognormal, and PLD generalized gamma fits to the same
observations. The controlled tests use 32 exact PLD and 32 finite-age-sampling
histories, plus 16 histories each for rising and constant natural mortality.
A separate ablation removes the encounter contribution while keeping the
PLD lognormal positive density fixed.

The whole-area SimSurvey experiment fits the eight original histories at 20+
and 14+. Both assessment resolutions supply the prescribed constant M.
Histories use a dome-selected fishery and one RV survey over 30 years, a
temporary fishing pulse, and a true age-20+ group. Abundance at ages 14 and
older during years 16–23 is the primary stock recovery measure.

## Results

All 288 controlled family fits and 64 encounter-ablation fits passed. Under
exact PLD sampling, PLD lognormal reduced mean absolute log error for age-14+
abundance by 0.0217 versus ordinary lognormal (95% paired interval 0.0127 to
0.0306 lower). Adding the encounter terms to the same positive likelihood
changed error by −0.0001 (−0.0016 to 0.0013). Across the four controlled
experiments, switching the PLD positive component to generalized gamma
changed mean error by less than 0.002 in either direction. In the SimSurvey
20+ pilot, ordinary lognormal, PLD lognormal and PLD generalized gamma
passed in 8/8, 5/8 and 7/8 histories. The two PLD families had only four
successful pairs with each other. The estimated effect of generalized
gamma on older-age error in those pairs was +0.091 (−0.477 to 0.658). The
small pilot and departures from prescribed M limit inference from this
comparison.

## Fit and render

Run from the package root with R, tinyAM development dependencies and Pandoc.
To render the saved comparison, run

```sh
Rscript analysis/family_expansion/simulate_pld.R report
```

To refit on the frozen observations, use a new results directory so the
original run and its source hashes are preserved. Supply the controlled and
spatial input directories after the output directory to override their
defaults.

```sh
Rscript analysis/family_expansion/simulate_pld.R comparison analysis/family_expansion/results/pld-reproduction
Rscript analysis/family_expansion/simulate_pld.R report analysis/family_expansion/results/pld-reproduction
```

Input observations default to `pld-simple-refined` and `pld-whole-area`
under the previous study's ignored `analysis/formula_expansion/results/`.
These directories also preserve the original fits, logs and HTML. The
original report source and runners are in Git at `b161dbcd`. The new run
writes to `analysis/family_expansion/results/pld-family-comparison/`.
Checkpoints record every attempt, warnings, errors and input and source hashes.
Rendering reads summary files. The comparison ran against a frozen package
checkout at `b161dbcd`. Eighteen representative fits matched the updated
package to within 1e-9 for objectives, recovery and family parameters.

To regenerate the samples, run `controlled`, `audit` and `spatial` with
the same fresh output directory, then render the report there. Use SimSurvey
at [29a1b03](https://github.com/PaulRegular/SimSurvey/commit/29a1b03350c274f88dac108fb95f9446eac79233).
The audit checks the first two spatial histories. The spatial run fits the
same eight-history design as the original pilot. Set
`TINYAM_PLD_PACKAGE_SOURCE` to run from a different tinyAM checkout.

```powershell
$study = './analysis/family_expansion'
$output = "$study/results/pld-regenerated"
& "$study/run_pld_study.ps1" -Rscript $rscript -Phase controlled -Results $output
& "$study/run_pld_study.ps1" -Rscript $rscript -SimSurveySource $simsurvey -Phase audit -Results $output
& "$study/run_pld_study.ps1" -Rscript $rscript -SimSurveySource $simsurvey -Phase spatial -Results $output
& $rscript "$study/simulate_pld.R" report $output
```

The Windows launcher runs one below-normal worker and logs CPU and RAM.
It pauses after two consecutive checks below 2 GiB available RAM. Run the
same command with the same sources and inputs to resume from checkpoints.

[Sampled-state quantile residual validation](sampled_quantile_validation.md)
uses a separate simulation to compare marginal quantile residuals evaluated
at fitted states and sampled states.
