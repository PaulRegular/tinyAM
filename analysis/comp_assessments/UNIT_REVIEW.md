# Dashboard and translation units

The previous converter used individual fish and kg. That was mathematically
valid, but differed from the packaged example and some comparison displays.
Translated catch, N and recruitment now use **thousand fish**, average weight
uses **kg**, and biomass, SSB and catch biomass use **tonnes**. Native database
values and source reconstructions retain their original units.

Count-based surveys use the same count scale. Rates retain their effort
denominator, for example thousand fish per tow; no area expansion is invented.
Relative survey indices retain their source scale. Generic or unresolved unit
descriptions remain in provenance but are omitted from plot labels. The
converter records final model units in its translation metadata.

For matching total-abundance counts, q is unitless. For a count-rate index it
has inverse-effort units, such as per tow. Neither the magnitude of q nor its
units make a relative/rate index an estimate of the whole population. A q of
one means numerical correspondence only for comparable abundance definitions,
ages and sampling times. Unresolved q units are not invented.

Reference trends and matrix views use the same displayed metric units.
Known natural-scale estimates and intervals are converted together; log-scale
SEs are unchanged. An unclassified native SE is not rescaled as though its
scale were known; its original value remains available in the native output.
Spawning-time definitions remain in Settings and comparison metadata, while
the axis title is simply SSB.

Only the Population patterns and Mortality patterns colour-ramp titles place
units on a second line. Other plot axes use a single-line label.

The recruitment ribbon defect was caused by missing terminal reference
intervals leaving invalid vertices at the polygon's closing edge. Ribbons now
close separately within valid runs, including interior gaps. Signed trend axes
include negative estimates and intervals.

Tests verify count-scale equivalence of objectives, gradients and population
surfaces with mapped starts, as well as unit labels, rate dimensions, reference
conversions and missing-interval polygons. Model specifications are unchanged;
source-based starting values are expressed on the new count scale.

The 35 stock scripts retain identical literal model/preparation calls compared
with `2cc0680`. The whiting refit differs by less than 3.4e-9 in log N after
conversion to a common physical scale; its objective differs by 1.8e-12.
Whiting and GOA pollock both passed numerical convergence and Hessian checks.

Fresh fits, dashboards, warnings and dashboard failures are cached locally.
The complete batch alone refreshes aggregate diagnostics/comparison tables and
the comparative figures. Non-converged dashboards are exploratory and do not
change aggregate comparison eligibility.

## Validation notes

The package's full suite passed 2,711 assertions without warnings, and
`R CMD check --no-manual` returned OK. Its single interactive-only browser
test was subsequently removed; HTML-generation assertions remain. The
updated dashboard tests then passed 196 assertions with browser and viewer
launches forbidden, without warnings or skips. The
colour-ramp-only formatting adjustment also passed 182 focused unit and
matrix assertions. Assigned examples and an HTML-file test avoid browser
launches. Ribbon tests cover missing intervals and stable animation traces.

The analysis-local suite retains six previously recorded failing files:
`test_age_subsets.R`, `test_bc_petrale_source.R`,
`test_dfo_herring_sgsl_translation.R`, `test_f_process_translations.R`,
`test_herring_translation.R` and `test_ices_norway_pout_translation.R`.
These concern existing optimization/source-validation issues or stale recipe
expectations; settled recipes were not changed to satisfy them. The general
database validator passed, and canonical database files were not modified.

The packaged cod data received only the average-weight unit-label change.
NCAM was unavailable, so its construction script was updated but full source
regeneration could not be exercised locally.

The first shared-worker batch was stopped after one worker consumed more
than six hours without returning the complete batch. A fresh run isolates
each stock in a separate R process, checkpoints every result and limits each
stock to 20 minutes. Scientific settings and the existing initialization
rules are unchanged; time limits, warnings and errors are retained locally.
The cause of the original stall is unresolved. Northern cod subsequently
converged in the isolated run in about five minutes. Its unit transformation
also preserved the joint objective and gradient to roughly 1e-12, reproduced
the previous marginal objective exactly, and preserved physical N, F and M.

The final isolated batch covers all 36 current assessment records: 21
converged, eight failed numerical checks, four have documented translation
blockers, two lack stock recipes, and BSAI Pacific ocean perch exceeded the
20-minute limit. The latter is recorded as a runner failure, not a fitted
comparison. No previous cached result substitutes for a failed attempt.

All 29 returned fits have newly rendered local dashboards, with no rendering
errors. The eight non-converged previews retain their explicit exploratory
warning. Saved whiting and GOA pollock Plotly data passed checks for finite
recruitment polygons, retained reference lines, negative recruitment-residual
axes and single-line ordinary axis labels. Missing reference projection rows
remain missing. All five comparative PNGs were regenerated and visually
inspected for legibility and clipping.

Fitting and package checks used R 4.4.1. Figure-only rendering used R 4.6.1 to
match the installed local plotting library; no models were refitted there.
Figure warnings (Windows font lookup, empty annotation groups and the existing
`geom_label(label.size)` deprecation) are retained in the ignored cache, as
are dashboard warnings and the isolated-fit progress log.

Previously committed aggregate results predate the current initialization
recipes and are not an isolated unit-conversion benchmark. The parity checks
above use identical specifications and starts mapped to the new count scale.
Numerical failures in the complete refresh are reported without tuning models
or changing their scientific specifications.
