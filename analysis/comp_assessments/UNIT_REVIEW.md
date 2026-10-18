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
