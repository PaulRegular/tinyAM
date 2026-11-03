# Herring in ICES divisions 6.a and 7.b-c (combined historical unit): tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_herring_west_scotland_ireland.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_herring_west_scotland_ireland.md).

## Expanded-library review

### tinyAM translation result

The translation combines the two catch-number fleets by year and age for the
tinyAM fit; the database retains the fleet-specific inputs. It fits one pooled
correlated F process, preserves the three surveys and their q/SD sharing, and
uses independent survey-age residuals because tinyAM does not represent the
source AR(1) correlation across survey ages. Source F-at-age is summed across
the two areas only for the common total-F comparison. The fit uses package-
default starting values; no retry was made because the default fit passed
numerical checks.

The default fit converged (optimizer code 0; maximum absolute gradient 0.000716;
positive-definite Hessian). The dashboard and fitted object are cached locally
under `results/cache/ices_herring_west_scotland_ireland_2020/`. Observation
predictions are close to their inputs on average, but the derived population
series disagree substantially with the official outputs. Mean absolute percent
differences are 276.7% for N-at-age, 217.0% for total abundance, 197.3% for
biomass, 184.6% for SSB, and 266.5% for recruitment. Their terminal differences
are especially large and trend correlations are 0.535, 0.099, 0.392, 0.457,
and 0.119, respectively. Approximate total F-at-age differs by 54.9% on average
(terminal difference -78.4%; correlation 0.430), while the common N-weighted
Fbar differs by 45.2% (terminal difference -76.4%; correlation 0.600). Fixed M
is reused from the input and is not an independent fit comparison. These
results indicate weak agreement despite numerical convergence; the age-specific
source F output remains available by area in the canonical database. The
source report provides no numerical uncertainty tables, so none are shown for
the Official series.
