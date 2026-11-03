# Gulf of Bothnia herring, 2026: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_herring_gulf_bothnia_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_herring_gulf_bothnia_2026.md).

## Expanded-library review

### Initial tinyAM fit

The package-default initialization converged (optimizer code 0, objective
561.3017, maximum absolute gradient 0.000285, positive-definite Hessian;
26 fixed and 1,365 random effects). The only numeric aggregate comparison is
approximate SSB: mean absolute percent difference 28.6%, terminal difference
+34.6%, and annual trend correlation 0.20 over 1980–2025. This is a sizable
disagreement under the documented differences above; no source observations
or starting values were changed to improve it. Recruitment and Fbar are marked
non-equivalent, and the M match is only the reused fixed input. The cached
dashboard and fitted object are available locally under
`results/cache/ices_herring_gulf_bothnia_2026/`.

## Translation limits

The initial tinyAM model uses years 1980–2025 and ages 1–15+, matching the
reported catch-composition age range. Recruitment is therefore represented
at age 1, not the source assessment's age 0, and recruitment is not treated
as a like-for-like comparison. The BIAS index is translated to thousand fish
by the shared unit converter; the script assigns within-year time 0.75 as an
approximation for the reported late-September/early-October survey.

tinyAM reconstructs catches by combining annual total biomass, age-number
proportions, and catch weights, then fits the resulting age-specific catches
with its lognormal likelihood. This does not reproduce SS3's total-catch and
age-composition likelihoods or their unknown precision settings. tinyAM uses
one shared estimated error for the acoustic index and age-specific
catchability; the assessment's exact index weighting is unavailable. SSB is
compared only as an approximate aggregate because tinyAM starts at age 1,
uses catch weight as its stock-weight proxy, and the report does not expose
the exact spawning-time settings. Source Fbar is unweighted, whereas tinyAM's
Fbar is abundance-weighted; it is not a direct comparison.
