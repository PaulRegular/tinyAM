# Northern cod (2J3KL) source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/dfo_cod_2j3kl.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/dfo_cod_2j3kl.md).

## Expanded-library review

### Smith Sound translation

The production report supplies biomass (Table 11) and sampled counts (Table 12)
as separate observations. Seven years have both: 1995 and 1998-2003. These can
be translated to an index using `N[a] = B * p[a] / sum(p * w)`, with B converted
to kg, p the sample proportions and w in kg/fish. This is a reconstruction,
not a published abundance estimate or a replacement for canonical counts.

The framework research document (2025/034, p. 10, equation 2.14) uses stock
body mass in the Smith Sound biomass equation. Beginning-of-year stock weights
(production Table 9) are therefore preferable to mid-year catch weights
(Table 10, used for predicted fishery landings). They still approximate mass
at survey time: neither Smith-specific weights nor a within-year growth
surface is supplied. The zero age-0 and age-14-16 sample categories can be
removed before conversion without changing proportions or reconstructed
biomass. The age-1 biomass contribution is retained in the denominator even
when the tinyAM fit starts at age 2.

The integration trial fits independent Smith catchability by informative age,
separate from RV catchability, and a separate observation SD. It does not
reproduce xteNCAM's latent local population fraction, age-year availability
process or normalization `max(q) = 1`. The existing 1995-2007, age-block RV
catchability adjustment is retained as a proxy for offshore availability
during the period with approximately more than 10 kt in Smith Sound. It is
not an exact translation of the source's 1995-2009 availability process.

## Common definitions

Compare N and mortality at matching ages. Sum official N, biomass and mature
biomass only over the tinyAM ages; do not compare official age-0+ totals with
tinyAM age-2+ totals. Table 17's native age-0 recruitment is not compared with
tinyAM recruitment. The official N-at-age surface supplies a separately
labelled age-2 abundance comparison with tinyAM recruitment. Common F/M means use the matching N and
age-specific mortality surfaces with tinyAM's population weighting. Preserve
the native Table 17/18 aggregates for dashboard context and label comparisons
derived from rounded age-specific tables separately. No official output time
series are digitized, and no standard errors are invented.

## Juvenile timing and integration

The production report (2026/026, pp. 12-13) and framework research document
(2025/034, pp. 14-15) specify shared, time-invariant catchability between
Fleming and Newman, with independent q at ages 0 and 1. Newman sampling spans
July-November. Both reports print the Fleming season in reversed order as
October-September. The cited Fleming survey report (2022/056, pp. 1, 3 and
Table 8.3) resolves this to September-October; it also gives 2020 dates of
September 30-October 29. That supporting report is cached locally.

The integration trial uses season midpoints 9.5/12 for Fleming and 9/12 for
Newman. These are explicit translation approximations, not recovered annual
xteNCAM timing parameters; canonical sampling_time remains missing. The trial
extends tinyAM to ages 0-14. The source fixes F at ages 0-1 to zero
(2025/034, p. 15); tinyAM's current process estimates F at every modeled age.
That constraint cannot be reproduced with the current interface and is not
silently relaxed in the main age-2+ recipe. Trial diagnostics determine whether
the extra indices can be fitted numerically, not whether this limitation is
scientifically negligible.

## Recruitment and juvenile limitations

Framework 2025/034, equation 2.26, relates age-zero recruitment to same-year
SSB. The age-2+ translation retains RW recruitment rather than treating a
survivor approximation as the same stock–recruit relationship. Fleming and
Newman observations remain canonical but are omitted from the retained fit.
Their timing approximations and the source zero-F constraint at ages 0–1
prevent a direct translation with the current age structure.
