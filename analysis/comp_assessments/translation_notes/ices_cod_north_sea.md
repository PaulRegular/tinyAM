# North Sea / Northern Shelf cod source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_cod_north_sea.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_cod_north_sea.md).

## Expanded-library review

### tinyAM translation of the 2025 assessment

The recipe fits 1983–2022 and ages 1–7+ as one population. Recruitment is a
separate RW; older-cohort survival departures are IID with a common SD and
free initial abundances. Source substock relationships are not represented.
F follows age-correlated random walks, rather than the source's complete
shared/scaled substock selectivity model. Its shared-selectivity code 4 combines
a cubic age effect with a stationary AR1 temporal level and is not identical
to `cor_rw`.

Only the mixed Q3+Q4 survey is fitted, retaining timing 0.75 and supplied
relative SDs. The six substock/recruitment index streams are omitted. Ordered
`mono(q_key)` imposes non-decreasing age catchability; this is a translation
assumption, not a recovered source constraint. Ordering must follow increasing
survey ages. The constraint reduces age-q flexibility without imposing a
specific logistic curve.

The supplied estimated M surface is treated as fixed, motivating the restricted
period with complete M. Weight and maturity are arithmetic means across available
Northwestern, Southern and Viking values, not abundance-weighted substock
aggregates. Biological GMRF uncertainty is omitted. Combined commercial catch
uses one independent lognormal SD; correlated catch errors, substock landings
and quarterly composition likelihoods are not reproduced.

For comparisons, source N, recruitment and biomass are summed across substocks;
F/M-at-age are N-weighted. Substock Fbar uses relative abundance at ages 2–4.
Intervals are not combined. These derived quantities support broad comparisons,
not exact equivalence with the original multistock model.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_cod_north_sea.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The revised translation fits one combined population over 1983–2022, with IID older-cohort survival deviations, correlated F random walks, fixed supplied M, and monotone non-decreasing survey catchability.

The `mono(q_key)` constraint is particularly important: it prevents the reduced observation model from estimating declining catchability at older ages, while allowing flexible increases and plateaus.

The comparison evaluates whether this simpler configuration can approximate broad official stock dynamics. It does not reproduce the official three-substock likelihood, selectivity model, or full uncertainty.
