# 4X5Y haddock: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/dfo_haddock_4x5y.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/dfo_haddock_4x5y.md).

## Expanded-library review

### tinyAM translation

The recipe uses 1985–2015, ages 1–11+, the source catch and survey index series, fixed source M, age-specific q, and the framework's per-recruit maturity and spawning weights as a proxy. tinyAM's catch likelihood and AR(1) F process are simplifications of ADAPT VPA, which treats catch-at-age as VPA input and estimates population states through its own VPA formulation. It also uses one midyear survey timing approximation. These differences are documented in the model background and mean that agreement is descriptive rather than a likelihood replication.

## Curated fit on SAIO `cecb45ad`

The 1985–2015 fit passed numerical checks from default starts. The curated recipe sets `F_settings$mean_ages = 6:10` to match the abundance-weighted age range in DFO Research Document 2017/006, Table 12. The fitted likelihood and SSB series are identical to the preliminary baseline. The corrected Fbar comparison has a median absolute difference of 8.5%, trend correlation 0.923 and a terminal-year difference of +49.9%.

Over 31 years, the approximate SSB comparison has a median absolute difference of 4.9%, trend correlation 0.988 and a terminal difference of +8.0%. Median differences for abundance and recruitment are 4.1% and 10.3%, with correlations of 0.942 and 0.784 respectively. The translated stock biology uses the source framework's static per-recruit weights and maturity.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/dfo_haddock_4x5y.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The fit recipe selects the official 1985–2015 assessment years and age ranges and pools source catch ages 11–16 into age 11+.

Table 15 of the framework report supplies maturity and spawning weight-at-age for per-recruit calculations. They are not represented as annual input surfaces used by the 2016 VPA. The fit comparison therefore treats SSB as approximate; the official model reports its native age-4+ biomass measure.
