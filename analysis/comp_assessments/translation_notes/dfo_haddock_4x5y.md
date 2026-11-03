# 4X5Y haddock: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/dfo_haddock_4x5y.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/dfo_haddock_4x5y.md).

## Expanded-library review

### tinyAM translation

The recipe uses 1985–2015, ages 1–11+, the source catch and survey index series, fixed source M, age-specific q, and the framework's per-recruit maturity and spawning weights as a proxy. tinyAM's catch likelihood and AR(1) F process are simplifications of ADAPT VPA, which treats catch-at-age as VPA input and estimates population states through its own VPA formulation. It also uses one midyear survey timing approximation. These differences are documented in the model background and mean that agreement is descriptive rather than a likelihood replication.
