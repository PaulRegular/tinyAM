# Icelandic Atlantic wolffish: preliminary 2026 translation

The expansion-batch recipe uses SAIO preliminary revision
`6f08c06c49e4b73591427a1495bfb4568caf49f0`. It fits 1998–2025 at ages 4–16+,
the contiguous recent catch-at-age period. Original survey values remain on
their unresolved native scales; the script records the label adaptation
without rescaling observations. Source facts and gaps remain in the
[SAIO review](https://github.com/PaulRegular/saio/blob/e2ef4ae682bf4e80fd4b85edd0c19f1c62780c08/discovery/source_reviews/ices_atlantic_wolffish_iceland_2026.md).

The first attempt converged with package defaults: optimizer code 0, objective
325.7048, maximum gradient 0.0013464 and a positive-definite Hessian. Standard
errors were obtained; no source-start retry was used. The model is preliminary,
with simplified N/F and observation processes, free initial abundance and fixed
source M. No package mathematics or source observations changed.

Agreement requires review. Total abundance has mean absolute percent difference
48.0%, terminal difference +53.1% and correlation 0.994. N-at-age has 89.7%,
+99.2% and 0.945; F-at-age has 40.1%, −43.8% and 0.642. The common-definition
mature-biomass comparison is approximate because native SSB weights, sex and
spawning conventions are unresolved; its mean absolute percent difference is
141.3%. Fixed-M agreement is by construction.

The local dashboard and fit are under
`results/preliminary/cache/ices_atlantic_wolffish_iceland_2026/`. They are not
committed. Convergence does not promote this case into the curated collection.
