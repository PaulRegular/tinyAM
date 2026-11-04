# Faroes ling: preliminary 2026 translation

The expansion-batch recipe uses SAIO preliminary revision
`6f08c06c49e4b73591427a1495bfb4568caf49f0`, retaining 1996–2025 and ages 3–12+.
Age-specific relative survey values are retained without numerical rescaling.
February/March and August survey fractions, 2/12 and 7.5/12, are explicit
approximations. Unavailable effort definitions and gear-change downweights
are not invented. Official source details remain in the
[SAIO review](https://github.com/PaulRegular/saio/blob/e2ef4ae682bf4e80fd4b85edd0c19f1c62780c08/discovery/source_reviews/ices_ling_faroe_5b_2026.md).

The package-default attempt failed with an NA/NaN gradient. One permitted
initialization-only retry used aligned official N and F. It returned objective
764.7211, maximum gradient 0.0000422, a positive-definite Hessian and standard
errors, but optimizer code 1 (`singular convergence (7)`). The case therefore
remains numerically non-converged; no successful comparison is claimed.

Both attempt records and the fitted object are retained under
`results/preliminary/cache/ices_ling_faroe_5b_2026/`. No further starting-value
retry, observation change or model rescue was attempted.
