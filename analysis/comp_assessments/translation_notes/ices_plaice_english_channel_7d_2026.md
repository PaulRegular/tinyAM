# Eastern English Channel plaice (`ple.27.7d`): tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_plaice_english_channel_7d_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_plaice_english_channel_7d_2026.md).

## Expanded-library review

### tinyAM fit and comparison

The package-default fit failed with `NA/NaN gradient evaluation` and returned no optimizer result. One retry used N, F and q starting values from the local refit of the pinned SAM model to the archived inputs; the same observations, settings and likelihood were retained, and tinyAM estimated all parameters freely. The retry converged with optimizer code 0 (`relative convergence (4)`), objective 592.745, maximum absolute gradient 0.000958, successful standard-error estimation and a positive-definite Hessian. The local SAM refit did not reproduce ICES SAG's aggregate recruitment and SSB, so none of its fitted states are treated as official assessment outputs.

The converged tinyAM fit has elevated parameter correlation and residual advisories: five fixed-parameter pairs have absolute correlation above 0.95, some standardized-residual groups fall outside the expected spread, some residuals exceed three standard deviations, and some adjacent-year residual correlations exceed 0.5. Across 1980–2025, age-1 recruitment differs by a mean absolute 66.2%, has a terminal difference of +68.1%, and trend correlation 0.608. Start-of-year SSB differs by a mean absolute 226.4%, has a terminal difference of +250.6%, and trend correlation 0.929. These are substantial scale disagreements; no model settings were changed to tune them.

A deterministic cohort-survival sensitivity, retaining the recruitment RW and
all other settings, failed package-default initialization and then converged
from the permitted local-SAM starting values (objective 629.905, gradient
0.001074 and positive-definite Hessian). It improved trend correlation to
0.866 for recruitment and 0.948 for SSB, but worsened median absolute percent
differences from 60.9% to 113.7% for recruitment and from 236.8% to 336.8%
for SSB. Terminal differences also increased. This alternative is rejected;
the existing N-RW translation remains the closer scale match, though neither
fit closely reproduces the published series.
