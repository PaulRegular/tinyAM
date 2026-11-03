# Summer flounder source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/nefsc_summer_flounder.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/nefsc_summer_flounder.md).

## Expanded-library review

### Translation status

The source inputs and outputs above have been added to the canonical database.
The initial tinyAM translation uses the published aggregate catch series once,
the Bigelow spring/fall age indices, and the static SSB-weight approximation.
Its 1982-2016 fit did not pass tinyAM's convergence check: the optimizer
returned code 0 with objective 757.41 and maximum gradient 0.00045, but the
Hessian was not positive definite (24 fixed and 552 random effects). The
standard-error calculation returned, but it is not treated as a successful
fit. No model-comparison dashboard was generated. No alternate settings were
tried to force convergence.

The official run retained more observations and assessment detail than this
tinyAM translation can represent. Any later comparison should be limited to
definitions and years that can be matched; in particular, the available
maturity schedule ends in 2016.
