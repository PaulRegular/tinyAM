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

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/nefsc_summer_flounder.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

- Catch-at-age values are transcribed from Tables A3, A6, A10, A23, A28, and
  A32. Fleet tables report ages 0-10; ages 7-10 are summed to the official
  model's 7+ group, and Maine-Virginia and North Carolina commercial landings
  are summed for the commercial-landings fleet. Table A32's published total is
  preserved as a separate aggregate series, not as a fifth fleet. The fleet
  age series and published aggregate differ by as much as 941 thousand fish in
  an age-year cell, although annual totals are close. The difference is
  documented rather than reconciled by altering source values. Sampling times 0.25 and
  0.75 are seasonal midpoint approximations; the source report does not give
  within-season timing for these series. The other indices used by
  `F2018_BASE_V2` are not yet represented.
- M is fixed at age in the official model. The age 0-7+ values are 0.26, 0.26,
  0.26, 0.25, 0.25, 0.25, 0.25, and 0.24 per year (Table A90). These values
  are stored as source inputs, not as estimated M outputs.
- Table A90 reports mean November SSB weights at age for 2013-2017. Those
  published values are stored as the source schedule. No 2017 row is reported because 2017 fall maturity data were
  unavailable. A fit that extends through 2017 must state how it handles that
  final year; the conservative comparison period using the available maturity
  data is 1982-2016.
- Tables A87-A89 provide annual SSB, recruitment, F-at-age, and January 1
  abundance-at-age through 2017. Recruitment is at true age 0. The report
  provides terminal-year 90% MCMC intervals for aggregate F and SSB, but no
  matched age-specific uncertainty table was recovered.
