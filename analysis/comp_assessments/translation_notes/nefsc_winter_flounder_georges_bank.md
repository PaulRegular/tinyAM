# Georges Bank winter flounder source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/nefsc_winter_flounder_georges_bank.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/nefsc_winter_flounder_georges_bank.md).

## Expanded-library review

### tinyAM translation

The recipe fits the detailed 2022 VPA input period, years 1982-2021 and ages
1-7, with age 7 as the plus group. It uses the source catch numbers, annual
weights and maturity, fixed M=0.3, and the three native survey series. tinyAM
uses an IID abundance process, an AR1 F process, survey- and age-specific q,
and a shared catch observation SD. These are documented approximations to the
VPA estimation procedure.

The native MFSPAWN = 0.25 0.25 pair is preserved but not interpreted as a
tinyAM spawning time. The SSB comparison is therefore approximate. Recruitment
uses the same age-1 calendar-year definition. FFull is compared only as an
approximate ages-4-6 metric because official age-specific F values are not
available to rebuild a common weighted average. Fixed M agreement checks the
source input mapping, not estimation.

The first fit used package default initialization with no `start_par`. It
returned optimizer code 1 (false convergence), objective 1618.301, maximum
absolute gradient 0.000215, a successful `sdreport`, and a positive-definite
Hessian; it is not classified as converged. One retry from the first fit's
estimates returned optimizer code 0 (relative convergence), but the Hessian
was not positive definite and fixed-parameter covariance was unavailable.
The default fit is retained for the exploratory dashboard; neither fit passes
the numerical checks. No source observation or model assumption was changed
to improve agreement.
