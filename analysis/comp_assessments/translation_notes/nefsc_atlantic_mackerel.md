# Atlantic mackerel source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/nefsc_atlantic_mackerel.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/nefsc_atlantic_mackerel.md).

## Expanded-library review

### tinyAM fit and comparison

The stock recipe fits 1968-2016, ages 1-10+, with an IID N process, AR1 F
deviations around an age-specific mean, fixed M=0.2, one lognormal catch
observation error, and separate catchability by survey-age. The official model's
time-constant fishery selectivity and its catch and age-composition likelihoods
are not represented directly. The aggregate egg SSB index is also omitted.
Accordingly this is a simplified comparison model, not a reproduction of ASAP.

The fit converged (optimizer code 0, relative convergence), with maximum
absolute gradient 0.00066, successful sdreport, and a positive-definite Hessian.
The terminal-year trajectory comparisons are strong for several quantities
(SSB trend correlation 0.95; N-at-age aggregate 0.88), but scale differences
remain material: terminal SSB is about 29% lower and terminal Fbar about 35%
higher in tinyAM. Recruitment comparisons use age 1 on both sides and have a
more modest trend correlation (0.61). M matches exactly by construction because
both models use fixed M=0.2; that agreement is not an independent validation.
These differences should be read alongside the structural approximations
above.
