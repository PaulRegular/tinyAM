# Celtic Sea common sole (ICES divisions 7.f-g): tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ck_ices_wgsce_solea_cs7f_7g.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ck_ices_wgsce_solea_cs7f_7g.md).

## Expanded-library review

### tinyAM translation

The recipe will fit 1971–2025 and ages 1–10. It retains source catch and stock weights and maturity, fixes M at the source value, uses an AR(1) F process, and uses the age-specific UK(E&W)-BTS-Q3 index at ages 1–5. The source's Q3 survey timing is represented by 0.625, the midpoint of the quarter; this is a translation assumption because exact within-quarter timing is not specified. The age-specific q and observation-error structures are simplified to forms tinyAM can express.

The five commercial biomass indices remain in the database but are omitted from the fit because tinyAM's current direct survey likelihood expects age-specific abundance information. Fitting only the age-specific BTS index is an explicit simplification, not a claim that the source model used no commercial indices. Source configuration also specifies a random walk on log age-1 recruitment, an AR(1) F state, and Fbar ages 3–8. tinyAM's abundance and observation-error structures do not reproduce the complete SAM likelihood or covariance specification.

Normal tinyAM starting values are tried first. If numerical checks fail, at most one retry may use the report's N-at-age and F-at-age values as optimizer starting values only. The observations, settings, and estimand remain unchanged; source estimates do not constrain tinyAM parameters.

### Numerical review

The default fit failed numerical checks; the permitted source-N/F-start retry
passed them (optimizer code 0, objective 429.111, maximum gradient 0.000099,
positive-definite Hessian, and successful uncertainty estimation). Despite that
pass, its absolute population scale is not credible. In 2025 the tinyAM age-1
estimate is 30,021,906,910 thousand fish versus 3,780 thousand in the report;
the age-10+ estimates are 485,819,524,786 versus 1,214 thousand. The model
matches recent catch-at-age predictions, but fitted F is only about
6e-11–2e-8 per year and age-specific q about 1e-13–4e-12. This is a strong
practical scale-confounding signal between N, low F and freely estimated q;
the numerical convergence checks do not certify a meaningful abundance scale.

A source-closer sensitivity turned off the IID cohort-survival process while
keeping the recruitment random walk and all observations and other settings.
It also passed the numerical checks (objective 437.337, gradient 0.000116 and
positive-definite Hessian), but its median percent differences remained
extreme: about 21 billion percent for N and total abundance and 51 billion
percent for SSB. Removing IID N deviations therefore does not resolve the
scale problem. Neither fit should be presented as a close assessment
replication or used to interpret tinyAM absolute N, F, q or derived biomass.
The current translation is retained only as a diagnostic experiment; the
five commercial biomass indices and their model likelihoods remain
unrepresented, and no unsupported q prior or fixed scaling was introduced.
