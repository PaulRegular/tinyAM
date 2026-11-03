# North Sea Norway pout: 2026 source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_norway_pout_north_sea.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_norway_pout_north_sea.md).

## Expanded-library review

### tinyAM translation revision (2026-10-14)

The revised translation retains the annualized observations and biological inputs while improving correspondence with the official population and fishing-mortality processes.

The principal changes are:

- Free initial abundance replaces exponential initialization.
- IID older-cohort survival deviations use age-grouped process SDs.
- Correlated F random walks replace stationary AR1 F.
- Catch observation SDs retain the translated native age groups.
- Survey catchability and observation SDs retain the source-derived sharing groups.

The revised configuration is reported to produce a working fit. Updated tinyAM convergence statistics and comparison metrics should be recorded separately from the official SESAM diagnostics.

#### Model period and age structure

The translation fits 1984–2024 and ages 0–3+, with age 3 as the plus group.

The model annualizes quarterly population dynamics and catch observations.

Recruitment is represented at age 0 in the annual model, while the source's Q3 recruitment timing is accounted for where possible in translating mortality exposure and survey observations.

This annualization is a substantive structural approximation.

#### Population process

The revised population settings are:

`N_settings = list(process = "iid", init = "free", sd_form = ~0 + N_sd_key)`

The model estimates initial older-age abundance freely.

Older-cohort survival deviations follow IID processes, with process SDs determined by the `N_sd_key` age groups.

The source's principal older-age variance distinction is represented by a shared SD for ages 1–2 and a separate SD for age 3+.

Recruitment at age 0 retains the separate tinyAM recruitment process. The source's age-0 recruitment/survival variance cannot be combined exactly with the older-age N-process variance under this parameterization.

Consequently, the age-specific variance structure is approximated rather than reproduced exactly.

The official quarterly process is represented through annual cohort survival.

#### Fishing-mortality process

The revised translation specifies:

`F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 1:2)`

Fishing mortality follows correlated random walks on the log scale.

For each age, annual increments are:

`delta_log_F[a,t] = log(F[a,t]) - log(F[a,t-1])`

The annual increment vector follows a multivariate Normal distribution with covariance:

`Sigma[a,b] = sigma_F^2 * rho^abs(a-b)`

Here, `rho` describes AR(1) dependence between neighbouring ages, and `sigma_F` is the common annual increment SD.

This reproduces the main temporal and cross-age features of the official SESAM F process more closely than the previous stationary AR1 model.

The official model's quarterly increments, separate age-grouped F innovation variances and large-increment adjustment are not reproduced.

The annual translation uses one shared F-process SD rather than the source's two variance groups.

The initial F state remains subject to tinyAM's random-walk initial-state convention.

Fbar is calculated over ages 1–2.

#### Natural mortality

The translation uses:

`M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)`

The original quarter-specific natural mortality values are converted to annual-equivalent mortality inputs.

No natural-mortality process is estimated.

For ages 1–3+, annual exposure includes all four quarters.

For age 0, annual exposure includes only Q3 and Q4 because recruitment enters in Q3.

This accounts for the source's differing age-specific exposure periods without reproducing the entire quarterly population model.

#### Catch observations

Quarterly catch numbers are summed into annual catch-at-age observations.

The translation retains the three source-derived error groups using:

`catch_settings = list(sd_form = ~0 + sd_block)`

The groups distinguish age 0, ages 1–2 and age 3+.

Observation SDs are estimated from the annualized data while preserving the source's age-group sharing.

The official quarterly catch likelihood is not reproduced: aggregating four quarterly observations changes both their temporal resolution and their error structure.

The translation also avoids treating incomplete 2025 catch as a complete annual observation.

#### Survey indices

The translation uses:

`index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)`

Catchability parameters are estimated according to source-derived `q_key` groups.

Observation-error SDs are estimated using the translated `sd_block` groups.

These groupings preserve the source's main catchability and variance-sharing patterns without requiring a separate catchability coefficient for every age-fleet combination.

Survey timing is represented by quarter midpoints:

- Q1: 0.125.
- Q3: 0.625.

The model retains individual-age survey observations on their native scales.

The grouped IBTS other-countries Q3 age-2–3+ observation remains excluded from the tinyAM likelihood but is retained in the canonical record.

No missing 2025 IBTS other-countries observations are inferred.

#### Biological inputs

The translation retains the official maturity schedule.

Annualized stock weights use:

- Q1 stock weights for ages 1–3+.
- Q3 stock weight for age 0.

This aligns the weights with the selected annual population-state convention while accounting for the later entry of age-0 recruits.

Natural mortality is converted using age-specific quarterly exposure.

The retained biological inputs are therefore source-based annual equivalents, not alternative estimated biological surfaces.

#### Spawning-stock biomass comparison

The official SESAM assessment calculates SSB at the beginning of each quarter.

tinyAM's common-period comparison uses Q1 spawning biomass calculated from beginning-of-year abundance, the translated Q1 stock weights and the official maturity vector.

For ages contributing to SSB, these definitions align at the start of the year.

Age-0 maturity is zero, so its different within-year recruitment timing does not directly contribute to SSB.

The common-period Q1 SSB comparison is therefore more directly comparable than other quarterly-to-annual outputs.

#### Abundance and recruitment comparisons

The official assessment's Q1 age-0 abundance is measured before current-year age-0 recruitment enters the population in Q3.

Consequently, that source age-0 abundance state does not correspond directly to tinyAM's annual recruitment state.

The official Q1 age-0 N row is excluded from the like-for-like age-specific abundance comparison.

Age-0 recruitment is compared separately using the source's recruitment output.

This avoids treating biologically different states as equivalent.

#### Starting parameters

The revised fit retains:

`start_par = start_par`

These starting values assist numerical optimization of the annualized model.

They do not fix the population or fishing-mortality states to the official SESAM estimates.

Any official N or F values supplied as initial parameter values are re-estimated freely during the tinyAM fit.

Retaining informed starting values is reasonable for this translation, given the substantial differences between quarterly SESAM dynamics and the annual tinyAM parameterization.

#### Interpretation of the revision

The revised model is a closer approximation to official SESAM in two important respects.

First, age-grouped IID survival-process SDs allow different survival variability between younger ages and the plus group.

Second, correlated F random walks reproduce the main random-walk and cross-age correlation structure used by SESAM.

These changes improve structural correspondence without reproducing the source model's quarterly dynamics, age-0 recruitment timing, grouped-age survey likelihood or special large-increment treatment.

The model should therefore be assessed primarily as an annual approximation of population dynamics and stock-status trajectories rather than as an exact reproduction of the quarterly likelihood.

Numerical success should be established using the current tinyAM convergence checks and comparison results, independently of the official benchmark fit.
