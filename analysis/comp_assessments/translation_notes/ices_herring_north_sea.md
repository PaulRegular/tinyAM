# North Sea autumn-spawning herring: 2026 source review: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_herring_north_sea.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_herring_north_sea.md).

## Expanded-library review

### tinyAM translation revision (2026-10-14)

The tinyAM translation was revised to more closely approximate the official assessment's population and fishing-mortality process structures.

The previous translation used random-walk cohort-survival deviations and independent age-specific F random walks. The revised translation uses independent cohort-survival deviations and correlated F random walks.

The translated model retains the same fitted period, population ages, biological inputs and selected survey streams.

#### Population process

The revised model specifies:

`N_settings = list(process = "iid", init = "free")`

Older-cohort survival deviations are independent rather than following a random walk.

Recruitment retains the separate random-walk process specified by the default tinyAM recruitment settings. Changing `N_settings$process` to `iid` does not change recruitment to IID.

Initial abundance for older ages is freely estimated.

The revised model uses one shared SD for older-cohort survival deviations. This broadly matches the separation between recruitment and older-age survival variation in the official assessment, but does not reproduce every aspect of FLSAM's state initialization or process specification.

The official N and F estimates are used as optimizer starting values only. They do not constrain the fitted tinyAM states.

#### Fishing-mortality process

The revised model specifies:

`F_settings = list(process = "cor_rw", mu_form = NULL, mean_ages = 2:6)`

This replaces independent F random walks with correlated random walks.

The `cor_rw` process estimates annual changes in log fishing mortality with AR(1) correlation across neighbouring ages. It therefore captures a central feature of the official FLSAM fishing-mortality process that was missing from the previous translation.

The distinction is important: correlation applies to annual F increments across ages, not simply to the absolute F levels.

The revised tinyAM model remains simpler than the official FLSAM configuration in two respects:

- It uses one shared F-process SD rather than the source's three innovation-variance groups for ages 0–1, 2–5 and 6–8.
- It retains separate F states at ages 7 and 8 rather than imposing the source's shared oldest-age F state.

Fbar is calculated over ages 2–6, consistent with the official assessment.

The revised specification improves structural agreement without reproducing the complete FLSAM covariance and age-sharing configuration.

#### Natural mortality

The translation specifies:

`M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)`

The translated numerical M surface is treated as fixed. The fitting period ends in 2025, when the published biological inputs remain available.

The official SMS-2023 surface and additional 0.02 adjustment must be represented consistently with the imported numerical values. No separate mortality process is estimated in tinyAM.

#### Commercial catch

The translation retains direct commercial catch-at-age observations for the common fitted period, retaining counts in thousand fish.

The catch observation model uses:

`catch_settings = list(sd_form = ~0 + sd_block)`

This permits separate estimated observation SDs for the translated catch groups.

The fishery closure in 1978–1979 is represented by missing catch observations, consistent with the official input construction.

Positive catch observations are fitted using a lognormal likelihood. Exact zeros are not fitted as positive lognormal observations.

The translation does not claim to reproduce every feature of the official catch observation model.

#### Survey indices

The tinyAM translation retains HERAS, IBTS-Q1, IBTS0 and IBTS-Q3.

The four larval abundance indices are omitted because their partial spawning-component likelihoods are not equivalent to ordinary abundance-at-age indices.

The revised survey model specifies:

`index_settings = list(q_form = ~0 + q_key, sd_form = ~0 + sd_block)`

This retains the translated survey-age catchability-sharing groups and observation-SD groups.

Catchability is freely estimated within those groups. No monotone or logistic age constraint is imposed.

The source's cross-age observation-error correlation for IBTS-Q3 is not represented.

Survey observations retain their published numerical scales. Because their calibrated units have not been independently established, estimated catchability absorbs differences in index scale.

Sampling times remain documented approximations based on survey seasons, because the official `fleet.txt` was unavailable.

#### Spawning-stock biomass

The official assessment reports spawning-stock biomass under its native spawning-time convention.

tinyAM calculates SSB using the translated stock weights, maturity and fitted abundance under tinyAM's spawning-time convention.

The official spawning fractions were not recovered. Consequently, native SSB and translated SSB definitions may not be exactly equivalent.

Comparisons should retain this caveat even when the estimated trajectories agree closely.

#### Interpretation of the revision

The revision makes the two principal population processes more consistent with the official FLSAM structure without adding stock-specific process complexity.

The improvement is attributable to scientifically motivated changes in the survival and fishing-mortality process assumptions, rather than changes to the official assessment data.

The translation remains an approximation. In particular, the four omitted larval indices, approximate survey timing, incomplete F age-sharing and variance groups, and differences in observation-error covariance may influence the fitted population states.

The revised model's numerical convergence and comparison statistics should be evaluated separately from the previous random-walk N and independent F random-walk configuration.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_herring_north_sea.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The revised translation is closer to the official FLSAM process structure, but it is not an exact reproduction.

The main remaining differences are:

- **Observation coverage:** Four regional larval indices and their partial spawning-component likelihoods are omitted.
- **Fishing mortality:** Correlated random-walk increments are represented, but the official age-specific variance sharing and common F state at ages 7–8 are not.
- **Observation errors:** The official IBTS-Q3 cross-age correlation is omitted.
- **Survey timing:** Sampling fractions for the four age-structured indices are approximations, not independently recovered native values.
- **SSB definition:** Official spawning-time fractions remain unavailable.
- **Native assessment recovery:** Complete current-run FLSAM inputs and fitted objects have not been recovered for a direct likelihood-level comparison.

These are documented limitations of the simplified translation, not errors in the official assessment.

These settings more closely approximate the official assessment's core population dynamics than the previous translation.

The case tests whether a comparatively simple age-structured state-space model can recover broad official stock dynamics without reproducing the complete FLSAM observation likelihood, state-sharing configuration or biological uncertainty.
