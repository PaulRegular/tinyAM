# ices cod irish sea: tinyAM translation notes

Inherited settings and diagnostics; no models were rerun. Check the registered recipe and review tier for current settings.

## Outputs and readiness

Table 6.14 reports annual SSB, age-0 recruitment and Fbar ages 2–4 with pointwise 95% confidence limits. The reported summary values were cross-checked against the official ICES SAG export for key 23721, within report rounding. The source report and SAG export do not provide recoverable age-specific N-at-age or F-at-age; M is fixed rather than estimated.

The source table contains one missing cell in the weight surface (2019 age 0). No authoritative value was found, so it remains missing in the database. The translation fits ages 0–6+ over 1968–2018, where weight and maturity are complete. It retains Q1 ages 1–4, Q4 age 0 and UK FSP Q1 ages 2–6+, with timing approximated at 0.125, 0.875 and 0.125. Catch years 2003–2005 are excluded as in the assessment. The 2020 weights were marked by the report as unavailable for model use; the fitted period ends before them.

Index scale multipliers and exact likelihood allocations are unresolved. The native Q1 CV and FSP effective sample size are not used because their current likelihood roles were not recovered; tinyAM estimates a separate index observation SD by survey. Catch-at-age remains in thousand fish. The official Beverton-Holt recruitment blocks are represented by a random-walk recruitment process, and the source selectivity process is approximated with an AR1 temporal F process and a shared mean for ages 3+.

The default-initialization fit converged (optimizer code 0, maximum gradient 0.000364, positive-definite Hessian). Over 1968–2018, approximate SSB agreement was mean absolute percent difference 36.5%, terminal difference -9.3%, and trend correlation 0.971. Age-0 recruitment agreement was weaker: mean absolute percent difference 61.0%, terminal difference +1568%, and trend correlation 0.731. The simpler recruitment process remains an important limitation; no starting values were supplied.

The comparison period does not include the assessment's terminal 2019–2025 estimates. Official aggregate SSB timing is not stated, while tinyAM calculates beginning-of-year mature biomass, so SSB is approximate. Recruitment is age-aligned but its published point-estimate convention is unspecified, so it is also approximate. Official age-specific N and F are unavailable. Official Fbar is ages 2–4, while tinyAM Fbar is abundance-weighted and the required official N-at-age is unavailable; therefore the Fbar values are not compared. M is fixed at the published input values and shown as an input rather than a fitted output. The maintained translation and structural checks are in [preliminary recipe](../scripts/translation/preliminary/ices_cod_irish_sea_2026.R) and [translation check](../tests/test_irish_sea_cod_translation.R); source values are validated by [SAIO source check](https://github.com/PaulRegular/saio/blob/main/tests/test_irish_sea_cod_source.R).

[Official source review](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_cod_irish_sea.md).
