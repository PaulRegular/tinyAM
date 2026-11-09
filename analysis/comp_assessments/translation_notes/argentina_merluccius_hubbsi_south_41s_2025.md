# Southern Argentine hake, south of 41 degrees S (2025)

This is a preliminary tinyAM translation of the INIDEP 2025 production assessment. The source review and numerical records remain in SAIO. The selected tinyAM source snapshot is e7605c572f55956925d019725acfb2a229312c98.

The recipe fits 1986-2024 and ages 1-6 plus 7+. It uses the official estimated catch numbers-at-age and the winter global-survey numbers-at-age, with the 1 July population weights, maturity schedule and fixed natural mortality. The single survey is assigned one shared q and one shared observation SD. The survey timing of 0.5 is an explicit midpoint assumption because the source report gives a season, not a precise date.

The source assessment includes multiple fishing fleets, additional CPUE and juvenile indices, and a source-specific recruitment structure. The tinyAM fit uses one aggregate catch series, only the global winter age index, IID recruitment, an RW process on log F with independent increments across ages, and no extra N process. It is an approximation, not a recreation of the official likelihood.

The tinyAM RW process has an unpenalized first row and independent increments across ages. The comparison maps the source output named spawning_biomass to tinyAM's SSB panel. It compares available N-at-age, biomass, SSB and recruitment. The source report's age-specific F series is not available in the curated output table, so F-at-age remains unavailable for comparison. Aggregate CVs remain in the SAIO source records and are not converted to standard errors in this recipe.

The fitting policy is package defaults first. If numerical checks fail, the single permitted retry initializes aligned N states from official N-at-age estimates; those values do not constrain the fit. Inspect the current attempt log and dashboard before interpreting estimates.

## Fit and comparison

The fit converged on the first attempt using package defaults. The optimizer
returned code 0, the maximum absolute gradient was 2.68e-6, sdreport succeeded,
and the Hessian was positive definite. The objective was 135.29; no source-based
starting-value retry was used.

Over 1986-2024, mean absolute percent differences were 17.96% for numbers at
age, 6.25% for total abundance, 10.84% for total biomass, 6.87% for recruitment,
and 16.07% for SSB. Terminal tinyAM estimates were lower by 31.01%, 20.93%,
24.23%, 22.07%, and 28.61%, respectively. Trend correlations were 0.880 for
numbers at age, 0.940 for total abundance, 0.890 for biomass, 0.938 for
recruitment, and 0.880 for SSB. These are preliminary comparisons; the
population model, recruitment process and biomass definitions remain
approximate. Age-specific F comparison is unavailable.

Fixed M reproduces the recorded 0.3 per-year source input exactly. This is an
input-consistency check, not evidence that tinyAM estimated M or independently
replicated an M output. The recipe now takes its fixed-M settings from the
translation helper so it uses the recorded source surface rather than tinyAM's
generic 0.2 default.

The local dashboard is cached under
results/preliminary/cache/e7605c572f55956925d019725acfb2a229312c98/argentina_merluccius_hubbsi_south_41s_2025/dashboard.html.