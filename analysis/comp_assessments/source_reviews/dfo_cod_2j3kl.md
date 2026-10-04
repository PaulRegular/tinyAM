# Northern cod (2J3KL) source review

For Northern cod, the database represents the 2025 production assessment under
the 2023 framework. DFO's detailed research document for that run was published
in April 2026; its model inputs extend through 2024 and selected estimates are
reported through 2025. DFO later published a 2026 assessment summary to 2026
(Science Advisory Report 2026/030), but it does not provide the full numerical
input and model-output series needed to link those data to that run. Following
the selection rule above, the 2026 summary is not entered as a separate
assessment and none of its values are mixed into the 2025 record.

The 2025 record contains 819 commercial catch numbers-at-age values for ages
2-14 in 1962-2024, plus all 71 reported 2J3KL landings values from 1954-2024.
The model uses landings as a bounded catch input and converts catch-age counts
to proportions for its composition likelihood. The fall RV data include 507
age-specific values and 40 total biomass index values used as the lagged cod-biomass
covariate for M; the total-abundance series in Table 6 is not included because it
is not an xteNCAM input. 2004 and 2021 are excluded from the age series, and
2022 was not surveyed.
The database also includes 239 age-specific sentinel values, 13 Smith Sound
biomass estimates and 119 sampled age counts, and 88 Fleming/Newman juvenile
index values. Table 9 and Table 10 each contribute 1,065 age-weight rows for
1954-2024. Table 8 reports 1,065 female maturity-at-age estimates as
calendar-year-by-age values for 1954-2024. The maturity model used a cohort
effect, but that does not change the table year-age indexing. Reported catch
counts remain counts in thousand fish; no
proportion-to-count conversion was needed. Fall timing is represented by a
0.75 season-level approximation; Smith Sound timing is derived from reported
months, while sentinel and juvenile survey timing remains unknown.

The assessment remains partial. The spring Capelin acoustic series is shown in
the detailed 2025 Capelin report but not tabulated; its exact spatial match to the
model's 3L covariate is unresolved. Detailed tagging observations and
reporting-rate likelihood data are not represented. The report also does not
explain how Smith Sound sample ages 15-16 map to the model's terminal age 14. Age-specific population and mortality
outputs are shown graphically but are not available as numerical tables. Table 16
does report fitted F and natural-mortality process parameters; the F correlations
and variance and the M-process correlations, variance, baseline M, and Capelin
effect are now described in assumptions.csv. M is estimated in the model, not
supplied as an input.
