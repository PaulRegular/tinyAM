# Atlantic striped bass: tinyAM translation notes

These inherited notes concern the separate 2018 assessment, not the current 2024 source record. No fit was rerun.

For analysis only, the two catch fleets are each converted from totals times
age proportions and then pooled. Aggregate relative indices are combined
with their available age compositions to create age-specific index values;
composition proportions are normalized across available ages, and zero cells
are omitted from tinyAM's lognormal likelihood. These transformed observations
do not replace canonical totals or compositions. tinyAM uses a simpler
age-specific lognormal likelihood and does not reproduce the source's
multinomial compositions, time-varying CV weights, or fleet-specific
selectivity. Female maturity is multiplied by the fixed female proportion for
tinyAM's single maturity surface. Source SSB and derived January-1 biomass
remain distinct comparisons; derived biomass has no uncertainty.

The 1982–2017 tinyAM fit passed the numerical checks after the one permitted
retry with source-based N/F starting values; the default-start attempt did not
pass. The successful fit had optimizer code 0, maximum absolute gradient
0.00283, a positive-definite Hessian and successful SD reporting (168 fixed
and 575 random effects). The main disagreement is fishing mortality: its
mean absolute percent difference is 157%, terminal-year difference is +49%,
and trend correlation is 0.775. By comparison, total abundance has 12.0% mean
absolute percent difference and 0.991 trend correlation; recruitment has
11.5% and 0.986. SSB differs by 36.7% on average with a -27.4% terminal-year
difference, while its trend correlation is 0.986; the source and tinyAM SSB
definitions remain approximate. Fixed M matches because it is supplied as an
input, not estimated. See the local dashboard for the year- and age-specific
comparisons.

[Official source review](https://github.com/PaulRegular/saio/blob/main/source_reviews/ck_nefsc_garmiii_saxatilis_gom_chatt.md).

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ck_nefsc_garmiii_saxatilis_gom_chatt.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

The published index tables provide numerical values for 1982–2017. The model
text also describes Maryland YOY and age-1 indices in 1970–1981, but those
values are not tabulated in the available numerical appendix. They are not
invented or used in the translation. Appendix B10 links YOY indices to age 1 and age-1 indices to age 2. Its “VA Age 1” label conflicts with the “NY Age-1” label in Table B5.6 and
assessment diagnostics. The translation retains the table label, links the
series to age 2, and records the discrepancy.
