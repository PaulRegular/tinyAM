# British Columbia Pacific ocean perch — translation notes

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/dfo_pop_bc.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

**Current comparison status: blocked (2026-10-17).** A complete matching
historical population comparison has not been recovered. A subsequent search
found regional biomass draws for 3CD and 5DE, but not the matching 5ABC series
(see the follow-up below).
The source record stays in the database; it is excluded from comparison
tables and figures. Earlier fitting results below are historical experiments.

A definition issue needs confirmation from final native posterior files or
the assessment authors. The public R21v3 control specifies fecundity option 2
with coefficients 1/0 and first mature age 9.49256964. Matching SS3 reads the
latter as integer 9 and calculates reproductive output separately from mature
biomass. Its code explicitly stores both `SSB_yr` and `SSB_B_yr`. The public
reporting code also handles `SpawnBio`, while the paper tables label B in
tonnes. Without the final posterior extraction, do not silently equate those
quantities or "repair" the official assessment. Published values are retained
as reported and SSB comparison is marked non-equivalent.
