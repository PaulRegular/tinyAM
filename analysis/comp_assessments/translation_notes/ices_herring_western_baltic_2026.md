# Western Baltic spring-spawning herring (her.27.20-24), 2026: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_herring_western_baltic_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_herring_western_baltic_2026.md).

## Expanded-library review

### Translation revision (2026-10-13)

The tinyAM translation was revised to address an inconsistency between the age structure of the official SAM population model and its survey indices.

The official assessment models ages 0–8+, whereas the HERAS index has a plus group at age 6 and the IBTS/BITS Q1 index has a plus group at age 5. The earlier tinyAM translation retained the population age range of 0–8+ but treated these survey plus-group observations as if they represented individual ages. This was inconsistent with the observation definitions in the official SAM model.

The revised translation fits ages **0–5+**, making the population plus group coincide with the earliest survey plus group. The translation helper aggregates catch and survey observations above the selected maximum age.

Specifically:

- **HERAS:** The age-5 index and age-6+ index are summed into a single age-5+ observation.
- **IBTS/BITS Q1:** The existing age-5+ observation is retained without additional age aggregation.
- **GERAS and N20:** Their observed age ranges are unchanged.
- **Commercial catch:** Catch numbers for ages 5–8+ are summed into the new age-5+ group.

This approximation preserves the aggregate observed quantities while avoiding the incorrect interpretation of plus-group survey indices as observations of individual ages.

The population and fishing mortality processes were also revised. The tinyAM N process now uses `iid` rather than `rw`, with `init = "free"`. The F process now uses `cor_rw` instead of independent random walks, allowing correlated temporal F changes across ages.

These changes bring the translation closer to the source SAM process structure, although the official hockey-stick recruitment relationship, complete F-state configuration, and correlated observation errors are not reproduced exactly.

#### Consequences of the reduced age structure

The revised plus group represents ages 5 and older, whereas the official population model estimates ages 5–7 separately and retains age 8 as its plus group.

The aggregation is exact for additive catch numbers and observed survey-index totals, subject to the source index definitions and unit consistency. It is not an exact reduction of the official population model.

In particular, the translation helper does not automatically calculate population-weighted stock weights, maturity, natural mortality, or fishing mortality for the new 5+ group. The translated model uses the source biological inputs at age 5 as an approximation for the combined older ages.

Combining survey observations also changes their statistical properties. The sum of correlated lognormal observations does not generally follow the same likelihood as its separate components. The revised translation therefore preserves the aggregate index values but not the official model's complete survey observation-error structure.

Comparisons involving abundance, SSB, and Fbar require particular care. Official age-specific abundance should be aggregated consistently to 5+ before comparison. The source age-5 F estimate is not directly equivalent to the effective F for the translated 5+ group, and the translated Fbar over ages 2–5 includes an aggregate terminal age rather than the official model's individual age-5 state.

The revised translation is intended as a simpler, internally consistent approximation of the official assessment, not an exact reproduction of its age structure or likelihood.

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_herring_western_baltic_2026.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

Remaining differences include the official hockey-stick recruitment relationship, cross-age survey observation correlations, the full age-specific F-state structure, and the use of source age-5 biology for the translated 5+ population group.

Inspect the current diagnostics and comparison definitions before interpreting this translation.
