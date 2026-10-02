# Comparative age structured assessments

This directory holds a small, inspectable record of production age structured
fish assessments. The aim is to see how widely the assessment approaches used
for real stocks differ, and later to study changing population processes with
consistent models.

The starting list comes from Julie Charbonneau and David Keith's [Age structured
marine fish database](https://github.com/JulieCharbonneau/Age-structured-marine-fish-database),
especially its [PLOS 2026 analysis data](https://github.com/JulieCharbonneau/Age-structured-marine-fish-database/tree/main/PLOS_2026_analysis).
That curated work helps identify candidate stocks and provides historical
metadata and values to cross-check. It is not a substitute for the original
assessment data: some entries are reconstructed, digitized, transformed, or
derived from model outputs. We prefer the latest accepted assessment and its
original data and model documentation. The curated database warns that its
quality checks cover stocks used in its analyses; consult the stock-specific
source material before relying on an extracted value. This work does not replace
or supersede the Charbonneau–Keith database. Downloaded assessment PDFs and useful public data/model files are kept locally by stock under the gitignored `source_cache/` folder; each manifest records the source URL, retrieval date, and checksum.

## What is recorded

Five ordinary CSV tables live in `database/`:

- `stocks.csv`: one row per management stock, with a stable local ID and the
  separate identifier used by the assessment agency.
- `assessments.csv`: the production run investigated, its data end year,
  model-defining framework, source links, and collection status.
- `assumptions.csv`: source assessment settings in long form.
- `inputs.csv`: observed catch-at-age counts, survey data, biological inputs,
  and natural mortality. `catch_at_age` stores numeric values, not proportions.
  When a source reports number proportions, multiply them by the matching total
  number of fish for that year, stock, fleet/area, and period. If only total
  landed biomass is available, the product is biomass at age; convert that to
  fish numbers only with compatible age-specific weights. Cite both source
  values and record the calculation and units. Keep total landings only when
  needed for a documented conversion; detailed gear splits are not needed for
  this database. `weight` and
  `catch_weight` stay separate when a source uses distinct stock and fishery
  weights. `maturity_cohort` preserves maturity ogives indexed by birth cohort.
  `samp_time` records survey timing as a fraction of the year from 0 to 1 for
  each survey observation. Seasonal approximations are explained in row notes.
- `outputs.csv`: reported population quantities and uncertainty where available.
  `age_group` records grouped estimates such as F for ages 5-8 or M for ages 9+.

Relevant assessment PDFs and public data/model files are cached locally by
stock in the gitignored `source_cache/` folder when they can be retrieved. The
manifest records each source and any retrieval limitation; cached files are
not committed.

The database records what each source assessment did. It does not store tinyAM
settings, proposed formulas, or judgements about whether tinyAM can reproduce a
feature. `R/audit_assumptions.R` compares documented source assumptions with
current tinyAM capabilities; its results belong in `results/audits/`.

Use official agency and assessment team sources where possible, in this order:
native fitted object, native model files, reproducible assessment repository,
official machine readable files or tables, then official reports. Digitized or
reconstructed values are labelled as such. `source_type` uses `native_model`,
`official_machine_readable`, `official_table`, `digitized`, `reconstructed`,
or `charbonneau_seed`. Source locations are stored in `assessments.csv`; each
assumption and value also identifies its source location. `official_table` covers
values copied from official report tables or figures; the precise table, figure,
page, or file appears in `source_reference`. Digitized or reconstructed values
are labelled accordingly. Unknown information is recorded as `unknown` with an
explanation, not guessed.

Select the newest production assessment whose public inputs and outputs can be
linked to the same model run and recorded clearly. A newer advisory report may
provide selected updates before detailed sources are available for that run; do
not mix those values with inputs from an earlier assessment. Use the detailed
research document for the selected run's inputs and results, even when it was
published after the advisory report, and use the framework research document to
describe model assumptions. Record the assessment year, last input year, and
last reported estimate year separately.

## Initial candidates

The initial sampling frame aims for prominent production assessments across
DFO, ICES, NOAA AFSC, and NOAA NEFSC, and across NCAM, SAM, ASAP, Stock Synthesis,
and other age structured models. The ten preferred candidates are Northern cod
(2J3KL), southern Gulf cod (4T-4VN), North Sea cod, Northeast Arctic cod, North
Sea haddock, North Sea herring, eastern Bering Sea pollock, Gulf of Alaska
pollock, Gulf of Alaska Pacific cod, and Georges Bank haddock. Their PLOS 2026
curated seed records are listed by `scripts/001_seed_charbonneau.R`; the seed's
year range and model label are not treated as current production assessment
facts. A stock may be replaced only if it is absent from the curated starting
frame or has no identifiable comparable production assessment. Any substitution
will be documented here.

## Using the records

Source `R/database_to_tiny_obs.R` and call
`database_to_tiny_obs(assessment_id, inputs)` to reshape compatible `catch`,
`index`, `weight`, and `maturity` rows into tinyAM observation tables. The
helper does not yet convert source `catch_at_age` rows automatically because
those tables can differ from the assessment likelihood inputs (for example,
reported landings-at-age versus the model full catch composition). It preserves
stored values, units, survey names, and observation timing. It does not choose
model settings or infer population processes.
Natural mortality stays in `inputs.csv` for a separate, informed model setup.
If a source structure cannot be represented without combining fleets, sexes,
or seasons, the helper reports that limitation instead of silently combining
the data.

From the repository root, run
`Rscript analysis/comp_assessments/scripts/002_validate_database.R` to check
required columns, identifiers, links, duplicate rows, value types, source labels,
and assessment status. To review the selected curated records, run
`Rscript analysis/comp_assessments/scripts/001_seed_charbonneau.R path/to/metadata_no_age_corection.csv dfo_cod_2j3kl`.
The final argument selects one candidate for review; the script does not add it
to the source-assessment tables.

Run `Rscript analysis/comp_assessments/scripts/003_fit_readiness.R` to export
`results/fit_readiness.csv`. It reports each assessment's modeled years and ages,
recruitment age, plus group, input counts, surveys and times, full-grid coverage,
model-output availability, converter result, and `tinyAM::check_obs()` result.
Grid checks use the source model's full year and age range, not merely the span
of the rows that happened to be extracted. A partial source record remains
partial even when some rows can be converted.

This first database pass documents source assessments only. It does not fit
tinyAM across the stock set, guarantee that every reported quantity can be
reproduced, or resolve differences among assessment methods. Input and output
availability varies by agency and stock; gaps remain explicit in the tables.

## Notes from the cod extraction review

For Northern cod, the database represents the 2025 production assessment under
the 2023 framework. DFO's detailed research document for that run was published
in April 2026; its model inputs extend through 2024 and selected estimates are
reported through 2025. DFO later published a 2026 assessment summary to 2026
(Science Advisory Report 2026/030), but it does not provide the full numerical
input and model-output series needed to link those data to that run. Following
the selection rule above, the 2026 summary is not entered as a separate
assessment and none of its values are mixed into the 2025 record.

The 2025 record contains 819 official commercial catch numbers-at-age values
for ages 2-14 in 1962-2024 and 507 age-specific fall RV survey values. The
published catch counts are retained directly; monthly and division-by-gear
landings tables are omitted because they are not needed for this catch-at-age
record. The fall survey omits 2004 and 2021 because of coverage problems and was not
conducted in 2022; season-only timing is stored as an explicit 0.75
approximation.
Table 9 contributes 825
beginning-of-year stock weight-at-age values for ages 0-14 in 1954-2008;
values for 2009-2024 are not yet transcribed. Table 10 contributes 1,065
mid-year catch weight-at-age values for ages 0-14 in 1954-2024. Both weight
surfaces come from a cohort-based growth model, and their pre-1983 cells are
hindcasts rather than direct survey measurements. Table 8 contributes 1,065
female maturity-at-age estimates for cohorts 1954-2024. Their year field is a
birth-cohort index, so an assessment-specific cohort-to-calendar-year mapping
is needed before using them as tinyAM maturity inputs. Numerical Sentinel,
Smith Sound, juvenile-survey, Capelin, tagging, and age-specific model-output
series remain untranscribed. M is estimated within the model, not supplied.
All catch-at-age values currently stored in the database are already counts in
thousand fish; no proportion-to-count conversion was needed.

For Southern Gulf cod, the 2024 assessment to 2023 reused the SCA model adopted
in 2012 and updated through 2019. Its record includes annual outputs: SSB for
1950-2023; recruitment of fish younger than age 4; and F for ages 5-8 and M for
ages 5-8 and 9+ for 1950-2023. These output tables are transcribed from DFO
2024 rebuilding-plan materials, which cite the 2024/026 assessment. The plan
2023 SSB (11.9 kt; 95% CI 7.8-16.5) differs from the assessment report direct
value (12 kt; 10.5-21.6); the database uses the report value for 2023 and the
plan table only through 2022. The 2024 record still lacks updated catch-at-age
proportions, survey observations, weight, maturity, and M values. The 2019
record contains 480 reported landed numbers-at-age values for 1971-2018 (ages
3-12+, in thousand fish). These are source-reported counts, not the full
catch-composition inputs used by the SCA model. It also retains 759 annual rows
derived by carrying source-listed maturity ogives forward to their next stated
change year; these do not complete the 2024 assessment inputs.
