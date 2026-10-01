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
or supersede the Charbonneau–Keith database.

## What is recorded

Five ordinary CSV tables live in `database/`:

- `stocks.csv`: one row per management stock, with a stable local ID and the
  separate identifier used by the assessment agency.
- `assessments.csv`: the production run investigated, its data end year,
  model-defining framework, source links, and collection status.
- `assumptions.csv`: source assessment settings in long form.
- `inputs.csv`: original catch-at-age, total landings, survey, biological, and
  natural mortality inputs. Aggregate landings have a blank age and are retained
  separately from catch-at-age observations.
- `outputs.csv`: reported population quantities and uncertainty where available.

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
assumption and value also identifies the table, page, model field, or file from
which it came. Unknown information is recorded as `unknown` with an explanation,
not guessed.

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
`database_to_tiny_obs(assessment_id, inputs)` to reshape one assessment's
canonical catch-at-age, survey, and biological input rows into tinyAM's `catch`,
`index`, `weight`, and `maturity` tables. Aggregate landings are retained in the
database but are not used as catch-at-age observations. The helper preserves reported values, units, survey names, and observation
timing. It does not choose model settings or infer population processes.
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

This first database pass documents source assessments only. It does not fit
tinyAM across the stock set, guarantee that every reported quantity can be
reproduced, or resolve differences among assessment methods. Input and output
availability varies by agency and stock; gaps remain explicit in the tables.
