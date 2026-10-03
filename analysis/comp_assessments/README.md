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
- `assessments.csv`: the detailed run represented, its data and estimate end
  years, model-defining framework, source links, and collection status.
- `assumptions.csv`: source assessment settings in long form.
- `inputs.csv`: catch numbers-at-age, survey data, biological inputs, and
  natural mortality. `type` identifies a broad input family, while `measure`,
  `basis`, and `unit` preserve the exact quantity and its scale.
  When a source reports number proportions, multiply them by the matching total
  number of fish for that year, stock, fleet/area, and period. If only total
  landed biomass is available, the product is biomass at age; convert that to
  fish numbers only with compatible age-specific weights. Cite both source
  values and record the calculation and units. Keep total landings when they
  are a model input or are needed for a documented conversion; detailed gear
  splits are not needed unless the fitted model uses them. `weight` and `catch_weight` stay separate when a source uses
  distinct stock and fishery weights. `year_basis = birth_cohort` preserves
  maturity ogives indexed by cohort without treating cohort as calendar year.
  `sampling_time` records survey timing as a fraction of the year from 0 to 1.
  Seasonal approximations are explained in row notes.
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
`official_machine_readable`, `official_table`, `official_document`, `digitized`,
`reconstructed_source_input`, or `charbonneau_seed`. Source locations are stored in `assessments.csv`; each
assumption and value also identifies its source location. `official_table`
covers values copied from official report tables, while `official_document`
covers values stated in report prose. The precise table, page, or file appears
in `source_reference`. Digitized or reconstructed values are labelled
accordingly. Unknown information is recorded as `unknown` with an explanation,
not guessed.

Represent the latest accepted production assessment, using the most recent
authoritative detailed source that documents that run. A detailed source may be
published after the advice or summary product, so distinguish the source's
publication date from the assessment year and data terminal year. If no detailed
source for the latest accepted assessment is available, use earlier detailed
material and the relevant framework only for context; do not substitute earlier
run inputs or outputs. Record unresolved gaps and keep completeness statuses
partial where needed. Record the assessment year, last input year, and last
reported estimate year separately.

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
`database_to_tiny_obs(assessment_id, inputs)` to reshape compatible catch,
index, weight, and calendar-year maturity rows into tinyAM observation tables.
The helper does not choose model settings, infer population processes, or map
cohort-indexed maturity onto calendar years. It preserves stored values, units,
survey names, and observation timing.
Use `R/read_committed_assessment.R` when analysis scripts must read only a
reviewed database snapshot; it returns the commit identifier alongside the
selected stock and assessment records.
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
of the rows that happened to be extracted. Missing catch cells are allowed by
tinyAM and are reported as source coverage gaps rather than treated as a
conversion failure. A partial source record remains partial even when some rows
can be converted.

This first database pass documents source assessments only. It does not fit
tinyAM across the stock set, guarantee that every reported quantity can be
reproduced, or resolve differences among assessment methods. Input and output
availability varies by agency and stock; gaps remain explicit in the tables.

## Stock source reviews

Stock-specific source checks, coverage notes and unresolved questions are kept
in `source_reviews/`, with one Markdown file per stock.

- [Northern cod (2J3KL)](source_reviews/dfo_cod_2j3kl.md)
- [Southern Gulf cod (4T–4VN)](source_reviews/dfo_cod_4t4vn.md)
- [North Sea / Northern Shelf cod](source_reviews/ices_cod_north_sea.md)
- [Northeast Arctic cod](source_reviews/ices_cod_northeast_arctic.md)
- [North Sea haddock](source_reviews/ices_haddock_north_sea.md)
- [North Sea herring](source_reviews/ices_herring_north_sea.md)
- [Eastern Bering Sea pollock](source_reviews/afsc_pollock_ebs.md)
- [Gulf of Alaska pollock](source_reviews/afsc_pollock_goa.md)
- [Gulf of Alaska Pacific cod](source_reviews/afsc_cod_goa.md)
- [Georges Bank haddock](source_reviews/nefsc_haddock_georges_bank.md)

