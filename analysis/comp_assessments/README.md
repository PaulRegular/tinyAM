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

## Notes from the cod extraction review

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
index values. Table 9 and Table 10 now both contribute 1,065 age-weight rows
for 1954-2024, and Table 8 contributes 1,065 female maturity values indexed by
birth cohort. Reported catch counts remain counts in thousand fish; no
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
For Southern Gulf cod, the 2024 Science Advisory Report is the latest accepted
assessment, using the SCA model through 2023. DFO identifies the 2019 run as the
last full assessment and says the same population model was used again in 2024.
The 2024 database record is current but partial: it contains the reported 2023
SSB estimate and does not borrow observation series or estimated outputs from
the 2019 run. Structural settings are cross-referenced to the 2019 detailed
report only where the 2024 summary confirms the same population model was used.
The detailed 2019 SCA assessment remains as a historical record under the 2012
framework. It contains 54 annual stock-catch values for 1965-2018
(tonnes), plus 480 source-reported landed numbers-at-age values for 1971-2018
(ages 3-12+, in thousand fish). These age counts do not recover the full fitted
catch-composition series, which the model describes as ages 2-12+. It also
retains 759 annual maturity rows carried forward between source-listed change
years. Its survey records keep the model's aggregate index separate from its
age-composition input. Both are explicitly marked reconstructions from
published age-specific source tables, and no derived age-specific abundance
series is entered. The reconstructed RV aggregate biomass series has 46 years because
weights-at-age are not tabulated for 1980 and 1985; the mobile-sentinel series
covers 2003-2018; and the longline index covers 1995-2017. These inputs remain
partial because the report does not publish all original composition samples,
all RV biomass-index years, or every input over the model's 1950-2018 span. The
recorded population outputs are maximum-likelihood estimates from Tables
21-23; the report's other population summaries are generally posterior
medians. The report also states terminal estimates for M at ages 5-8 and 9+
and fully recruited q for the RV and mobile-sentinel surveys; these are stored
as grouped or time-invariant outputs with their source identified as report
text.

### Northern Shelf cod (2025)

The North Sea cod catalogue entry is represented by the accepted three-substock
Northern Shelf assessment after the 2023 benchmark. Original biological
observations remain separate from fitted biological surfaces. All seven index
streams retain their supplied log-scale SDs; model timing is 0.125 for Q1, 0.75
for Q3+Q4 and 0 for forward-shifted recruitment indices. Fitted inputs extend
through the 2025 Q1 survey, while catches and reported F end in 2024.

See [the source review](source_reviews/ices_cod_north_sea.md) for coverage, pinned
workflow revisions and remaining gaps. Cached sources remain gitignored. The
source-specific importer requires Python with pdfplumber and base R:

```sh
python analysis/comp_assessments/scripts/004_import_north_sea_cod.py . --rscript Rscript
Rscript analysis/comp_assessments/scripts/002_validate_database.R
Rscript analysis/comp_assessments/scripts/005_validate_north_sea_cod.R
```

The importer checks all 2,688 native N/F/M values against the detailed report
and can be rerun without changing other assessments. Completeness statuses
remain partial; this record does not borrow observations from older runs.
