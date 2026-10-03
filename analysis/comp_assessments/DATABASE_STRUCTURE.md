# Assessment Database Structure

## Purpose

This document defines the schema for the stock-assessment database in `analysis/comp_assessments/database/`.

The database preserves, in a tidy and auditable form, the inputs, assumptions, and outputs of accepted age-structured stock assessments.

The canonical database contains five CSV files:

``` text
analysis/comp_assessments/
└── database/
    ├── stocks.csv
    ├── assessments.csv
    ├── assumptions.csv
    ├── inputs.csv
    └── outputs.csv
```

All canonical tables should remain ordinary CSV files so they can be inspected directly in R, Excel, or a text editor.

This document defines source records, not the tinyAM analysis pipeline.
Use `PROTOCOL.md` to obtain, verify, or repair those records, and
`TINYAM_TRANSLATION.md` to convert them into observations, fit models, and
compare results. A source-data gap discovered during translation may require
a validated database revision; downstream model choices remain separate.

## Core principles

1.  **Represent the accepted assessment faithfully.** Preserve the dimensions and meaning of the accepted model, including fleets, surveys, sexes, regions, and seasons where applicable.
2.  **Preserve provenance.** Material data and assumptions should be traceable to an authoritative file, object, table, document, or documented reconstruction.
3.  **Keep the canonical database source-focused.** Tidy reshaping is allowed, but downstream-analysis transformations and settings do not belong in the canonical tables.

------------------------------------------------------------------------

# Relationships

``` text
stocks.csv
    1
    |
    +----< assessments.csv
               1
               |
               +----< assumptions.csv
               |
               +----< inputs.csv
               |
               +----< outputs.csv
```

`stock_id` identifies the biological or management stock.

`assessment_id` identifies a particular assessment event or accepted assessment run.

------------------------------------------------------------------------

# 1. `stocks.csv`

One row per biological or management stock.

## Required columns

| Column | Description |
|------------------------------------|------------------------------------|
| `stock_id` | Stable internal identifier for the stock |
| `charbonneau_id` | Identifier in the Charbonneau–Keith database, when applicable |
| `authority` | Assessment authority, e.g. DFO, ICES, NOAA-AFSC |
| `authority_stock_id` | Current official stock identifier/code |
| `scientific_name` | Scientific name |
| `common_name` | Common name |
| `area` | Assessment or management area |
| `region` | Broader geographic region |
| `ocean` | Ocean basin |
| `notes` | Stock-identity, boundary, split/merge, or other notes |

## Identifier rules

`stock_id` should be readable, unique, stable through time, and independent of assessment year.

Examples:

``` text
dfo_cod_2j3kl
ices_cod_north_sea
afsc_pollock_ebs
```

Do not embed an assessment year in `stock_id`.

------------------------------------------------------------------------

# 2. `assessments.csv`

One row per assessment event or model run represented in the database.

## Required columns

| Column | Description |
|------------------------------------|------------------------------------|
| `assessment_id` | Stable identifier for this assessment |
| `stock_id` | Foreign key to `stocks.csv` |
| `assessment_year` | Year of the assessment/advice |
| `terminal_year` | Last year of data fitted by the model |
| `estimate_terminal_year` | Last year with a population estimate reported; may be later than `terminal_year` |
| `assessment_type` | Annual, update, benchmark, framework, etc. |
| `model_family` | SAM, Stock Synthesis, WHAM, NCAM, VPA/ADAPT, etc. |
| `model_version` | Model/package/version identifier where known |
| `is_current` | Whether this is the current represented assessment for the stock |
| `is_applied` | Whether this is the accepted assessment model used for advice |
| `framework_year` | Benchmark/framework year defining the current model, if distinct |
| `assessment_url` | Main authoritative assessment/report URL |
| `framework_url` | Benchmark/framework/methodology URL |
| `data_url` | Official input-data URL where distinct |
| `model_url` | Native model object/run URL where distinct |
| `repository_url` | Reproducible assessment repository where available |
| `assumptions_status` | Completeness of assumptions extraction |
| `inputs_status` | Completeness of accepted-model input extraction |
| `outputs_status` | Completeness of accepted-model output extraction |
| `notes` | Assessment-level caveats or provenance notes |

## Status values

Use:

``` text
not_started
partial
complete
not_applicable
```

Use `not_applicable` sparingly.

Status fields describe completeness relative to the accepted assessment, not relative to a downstream analysis.

Keep translation readiness, fitting failures, convergence results, and
comparison completion in analysis diagnostics. A partial source record can
support a documented limited analysis, and a complete source record can still
require modeling approximations. Neither case changes what these status
fields mean.

------------------------------------------------------------------------

# 3. `assumptions.csv`

Long-format description of biological and statistical assumptions.

## Required columns

| Column | Description |
|------------------------------------|------------------------------------|
| `assessment_id` | Foreign key |
| `component` | Population, recruitment, N, F, M, catch, index, q, biology, etc. |
| `fleet` | Fleet to which the assumption applies |
| `survey` | Survey/index to which the assumption applies |
| `sex` | Sex to which the assumption applies |
| `region` | Region to which the assumption applies |
| `season` | Season to which the assumption applies |
| `setting` | Name of the assumption |
| `value` | Value or concise description |
| `source_reference` | Specific file, object field, table, page, or section |
| `notes` | Clarification, interpretation, or unresolved ambiguity |

Blank `fleet`, `survey`, `sex`, `region`, or `season` fields mean the assumption is not specific to that dimension.

If an important assumption cannot be resolved, use:

``` text
value = unknown
```

and explain the unresolved issue in `notes`.

------------------------------------------------------------------------

# 4. `inputs.csv`

Long-format representation of numerical inputs supplied to the accepted assessment.

## Required columns

| Column | Description |
|------------------------------------|------------------------------------|
| `assessment_id` | Foreign key |
| `type` | Broad input family |
| `measure` | Exact quantity represented |
| `basis` | Numbers, biomass, number proportion, biomass proportion, etc. |
| `fleet` | Fleet identity |
| `survey` | Survey/index identity |
| `sex` | Sex |
| `region` | Spatial region used by the fitted model |
| `season` | Season |
| `year` | Calendar/model year, or cohort identifier when year_basis says so; blank for time-invariant inputs |
| `year_basis` | Optional meaning of `year`; use `calendar_year` or `birth_cohort` |
| `age` | Age |
| `value` | Numerical value |
| `unit` | Explicit unit |
| `sampling_time` | Timing within year as a fraction from 0 to 1 when part of the source/model definition |
| `source_type` | Provenance category |
| `source_reference` | Specific file/table/object/page |
| `transformation` | Transformation used only to recover the accepted-model input |
| `notes` | Additional clarification |

## Initial `type` vocabulary

``` text
catch
index
weight
catch_weight
maturity
M
```

Add new types only when a real assessment requires them.

## Initial `measure` vocabulary

### Catch/removals and indices

``` text
numbers_at_age
biomass_at_age
total_numbers
total_biomass
proportion_at_age
```

### Biology

``` text
weight_at_age
maturity_at_age
natural_mortality_at_age
landings_numbers_at_age
landings_fraction_at_age
landings_weight_at_age
discard_weight_at_age
log_index_sd
landings_proportion
```

`landings_numbers_at_age` and `landings_fraction_at_age` retain removal components
used by a source model alongside its combined catch. Landings and discard weights
are distinct from combined catch weights. `log_index_sd` stores a supplied
standard deviation of the log survey index as a separate year record, with age
when the index is age-specific, and
the same survey identity as the corresponding observation; use `type = index`,
`basis = log_scale`, and document whether it fixes variance or supplies a relative
weight. `landings_proportion` retains a source landings-weight composition by substock
or quarter (`basis = proportion_biomass`); these records have no age. Reference
components are retained, and report rounding is not silently renormalized.
These measures are required by the 2025 Northern Shelf cod assessment.

`larval_abundance_index` preserves the partial spawning-component larval indices used by North Sea herring. Use `type = index`, `basis = numbers`, blank fish age, and explicit spawning component in `region` and survey time window in `season`. Native units and configured model timing must be documented separately; time-window column numbers are not fish ages.

Additional measures may be added when necessary, but should remain explicit and biologically interpretable.

## Initial `basis` vocabulary

``` text
numbers
biomass
proportion_numbers
proportion_biomass
kg_per_fish
proportion
per_year
log_scale
```

For example:

``` text
type = index
measure = proportion_at_age
basis = proportion_numbers
```

is distinct from:

``` text
type = index
measure = proportion_at_age
basis = proportion_biomass
```

## `source_type`

Use:

``` text
native_model
official_machine_readable
official_table
official_document
digitized
reconstructed_source_input
charbonneau_seed
```

Use `official_document` for values stated in an official report's prose when
they are not presented in a numerical table. Use `official_table` for values
transcribed from tables and `digitized` for values recovered from figures.

`reconstructed_source_input` means that the actual accepted-model input was not directly available but was reconstructed from authoritative source quantities. Record the reconstruction in `transformation`.

Do not use `reconstructed_source_input` for values derived solely for a downstream analysis.

## Reshaping

Changing a source table from wide to long form is allowed when scientific meaning is unchanged.

Catch-at-age is the working quantity required by the tinyAM translation.
Distinguish two source situations:

- If proportions describe an underlying numbers-at-age input that is not
  otherwise available, recover that accepted-model input using matching total
  removals. Record the calculation, year/fleet/area, source references, and
  units in `transformation`, with `source_type = reconstructed_source_input`.
- If the accepted model fits a total and an age composition separately, retain
  both native inputs. Their conversion to tinyAM numbers-at-age is a derived
  analysis product and belongs outside the canonical tables.

Number proportions multiplied by total numbers give numbers-at-age; multiplied
by total biomass, they do not directly give fish counts. Conversions involving
biomass require compatible catch weights and the documented formulas in
`TINYAM_TRANSLATION.md`. Retain the source quantities needed to reproduce a
reconstruction, and do not sum overlapping direct and reconstructed values.
Total landings and gear detail are otherwise optional unless material to the
accepted model or the conversion.

Example:

``` text
year | age1 | age2 | age3
```

may be represented as:

``` text
year | age | value
```

------------------------------------------------------------------------

# 5. `outputs.csv`

Long-format outputs from the accepted assessment.

## Required columns

| Column             | Description                     |
|--------------------|---------------------------------|
| `assessment_id`    | Foreign key                     |
| `type`             | Broad output family             |
| `measure`          | Exact output represented        |
| `fleet`            | Fleet where applicable          |
| `survey`           | Survey where applicable         |
| `sex`              | Sex where applicable            |
| `region`           | Region where applicable         |
| `season`           | Season where applicable         |
| `year`             | Year; blank for time-invariant estimates such as reported catchability parameters |
| `age`              | Age where applicable            |
| `age_group`        | Optional age range for grouped measures such as Fbar or Mbar |
| `value`            | Point estimate                  |
| `se`               | Standard error if available     |
| `lwr`              | Lower interval if available     |
| `upr`              | Upper interval if available     |
| `unit`             | Unit                            |
| `source_type`      | Provenance category             |
| `source_reference` | Specific file/table/object/page |
| `notes`            | Additional clarification        |

## High-priority output measures

``` text
numbers_at_age
fishing_mortality_at_age
natural_mortality_at_age
SSB
recruitment
```

Use `natural_mortality_at_age` when M is estimated.

## Secondary output measures

``` text
total_biomass
predicted_catch
predicted_index
Fbar
Mbar
q
Fmsy
Bmsy
msy
```

Use measures only when they are defined by the accepted assessment.

Retain uncertainty where readily available.

## Uncertainty and estimate definitions

The `se` column records a standard error, not a coefficient of variation (CV).
Document its scale in `notes`: normally the natural scale of `value`, or an
explicitly reported log scale. Identify delta-method natural-scale SEs as
approximations. Do not silently put a CV in `se` or infer a scale from the size
of the number. If only a CV is available, preserve it with its source meaning
in notes; any conversion to an SE must be justified and documented.

For `lwr` and `upr`, record the interval level and whether it is a confidence
or credible interval. Identify whether `value` is a mean, median, or another
reported estimate when that distinction matters. Unknown uncertainty or
interval construction remains unknown; missing values must not be interpreted
as zero uncertainty or used to invent new intervals.

Keep native definitions, such as recruitment age, SSB timing and sex convention,
Fbar age range and weighting, and terminal age groups, in the associated
assumptions and row notes. Common-definition comparisons calculated for tinyAM
belong in analysis outputs rather than replacing source values here.

`type` identifies the broad output family; `measure` identifies the exact
quantity. For example, an age-specific fishing mortality row uses
`type = mortality` and `measure = fishing_mortality_at_age`.

------------------------------------------------------------------------

# Canonical-data integrity

Validation should check:

- required columns;
- unique `stock_id`;
- unique `assessment_id`;
- valid foreign keys;
- controlled `source_type` values;
- sensible year, age, and numerical fields;
- maturity proportions in `[0,1]`;
- `sampling_time` in `[0,1]` where represented numerically;
- duplicate canonical rows;
- missing fleet or survey labels where the model distinguishes them;
- unexplained multiple current assessments for the same stock.

Per-assessment summaries should include:

- input rows by `type`, `measure`, fleet, and survey;
- fleets represented;
- surveys represented;
- year and age coverage;
- outputs represented;
- assumptions represented.

Downstream-analysis-specific transformations, settings, compatibility judgments, and fitted results do not belong in the canonical database.

`spawning_weight_at_age` distinguishes supplied spawning-stock weights from
population weights when the source uses different biological surfaces. Use
`type = weight`, `basis = kg_per_fish`; do not label biological purposes as spatial regions.

## Native length and conditional age compositions

`proportion_at_length` and `conditional_proportion_at_age` retain the composition
representation consumed by an assessment. Use `basis = proportion_numbers` and
`unit = proportion`. Do not replace these with an unconditional age distribution.

Optional input columns preserve their dimensions:

- `observation_id`: stable source record identifier, shared by its composition cells.
- `length_bin`: native length-bin label for a length composition, in cm.
- `length_bin_lower`, `length_bin_upper`: native conditioning-bin codes/bounds for
  conditional age composition; specify the coding method and units in notes.
- `sample_size`: supplied input sample size, distinct from an estimated effective size.
- `age_error`: native ageing-error definition code; interpretation belongs in assumptions.
- `partition`: native retained/discard/combined observation code, documented in notes.

These columns are blank for inputs to which they do not apply. Preserve the
source month in notes and record survey timing separately as a fraction of year.

## Environmental covariates

Use `type = covariate`, `measure = environmental_covariate` and
`basis = native_covariate` for supplied environmental series used by the model.
Signed finite values are permitted. Preserve the native scale and identify the
variable and consuming model component in provenance/notes; do not standardize
it again or treat it as an abundance index.

`type = biology` with `fraction_F_before_spawning` or
`fraction_M_before_spawning` records supplied year-age mortality fractions used
in spawning biomass, with `basis = proportion`. Preserve fleet dimensions for F
where present. These are biological inputs, separate from survey sampling timing.

Time-invariant `q_power` outputs may have a blank year, like `q`. Record the
power separately from the catchability coefficient and document the observation
equation; a coefficient in a power relationship is not a simple catch fraction.
Natural-scale SEs obtained by the delta method must be identified as such.
