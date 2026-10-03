# Translating the Assessment Database to tinyAM

## Purpose

This document defines how a curated source assessment is transformed into data and model settings suitable for a tinyAM analysis.

The canonical assessment database is intentionally richer than tinyAM. It preserves the fleets, surveys, biological inputs, and statistical assumptions of the accepted source assessment.

Translation therefore has two distinct stages:

1.  **data processing** — convert curated source inputs into the data structures required by tinyAM;
2.  **model specification** — audit source assumptions and choose the closest defensible tinyAM representation.

These steps are analysis choices and must not alter the canonical source database.

------------------------------------------------------------------------

# 1. Translation workflow

For a selected `assessment_id`:

``` text
canonical database
        ↓
select source inputs
        ↓
prepare catch-at-age
        ↓
prepare survey indices-at-age
        ↓
prepare stock weight and maturity
        ↓
prepare M
        ↓
construct tinyAM obs
        ↓
audit source assumptions
        ↓
choose tinyAM settings
        ↓
make_dat()
        ↓
fit_tam()
        ↓
compare with source-assessment outputs
```

A useful analysis-local structure is:

``` text
analysis/comp_assessments/
├── R/
│   ├── prepare_tiny_obs.R
│   ├── audit_assumptions.R
│   └── helpers_*.R
├── results/
│   ├── processed_inputs/
│   ├── audits/
│   └── fits/
└── ...
```

The exact file names may evolve. The important boundary is that translation code remains outside the canonical database.

------------------------------------------------------------------------

# 2. Current tinyAM observation structure

tinyAM expects:

``` r
obs <- list(
  catch = ...,
  index = ...,
  weight = ...,
  maturity = ...
)
```

Every component contains:

``` text
year
age
obs
```

The index table additionally requires:

``` text
survey
samp_time
```

Formula covariates such as `q_block` are also needed when referenced by model
settings. The analysis-local database translator adds `q_block` for age and
`q_key` for observed survey-by-age combinations, allowing the settings to show
whether catchability is shared across surveys or estimated separately.

Natural mortality is supplied through model settings rather than as a fifth `obs` table.

------------------------------------------------------------------------

# 3. Separate deterministic data processing from model choices

Data processing should be as mechanical and reproducible as possible.

Examples:

- unit conversion;
- reconstructing numbers-at-age from totals and compositions;
- aggregating fishing fleets;
- expanding constant biological vectors over years;
- building complete year-age grids;
- renaming database `sampling_time` to tinyAM `samp_time`.

Model specification is different. It includes decisions such as:

- whether N process deviations are IID, RW, AR1, or off;
- how F is represented;
- how M is represented;
- how q is parameterized;
- how observation SDs are shared.

Do not hide model choices inside basic data-wrangling helpers.

------------------------------------------------------------------------

# 4. Preparing catch-at-age

tinyAM currently uses one aggregate catch-at-age observation table.

The source assessment may instead contain:

- one fleet;
- multiple fleets;
- landings and discards separately;
- totals plus age compositions;
- biomass and composition rather than direct numbers-at-age.

The translation step should first reconstruct a numbers-at-age series for each source catch/removal stream, then aggregate those streams as appropriate.

## 4.1 Direct numbers-at-age

If the source database contains direct catch numbers-at-age, harmonize units only as needed.

## 4.2 Total numbers + number proportions

If total catch in numbers $C_{f,t}$ and number proportions $p^N_{f,t,a}$ are available:

$$C_{f,t,a} = C_{f,t} p^N_{f,t,a}.$$

`database_to_tiny_obs()` applies this calculation only when the number-
proportion rows match one annual `total_numbers` row on year and catch-stream
identity. It retains the source proportions as reported rather than
renormalizing them. Missing matches, duplicate totals, and overlapping direct
and reconstructed year-age values are errors; the helper never silently drops
or double counts catch input.

## 4.3 Total biomass + biomass proportions

If total catch biomass $B_{f,t}$ and biomass proportions $p^B_{f,t,a}$ are available:

$$C_{f,t,a}
=
\frac{B_{f,t}p^B_{f,t,a}}
{w^C_{f,t,a}},$$

where $w^C_{f,t,a}$ is catch weight-at-age.

## 4.4 Total biomass + number proportions

If total catch biomass $B_{f,t}$ and number proportions $p^N_{f,t,a}$ are available:

$$C^{N}_{f,t}
=
\frac{B_{f,t}}
{\sum_a p^N_{f,t,a} w^C_{f,t,a}},$$

then:

$$C_{f,t,a}
=
C^N_{f,t}p^N_{f,t,a}.$$

## 4.5 Total numbers + biomass proportions

If total catch in numbers $C^N_{f,t}$ and biomass proportions $p^B_{f,t,a}$ are available:

$$C_{f,t,a}
=
C^N_{f,t}
\frac{p^B_{f,t,a}/w^C_{f,t,a}}
{\sum_j p^B_{f,t,j}/w^C_{f,t,j}}.$$

All four total/composition combinations above are supported when the required
matching total and composition rows are present. Biomass-based catch
reconstruction requires compatible catch-weight rows and never substitutes
stock weights. It converts the full source age composition before applying a
requested age range, so biomass-to-number calculations retain their correct
denominators.

## 4.6 Aggregate fleets

When the tinyAM analysis intentionally represents the fishery as one aggregate catch process:

$$C^{tinyAM}_{t,a}
=
\sum_f C_{f,t,a}.$$

This is a deliberate simplification relative to assessments that estimate fleet-specific F/selectivity.

The source database remains fleet-specific.

The translation/audit should record that the tinyAM model is using an aggregate fishery.

## 4.7 Landings and discards

If landings and discards are both true removals in the accepted assessment and can be converted to compatible numbers-at-age:

$$C_{t,a}^{total}
=
C_{t,a}^{landings}
+
C_{t,a}^{discards}.$$

Do not combine removal streams without understanding how the source assessment treats them.

------------------------------------------------------------------------

# 5. Preparing survey/index-at-age data

Every survey used in the accepted assessment should be considered.

Unlike fishing fleets, distinct surveys should generally remain distinct because tinyAM supports multiple survey series with separate catchability structures.

For each source survey, create an age-specific abundance index where scientifically defensible.

## 5.1 Direct numbers-at-age

If the survey provides direct abundance-at-age:

$$I^N_{s,t,a}$$

retain the survey identity and harmonize units only as needed. If the source
reports values on a native survey-index scale, keep those values unchanged;
do not treat an index as a count of fish.

## 5.2 Total abundance + number proportions

$$I^N_{s,t,a}
=
I^N_{s,t} p^N_{s,t,a}.$$

## 5.3 Total biomass + biomass proportions

$$I^N_{s,t,a}
=
\frac{I^B_{s,t}p^B_{s,t,a}}
{w_{s,t,a}}.$$

## 5.4 Total biomass + number proportions

First estimate total abundance:

$$I^N_{s,t}
=
\frac{I^B_{s,t}}
{\sum_a p^N_{s,t,a}w_{s,t,a}},$$

then:

$$I^N_{s,t,a}
=
I^N_{s,t}p^N_{s,t,a}.$$

## 5.5 Total abundance + biomass proportions

$$I^N_{s,t,a}
=
I^N_{s,t}
\frac{p^B_{s,t,a}/w_{s,t,a}}
{\sum_j p^B_{s,t,j}/w_{s,t,j}}.$$

## 5.6 Choice of weight for survey conversion

Use, in order of preference:

1.  survey-specific weight-at-age used by the source assessment;
2.  another weight series explicitly associated with that index;
3.  stock weight-at-age as a documented approximation.

Do not silently use stock weight-at-age when the source assessment defines a different survey weight.

## 5.7 Survey timing

Map the database field:

``` text
sampling_time
```

to tinyAM:

``` text
samp_time
```

`samp_time` must be numeric in `[0,1]`.

When the exact fractional timing is not available:

- derive it from documented survey timing only when defensible;
- record the approximation;
- do not default all surveys to `0.5`.

## 5.8 Spatial scale

Use the survey series corresponding to the spatial extent actually fitted by the accepted assessment.

Do not use finer-scale subarea series when the assessment fits a combined index unless the analysis explicitly intends to deviate from the accepted model.

The observation translator stops when a selected survey has a measure that
cannot be interpreted as an age-specific abundance index. Do not force a
special likelihood, such as a larval or spawning-component index, into the
standard abundance-index table. Exclude it explicitly for a limited
approximation or define and audit a scientifically defensible mapping first.

For North Sea herring, retain `HERAS`, `IBTS0`, `IBTS-Q1`, and `IBTS-Q3`
as the selected surveys. The four `LAI-*` spawning-component series are
excluded pending mapping, as recorded in
`results/audits/ices_herring_north_sea_2026_translation_decisions.csv`.
Native survey-index values are retained even when their numerical units are
unresolved; their scale is absorbed by survey catchability. This does not
resolve sampling timing, which must still be documented before fitting.
The observation translation records selected and excluded survey names in
its provenance attribute.

------------------------------------------------------------------------

# 6. Preparing stock weight-at-age

tinyAM requires a complete modeled year × age stock-weight table with no missing values.

The source database may contain:

- an annual age-year matrix;
- one constant age vector;
- values extending above the modeled plus age;
- sex-specific values.

## 6.1 Annual values

Use the source annual values directly after unit harmonization.

Prefer kg per fish.

## 6.2 Constant vector

If the source assessment assumes a constant vector $w_a$, expand it across modeled years during translation:

$$w_{t,a}=w_a.$$

Do not duplicate these rows in the canonical database merely to satisfy tinyAM.

## 6.3 Missing values

Do not silently interpolate missing source values.

If a complete weight surface cannot be constructed from defensible source information, flag the assessment as not yet fit-ready.

------------------------------------------------------------------------

# 7. Preparing maturity-at-age

tinyAM requires a complete modeled year × age maturity matrix with no missing values.

If the source assessment uses a constant age vector $m_a$:

$$m_{t,a}=m_a$$

for all modeled years.

If maturity varies annually, use the annual source values.

Preserve the source convention for sex whenever possible.

Do not average male/female maturity without a documented reason.

All translated maturity values must lie in `[0,1]`.

------------------------------------------------------------------------

# 8. Preparing natural mortality

Natural mortality is not part of `obs`, but a baseline or mean M representation is required to fit tinyAM.

The source assessment may use:

- constant M;
- age-specific M;
- year-specific M;
- age-year M;
- estimated M;
- time-varying latent M.

## 8.1 Fixed or externally supplied M

Translate the numerical source M into a form usable in `M_settings$mu_supplied`.

For example, a constant source value may become:

``` r
M_settings = list(
  process = "off",
  mu_supplied = ~ I(0.2)
)
```

An age-year source surface may be joined to the weight table and referenced through a column such as:

``` r
~ M_assumption
```

## 8.2 Estimated or time-varying M

Do not automatically reproduce source-model M dynamics.

Use the assumption audit to determine:

- whether tinyAM can represent the structure directly;
- whether a simpler baseline M should be used;
- whether an M-process sensitivity analysis is warranted.

This is a model-specification decision, not data wrangling.

------------------------------------------------------------------------

# 9. Build the tinyAM observation object

After processing:

``` r
obs <- list(
  catch = catch,
  index = index,
  weight = weight,
  maturity = maturity
)
```

## `catch`

Required:

``` text
year
age
obs
```

The current tinyAM validator requires exactly one row per modeled year × age combination.

NA catch observations are allowed.

Zero catch is structurally accepted but is treated as missing by the log-scale observation model.

## `index`

Required:

``` text
year
age
obs
survey
samp_time
```

Survey tables may be sparse.

`samp_time` must contain no missing values and must lie in `[0,1]`.

## `weight`

Required:

``` text
year
age
obs
```

Must form a complete modeled year × age grid.

No missing `obs`.

## `maturity`

Required:

``` text
year
age
obs
```

Must form a complete modeled year × age grid.

No missing `obs`.

Then run:

``` r
tinyAM::check_obs(obs)
```

A passing result means the translated observation object satisfies current tinyAM structural requirements.

It does not establish that the chosen model settings faithfully reproduce the source assessment.

------------------------------------------------------------------------

# 10. Plus groups

The source assessment's modeled plus group should be identified during curation.

tinyAM may model the same plus age or a deliberately different age range.

When biological weight/maturity data extend above the modeled plus age, retain those source values where possible because tinyAM can use hidden older ages in its internal plus-group biology.

Do not silently collapse source data before deciding on the modeled age range.

If a source survey reports one terminal group such as 12+ while the assessment
models population ages through 15, keep the survey observation as its original
12+ group (stored at age 12) and retain the model's 13–15 population ages. A
fitted model object may have age-specific population outputs at 13–15 without
having separate survey observations for those ages. Do not copy fitted values
or split the 12+ survey observation into invented age-specific inputs.

------------------------------------------------------------------------

# 11. Audit source assumptions before choosing settings

The translation should produce an analysis-local audit comparing source assumptions with current tinyAM capabilities.

Suggested output fields:

``` text
component
fleet
survey
source_setting
source_value
tinyam_support
proposed_representation
notes
```

Suggested support classes:

``` text
supported
partially_supported
unsupported
not_checked
```

The audit is not part of the canonical source database.

------------------------------------------------------------------------

# 12. Choosing `N_settings`

Review source assumptions concerning:

- recruitment treatment;
- process deviations;
- cohort survival;
- correlation through age/year;
- first-year abundance treatment.

Possible current tinyAM process options:

``` text
off
iid
rw
ar1
```

Initialization options:

``` text
exp
free
random
```

Do not assume a source model's process is equivalent merely because both use the words "random walk" or "AR1".

Compare the actual mathematical structure.

For broad comparative analyses, a standardized tinyAM N process may intentionally differ from the source assessment. Record that explicitly.

------------------------------------------------------------------------

# 13. Choosing `F_settings`

The source assessment may contain several fleets and fleet-specific selectivities/F processes.

The standard tinyAM translation may instead aggregate catch across fleets and estimate one flexible F-at-age surface.

This is a deliberate modeling simplification.

Review:

- fleet-specific F structure;
- selectivity;
- time processes;
- age correlations;
- F-state sharing.

Then choose among current tinyAM:

``` text
iid
rw
ar1
```

and an optional mean-F formula.

Do not claim fleet-level replication when the translated model uses aggregate catch.

------------------------------------------------------------------------

# 14. Choosing `M_settings`

Review:

- source M level/surface;
- whether M is fixed or estimated;
- age/time variation;
- correlation/process structure.

Then determine:

- supplied baseline/mean M;
- whether `process = "off"` is most faithful;
- whether an IID/RW/AR1 M process is warranted;
- age blocks;
- first deviation year.

For standardized cross-stock analyses, document any deliberate departure from source-model M treatment.

------------------------------------------------------------------------

# 15. Choosing catch observation settings

Review the source catch likelihood and observation variance structure.

Current tinyAM allows:

``` r
catch_settings = list(
  sd_form = ...,
  sd_supplied = ...,
  fill_missing = ...
)
```

When multiple source fleets have been aggregated, source fleet-specific observation-error assumptions may no longer map directly.

The audit should explicitly identify this loss of structure.

------------------------------------------------------------------------

# 16. Choosing survey q and index settings

For every survey retained in the translated data, review:

- q age sharing;
- q time variation;
- q constraints;
- observation SD structure;
- observation correlation.

Current tinyAM can represent q using formula-based structures and optional monotone age effects.

A source assessment may have survey-specific parameter blocks that can be represented using covariates and formulas.

However, distinguish carefully between:

- sharing mean q parameters;
- sharing latent states;
- correlated observation/process errors.

Formula equivalence does not imply stochastic-process equivalence.

------------------------------------------------------------------------

# 17. Replication-oriented versus standardized models

The same curated assessment may support at least two different tinyAM objectives.

## 17.1 Replication-oriented model

Goal:

> approximate the accepted assessment as closely as current tinyAM reasonably allows.

This model may use source-specific q blocks, M assumptions, observation-error structures, and process choices.

## 17.2 Standardized comparative model

Goal:

> fit a common tinyAM formulation across many stocks so that estimated population-process deviations are comparable.

This model may intentionally simplify:

- fleet structure;
- selectivity;
- q structure;
- observation-error structure;
- process assumptions.

These two objectives should not be conflated.

A successful replication-oriented model does not automatically define the standardized comparative model.

------------------------------------------------------------------------

# 18. Save processed translation outputs separately

Useful derived artifacts may be saved under:

``` text
results/processed_inputs/
```

Examples:

- fleet-specific reconstructed catch-at-age;
- aggregate catch-at-age;
- survey numbers-at-age derived from totals/compositions;
- expanded weight/maturity matrices;
- M surfaces;
- translation diagnostics.

These are analysis products.

They should not be written back into canonical `inputs.csv` unless they are themselves verified accepted-model inputs.

------------------------------------------------------------------------

# 19. Translation provenance

Every nontrivial transformation should be auditable.

For each derived series, record enough information to reproduce:

- source rows used;
- formula;
- units;
- weights used;
- aggregation rule;
- any approximation;
- reason for the transformation.

A practical processed-data table may include fields such as:

``` text
assessment_id
source_type
source_measure
target_measure
formula
notes
```

or equivalent analysis-local metadata.

------------------------------------------------------------------------

# 20. Fit-readiness checklist

A curated assessment is ready for an initial tinyAM fit when the translation can produce:

## Catch

- numerical catch-at-age;
- one complete modeled year × age grid;
- common units across fleets before aggregation.

## Surveys

- age-specific abundance index for each retained survey;
- survey identity;
- valid `samp_time`.

## Weight

- complete modeled year × age stock-weight matrix.

## Maturity

- complete modeled year × age maturity matrix.

## Natural mortality

- a defensible baseline/mean M representation.

## Validation

- `tinyAM::check_obs(obs)` passes.

Fit-readiness does not imply source-assessment equivalence.

------------------------------------------------------------------------

# 21. Compare fitted outputs with the accepted assessment

Use source `outputs.csv` to compare, where available:

- N-at-age;
- F-at-age;
- SSB;
- recruitment;
- Fbar;
- predicted catch/index;
- uncertainty.

Differences should be interpreted in light of documented translation choices, such as:

- aggregated fleets;
- simplified selectivity;
- different process structure;
- simplified q;
- simplified observation errors.

Do not interpret a difference as a model failure before checking whether it follows from an intentional translation choice.

In addition to basic percent difference and bias summary statistics for the abovmentioned outputs, build a tam_list using available output data from the accepted assessment and use it in `vis_tam()` to produce an interactive html to ease visual comparison.

------------------------------------------------------------------------

# 22. Scope discipline

Do not preemptively build exported package converters.

Prefer small analysis-local helpers that emerge from repeated transformations.

Keep the source database independent of tinyAM.

Keep tinyAM translation transparent, reproducible, and scientifically explicit.
