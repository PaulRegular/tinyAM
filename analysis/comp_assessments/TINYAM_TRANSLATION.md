# Translating the Assessment Database to tinyAM

## Purpose

This document defines how a curated source assessment is transformed into data and model settings suitable for a tinyAM analysis.
Read the [scope and provenance statement](README.md#scope-and-provenance)
before interpreting these illustrative model comparisons.

The canonical assessment database is intentionally richer than tinyAM. It preserves the fleets, surveys, biological inputs, and statistical assumptions of the accepted source assessment.

Within an assessment, preserve the finest authoritative data available even
when the accepted model pools ages or years in its likelihood. A translation
may use a more detailed reported series when tinyAM can represent it, while
keeping the source likelihood's pooling and any report rounding explicit.
Aggregate reported ages only as required by tinyAM's model age grid; for
example, a model plus group at age 10 receives the sum of reported ages 10 and
older. Do not fit both overlapping detailed and pooled versions of the same
observations.

The expected result for each selected assessment is a documented tinyAM model,
a record of its convergence, and, when it converges, a comparison with the
accepted assessment. Fits and dashboards can be cached locally when useful;
they are not routine committed outputs. Preparing observation tables alone does
not complete a translation task.

Data conversion and model specification are separate: the shared converters
reshape recorded information, while the stock script makes the scientific
choices. Missing or incorrect source records must be repaired under
`PROTOCOL.md` before translation resumes. Such repairs improve the source
database; tinyAM-specific transformations and simplifications remain separate
analysis products.

`DATABASE_STRUCTURE.md` defines the canonical tables, `PROTOCOL.md` defines
source curation, and this document defines the workflow from those tables to
fitted models and comparisons. Use R for new extraction, validation,
translation, and analysis code.

------------------------------------------------------------------------

# 1. Translation workflow

Use `run_assessment()` for one stock while developing a translation and
`run_assessments()` for a reproducible batch. Normally select the record marked
`is_current`: the most recent accepted assessment for which detailed inputs,
assumptions, and outputs can be recovered for the canonical database. This is
the current detailed assessment and may predate the newest advice, FSAR, or
publication when later material is summary-only. Use another record only when
that assessment is explicitly part of the task. Keep all inputs, assumptions,
and outputs tied to the same `assessment_id`.

``` text
load canonical database
        ↓
select assessment and its stock specification
        ↓
database_to_tam_obs(): catch, index, weights, maturity, and M assumption
        ↓
review source assumptions and specify tinyAM settings
        ↓
check_obs() and prepare_tam()
        ↓
fit_tam()
        ↓
check convergence and record diagnostics
        ↓
if converged: construct source reference with database_to_tam_ref()
        ↓
return results; optionally render a dashboard or cache the fit locally
        ↓
for a batch, write shared diagnostics and comparison summaries once
```

When a source-data gap prevents translation or a requested comparison:

1. identify the missing quantity and whether it is an input, assumption, or output;
2. review `PROTOCOL.md`, the stock source review, and cached authoritative files;
3. obtain or correct the source records in the canonical database, retaining provenance;
4. validate the revised records and reload the database before retrying translation.

Record the database revision used for every fit. A reader restricted to committed
records must read the newly validated and committed revision after a repair,
rather than continue using its old snapshot. If authoritative information cannot
be recovered, report the blocker and continue with other requested stocks.

An unsupported model feature is not automatically a missing-data problem.
For example, a source model may estimate M without a fixed M input. That case
requires a documented tinyAM model choice, not a fabricated input table.

## Shared code and stock scripts

Keep the routine workflow small:

- one R driver loads the database and loops over stock specifications;
- `database_to_tam_obs.R` creates the observation list, adds numerical M to
  `obs$weight$M_assumption`, and records the source M treatment;
- `database_to_tam_ref.R` creates an assessment reference from recorded outputs;
- one small R script per stock prepares observations, shows a literal
  `fit_tam()` call, and defines background and comparison metadata at top level.

Shared code handles validation, diagnostics, references and exports. Stock
scripts should not duplicate converters or parse native source files. Keep
source-import scripts separate; they are used when repairing the database,
not on every model run. Keep these database-specific helpers analysis-local.

Do not return a settings list from `translate_stock()`. Keep years, ages,
survey selections, formulas, processes and initialization visible in the stock
script. Guard all fitting with `if (do_fit)`, including any explicit warm-start
call. The runner supplies `source`, `do_fit` and `silent` in an isolated
environment; `fit = FALSE` must not run a preliminary fit. Fitted settings live
in `fit$dat`, rather than a duplicate `result$settings` object.

A useful directory structure is:

``` text
analysis/comp_assessments/
├── R/
│   ├── read_database.R
│   ├── database_to_tam_obs.R
│   ├── database_to_tam_ref.R
│   ├── audit_assumptions.R
│   └── run_assessment.R
├── scripts/
│   ├── database/
│   └── translation/
│       ├── run_stock.R
│       ├── run_all.R
│       └── stocks/
│           └── <assessment_id>.R
├── tests/
├── results/
│   ├── observation_readiness.csv
│   ├── fit_diagnostics.csv
│   ├── comparison_summary.csv
│   ├── sensitivity_summary.csv
│   └── cache/                    # gitignored local fits and dashboards
└── source_cache/                 # gitignored authoritative source files
```

`run_stock.R` reads the working-tree database, sets the execution controls and
sources the selected stock script, leaving observations, fit and background
in the RStudio workspace. Use `run_assessment()` for shared comparisons and
optional dashboards/cache. `run_all.R` reads the committed snapshot and
uses `future::multisession` through `furrr`; workers load tinyAM from the current
repository checkout. Workers return objects only. The parent process writes
aggregate CSV files after all workers finish, and restores the previous future
plan. Set `parallel = FALSE` when debugging.

## Stock background

Each stock script supplies a short Markdown table aligned with `prepare_tam()`:

| Component | Accepted assessment | tinyAM representation | Reason for difference |
|---|------|------|------|
| Years | Fitted historical years | Years used in this fit | Explain any restricted period |
| Ages | Recruitment age and terminal age group | Model ages and plus group | Explain any age restriction or grouping |
| N | Recruitment, survival, process variation, and initial abundance | `N_settings` and recruitment treatment | Explain omitted or changed population assumptions |
| F | Fleets, selectivity, and changes through time | `F_settings` and catch aggregation | Explain the simplified fishery structure |
| M | Fixed or estimated mortality, age groups, and time variation | M assumption and `M_settings` | Explain the baseline and any estimated process |
| Catch | Removal streams and observation-error model | Catch observations and `catch_settings` | Explain reconstruction and error-model differences |
| Index | Surveys, timing, catchability, and observation errors | Retained surveys and `index_settings` | Explain exclusions, q sharing, and error-model differences |
| Weights | Stock, catch, survey, and spawning weights | Weight series used and conversions | Explain any substituted weights |
| Maturity | Age/year structure, sex convention, and spawning timing | Maturity and biomass convention | Explain differences in SSB definition |

Every stock background should include an `Assessment documentation` section
generated from the canonical assessment metadata, preferably with
`print_sources(source$assessment)`. The stock script controls where this section
appears in its background; `assessments.csv` remains the source of truth for
the URLs, and dashboard rendering does not append links automatically.

Define consecutive age and year blocks with `cut_ages()` and `cut_years()`.
Check boundary memberships and preserve the intended parameter-sharing labels.
Non-consecutive sharing, such as one age versus all other ages, can remain an
explicit conditional assignment.

Use short biological explanations rather than only process names or formulas.
Identify whether each choice preserves, approximates, or omits the source
assumption. Honor previously agreed stock-specific choices and record them
here. Include source references where needed, without repeating the full
curation history from `source_reviews/`.

Pass this Markdown to `vis_tam(background = background)`. The dashboard should
display a Background page only when text is supplied; ordinary dashboards
without background text retain their existing layout.

The Background option accepts general Markdown and imposes no table-specific
styling. For these stock tables, the separator `|---|------|------|------|`
sets relative column widths of 1:2:2:2 in Pandoc. Keep layout choices in the
stock scripts.

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

`database_to_tam_obs()` incorporates a numerical M assumption in
`obs$weight$M_assumption`; its M translation helper is defined in that same
file, so there is no separate M converter to source. When the database has
relative index-SD inputs, the converter joins them to the matching survey,
year, and age rows as `relative_sd` for an explicit `sd_supplied` setting.
`M_settings$mu_supplied` can then reference `~ M_assumption`.
`M_settings` determines whether mortality is fixed or has an estimated process;
there is no fifth M observation table or new observation class.

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

Use the accepted model's year and age range as the starting point. Do not
silently shorten the period, change the plus group, repeat an annual biological
series, or drop a survey merely to make validation pass. First investigate
whether the missing source information can be recovered. An explicitly chosen
restricted analysis must explain its limits in the stock background and use
the corresponding common period and age groups for comparisons.

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

`database_to_tam_obs()` applies this calculation only when the number-
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

Do not add overlapping direct catch-at-age and reconstructed totals/compositions
twice. Retain the source representation needed to understand the accepted
likelihood; follow the canonical catch rules in `DATABASE_STRUCTURE.md` and
`PROTOCOL.md` when deciding whether a reconstructed series is a source input
or a derived translation product.

## 4.7 Landings and discards

If landings and discards are both true removals in the accepted assessment and can be converted to compatible numbers-at-age:

$$C_{t,a}^{total}
=
C_{t,a}^{landings}
+
C_{t,a}^{discards}.$$

Do not combine removal streams without understanding how the source assessment treats them.

------------------------------------------------------------------------

## 4.8 Reporting total catch biomass

Stock recipes may optionally provide `catch_reporting`, with `weights`
(`year`, `age`, `weight` in kg/fish) and `totals` (`year`, `yield` in kg).
The runner then reports the original total catch biomass and sums predicted
catch numbers times catch weights over all model ages. This avoids comparing
an incomplete composition with an all-age prediction, or substituting stock
weights for catch weights. It changes the yield reporting tables only; it does
not add a total-catch likelihood or alter fitted predictions.

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

`database_to_tam_obs()` uses a matching survey-weight series by default and
falls back to the selected biological weight series when no matching series is
available. Set `index_weight_source = "stock"` when a documented approximation
requires stock weights for every survey, even when survey-specific weights are
also present. Record that choice in the stock background and translation
provenance. Do not silently use stock weight-at-age when the source assessment
defines a different survey weight.

For tinyAM's biological weight surface, an unlabeled stock-weight series is
selected by default when present. If the database contains only one labeled
series, it is used; if multiple labeled series are present, select one with
`weight_survey`. Use `weight_survey = ""` to explicitly select unlabeled stock
weights.

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
- consult the framework/benchmark research document and cited survey reports
  when the production report does not give sufficient timing detail;
- record the approximation;
- do not default all surveys to `0.5`.

A documented sampling season is not an exact annual index time. If a season
midpoint is used, pass it explicitly through `sampling_times` and record it in
the stock background; do not label it as a recovered source-model parameter.

## 5.8 Spatial scale

Use the survey series corresponding to the spatial extent actually fitted by the accepted assessment.

Do not use finer-scale subarea series when the assessment fits a combined index unless the analysis explicitly intends to deviate from the accepted model.

The observation translator stops when a selected survey has a measure that
cannot be interpreted as an age-specific abundance index. Do not force a
special likelihood, such as a larval or spawning-component index, into the
standard abundance-index table. Exclude it explicitly for a limited
approximation or define and audit a scientifically defensible mapping first.
Record stock-specific survey selection, exclusions, unresolved units, and
timing choices in that stock's background or translation-decisions file. The
observation translation records selected and excluded survey names in its
provenance attribute.

Deriving age-specific observations from one total and an age composition changes
the observation model. The derived ages share information from the same total;
tinyAM's age-specific lognormal likelihood does not automatically reproduce
that dependence or the source composition likelihood. Explain this approximation
in the stock background rather than claiming an exact likelihood mapping.

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

Follow the database-repair loop in Section 1 before accepting a restricted
period or age range. Expand a constant vector only when the source assumption
actually makes it constant.

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

A defensible baseline or mean M representation is required to fit tinyAM.
`database_to_tam_obs()` includes its numerical values in
`obs$weight$M_assumption`; `M_settings` specifies how those values enter the
model. The converter must distinguish a source-supplied M surface from an
explicitly chosen baseline for a tinyAM approximation.

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

The shared observation converter joins an age-year source surface to the weight
table by year and age, checks coverage and uniqueness, and makes it available as:

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

If numerical baseline values are chosen for the approximation, supply them
explicitly to the observation converter and record their origin and purpose.
Do not silently substitute source-estimated M outputs as fixed inputs, choose
an arbitrary default, or describe starting values as a fixed mortality
assumption. If M is estimated, the stock script must explain what is estimated
and how the supplied values or mean structure are used.

The absence of fixed M inputs in an assessment that estimates M is not, by
itself, a curation gap. Recover the source assumptions and available estimated
outputs, then specify a defensible tinyAM representation.

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

Also run `prepare_tam()` with the stock's proposed settings to check formula
covariates, supplied M, blocking, and model dimensions before fitting. Valid
observations alone do not establish that a model is ready to estimate.

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

Preserving the 12+ label does not make it a single-age-12 observation. Its
prediction must represent the matching age group, such as a sum over ages
12–15. If the current tinyAM observation model cannot represent that group,
use a documented compatible alternative or exclude the grouped observation
explicitly. Do not fit it against age-12 abundance alone. Apply the same rule
when comparing grouped population or mortality outputs.

For fitted-output comparisons, a report's display groups are separate from the
model's plus group. A stock recipe can declare an explicit
comparison_age_groups mapping, such as N = list("10+" = 10:15), so a source
10+ estimate is compared with the sum of tinyAM ages 10-15 without changing the
tinyAM model plus age. Never compare a reported 10+ value with tinyAM age 10
alone. In database_to_tam_ref(), additive N and biomass values are summed
when multiple source ages are collapsed; F and M are N-weighted only when
multiple source ages are actually combined. A single source rate already at
the requested terminal age is retained as reported.
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

Recruitment is separate from older-cohort deviations: `N_settings$rec_form`
defaults to `~ rw(year)`, while `N_settings$process` controls survival variation.
Supported alternatives include `~ iid(year)`, `~ bh(ssb) + iid(year)` and
`~ ricker(ssb) + ar1(year)`, with fixed annual covariates when justified.
Recruitment means abundance at the youngest modeled age. Parent SSB is
start-of-year mature biomass; the default lag is that recruitment age. Early
states with no modeled parent SSB are fixed boundary parameters, not fabricated
historical SSB. The curves describe median recruitment, without an automatic
lognormal correction. Check the source equation, age, lag, spawning convention,
bias correction and fixed versus estimated parameters before calling these
options equivalent. Fixed BH steepness one does not justify estimating a free
two-parameter BH curve.

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

Mean formulas can include `iid()`, `rw()` or `ar1()` effects, for example
`mu_form = ~ factor(age) + rw(year)` with an IID residual F process. Follow the
package's safeguards for overlapping temporal processes. A shared temporal
mean plus independent residuals does not reproduce correlated RW innovations
across ages. `logistic()` is a catchability curve, not a fishery-selectivity term.

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

Structured IID/RW/AR1 mean effects are also available for M. Extra M mean and
residual variation may be weakly separable; inspect uncertainty and advisories.
Source priors or fixed latent F/M process SDs are not reproduced merely by
using their values as starting parameters.

A standalone stochastic M mean can use a known increment SD, for example
`process = "off", mu_form = ~ 0 + age_group + rw(year, by = age_group, sd = 0.075)`.
Here M varies through the mean formula even though the separate residual
process is off. Its first effect is anchored and the group coefficients set
the starting levels; this does not supply a prior on those levels. Check the
source first-increment year and parameter-integration convention explicitly.

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

Catch and index `sd_form` can contain fixed age trends or Gaussian formula
effects. Age effects may help represent noisier young and old observations,
but mean effects, variance effects and observation correlation are different
models. Temporal variance effects remain experimental; follow the package's
support warnings and compare simpler alternatives. These formulas do not
change the scalar SDs of the latent F, M or N residual processes.

------------------------------------------------------------------------

# 16. Choosing survey q and index settings

For every survey retained in the translated data, review:

- q age sharing;
- q time variation;
- q constraints;
- observation SD structure;
- observation correlation.

Current tinyAM supports fixed q formulas, `mono()` and rising `logistic()` age
curves, IID/RW/AR1 Gaussian effects and random intercepts `(1 | group)`.
Gaussian effects allow numeric multipliers or categorical `by` groups;
logistic curves allow categorical groups. Preserve the retained `q_link`
(log or logit) unless source evidence supports changing it. Effects act on
that link scale; the logistic age curve multiplies q after applying the link.
It does not represent a dome or double-logistic selectivity curve.

A source assessment may have survey-specific parameter blocks that can be represented using covariates and formulas.

However, distinguish carefully between:

- sharing mean q parameters;
- sharing latent states;
- correlated observation/process errors.

Formula equivalence does not imply stochastic-process equivalence.

Preserve source q and observation-SD sharing through formulas wherever the
mapping is equivalent. A nonlinear or power catchability relationship is not
the same as a simple proportional q coefficient; identify any unsupported
structure and explain the chosen approximation. Do not adjust source
observations using fitted catchability or other fitted parameters to make a
simpler observation model appear equivalent.

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

## Review and retention policy

Before fitting candidates, review source equations and which parameters are
active, fixed or penalized. Specify a small candidate set, test individual
changes first, and combine them only when those results support it. Keep
observations, historical periods, ages, units and comparison definitions
unchanged. Retain every warning, failure and attempted model in the local
review cache; track reproducible candidate code and one compact results table.

Use `check_tam()` for optimizer status, gradients, Hessian and numerical checks.
Review structural checks and statistical-support advisories separately.
Inspect residuals and parameter uncertainty, not just a convergence flag.
Evaluate scale differences and `trend_correlation` separately for matched
SSB, recruitment, N and F. Label approximate comparisons, exclude unavailable
or non-equivalent ones, and do not count fixed-M agreement as improvement.
Keep the baseline when results are mixed or weakly supported; flag the
trade-off for review. Raw objectives across different process formulations
and likelihood equality are not model-selection criteria.

Stock backgrounds describe only the retained model. Detailed development
history stays in review code/results and Git. Refresh aggregate diagnostics
and comparisons only with a complete final batch, not a subset of stocks.

------------------------------------------------------------------------

# 18. Keep generated outputs small

The committed results are the assessment-level observation/data-readiness
table (`observation_readiness.csv`), one
diagnostics table, one comparison-summary table, and the small sensitivity
summary retained from prior focused checks. Do not commit duplicate per-stock
fits, translated observations, audits, settings, or dashboards; those can be
recreated from the database and stock translation script.

`run_assessment(..., cache = TRUE)` optionally saves a fitted object under
`results/cache/<assessment_id>/`. A dashboard is optional through
`dashboard = TRUE` and is cached there only when caching is enabled. Both
`run_assessment()` and `run_assessments()` leave caching off by default.
`run_assessments(..., save_results = TRUE)` writes the shared diagnostic and
comparison-summary tables once in the parent process.

These products must not be written back into canonical `inputs.csv` or
`outputs.csv`. Source-data corrections identified during translation follow
the separate curation and validation process.

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

# 20. Model-readiness checklist

A selected assessment is ready for an initial tinyAM fit when its stock recipe,
translated observations, and proposed settings satisfy the checks below.
Canonical completeness statuses and model readiness are different: a partial
assessment may support a clearly limited fit, while a complete source model may
contain features that tinyAM cannot represent.

The repository's `results/observation_readiness.csv` is an observation/data-readiness
screen. It checks source coverage, observation conversion, and `check_obs()`;
it does not load stock recipes or call `prepare_tam()`. The model-readiness checks
below apply after the recipe and settings are available.

## Catch

- numerical catch-at-age;
- one complete modeled year × age grid, with unavailable catch recorded as NA;
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

- a defensible baseline/mean M representation;
- numerical supplied values included in `obs$weight$M_assumption` when used;
- explicit treatment of source-estimated M and any missing source assumptions.

## Validation

- `tinyAM::check_obs(obs)` passes;
- `prepare_tam()` succeeds with the proposed years, ages, and settings;
- all retained observations have compatible age groups and prediction definitions;
- the stock background explains the accepted assumptions and chosen differences.

Passing these model-readiness checks does not imply source-assessment equivalence.

## Fit and assess convergence

Once the checks pass, attempt the stock's documented `fit_tam()` model.
Do not stop at exporting observations or a settings template. Ensure the R
environment can load tinyAM and the packages needed for fitting and rendering;
missing software is an execution blocker, not a source-data gap.

Use `is_converged(fit)` together with the optimizer result and uncertainty
diagnostics. Record at least:

- assessment ID, database revision, and attempt identifier;
- convergence result;
- optimizer convergence code and message;
- objective value and maximum absolute gradient;
- whether standard errors were obtained and whether the Hessian was positive
  definite, meaning the local curvature supports uncertainty estimation;
- errors or warnings that affect interpretation.

Catch failures per stock so one error does not stop the whole loop. Retain
diagnostics for failed or non-converged attempts, identify the reason, and
continue with other requested stocks. Report numerical convergence separately
from biological plausibility and source-model agreement.

An additional attempt may investigate a clear numerical issue, such as poor
starting values, but must record what changed. Do not repeatedly change data,
processes, or q structures merely to obtain convergence. Source estimates may
provide starting values; they must not silently constrain tinyAM estimates.

For a converged fit, `run_assessment()` returns the fitted object, source
reference, detailed differences, and summary in memory. The fit is cached only
when `cache = TRUE`. A dashboard is rendered only when `dashboard = TRUE` and
is stored in the local cache only when caching is also enabled; otherwise it is
written to a temporary file. A successful run does not require a per-stock fit
or detailed-difference file on disk. The batch runner writes aggregate
diagnostics and comparison summaries only when `save_results = TRUE`. Report
investigated blockers, failed attempts, and non-converged fits as separate
outcomes, rather than counting them as successful models. Preparing a usable
`obs` list alone is not completion.

------------------------------------------------------------------------

# 21. Compare fitted outputs with the accepted assessment

## Build the source reporting reference

Use `database_to_tam_ref()` to translate available canonical `outputs.csv`
records for the same assessment into a `tam_ref`. Pass the converged tinyAM fit
as `template` so the comparison list keeps the same output tables, years, ages,
and groups. The helper blanks reported values first, then fills values available
from the assessment; missing values remain `NA`. This `tam_ref` object has
reporting tables arranged like a `tam_fit`; it cannot be fitted, updated,
simulated, projected, or used for retrospective estimation.

Include year-age population tables and available aggregate trends. Supply
translated source observations where needed to show the input data, but leave
unavailable predictions, q estimates, parameters, residuals, and uncertainty
missing. Do not construct an optimizer or invent random effects to satisfy the
dashboard. Known fixed M may be shown from the documented source input, clearly
labelled as fixed rather than estimated.

`tidy_tam()` and `vis_tam()` must accept these reporting tables without assuming
that every `tam_ref` contains an optimizer or full assessment output. Tables
show unavailable entries as blank `NA` values, and plots show only the available
values.

Use source outputs to compare, where available:

- N-at-age;
- F-at-age;
- M-at-age;
- SSB;
- recruitment;
- other reported quantities, predictions, and uncertainty when useful.

## Match the definitions

Match units, fitted historical years, recruitment age, age groups, fleet/sex/
region dimensions, and whether estimates refer to the beginning of the year
or spawning time. Keep projection years separate from fitted historical years.
Do not split grouped outputs into invented single-age estimates or pool
different assessment runs.

Retain native source outputs. When definitions differ, identify the difference
in the background and dashboard. A separately labelled common-definition
quantity may be calculated only when the necessary source states and biology
are available. Do not present percent differences between incompatible
definitions as a like-for-like comparison.

Reference translation defines comparable recruitment at tinyAM's first modeled
age. If the accepted recruitment age differs, it uses accepted beginning-of-year
N at that age in the **same calendar year**. It never shifts years to match a
cohort. Missing accepted N remains missing. Native recruitment, original units
and provenance are retained in the database and the reference's `native_pop`
attribute. The aligned series supplies both the dashboard and numerical
comparison.

Reference translation fills the requested year–age grid, including missing
cells, and constructs its `comparisons` table before downstream reporting.
Numerical summaries use this table; they do not repeat age/year matching.
Annual common-definition outputs appear in the reference's reporting tables,
while native aggregates remain in `native_pop`. Derived aggregates have no
invented uncertainty.

It sums accepted N over exactly the modeled ages only when every required age is available in that
year. Partial age coverage does not become a partial annual total. Mortality
means are recalculated using matching ages and accepted N as weights; an
arithmetic source Fbar is not compared directly with tinyAM's population-weighted
Fbar.

For biomass and SSB, use matching reported biomass or mature-biomass-at-age
tables when available. Otherwise, when accepted N is available, a labelled
common-definition calculation can use the same translated weights and maturity
as tinyAM. This latter calculation compares population estimates under shared
biology, not the source model's native biomass or spawning-time SSB. The native
values stay in `native_pop` for separate source context. A tinyAM plug-in M baseline
must not be presented as accepted M; only reported estimated M or a documented
fixed numerical source input can fill that reference surface.

Unreconciled definitions are marked `non_equivalent`; missing numerical outputs
are marked `unavailable`. Both receive missing differences and an explanatory
reason. Grouped outputs are compared only after a documented age mapping, and
dimensions must map uniquely. Derived sums and mortality means do not receive
invented uncertainty from individual SEs without the required covariance.

Respect the reported uncertainty scale and interval meaning documented in
`DATABASE_STRUCTURE.md`. Keep source intervals when their confidence level or
construction cannot be reconstructed. Missing uncertainty is not zero
uncertainty; do not manufacture SEs or intervals from point estimates alone.

## Review the dashboard and comparison summary

For a converged fit, optionally use:

``` r
vis_tam(
  model_list = list(Assessment = assessment_reference, tinyAM = tam_fit),
  background = background,
  output_file = dashboard_file,
  open_file = FALSE
)
```

This example describes the intended Background interface noted in Section 1.
The dashboard is the primary way to inspect differences in trajectories and
age patterns. Include available SSB, recruitment, N, F, and M, and show other
panels only when their data are available.

For comparable metric values, the runner returns detailed differences in
`x$differences` and summary statistics in `x$summary`. The batch runner writes
the aggregate comparison summary when requested; detailed per-stock
differences remain in memory unless deliberately exported for a specific
review. Avoid routine per-stock exports that duplicate the dashboard.

For comparable metric values, calculate:

$$\text{percent difference}
= 100\frac{\text{tinyAM estimate} - \text{source estimate}}
{\text{source estimate}}.$$

Zero or missing source denominators give an unavailable percent difference
with an explanation. Include available uncertainty without forcing identical
interval definitions. Additional summary statistics are optional; avoid
exporting many tables that repeat the dashboard. Impose no arbitrary agreement
threshold or likelihood-equality requirement.

`mean_absolute_difference`, mean/median absolute percent differences, and the
terminal percent difference describe scale agreement. Absolute differences
retain the reported comparison unit and include valid zero source values;
percent differences exclude zero denominators. `trend_correlation` describes
trajectory agreement: annual metrics use correlation across common years;
age-specific metrics use the mean of within-age correlations across years.
Constant series have no defined correlation and are omitted from that mean.
This avoids treating the age gradient itself as evidence of temporal agreement.

Differences should be interpreted in light of documented translation choices, such as:

- aggregated fleets;
- simplified selectivity;
- different process structure;
- simplified q;
- simplified observation errors.

Do not interpret a difference as a model failure before checking whether it follows from an intentional translation choice.

------------------------------------------------------------------------

# 22. Scope discipline

Do not preemptively build exported package converters.

Prefer small analysis-local helpers that emerge from repeated transformations.

The shared database converters described here are justified by this multi-stock
workflow. Keep them analysis-local; the reusable Background argument belongs
in the package dashboard interface. No broader model mathematics or package
API redesign is required by this task.

Keep source curation, translation, and fitting responsibilities clear. New
code uses R; existing source importers may be retained separately for provenance
without becoming dependencies of the model-running loop.

Keep the source database independent of tinyAM.

Keep tinyAM translation transparent, reproducible, and scientifically explicit.
