# Assessment Database Curation Protocol

## Purpose

This document defines the standing protocol for adding, reviewing, and revising stock assessments in `analysis/comp_assessments/`.

The governing question is:

> **What data and assumptions entered the accepted assessment, and what did that assessment estimate?**

The objective is to represent the accepted model faithfully, not to collect every table in every report.

Follow the schema and controlled vocabularies in `DATABASE_STRUCTURE.md`.

Use this protocol to curate or repair source records. Follow
`TINYAM_TRANSLATION.md` for translating those records, fitting tinyAM models,
and comparing results. New extraction, validation, translation, and analysis
code should use R. Keep source-import scripts separate from the routine
database-to-model workflow; fitting must not depend on rerunning an importer
or reading source-cache files directly.

------------------------------------------------------------------------

# 1. Starting point: Charbonneau–Keith repository

The primary starting point for identifying age-structured assessments is the curated:

**Julie Charbonneau / David Keith — Age-structured-marine-fish-database**

Repository:

<https://github.com/JulieCharbonneau/Age-structured-marine-fish-database>

For the current curated analysis set, begin with:

<https://github.com/JulieCharbonneau/Age-structured-marine-fish-database/tree/main/PLOS_2026_analysis>

Use this repository as the initial sampling frame and stock catalogue.

For each candidate stock:

1.  locate the stock in the Charbonneau–Keith repository;
2.  record its Charbonneau identifier;
3.  inspect the metadata, notes, source assessment, and cited documents;
4.  trace the stock to the relevant current accepted assessment;
5.  retrieve the accepted assessment's authoritative inputs, outputs, and assumptions whenever possible.

Charbonneau–Keith is the default place to **start**, but not automatically the preferred source for every numerical value.

Its contents may include:

- extracted assessment data;
- fitted outputs;
- reconstructed or digitized quantities;
- transformed values;
- collapsed age classes;
- quantities derived from model results.

Do not assume a Charbonneau–Keith value is an original assessment input.

When a stock has received a newer accepted assessment than the one represented in Charbonneau–Keith, retain the Charbonneau stock identity as the starting point and curate the current accepted assessment unless the active task requests a historical assessment.

If a requested stock is not represented in the curated Charbonneau–Keith set, document that fact before adding it from another authoritative source.

------------------------------------------------------------------------

# 2. Source hierarchy

After identifying the stock, seek evidence approximately in this order:

1. native fitted assessment object;
2. native assessment input and output files;
3. official reproducible assessment repository;
4. official machine-readable assessment tables;
5. official assessment or framework tables and appendices;
6. documented reconstruction from authoritative source information;
7. Charbonneau–Keith extracted values when they are the best available source or a useful cross-check.

Prefer the original input structure actually consumed by the accepted model.

When native model files or fitted objects are publicly available, inspect them before relying on transcribed report tables.

## 2.1 Where to look for model objects, native files, and reproducible assessments

Before concluding that native model files are unavailable, search the main repositories and portals associated with the assessment authority and model family.

### ICES

Check the following resources:

- ICES Transparent Assessment Framework (TAF) GitHub organization:
  https://github.com/ices-taf

  TAF contains repositories for reproducible ICES assessments. Search using the ICES stock code, species, assessment year, and relevant expert-group acronym.

- ICES TAF portal:
  https://taf.ices.dk

- stockassessment.org:
  https://www.stockassessment.org

  This is particularly important for assessments fitted with SAM. When the short stock/case name is known, the R package `stockassessment` provides `fitfromweb()` for retrieving fitted SAM objects from stockassessment.org.

- `stockassessment` package documentation:
  https://fishfollower.r-universe.dev/stockassessment

  Use this to understand the contents of SAM objects and functions such as `fitfromweb()`, `read.data.files()`, `loadConf()`, and related extraction tools.

- ICES Scientific Reports / ICES Library:
  https://ices-library.figshare.com

  Search for the relevant assessment working group, benchmark workshop, stock code, species, and year. Working-group reports often contain detailed assessment methods, tables, annexes, and links or references to model files.

- ICES Stock Assessment Graphs:
  https://standardgraphs.ices.dk/stockList.aspx

  Useful for published stock-level graphs, tables, and downloadable assessment time series.

- ICES Advice and Scenarios Database:
  https://asd.ices.dk

  Useful for confirming current advice, stock codes, assessment year, advice category, and links to published advice.

Do not assume that the newest repository or stockassessment.org case is the accepted assessment. TAF and stockassessment.org may contain historical, test, benchmark, sensitivity, or development runs. Confirm the accepted assessment against ICES advice and working-group or benchmark documentation.

### Fisheries and Oceans Canada (DFO)

Check:

- Canadian Science Advisory Secretariat (CSAS):
  https://www.dfo-mpo.gc.ca/csas-sccs/index-eng.htm

- DFO science reports and publications:
  https://www.dfo-mpo.gc.ca/science/publications/reports-rapports-eng.html

Search for the species, stock area, assessment year, Research Document, Science Advisory Report, Science Response, and Proceedings.

CSAS Research Documents and Proceedings are particularly important because detailed assessment inputs, model descriptions, diagnostics, and appendices may appear there even when the Science Advisory Report contains only summary results.

Also follow any repositories, supplementary files, R packages, data archives, or model-object links cited in CSAS documents. DFO does not have a single public repository containing all assessment model objects, so stock-specific searches may be necessary.

### NAFO

Check:

- NAFO Library and Archives:
  https://www.nafo.int/Library/NAFOdocuments

- NAFO Scientific Council documents:
  https://www.nafo.int/Library/Science/SC-Documents/ev/1

- NAFO Scientific Council Reports:
  https://www.nafo.int/Library/Science-Council/Annual-Scientific-Council-Report-Compilations/ev/1

- NAFO stock advice:
  https://www.nafo.int/Science/Science-Advice/Stock-advice

Search by stock, NAFO Division/Subarea, assessment year, SCR document number, SCS document number, and Scientific Council meeting.

Scientific Council Research Documents (SCRs) can contain substantially more assessment detail than summary advice sheets, including model formulations, input tables, sensitivity analyses, and references to data or code.

Where a stock is assessed jointly or closely connected with ICES, also search the ICES Library and relevant joint working-group reports.

### NOAA Fisheries

Start with:

- NOAA Population Assessments resources:
  https://www.fisheries.noaa.gov/topic/population-assessments/resources

- Stock SMART:
  https://www.fisheries.noaa.gov/resource/tool-app/stock-smart

Stock SMART is useful for confirming assessment identity, year, stock status, reference points, and published output time series, but it will not usually contain the full native age-structured model input files.

For Alaska and North Pacific assessments, also search:

- Alaska Fisheries Science Center GitHub:
  https://github.com/noaa-afsc

Many AFSC stock-specific repositories contain assessment code, Stock Synthesis configurations, input files, output files, SAFE documents, and assessment run material.

For Northeast U.S. assessments, check:

- Northeast Stock Assessment Documents Search Tool:
  https://www.fisheries.noaa.gov/resource/publications-database/northeast-stock-assessment-documents-search-tool

- Stock Assessment Review Index (SARI):
  https://apps-nefsc.fisheries.noaa.gov/saw/sari.php

For Stock Synthesis assessments, search stock-specific repositories for files such as:

- `data.ss`;
- `control.ss`;
- `starter.ss`;
- `forecast.ss`;
- `Report.sso`;
- model run directories;
- associated R scripts and input tables.

NOAA GitHub repositories are distributed across programs and science centers rather than one universal stock-assessment repository. Search using the official stock name, assessment area, assessment year, model family, and science-center acronym.

## 2.2 Search strategy when the model object is not immediately obvious

For each assessment, search combinations of:

- official stock code;
- species name;
- stock area or management unit;
- assessment year;
- benchmark/framework year;
- expert-group or working-group acronym;
- model family, e.g. SAM, Stock Synthesis, WHAM, ASAP, VPA/ADAPT;
- phrases such as `assessment`, `benchmark`, `model`, `input`, `data`, `run`, `repository`, and `GitHub`.

When a report cites a model repository, supplemental dataset, package, or file archive, follow that trail before transcribing tables manually.

If multiple candidate model runs are found, determine which one corresponds to the accepted assessment by comparing:

- terminal data year;
- assessment year;
- model version;
- fleet and survey structure;
- age range;
- reported SSB, recruitment, and F trajectories;
- benchmark or working-group documentation.

Do not select a model object solely because it is the newest file or repository.

Use the most recent authoritative detailed source that documents the accepted production model. A detailed Research Document may be published a year or more after its Science Advisory Report or other summary product, so publication year and assessment year must not be treated as the same thing. Check the model version, terminal data year, and accepted run described by the detailed source, and record the assessment year and data terminal year separately from the source's publication date.

For example, if a later-published Research Document provides the inputs, assumptions, and outputs for the assessment summarized in an earlier advice document, use that Research Document for the detailed extraction. If no detailed source for the latest accepted assessment can be found, consult the most recent earlier detailed assessment and the relevant framework or benchmark document for context, but do not substitute earlier-run inputs or outputs for the accepted assessment or combine material from different runs. Record any unresolved gaps and keep completeness statuses partial where appropriate.

## 2.3 Source discipline

Do not:

- reconstruct an input from fitted outputs when the original input is available;
- manufacture missing values;
- assume a prominently reported quantity is necessarily a fitted-model input;
- assume that a public model object is the accepted assessment without verification;
- mix inputs or outputs from different assessment runs.

Keep raw observations, model inputs, transformed inputs, fitted values, outputs, and diagnostics conceptually distinct.

When the native assessment object or input files are found, preserve their provenance and use them as the preferred basis for populating the database. Report tables should generally be used as a cross-check or fallback when native inputs are available.

------------------------------------------------------------------------

# 3. Identify the accepted assessment and inventory the model

Before extracting numerical data, establish:

- stock identity and assessment authority;
- latest accepted assessment;
- terminal fitted data year;
- model family/version where known;
- model-defining benchmark or framework when distinct;
- accepted/base run when several runs are presented.

Do not infer accepted status from filenames alone.

Avoid test, sensitivity, preliminary, forecast-only, development, experimental, clone, and alternate runs unless official documentation identifies one as accepted.

If an annual assessment continues to use a model defined at an earlier benchmark/framework, record both assessment and framework years.

Then inventory the structure actually used by the accepted model:

- modeled years;
- modeled ages;
- recruitment age;
- terminal plus group;
- fishing fleets;
- survey/index series;
- sexes;
- spatial regions;
- seasons;
- landings, discards, and other removal components.

This inventory defines what must be curated.

A source table is material when it is:

1.  an input to the accepted model;
2.  an output from the accepted model; or
3.  needed to reconstruct an accepted-model input that is not directly available.

------------------------------------------------------------------------

# 4. Curate numerical inputs

## 4.1 Catch and removals

Capture every fitted catch/removal fleet separately.

Preserve the representation actually used by the model, including as applicable:

- numbers-at-age;
- biomass-at-age;
- total numbers;
- total biomass;
- proportions-at-age by number;
- proportions-at-age by biomass;
- landings;
- discards;
- other removals.

If the model uses a total plus an age composition, retain both.

Catch-at-age is the working quantity required for the tinyAM translation.
If a report supplies proportions describing an underlying catch-at-age input,
recover that input using matching totals and, when biomass is involved,
compatible catch weights; document the calculation. If the accepted model
instead fits totals and compositions separately, preserve those native inputs
and derive tinyAM catch-at-age outside the canonical database. Do not count
the native and reconstructed representations twice.

Total landings are otherwise optional unless the accepted model uses them or
they are needed for a documented conversion. Detailed gear splits are not
required unless they correspond to separate fitted fleets.

If the model combines several real-world fisheries into one modeled fleet, represent the modeled fleet rather than inventing finer fleet distinctions.

Capture catch weight-at-age when it is an explicit model input or needed to interpret a biomass-based catch input.

For catch compositions, determine whether values are:

- number proportions;
- biomass proportions;
- normalized frequencies;
- sample counts;
- another composition measure.

Do not infer the basis merely because values sum to one.

Capture composition weighting or effective sample size when it is a material model input.

## 4.2 Surveys and indices

Capture **every survey/index used by the accepted model**.

For each fitted survey, determine as applicable:

- survey identity;
- fitted spatial extent;
- years and ages used;
- native index measure;
- total abundance or biomass;
- direct abundance- or biomass-at-age;
- age composition and its basis;
- sex and season;
- sampling timing.

If the model uses an aggregate index plus an age composition, retain both.

Do not replace source totals/compositions with a derived age-specific abundance series unless that derived series was itself the accepted-model input.

Match the spatial scale used by the fitted model. If published material contains finer subareas than the model fits, do not treat them as separate fitted indices.

If the accepted model uses a combined index and only component values are published, reconstruct the combined input only when the aggregation is documented or defensible from authoritative information, and record the transformation.

For survey timing, prefer model configuration or documented timing. If only approximate timing can be recovered, document the limitation rather than inventing precision.

## 4.3 Biological inputs

Capture as applicable:

- stock weight-at-age;
- catch weight-at-age;
- survey-specific weight-at-age;
- maturity-at-age;
- natural mortality;
- spawning timing;
- fraction of fishing mortality before spawning;
- fraction of natural mortality before spawning.

Preserve annual variation when used by the accepted model.

If a model uses a constant vector, store the source vector rather than manufacturing repeated annual rows unless the source itself represents a complete matrix.

Keep stock weight and catch weight distinct.

For maturity, preserve the accepted sex convention rather than averaging sexes arbitrarily.

For natural mortality, distinguish:

- constant M;
- age-specific M;
- year-specific M;
- age-year M;
- internally estimated M.

Capture numerical M whenever it is a fixed or externally supplied input. If M is estimated internally, document the assumption and capture estimated M in `outputs.csv` where available.

------------------------------------------------------------------------

# 5. Curate model assumptions

Use `assumptions.csv` for the biological and statistical structure of the accepted model.

Investigate at least:

## Population and recruitment

- modeled years and ages;
- recruitment age;
- plus group;
- sex/spatial/seasonal structure;
- recruitment treatment;
- stock-recruit relationship where applicable;
- process-error structure;
- abundance/recruitment deviations and correlations.

## Fishing mortality

For each fitted fleet: - selectivity or F structure; - state sharing; - temporal process; - age correlation; - variance structure; - catch likelihood and observation variance.

## Natural mortality

- fixed versus estimated;
- age/time structure;
- process or estimation structure.

## Surveys and catchability

For every fitted survey: - sampling timing; - q/catchability structure; - q age blocks/sharing; - q time variation; - observation likelihood; - observation variance; - observation correlation.

## Biology

- stock and catch weights;
- maturity;
- spawning timing assumptions.

If an assumption cannot be resolved, record:

``` text
value = unknown
```

and explain the uncertainty.

Never infer statistical meaning from a field name, integer key, shorthand code, or parameter index without checking documentation or implementation.

------------------------------------------------------------------------

# 6. Curate outputs

Recover outputs needed to characterize the accepted assessment.

## High priority

- numbers-at-age;
- fishing mortality-at-age;
- natural mortality-at-age when estimated;
- spawning-stock biomass;
- recruitment.

## Secondary, where available

- total biomass;
- Fbar;
- Mbar when estimated;
- catchability estimates;
- predicted catch;
- predicted survey indices;
- Fmsy;
- Bmsy;
- MSY.

Preserve fleet, survey, sex, region, and season dimensions where the output is specific to them.

Retain uncertainty when practical.

Document the SE scale, interval level and type, and point-estimate convention
using the rules in `DATABASE_STRUCTURE.md`. A reported CV or log-scale SE
must not be presented as a natural-scale SE without a documented conversion.

Do not use fitted outputs as substitutes for source inputs.

------------------------------------------------------------------------

# 7. Transformations and reconstruction

The canonical database may standardize representation without changing scientific meaning.

Allowed examples include:

- wide-to-long reshaping;
- standardized variable names;
- explicit fleet/survey labels;
- standardized missing-value representation.

Do not silently:

- convert biomass to numbers;
- combine fleets or surveys;
- combine regions;
- change plus-group definitions;
- interpolate missing biological values;
- replace total-plus-composition inputs with derived age-specific values.

A scientific transformation is allowed only when needed to recover the **actual input used by the accepted assessment** and the model-ready input is not directly available.

In that case:

- use `source_type = reconstructed_source_input`;
- document the formula or aggregation in `transformation`;
- cite the source quantities;
- retain component source quantities when useful;
- do not present the reconstructed value as directly published.

This exception concerns reconstruction of the accepted assessment, not preparation for a downstream analysis.

These restrictions apply to canonical source records. Derived tinyAM
observations, settings, comparisons, and fitted results belong in analysis
outputs, as described in `TINYAM_TRANSLATION.md`.

------------------------------------------------------------------------

# 8. Completeness and source cache

An assessment is input-complete when all material input streams used by the accepted model have been represented or transparently reconstructed.

This includes, where applicable:

- all fitted fleets;
- all fitted surveys;
- biological inputs;
- natural mortality inputs;
- material sex, spatial, and seasonal structure.

Do not declare an assessment input-complete merely because one usable catch series and one survey have been found.

`outputs_status` and `assumptions_status` likewise refer to completeness relative to the accepted assessment.

## Stock source reviews

Keep stock-specific review notes in one Markdown file per stock under
`analysis/comp_assessments/source_reviews/`, named with the stable stock ID.
Use that file for accepted-run verification, source trails, input and output
coverage, transformations, validation findings, and unresolved gaps. When a
stock has several assessment records, identify the assessment IDs and keep
run-specific notes clearly separated within the stock's file.

Keep the README focused on the database's purpose, structure and usage, with
a compact index linking to the stock review files. Do not put detailed
stock-specific review notes in the README or duplicate them there. Move any
existing stock notes into the appropriate review file when revisiting a stock.
The database tables remain the source of numerical values, assumptions,
provenance and completeness statuses; review files explain the supporting
checks and limitations.

## Source cache

For every stock processed or revisited, cache the authoritative source files
used for curation locally under:

``` text
analysis/comp_assessments/source_cache/
```

This directory must remain gitignored.

Useful cached artifacts include assessment/framework PDFs, native input/output files, fitted objects, supplementary CSVs, and configuration files.

Organize files by stock and assessment run. Maintain a manifest identifying the
source URL, retrieval date, file checksum, and represented assessment. Record
failed or restricted downloads and any alternative source used. Keep useful
files already cached; do not download every document or unrelated material.

Use the cache to repair canonical records. The routine tinyAM translation and
fitting workflow reads the database, not native files from the cache.

Keep authoritative URLs in `assessments.csv`.

Do not mix artifacts from different assessment runs without preserving provenance.

------------------------------------------------------------------------

# 9. Per-assessment workflow and completion criteria

For every assessment:

1.  locate the stock in Charbonneau–Keith where applicable and record its identifier;
2.  verify the current accepted assessment and accepted/base run;
3.  identify the model-defining benchmark/framework;
4.  inventory modeled years, ages, fleets, surveys, sexes, regions, seasons, and plus group;
5.  enumerate all material input streams before extraction;
6.  curate catch/removal inputs for every fitted fleet;
7.  curate every survey/index used by the accepted model;
8.  curate biological and natural-mortality inputs;
9.  curate source-model assumptions;
10. curate high-priority outputs;
11. run structural validation;
12. compare represented fleets and surveys against the accepted-model inventory;
13. assign completeness statuses and record missing items;
14. commit the assessment as one coherent unit when the active task requires stock-by-stock commits.

Before declaring the assessment curated, report compactly:

``` text
Identity:
  stock
  Charbonneau ID
  assessment year
  terminal year
  model/framework/base run

Model dimensions:
  years
  ages
  recruitment age
  plus group
  fleets
  surveys
  sexes/regions/seasons

Inputs:
  fleet coverage
  survey coverage
  biology
  M

Outputs:
  N-at-age
  F-at-age
  M-at-age if estimated
  SSB
  recruitment
  other important outputs

Status:
  inputs_status
  outputs_status
  assumptions_status
  missing/unresolved items
```

Do not move on merely because one or two useful data streams have been recovered.

## Gaps discovered during translation

When a translator identifies a missing or inconsistent source input,
assumption, or output, return to this protocol and revisit that stock's
authoritative sources. Correct or complete the canonical records, update the
source review and relevant completeness statuses, run validation, and reload
the revised database before retrying translation. Record the database revision
used in the resulting analysis.

Do not change source values, drop required ages, or shorten the fitted period
simply to make tinyAM accept the data. Document a genuine unresolved source gap
and continue with other requested stocks when it cannot be filled.

Distinguish these gaps from unsupported source-model features. For example,
an assessment that estimates M may have no fixed numerical M input to recover.
Curate its actual assumptions and estimated outputs; choose any simplified
tinyAM M treatment outside the canonical tables.

Curation completion does not mean a tinyAM model has been fitted. When the
active task includes translation, follow `TINYAM_TRANSLATION.md` through
validation, a fitting attempt, convergence diagnostics, and comparison exports
for converged fits.

------------------------------------------------------------------------

# 10. Validation and final review

Run the structural validation defined in `DATABASE_STRUCTURE.md`.

Structural validation does not establish scientific completeness.

Also verify that:

- every expected fitted fleet is represented;
- every expected fitted survey is represented;
- spatial/sex/seasonal dimensions match the accepted model;
- outputs and assumptions are not being substituted for inputs;
- reconstructed inputs are explicitly identified and documented.

Prefer:

> incomplete but correct

over:

> complete but guessed.

When information is unclear, inspect native files, fitted objects, framework documents, appendices, official repositories, and Charbonneau–Keith source notes before recording `unknown`.

------------------------------------------------------------------------

# 11. Scope discipline

Do not:

- turn the database into a new R package;
- build database-server infrastructure;
- preemptively build generic model-format converters;
- archive every document;
- collect unrelated descriptive tables;
- aggregate away fitted fleets or surveys in canonical inputs;
- create derived series solely for a downstream analysis;
- alter source-assessment values to make them easier to use elsewhere.

Reusable source-specific helpers may be added later when repeated real use justifies them.

------------------------------------------------------------------------

# 12. Agent task completion

The protocol does not define how many assessments to add or review. The active task does.

For tasks that also require tinyAM translation, report curation progress and
modeling progress separately: translated, blocked, fit failed or did not
converge, and converged with comparison outputs. Keep these analysis statuses
outside the canonical completeness fields. An exported observation table or
settings template alone does not complete a requested fitting task.

If the goal says:

> Add 10 assessments

the task is complete only after 10 assessments have been added according to this protocol, unless an external blocker requires user intervention.

Completing a representative subset does not satisfy a numeric target.

For long-running tasks, report progress explicitly:

``` text
Progress: 4 / 10 assessments completed.
```

If one requested assessment is blocked by unavailable authoritative source material, document the blocker and continue with the remaining requested assessments unless user guidance is required before substitution.

A final multi-assessment report should state:

- target number;
- number complete;
- number partial;
- number blocked;
- assessment IDs;
- commits created;
- major missing data or unresolved assumptions.

Do not stop merely because the workflow has been demonstrated successfully on a subset.
