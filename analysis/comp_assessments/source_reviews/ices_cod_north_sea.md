# North Sea / Northern Shelf cod source review

Source verification: 2026-10-03.  
Translation revision: 2026-10-14.

The source-reviewed record remains partial; this is not an input-completeness declaration.

## Identity and official run

Charbonneau identifier: `ICES-WGNSSK_NS 4-7d,20_Gadus_morhua`.

The catalogue uses the 2020 report, with data ending in 2019. Its historical stock boundary must not be carried forward unchanged: the 2023 benchmark combined North Sea and West of Scotland cod in a three-substock Northern Shelf assessment (`cod.27.46a7d20`).

As checked on 2026-10-03, the ICES Advice and Scenarios Database lists `cod.27.46a7d20` with assessment year 2025 and status Active. Its advice is DOI `10.17895/ices.advice.27202566`. Assessment components are the combined stock (key 19639), Northwestern (19661), Southern (19662), and Viking (19663). The advice applies to 2026 catches.

The June 2026 WGNSSK report states that no cod advice was issued that spring. No later active cod assessment was listed in the ASD at the source-review date.

Detailed official-run source: [WGNSSK 2025 autumn cod chapter](https://ndownloader.figshare.com/files/59358443), ICES Library article 29085995.

Framework: [WKBCOD 2023](https://doi.org/10.17895/ices.pub.22591423).

[Stock annex revised in 2024](https://ndownloader.figshare.com/files/47179495).

## Official assessment: model structure and inputs

The official assessment uses `stockassessment` 0.12.4 and `multiStockassessment` 0.4.0.

The population model covers 1983–2025, ages 1–7+, for Northwestern, Southern, and Viking substocks. Commercial catch data end in 2024.

Seven fitted index streams are defined in Table 4.7c: one mixed Q3+Q4 index, three substock Q1 indices, and three substock recruitment indices.

The recruitment indices shift age-0 Q3+Q4 observations to age 1 at the following 1 January. This model-year convention must be preserved. The seven underlying field surveys should not be interpreted as seven independent fitted index streams.

### Material inputs

- Combined commercial catch numbers-at-age, Table 4.2c; catch mean weights, Table 4.4c; landings fractions and component weights where required by auxiliary likelihoods.
- Q1 indices and supplied uncertainties, Table 4.6a, separately for each substock.
- Mixed Q3+Q4 indices and uncertainties, Table 4.6b.
- Recruitment indices and uncertainties, Table 4.6c, separately for each substock.
- Q1 substock landings proportions, section 4.2.1.3 (1995–2024).
- Quarterly total landings proportions, section 4.2.1.3 (1995–2024).
- Substock stock-weight and maturity observations, Tables 4.5a–b.
- Supplied natural-mortality observations, Table 4.5c, common to the substocks, with missing recent years retained.

Spawning fractions are zero in the official setup (`prop.f = NULL`, `prop.m = NULL`). Model sampling times are 0.125 for Q1, 0.75 for Q3+Q4, and 0 for forward-shifted recruitment indices.

Recruitment and fishing mortality use temporally structured processes. The F specification includes cross-age dependence and shared, scaled selectivity among substocks. Catch and age-based survey observation errors have AR(1) correlation across ages.

Substock and quarterly composition observations use additive logistic-normal likelihoods. Stock weights, maturity, and natural mortality are fitted GMRF processes. Consequently, Tables 4.11–13 contain biological estimates rather than original inputs.

Recreational removals are excluded from the fitted assessment.

## tinyAM translation revision (2026-10-14)

The translation was revised to better represent the population and fishing-mortality processes while constraining survey catchability to a biologically plausible, non-decreasing age pattern.

The translation remains a single-stock approximation of the official three-substock assessment.

### Population structure

The model fits years 1983–2022 and ages 1–7+, preserving the source population age range.

The earlier translation used random-walk cohort-survival deviations. The revised translation specifies:

`N_settings = list(process = "iid", init = "free")`

This retains a separate random-walk recruitment process under the default recruitment settings while using independent cohort-survival deviations for older ages.

Initial abundance at older ages is freely estimated.

Unlike the official assessment, tinyAM uses one combined population and one shared SD for older-cohort survival deviations. Substock-specific population processes and relationships among substocks are not represented.

### Fishing mortality

The previous translation used independent fishing-mortality random walks. The revised configuration is:

`F_settings = list(process = "cor_rw", mu_form = NULL)`

This allows AR(1) correlation across ages in annual log-F increments, rather than treating changes in F at neighbouring ages as independent.

The revised process is closer to the correlated age structure of the official model. However, it does not reproduce the complete multistock F parameterization or the shared, scaled selectivity relationships among Northwestern, Southern, and Viking substocks.

The official implementation's shared selectivity code 4 combines a cubic age effect with a scalar stationary AR1 temporal level. That structure should not be described as identical to tinyAM's `cor_rw` process.

### Survey selection

Only the mixed Q3+Q4 survey (`Survey_Q34`) is fitted in tinyAM.

This index represents the combined population and is therefore more directly applicable to the single-stock approximation than the substock-specific Q1 and recruitment indices.

Its sampling time of 0.75 and supplied relative observation SDs are retained.

The other six fitted index streams are omitted. The resulting observation likelihood contains less information than the official assessment, particularly about recruitment and differences among substocks.

### Monotone survey catchability

The revised survey catchability specification is:

`index_settings = list(q_form = ~mono(q_key), sd_form = ~1, sd_supplied = ~relative_sd)`

The `mono(q_key)` term imposes **non-decreasing catchability across ordered survey age groups**.

tinyAM estimates a baseline catchability and cumulative non-negative increments between successive age groups on the log-catchability scale. Increments can equal zero, allowing adjacent ages to share catchability.

This specification permits increasing or flat catchability but excludes dome-shaped relationships in which catchability declines at older ages.

The constraint is intentional. With only the mixed Q3+Q4 survey and a combined population, unconstrained age-specific catchability may be weakly identified and could produce a dome-shaped pattern that is difficult to distinguish from changes in estimated abundance-at-age.

The monotone curve reduces this flexibility without requiring a specific logistic function or assuming that catchability must strictly increase with age.

**This is a translation assumption, not a recovered constraint from the official assessment.** The official multistock model uses a different observation and population structure, and its numerical catchability estimates have not been independently recovered.

The ordering of `q_key` factor levels must correspond to increasing survey age groups, because `mono()` follows the declared factor-level order.

### Natural mortality

The translation specifies:

`M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~M_assumption)`

The supplied year-age natural-mortality surface is treated as fixed. This motivates restricting the fit to 1983–2022, when the translated M surface is complete.

The official model estimates M through a biological GMRF. That estimated variation and its uncertainty are not reproduced by tinyAM.

### Commercial catch

The combined commercial catch numbers-at-age are retained.

The translation uses independent lognormal catch errors with a single estimated observation SD:

`catch_settings = list(sd_form = ~1)`

The source's age-correlated catch observation errors, substock landings compositions, and quarterly composition likelihoods are omitted.

### Biological inputs

Official weight and maturity observations are substock-specific.

The translation uses the arithmetic mean of the available Northwestern, Southern, and Viking values for each year and age.

This provides one biological surface for the combined population, but the resulting weights and maturity proportions are not weighted by the relative abundances of the substocks.

The method is therefore a practical approximation, rather than an exact aggregation of the official biological inputs.

### Comparison quantities

The official assessment reports outputs separately for each substock. These are combined for comparison with the single-population tinyAM fit.

Population numbers, recruitment, and biomass are summed across substocks.

Age-specific fishing and natural mortality are averaged across substocks using estimated abundance-at-age as weights.

Substock Fbar estimates are combined using their relative abundance over ages 2–4.

Uncertainty intervals are not combined.

These derived quantities allow comparisons of broad population dynamics but should not be interpreted as exact equivalents of the quantities estimated within the original three-substock model.

## Numerical output source and verification

The official [ICES technical service](https://doi.org/10.17895/ices.advice.32541738) releases FLStock objects from official assessments. The [updated archive](https://ndownloader.figshare.com/files/67087409) contains `Selectivity Indicators/FLStocks.RData`.

The stock `cod.27.46a7d20` is represented as three FLStock objects, labelled Northwest, South, and Viking, covering 1983–2025 and ages 1–7+.

Provenance check: Northwestern 1983 N-at-age begins 485859.551, 122291.133, and 22759.112 thousand fish, matching the rounded Table 4.9 values of 485860, 122291, and 22759.

The complete N, F, and M surfaces have been checked against Tables 4.8, 4.9, and 4.13 before importing.

All 2,688 tabulated N/F/M values across the three substocks match the official FLStock objects within report rounding:

- 882 F values for 1983–2024.
- 903 N values for 1983–2025.
- 903 M values for 1983–2025.

Zero mismatches were identified.

These objects contain official assessment outputs. Their biological or catch slots must not be treated as original observations without independent verification.

Uncertainty should be obtained from the official SAG XML or Table 4.14 rather than inferred from rounded estimates.

Technical-service DOI: `10.17895/ices.advice.32541738.v2`.

## Source search and remaining gaps

The official ICES-advice repository is [2025_cod.27.46a7d20](https://github.com/ices-advice/2025_cod.27.46a7d20), pinned at revision `390ff07e3d1de2e14b3a50b1ac27260566d73b62`.

It provides assessment and data-preparation scripts but not the complete raw input files, native fitted object, configuration files, or release assets.

`DATA.bib` refers to local files despite labelling access Public. The cited mixed-fisheries repository link returned HTTP 404 during the source review.

The TAF catalogue, stockassessment.org directory, and stock-annex references were checked.

The public user3 SAM archive contains the historical single-stock `WKCOD_combined_99` benchmark case, which is not the official multistock assessment.

Published input tables are retained with their report-rounding provenance. Auxiliary composition streams and supplied survey uncertainties are included in the canonical database.

Unrounded auxiliary values could not be recovered from the public files. Missing values are not filled and rounded proportions are not renormalized.

The official configuration's version 0.4.0 implementation, revision `0465d228884e0e1fe394276ef320e2d8f20ac16c` dated 2025-04-22, resolves `initN = 2` as separate recruitment-level parameters and initial age recursion under first-year F+M, with log SD 0.01.

Shared selectivity code 4 combines a cubic age effect and a scalar stationary AR1 temporal level. Source comments describing that temporal level as a random walk do not match the implementation's likelihood equations.

The native installed code revision, complete input data, numerical catchability estimates, and detailed prediction uncertainty remain unavailable. The biological-process mean and variance specifications are not fully recorded.

No input values are borrowed from the historical `WKCOD_combined_99` fit or the 2023 assessment.

## Recorded coverage and validation

Assessment ID: `ices_cod_north_sea_2025`.

The last fitted input year is 2025, corresponding to the Q1 survey. Catch and reported F end in 2024; N and M estimates extend through 2025.

The official record includes three substocks, one combined commercial catch fleet, seven fitted index streams, and two auxiliary composition likelihoods.

The canonical database contains 6,515 inputs, 3,327 outputs, and 75 assumptions.

The inventory includes:

- 294 combined catch-at-age values.
- 294 landings component counts and 294 reconstructed landings fractions.
- 882 catch, landings, and discard weights.
- 894 stock weights, with nine source missing cells omitted.
- 903 maturity values.
- 280 supplied M values, with 2023–2025 missing.
- 1,232 index observations with corresponding log-SD records, including a missing Q34 age-7 observation in 2004.
- 210 quarterly/substock landings-weight proportions.

Source-specific structural validation and output coverage checks pass. A rendered Q1 table was reviewed to confirm age-value and SD column alignment.

The assessment assumptions record parameter sharing, abundance innovations, recruitment and F processes, observation correlations, relative index weights, and biological GMRF treatments.

The statuses remain partial because the complete native inputs and fitted model have not been recovered.

## Record summary and reproduction

The North Sea cod catalogue entry now represents the official three-substock Northern Shelf assessment following the 2023 benchmark.

The canonical database preserves the official substock inputs, fitted outputs, survey age structures, observation uncertainties, and biological-process information without substituting inputs from earlier assessment models.

The revised tinyAM translation fits one combined population over 1983–2022, with IID older-cohort survival deviations, correlated F random walks, fixed supplied M, and monotone non-decreasing survey catchability.

The `mono(q_key)` constraint is particularly important: it prevents the reduced observation model from estimating declining catchability at older ages, while allowing flexible increases and plateaus.

The comparison evaluates whether this simpler configuration can approximate broad official stock dynamics. It does not reproduce the official three-substock likelihood, selectivity model, or full uncertainty.

Cached sources remain gitignored. The source-specific importer requires Python with pdfplumber and base R:

```sh
python analysis/comp_assessments/scripts/database/004_import_north_sea_cod.py . --rscript Rscript
Rscript analysis/comp_assessments/scripts/database/002_validate_database.R
Rscript analysis/comp_assessments/scripts/database/005_validate_north_sea_cod.R
```

The importer validates all 2,688 native N/F/M values against the detailed report and can be rerun independently of other assessments. Input completeness remains partial.
