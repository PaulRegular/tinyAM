# North Sea / Northern Shelf cod source review

Review date: 2026-10-02. A source-reviewed partial record has been added; this is not an input-completeness declaration.

## Identity and accepted run

Charbonneau identifier: `ICES-WGNSSK_NS 4-7d,20_Gadus_morhua`. The catalogue uses the 2020 report, with data ending in 2019. Its old stock boundary must not be carried forward unchanged: the 2023 benchmark combined North Sea and West of Scotland cod in a three-substock Northern Shelf assessment (`cod.27.46a7d20`).

The latest cod advice found in the ICES Library search is the September 2025 advice, DOI `10.17895/ices.advice.27202566`. The ASD lists that assessment as active (combined key 19639; Northwestern 19661; Southern 19662; Viking 19663). The June 2026 WGNSSK report states that no cod advice was issued in spring. Search for an autumn 2026 update again before finalizing current status.

Detailed accepted-run source: WGNSSK 2025 autumn cod chapter, https://ndownloader.figshare.com/files/59358443 (ICES Library article 29085995). Framework: WKBCOD 2023, https://doi.org/10.17895/ices.pub.22591423. Stock annex revised in 2024: https://ndownloader.figshare.com/files/47179495.

## Model inventory

The final assessment uses stockassessment 0.12.4 and multiStockassessment 0.4.0. Population surfaces cover 1983–2025, ages 1–7+, for Northwestern, Southern and Viking substocks; catch data end in 2024. Seven fitted index streams are defined in Table 4.7c: a mixed Q3+Q4 index, three substock Q1 indices, and three recruitment indices. The latter shift age-0 Q3+Q4 observations to age 1 at the following 1 January. Preserve that model-year meaning and do not treat the seven underlying field surveys as seven independent fitted indices.

Material input streams:

- Combined commercial catch numbers-at-age, Table 4.2c; catch mean weights, Table 4.4c; landings fractions and component weights where required by auxiliary likelihoods.
- Q1 indices and their supplied uncertainty, Table 4.6a, separately for each substock.
- Mixed Q3+Q4 indices and uncertainty, Table 4.6b.
- Recruitment indices and uncertainty, Table 4.6c, separately for each substock.
- Q1 substock landings proportions, section 4.2.1.3 (1995–2024).
- Quarterly total landings proportions, section 4.2.1.3 (1995–2024).
- Substock stock-weight and maturity observations, Tables 4.5a–b.
- Supplied natural-mortality observations, Table 4.5c, common to the substocks, with missing recent years retained.
- Spawning fractions are zero in the official setup (prop.f=NULL and prop.m=NULL). Exact model timing from the official index export is 0.125 for Q1, 0.75 for Q3+Q4 and 0 for the forward-shifted recruitment indices.

Recruitment and F follow random walks. F increments have AR1 correlation across ages; substock F is linked through scaled selectivity with flexible Fbar. Catch and age-based survey observations have AR1 correlation across ages. Substock and quarterly composition observations use additive logistic-normal likelihoods. Stock weights, maturity and M are fitted GMRF processes. Consequently Tables 4.11–13 contain biological estimates, not original inputs. Recreational removals are explicitly excluded from the fitted assessment.

## Numerical output source found

The official ICES technical service https://doi.org/10.17895/ices.advice.32541738 releases FLStock objects from accepted 2025/2026 assessments. The updated archive is https://ndownloader.figshare.com/files/67087409. `Selectivity Indicators/FLStocks.RData` contains `cod.27.46a7d20` as a list of three FLStock objects, labelled Northwest, South and Viking, covering 1983–2025 and ages 1–7+.

Provenance check: Northwestern 1983 N at age starts 485859.551, 122291.133, 22759.112 thousand fish, matching the rounded Table 4.9 values 485860, 122291, 22759. The complete N/F/M surfaces have now been checked against Tables 4.8, 4.9 and 4.13 before importing. These objects explicitly contain assessment outputs; their catch and biological slots must not be used as observations without independent verification. Obtain uncertainty from the official SAG XML or Table 4.14 rather than inventing standard errors.

## Source search and remaining gaps

The official ices-advice organization contains https://github.com/ices-advice/2025_cod.27.46a7d20, pinned at 390ff07e3d1de2e14b3a50b1ac27260566d73b62. It supplies the assessment/data-preparation scripts but no raw input files, native fitted object, configuration file or release assets. DATA.bib refers to local files, despite labelling access Public. The cited mixed-fisheries repository link returned HTTP 404. The TAF organization catalogue, stockassessment.org directory and stock-annex references were checked. The public user3 SAM archive contains the historical single-stock `WKCOD_combined_99` benchmark case; it is not this accepted multistock assessment. Published input tables are retained with explicit report-rounding provenance. Both auxiliary composition streams and supplied survey uncertainty are included. Unrounded auxiliary values cannot be recovered from the public files; no missing values are filled or rounded proportions renormalized.

Catalogue correction to keep in mind: `NEFSC-GARMIII_5Y_Melanogrammus_Aeglefinus` identifies Gulf of Maine haddock, not Georges Bank haddock. Do not use the current candidate-map label as evidence of stock identity. The requested NEFSC choice remains Georges Bank haddock or Gulf of Maine cod, pending source review.

Numerical verification completed: all 2,688 tabulated N/F/M values across the three substocks match the official FLStock objects within the report's rounding (882 F values for 1983-2024; 903 N values and 903 M values for 1983-2025). Zero mismatches. This supports output provenance but does not establish input completeness. The technical-service DOI is 10.17895/ices.advice.32541738.v2.

## Recorded coverage and validation

Assessment ID: `ices_cod_north_sea_2025`. Last fitted input year is 2025 (Q1 survey); catch and reported F end in 2024. N and M estimates extend through 2025. The compact initial stock inventory is three substocks, one combined commercial catch fleet, seven fitted index streams and two auxiliary composition likelihoods.

The record has 6,515 inputs, 3,327 outputs and 75 assumptions. The input inventory includes 294 combined catch age values; 294 landings component counts and 294 reconstructed landings fractions; 882 catch/landings/discard weights; 894 stock weights (nine source missing cells remain omitted); 903 maturity values; 280 supplied M values (2023-2025 remain missing); 1,232 indices with matching log-SD records (Q34 age 7 in 2004 remains missing); and 210 quarterly/substock landings-weight proportions.

All three N/F/M output surfaces match the report within rounding. Recruitment matches age-1 N; reported Fbar matches the mean of ages 2-4. Structural validation and the source-specific coverage checks pass. A rendered Q1 table was reviewed to confirm age-value and SD column alignment.

The assessment assumptions now record parameter-sharing keys, independent abundance innovations, recruitment and F random walks, AR1 observation correlations, log-index relative weights, and GMRF biological observation treatment. The cited version 0.4.0 implementation (revision 0465d228884e0e1fe394276ef320e2d8f20ac16c dated 2025-04-22) resolves initN=2 as separate recruitment-level parameters and initial age recursion under first-year F+M, with log SD 0.01. Shared selectivity code 4 combines a cubic age effect and a scalar stationary AR1 temporal level; source comments calling that level RW do not match the likelihood equations.

Statuses remain partial: exact native inputs and installed code revision are not publicly supplied, numerical catchability estimates and detailed prediction uncertainty have not been recovered, and the biological-process mean/variance specification is not fully recorded. No values from the historical WKCOD_combined_99 fit or 2023 assessment are transferred into this record.

## Record summary and reproduction

The North Sea cod catalogue entry is represented by the accepted three-substock
Northern Shelf assessment after the 2023 benchmark. Original biological
observations remain separate from fitted biological surfaces. All seven index
streams retain their supplied log-scale SDs; model timing is 0.125 for Q1, 0.75
for Q3+Q4 and 0 for forward-shifted recruitment indices. Fitted inputs extend
through the 2025 Q1 survey, while catches and reported F end in 2024.

Cached sources remain gitignored. The
source-specific importer requires Python with pdfplumber and base R:

```sh
python analysis/comp_assessments/scripts/database/004_import_north_sea_cod.py . --rscript Rscript
Rscript analysis/comp_assessments/scripts/database/002_validate_database.R
Rscript analysis/comp_assessments/scripts/database/005_validate_north_sea_cod.R
```

The importer checks all 2,688 native N/F/M values against the detailed report
and can be rerun without changing other assessments. Completeness statuses
remain partial; this record does not borrow observations from older runs.
