# Irish Sea cod source review

Source review: 2026-10-19.
Assessment record: `ices_cod_irish_sea_2026`.

The source record is partial and is not fit-ready for a tinyAM translation.

## Identity and official run

Charbonneau–Keith identifier: `ICES-WGCSE_morhua_IS7a`. The current ICES stock code is `cod.27.7a`, Atlantic cod in Division 7.a (Irish Sea).

The current applied assessment is the 2026 WGCSE run, assessment key 23721, with data through 2025 and advice for 2027–2028. The report says the final assessment used Stock Synthesis 3 (SS3). The latest stock annex says the stock was last benchmarked in February 2022 and was updated in May 2025. Native 2026 SS3 data and control files were not found.

The working-group chapter is divided into separate PDFs on Figshare. The cod chapter and the 2025 stock annex are cached locally with retrieval URLs, dates and checksums under `source_cache/ck_ices_wgcse_morhua_is7a/assessment_2026/`.

## Model dimensions and inputs

The assessment spans 1968–2025, ages 0–6+, with recruitment at age 0 and Fbar over ages 2–4. There is one Irish Sea population and a single commercial fleet; recreational removals are excluded.

Table 6.9 reports total catch numbers-at-age in thousands of fish. It marks 2003–2005 and 2020 as excluded from the assessment because sampling was limited. Those source values remain in `inputs.csv` with the exclusion noted; the translation must omit them from the catch likelihood. Table 6.11 separately reports discards-at-age in thousands. Discards before 2007 are reconstructed in the source, and no 2020 age-discard values are tabulated.

Table 6.10 reports mean weights-at-age in landings, used for the whole stock and catch. The age-0 value is missing for 2019. In 2020 the report marks the weights as available only from Q1 sampling and not included; 2023–2025 have low or incomplete age sampling, and the 2025 ages 4+ values use the preceding three-year average. These source limitations are retained rather than filled.

Table 6.12 reports annual maturity from 1996 through 2025. Ages below 2 are immature and ages 3+ are mature. The 2025 stock annex reports a constant pre-1996 age-2 value of 0.27; the 2026 table reports 0.292 for 1996. The ogive is held constant from 2025 onward. The stock annex reports fixed, time-constant natural mortality of 1.7, 0.714, 0.68 and 0.65 for ages 0, 1, 2 and 3+, respectively.

The report provides age-specific index series for NIGFS-WIBTS-Q1 (ages 1–4), NIGFS-WIBTS-Q4 (age 0) and UK FSP (ages 2–6+). The working-group text says FSP spans 2004–2025, while Table 6.13 values begin in 2005; the database follows the numeric table. Table 6.13 gives one annual CV beside the Q1 age-index vector; the canonical long table repeats that value across the four ages as a source-table expansion, not as evidence of separately estimated age-specific CVs. Q1 CPUE and its standard deviation are also reported separately, but their native units are not specified. The age-index scale multipliers are not stated either.

The 2026 report and stock annex summarize model structure, selectivity, survey ages and some effective-sample-size settings. The native SS3 data/control files were not recovered, so exact likelihood definitions, observation weights and current parameter activation remain uncertain. These assumptions are explicitly marked partial.

## Outputs and readiness

Table 6.14 reports annual SSB, age-0 recruitment and Fbar ages 2–4 with pointwise 95% confidence limits. The reported summary values were cross-checked against the official ICES SAG export for key 23721, within report rounding. The source report and SAG export do not provide recoverable age-specific N-at-age or F-at-age; M is fixed rather than estimated.

The source table contains one missing cell in the weight surface (2019 age 0). No authoritative value for that cell was found, and it remains missing in the database. Rather than interpolate it, the translation uses ages 1–6+, which have a complete annual weight and maturity surface. This is a limited approximation: it drops age-0 catch and the Q4 age-0 index, and tinyAM's age-1 recruitment is not comparable with the official age-0 recruitment.

The translation retains catch-at-age in thousand fish except for the assessment-excluded years 2003–2005 and 2020. It also uses the Q1 index at ages 1–4 and UK FSP Q1 at ages 2–6+, with both sampled at 0.125 as a first-quarter timing approximation. The age-0-only Q4 index is omitted. Unknown index scale multipliers are left unchanged and absorbed by freely estimated q; the source's Q1 CV and FSP effective sample size are not imposed because their current likelihood allocation is unverified. Observation SDs are estimated separately by survey.

The available outputs support an approximate aggregate SSB comparison only. Source SSB timing is not stated, while tinyAM calculates beginning-of-year mature biomass. Age-specific N and F are not available. Source Fbar is ages 2–4, but tinyAM reports an abundance-weighted mean and source N-at-age is unavailable, so the two Fbar values are not treated as comparable. Natural mortality is fixed at the published values and can be displayed as an input, not as an estimated output.

The initial age-1+ fit converged from package defaults (optimizer code 0, maximum gradient 0.000335, positive-definite Hessian). Its approximate SSB comparison was poor: mean absolute percent difference 55.5%, terminal difference -99.3%, and trend correlation 0.92. The model omitted the official age-0 index and did not represent age-0 recruitment; the resulting recruitment mismatch makes this full-period age-1+ candidate unsuitable as the retained comparison.

A complete age-0–6+ weight surface is available through 2018. The next candidate will retain the official age-0 recruitment definition and all three age-specific index series, but limit the fitted period to 1968–2018 rather than fill the missing 2019 age-0 weight or use unassessed 2020 weights. This is a shorter, more definition-aligned comparison; it will not represent the terminal 2019–2025 assessment period. The current candidate script and structural test are in `scripts/translation/stocks/ices_cod_irish_sea_2026.R` and `tests/test_irish_sea_cod_translation.R`; the source tables remain validated by `tests/test_irish_sea_cod_source.R`.
