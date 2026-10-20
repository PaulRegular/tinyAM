# Gulf of Maine haddock

## Stock identity and assessment selection

The metadata record is `ck_nefsc_garmiii_aeglefinus_5y` (Charbonneau–Keith identifier `NEFSC-GARMIII_aeglefinus_5Y`), Gulf of Maine haddock (*Melanogrammus aeglefinus*). It is distinct from Georges Bank haddock, which has its own stock record and assessment area.

The 2024 NOAA Management Track assessment updated data through 2023 and was the first Gulf of Maine haddock Management Track run using WHAM. The September 2024 peer-review panel concluded that the assessment was technically sufficient and represented best scientific information available for management. NOAA held a later Gulf of Maine haddock Management Track review on September 22–23, 2026. The 2026 source bundle currently cached for this review contains a stock assessment report, WHAM diagnostic PDFs and figures, and a presentation, but no numerical model inputs or fitted age-specific result tables.

The 2024 run is therefore the latest assessment found with detailed assessment materials suitable for source review. This does not mean its input series are fully recovered: the public 2024 bundle contains graphical input summaries rather than the numerical observations needed to reproduce the fitted model. Keep any future 2026 summary-only record separate from the 2024 run.

## 2024 materials and recoverable values

The official 2024 bundle includes the assessment report, `Data.zip`, and `Models.zip`. `Data.zip` contains 36 PNGs, three GIFs, and two JPGs; these show commercial and recreational catch-at-age, spring and fall NEFSC bottom-trawl survey information, bottom-longline survey information, age-specific weights, and maturity. `Models.zip` contains 243 PNG diagnostic and output plots. Neither archive contains machine-readable age-composition, index, biological-input, or fitted-model tables. The figures identify model streams and coverage, but digitizing their marks would create approximate values and would not recover the original observation weights or full precision. These plotted values are not entered as exact inputs.

The 2024 report does provide a numerical summary table for 2014–2023. It reports total fishery catch (metric tonnes), spawning-stock biomass (metric tonnes), full F, and age-1 recruitment (thousand fish):

| Year | Catch (t) | SSB (t) | Full F | Age-1 recruitment (thousand fish) |
|---:|---:|---:|---:|---:|
| 2014 | 1,410 | 37,765 | 0.06 | 55,825 |
| 2015 | 1,512 | 65,183 | 0.04 | 11,244 |
| 2016 | 3,406 | 86,275 | 0.06 | 12,355 |
| 2017 | 4,894 | 87,550 | 0.08 | 13,367 |
| 2018 | 3,710 | 107,367 | 0.05 | 8,639 |
| 2019 | 4,542 | 81,012 | 0.06 | 6,067 |
| 2020 | 4,810 | 51,247 | 0.11 | 7,559 |
| 2021 | 4,618 | 22,797 | 0.24 | 27,946 |
| 2022 | 3,353 | 25,766 | 0.17 | 25,791 |
| 2023 | 2,734 | 17,836 | 0.23 | 4,974 |

These are assessment summary outputs, not substitutes for the missing fleet-specific removals, catch and survey age compositions, biological matrices, or N-at-age and F-at-age surfaces. The report also gives short-term projection summaries for 2024–2027; they are projections and must remain distinct from fitted historical outputs.

## Source search and local cache

The 2024 and 2026 NOAA Stock Assessment Information System (SASI) listings and their downloadable bundles are cached under `source_cache/nefsc_haddock_gulf_of_maine_2026/`. The local `manifest.csv` records URLs, verification date, file names, and SHA-256 checksums. The 2024 archive's extracted report, data figures, and model figures are retained beside the original archive. The cached 2026 bundle contains four PDFs and two PNGs, with no CSV, R data object, native WHAM input, or numeric age-specific result file found.

The 2024 report and review record establish the assessment year, terminal year, model family, fleets, surveys, and summary values. The graphical age compositions and index plots do not expose the numerical sample sizes, observation weights, or original index scale needed for a faithful tinyAM translation. Do not combine the 2022 Research Track framework data with the 2024 Management Track production run: the framework's historical inputs and estimates are not the updated 2024 observations or fit.

## Translation status

No assessment, input, assumption, or output rows have yet been added for this stock. A limited output-only import could preserve the 2024 summary table, but it would not make the case ready for `fit_tam()`. Continue searching NOAA's SASI portal, assessment repositories, and linked data systems for the numeric 2024 inputs and fitted outputs. If no numerical source becomes available, keep the source limitation explicit and do not manufacture observations by digitizing the charts.

Sources:

- [2024 Gulf of Maine haddock assessment report](https://d23h0vhsm26o6d.cloudfront.net/a.iv_Gulf_of_Maine_haddock_Update_2024_09_03_193406.805521.pdf)
- [2024 September Management Track peer-review panel report](https://www.fisheries.noaa.gov/s3/2024-10/2024-Management-Track-Peer-Review-Panel-Report_10_8_24_508.pdf)
- [2026 September Management Track peer-review event](https://www.fisheries.noaa.gov/event/peer-review-2026-september-management-track-assessments)
- [NOAA Stock Assessment Information System](https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php)
