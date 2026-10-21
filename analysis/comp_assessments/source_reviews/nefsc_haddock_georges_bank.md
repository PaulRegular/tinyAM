# Georges Bank haddock source review

## Current official assessment

NOAA's June 2026 Management Track Peer Review Panel accepted the Georges Bank haddock WHAM assessment as Best Scientific Information Available. It updates the 2024 management track assessment under the 2021–2022 research track framework, with data through 2025. This is the full Georges Bank assessment, distinct from DFO's Eastern Georges Bank assessment.

The official input inventory includes US commercial landings/discards and Canadian catches, fishery catch-at-age and weights, three survey series (NEFSC fall, NEFSC spring and DFO), weight-at-age and maturity. The panel explicitly accepts excluding the spring 2023 survey and fall 2025 survey. A downloaded pre-review run must be checked against those exclusions before importing numerical results.

Official acceptance source:
https://www.fisheries.noaa.gov/s3/2026-07/june-2026-management-track-peer-review-panel-report_508_20260710.pdf

Official data portal:
https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php

Public search parameters: species_id=5, stock_id=1, year=2026, review_type_id=6. Stock 8 is Eastern Georges Bank and must not be substituted.

## Cached source files

The gitignored source_cache/nefsc_haddock_georges_bank_2026 directory contains the 2026 assessment plan, peer-review report, detailed nine-page assessment report, portal HTML/search results and the official 2026_HAD_GB.zip download. The latter is a document bundle, not a native fitted-model archive. Its nested diagnostic archives contain 62 files in the principal diagnostics bundle but no RDS/RData/R/CSV/DAT files were found by filename extension. Further native-source searches are required before concluding that model files are unavailable.

Detailed report URL:
https://apps-nefsc.fisheries.noaa.gov/saw/sasi_files.php?year=2026&species_id=5&stock_id=1&review_type_id=6&info_type_id=-1&map_type_id=&filename=Georges_Bank_haddock_Update_2026_06_17_100805.801712.pdf

Bundle URL:
https://apps-nefsc.fisheries.noaa.gov/saw/sasi_files.php?year=2026&species_id=5&stock_id=1&review_type_id=6&filename=2026_HAD_GB.zip

## Numerical cross-check targets

The report gives 2025 SSB 29,037 t, Fbar (ages 5–7) 0.21 and age-1 recruitment 2,284 thousand fish, with no retrospective adjustment. These values agree with the canonical imported outputs and are checked by the validator below.

The existing candidate association to NEFSC-GARMIII_5Y_Melanogrammus_Aeglefinus denotes Gulf of Maine haddock, not Georges Bank haddock; do not carry that identity into this stock.

The scientist's official NOAA profile links to github.com/liz-brooks. Her public repository inventory and wham_devel recursive tree were cached and inspected; stock-specific 2026 Georges Bank files were not found there. Example/vignette RData objects are not official assessment objects and were not used. Search continues through the framework and other authority repositories.

## Charbonneau identity verified

The checked PLOS_2026_analysis/metadata_no_age_corection.csv has one NEFSC haddock record: NEFSC-GARMIII_5Y_Melanogrammus_Aeglefinus. Its notes explicitly identify Gulf of Maine haddock. No Georges Bank haddock record was found. The erroneous Georges Bank association has therefore been removed from the candidate seed script; this stock will be added from authoritative NOAA sources with the absence from the curated catalogue documented, rather than borrowing another stock's ID.

## Detailed framework source

The NOAA portal's 2022 Research Track search (species 5, stock 1, review type 5) provides the full assessment text and GBHaddock_Benchmark21_WHAM_BASE_Rcode.zip, both cached. The framework script reads BASE2_FORWHAM_DROPY41_2019.DAT and configures NAA_re=list(cor='2dar1', sigma='rec+1') in the base run. This establishes a framework source trail, but neither its 2019 numerical data nor its example fitted RDS references replace the 2026 official inputs. The current annual control settings still require verification.

## Separate 2022 Management Track source bundle

The 2022 Management Track is a separate run from the 2022 Research Track framework. NOAA's SASI portal identifies it with `review_type_id=2`; the framework search is `review_type_id=5`. The Management Track bundle [`2022_HAD_GB.zip`](https://apps-nefsc.fisheries.noaa.gov/saw/sasi_files.php?year=2022&species_id=5&stock_id=1&review_type_id=2&filename=2022_HAD_GB.zip) was downloaded on 2026-10-20 (SHA-256 `7519E03B78AD5F2D3E1A6577B5C98FA9D2407D2DFF289067FC99AF77300C4004`). It contains the 2022 update report, a 14-page input-figure PDF, an 18-page results-figure PDF, survey figures, and diagnostics. No native fitted object, input data table, WHAM data file, or machine-readable age-specific result table was found in the archive.

The input figures show total catch and catch age composition; three rescaled survey indices and their age-composition plots; annual weight-at-age for spawning biomass, January 1 biomass, catch, the fishery, and the indices; and maturity curves for selected years. These are useful checks on the model's streams and coverage, but they are not the original numerical series. In particular, the age-composition bubble plots provide only reference bubble sizes, and the index series are explicitly rescaled without a recoverable conversion back to the fitted input values. Digitizing them would introduce unquantified measurement and scale error, so they are not imported as exact official inputs.

The results PDF provides plotted historical trends and model diagnostics, not a fitted object or full numerical N-at-age and F-at-age matrices suitable for a complete comparison. It therefore does not supply a source-complete historical assessment for a tinyAM translation. The figure PDFs, report, original bundle, and portal search response are retained in the gitignored source cache with provenance in its local manifest.

The six multipart files `2022_HAD_GB_1.zip` through `_6.zip` are the distinct 2022 Research Track archive. Its code and outputs end in 2019 and do not represent the 2022 Management Track production run. Neither 2022 track is mixed into the 2026 record. The 2026 assessment remains the latest official run; its missing age-structured inputs and outputs continue to prevent a faithful tinyAM fit.

## Initial canonical import

The 2026 assessment record now contains eight published Catch for Assessment totals (2018–2025), the constant nine-age M vector stored once, 24 historical outputs (SSB, Fbar and age-1 recruitment for 2018–2025), and 16 documented assumptions. Terminal SSB/F match the final peer-review report. Published terminal 95% intervals are retained; missing SEs and other intervals are not invented. The official process is 2DAR1 and M is fixed at 0.2, according to final TOR 3.

All three completeness statuses remain partial. This initial summary import is not a completed curation: numerical catch-at-age, survey inputs/timing, historical biological matrices, full N/F-at-age outputs and current native configuration are still outstanding. No older framework numerical data have been substituted.

## Follow-up source checks

A renewed search for the 2026 WHAM object found no current native input archive. The WHAM comparison vignette explicitly labels its 2019 Georges Bank example as very preliminary; it is not a substitute for the official 2026 run (https://timjmiller.github.io/wham/articles/ex08_compare.html). The 2026 numerical catch, SSB, recruitment and Fbar rows are now checked directly against cached report text by `scripts/database/022_validate_georges_bank_haddock.R`, together with model identity, terminal year, fixed-M representation and continued partial status. Missing survey and age-composition coverage is tested explicitly rather than concealed by historical example data.

On 2026-10-20, the 2026 portal listing and `2026_HAD_GB.zip` were rechecked. The portal lists the assessment bundle, the report, one OSA diagnostic archive, one HTML comparison of GB and EGB results, two peer-review presentations, and an HTML comparison of runs with and without the fall 2025 index. The bundle contains six top-level files; its nested diagnostic archive contains 61 PNGs and a directory entry, while each comparison archive contains one HTML file. No current age-structured inputs, source data tables, or native fitted-model object are included. The comparison files discuss configurations and historical GB/EGB outputs but do not provide the missing 2026 input series or native fitted age-specific output tables. The public AFSC [`flatfish_assessment_code` repository](https://github.com/noaa-afsc/flatfish_assessment_code/tree/9709548fafc7d69a2e3f88f80f60c74b7d3cfb0e/examples/atf/bsai/m0) contains a BSAI arrowtooth example explicitly identified as a 2019 model, so it is not substituted for the 2022 production inputs. These checks do not establish that files are unavailable through every NOAA system; the public material reviewed still does not support a faithful tinyAM fit for this record.

The report's exact Table 3 projections for 2026–2031 are now retained in `outputs.csv` for total fishery catch, SSB, and Fbar. The 2026 catch is labelled as a supplied planning assumption; 2027–2031 catch values retain their reported intervals. SSB limits and the 2026 Fbar interval are retained, while scenario Fbar is recorded without invented uncertainty. This adds the published projection series but does not resolve the missing historical age-structured inputs or allow a tinyAM fit. The fitted terminal year remains 2025 and the reported projection horizon is 2031.

On 2026-10-21, the 2024 NOAA bundle was rechecked. Its five files are PDFs only: the 2024 update report, two OSA diagnostic reports, a WHAM methods paper, and a draft assessment report. It contains no machine-readable observations or native WHAM object.

The newly published DFO Research Document 2026/048 was also checked as a possible source for the transboundary stock. It is a separate Eastern Georges Bank assessment, not the NOAA Georges Bank assessment in this record. The report documents the DFO model's 2025 run (data through 2024), including the addition of DFO survey stratum 5Z9 and an estimated M for 2010–2024. Its appendix lists input coverage and model settings; the numerical tables are catch totals, projections, sensitivity diagnostics, and projection inputs. Catch-at-age and survey age compositions are shown only as figures, and there are no tabulated numerical input matrices or N/F-at-age outputs. The DFO publication record has no data/code links. The report was cached locally with its retrieval date and checksum, but it cannot supply or validate the missing NOAA 2026 inputs. This check does not change the NOAA record or its partial completeness statuses.

On 2026-10-21, the official NOAA SASINF portal was queried again for the 2026 Management Track record, separately filtering the `Tables` and `Models` categories. Each result contained only the assessment report; neither category exposed numerical input tables or a native WHAM file. The exact HTML responses and SHA-256 hashes are retained in the ignored source cache manifest. Publicly accessible portal materials still do not support a faithful tinyAM translation of the NOAA 2026 run.
