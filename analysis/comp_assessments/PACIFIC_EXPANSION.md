# Pacific expansion

The expansion adds twelve stock cases, not separate sexes, model variants or
regional components. The eight primary catalogue identities below come from
Charbonneau–Keith's PLOS analysis set. Pacific ocean perch is in the broader
catalogue. Coastal hake and the two Canadian cases are explicit exceptions.

| Stock | Catalogue identity | Detailed production source being reviewed |
|---|---|---|
| Eastern Bering Sea Pacific cod | AFSC_BSAI_Gadus_macrocephalus | 2024 Model 24.1, native SS files and Appendix 2.4 workbook |
| Aleutian Islands Pacific cod | AFSC_AI_Gadus_macrocephalus | 2024 Model 24.1, final SAFE report and public native files |
| Aleutian Islands pollock | AFSC_AI_Gadus_chalcogrammus | 2024 Model 15.1, SAFE tables and AMAK framework |
| BSAI Atka mackerel | AFSC_BSAI_Pleurogrammus_monopterygius | 2024 Model 16.0b, SAFE tables and AMAK framework |
| Eastern Bering Sea yellowfin sole | AFSC_EEBSAI_Limanda_aspera | 2024 Model 23.0, native ADMB files |
| Eastern Bering Sea northern rock sole | AFSC_EEBSAI_Lepidopsetta_polyxystra | 2024 Model 24.2, native ADMB inputs and SAFE outputs |
| BSAI Alaska plaice | AFSC_BSAI_Pleuronectes_quadrituberculatus | 2024 Model 24.1b, native SS input and output files |
| BSAI arrowtooth flounder | AFSC_BSAI_Atheresthes_stomais | 2022 Model 18.9; 2024 document is a projection update |
| BSAI Pacific ocean perch | AFSC_alutus_BSAI | 2024 Model 24, native ADMB input and output files |
| Joint US–Canada coastal hake | Catalogue exception | 2026 catch-only update, 2025 detailed assessment and native SS files |
| Strait of Georgia herring | Canadian catalogue exception | 2025 update, SISCAH ensemble and interim estimation model |
| British Columbia Pacific ocean perch | Canadian catalogue exception | 2023 three-area assessment, Research Document 2025/004 and native SS inputs |

Historical catalogue labels do not establish current eligibility. In particular,
Aleutian Islands cod now has an age-structured Tier 3 assessment. The historical
arrowtooth catalogue identifier retains its original spelling; the species is
*Atheresthes stomias*.

The final December 2024 SSC report identifies the selected Alaska models.
NOAA's 2026 harvest-specifications review explains why the 2024 fitted models
remain the relevant detailed production runs after the 2025 funding lapse.
Those documents are cached in `source_cache/pacific_expansion/`, with manifests.
Current-year development directories are not evidence of an official new fit.

Each stock receives a source review and a database checkpoint before its
translation checkpoint. Cache manifests record retrieval dates, checksums and
pinned repository revisions. Numerical source blockers and non-equivalent
comparisons remain explicit; they do not count as completed translations.
