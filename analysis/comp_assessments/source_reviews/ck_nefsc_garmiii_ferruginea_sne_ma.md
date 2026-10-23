# Southern New England/Mid-Atlantic yellowtail flounder

## Source identification

The stock is `ck_nefsc_garmiii_ferruginea_sne_ma` (*Limanda ferruginea*), the SNE/MA component in the NOAA NEFSC assessment. The database currently contains only its Charbonneau–Keith catalogue metadata; there is no assessment event or numerical input/output record yet.

The 2024 Yellowtail Flounder Research Track Working Group report is the most detailed public source found for the SNE/MA model. Its numerical tables include total commercial catch-at-age in thousands of fish for 1973–2022 (Table 2.8.3) and age-specific abundance from the NEFSC spring, fall and winter surveys (Tables 3.6.2–3.6.4). The report also describes the final candidate model: ages 1–6+, a combined fleet, a 0.05 fleet CV, fixed M=0.5, and three NEFSC surveys. However, its SASI run bundle contains plots and PDFs rather than native numerical input arrays or a fitted object. The report does not provide all exact age-composition likelihood inputs, including sample sizes/weights, or complete annual weight and maturity matrices needed for the corresponding biomass calculations. Its fitted annual N/F/SSB surfaces are shown as figures; only terminal 2022 values are stated numerically (F=0.11, SSB=44 mt, recruitment=495,000). The published tables therefore support source data recovery, but not a complete verified reconstruction of the likelihood inputs or age-specific outputs.

The separate 2025 Management Track report updates the WHAM model through 2024, but the cached report is explicitly marked “draft working paper for peer review only.” It reports summary results and plots, not detailed age-structured fit files. The valid comparison workbook bundled with the 2025 archive compares model settings only; it does not contain fitted time series. A separately retrieved file bearing the same workbook name was an HTML response, not an Excel workbook, and is not treated as source data. Do not combine the 2025 update with the 2024 Research Track tables.

## Coverage and next step

The source cache is under `source_cache/ck_nefsc_garmiii_ferruginea_sne_ma/`; its manifest records the NOAA SASI archives, retrieval date, checksums, and represented run. Available report tables could be added as clearly identified source observations with partial status, but they do not resolve the missing likelihood weights/sample sizes or complete biology. Full annual and age-specific fitted outputs are not available as numerical values in the public run bundle. The present decision is to keep this stock out of the translation queue until a complete verified input set and defensible official-output reference are available. Do not fill gaps with tutorial data, digitized figures, or unverified biological reconstructions.

## Sources

- [2024 Yellowtail Flounder Research Track report and SASI files](https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php?year=2024&species_id=13&stock_id=11&review_type_id=5)
- [2025 SNE/MA Management Track SASI files](https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php?year=2025&species_id=13&stock_id=5&review_type_id=2)
- [2025 Management Track report](https://apps-st.fisheries.noaa.gov/sis/docServlet?fileAction=download&fileId=10439)
