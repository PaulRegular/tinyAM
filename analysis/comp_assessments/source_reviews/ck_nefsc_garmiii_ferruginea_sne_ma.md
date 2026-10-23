# Southern New England/Mid-Atlantic yellowtail flounder

## Source identification

The stock is `ck_nefsc_garmiii_ferruginea_sne_ma` (*Limanda ferruginea*), the SNE/MA component in the NOAA NEFSC assessment. The database currently contains only its Charbonneau–Keith catalogue metadata; there is no assessment event or numerical input/output record yet.

The 2024 Yellowtail Flounder Research Track Working Group report is the most detailed public source found for the SNE/MA model. Its SASI archive includes the report, model-selection material, and a `SNEMA_Runs.zip` bundle. The nested run archive contains diagnostic and result figures plus PDFs, but no native numerical input arrays, fitted object, or tabular time-by-age output. The report contains source tables for commercial data and survey abundance-at-age, and a terminal 2022 SSB/F summary, but these have not yet been reconciled to the exact WHAM likelihood inputs and output definitions. They are not imported as if they were verified model arrays.

The 2025 Management Track archive is a distinct later applied update with data through 2024. Its report and comparison workbook are cached separately. Do not combine its series with the 2024 Research Track run. The 2025 archive does not provide the detailed age-structured fit files needed to fill the missing source record.

## Coverage and next step

The source cache is under `source_cache/ck_nefsc_garmiii_ferruginea_sne_ma/`; its manifest records the NOAA SASI archives, retrieval date, checksums, and represented run. The 2024 report can support further extraction of catch, survey, and biological tables, but input completeness cannot be established until each table is matched to the fitted model's actual data structure. Full annual and age-specific fitted outputs are not available as machine-readable values in the public run bundle.

Keep this stock out of the translation queue until the source tables are reconciled to the production likelihood and a defensible official-output reference can be assembled. If those checks remain unresolved, report the case as source-limited and continue with another stock rather than estimating missing observations or digitizing plots without an explicit, validated method.

## Sources

- [2024 Yellowtail Flounder Research Track report and SASI files](https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php?year=2024&species_id=13&stock_id=11&review_type_id=5)
- [2025 SNE/MA Management Track SASI files](https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php?year=2025&species_id=13&stock_id=5&review_type_id=2)
- [2025 Management Track report](https://apps-st.fisheries.noaa.gov/sis/docServlet?fileAction=download&fileId=10439)