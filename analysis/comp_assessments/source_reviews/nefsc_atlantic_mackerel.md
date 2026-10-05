# Atlantic mackerel source review

## Assessment and available files

The current U.S. Northwest Atlantic (*Scomber scombrus*) assessment is the
September 2025 NOAA management-track update, using data through 2024. NOAA still
reports stock status from the 2025 assessment. The assessment was initially
planned as a WHAM implementation, but delayed availability of the egg index
meant the final update retained the ASAP approach. Its sources include fishery
catch and age data, the spring bottom-trawl survey, and a range-wide egg/SSB
index.

The 2025 NOAA SASINF package is cached locally under
`source_cache/nefsc_atlantic_mackerel_2025/`. It contains an 11-page report and
three diagnostic-plot PDFs (base model, MCMC, and retrospective plots). The
package has no age-specific input tables, model input files, or fitted model
object. The report and plots therefore do not provide the complete input and
fitted age series needed for a reproducible tinyAM fit and comparison.

## Translation status

No Atlantic mackerel observations or assessment outputs have been added to the
canonical database, and no tinyAM fit has been attempted. The age-specific data
cannot be recovered from the files available in the assessment package without
reconstructing values from figures. This stock is blocked pending recovery of
the assessment's age-specific catch, survey, biology, and fitted-output data.

## Sources

- [NOAA Atlantic mackerel science and current status](https://www.fisheries.noaa.gov/species/atlantic-mackerel/science)
- [2025 Atlantic Mackerel Management Track Assessment materials, NOAA SASINF](https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php)
- [October 2025 Scientific and Statistical Committee meeting materials](https://www.mafmc.org/ssc-meetings/october-2025)
- [2025 Atlantic Mackerel assessment summary](https://static1.squarespace.com/static/511cdc7fe4b00307a2628ac6/t/68ded5113fae5a12be117b6d/1759434002250/2025_Atlantic_Mackerel_Assessment.pdf)
- [February 2025 Assessment Oversight Panel report](https://www.fisheries.noaa.gov/s3/2025-03/February-2025-AOP-Report.pdf)
