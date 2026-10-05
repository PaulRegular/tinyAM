# Summer flounder source review

## Assessment and available files

The latest accepted assessment is the 2025 NOAA Fisheries management-track
assessment, based on data through 2024. The June 2025 peer-review panel found it
met the terms of reference and represented the best scientific information
available. It uses ASAP and updates catch, survey indices, weights, and reference
points. The assessment combines sexes, estimates recruitment at age 0, and
reports fully selected fishing mortality at age 4. Its four fishery fleets are
commercial and recreational landings and discards.

The 2025 NOAA SASINF download is cached locally under
`source_cache/nefsc_summer_flounder_2025/`. It contains a nine-page assessment
summary, a 323-page plots file, bridge-run material, and a presentation. The
portal's separate **Tables** and **Models** searches return only the nine-page
summary. The package contains no age-specific input tables, model input files,
or fitted ASAP object. The summary has recent total catch, SSB, fully selected
F, and recruitment values, but those do not supply the full age-specific series
needed for a tinyAM fit and comparison.

The 2023 SASINF package is also cached under
`source_cache/nefsc_summer_flounder_2023/`. It likewise contains a nine-page
summary, plots, and a presentation, but no machine-readable model inputs or
fitted object. It therefore does not close the input gap for a reproducible
age-structured fit.

## Translation status

No summer flounder observations or assessment outputs have been added to the
canonical database, and no tinyAM fit has been attempted. Reconstructing the
missing age-specific series from plots or substituting an older assessment's
inputs would create an undocumented approximation. The stock remains blocked
until age-specific catch, survey, biology, and fitted-output data for one
accepted assessment can be recovered.

## Sources

- [2025 June Management Track Peer Review Panel Report](https://www.fisheries.noaa.gov/s3/2025-07/2025-June-Management-Track-Peer-Review-Panel-Report-508.pdf)
- [NOAA SASINF assessment-data search](https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php)
- [2025 Summer Flounder Management Track Assessment Report](https://asmfc.org/wp-content/uploads/2025/08/SF_Management_Track_Assessment_2025.pdf)
- [NOAA Summer Flounder assessment status](https://www.fisheries.noaa.gov/species/summer-flounder/science)
