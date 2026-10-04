# Southern Gulf cod (4T–4VN) source review

For Southern Gulf cod, the 2024 Science Advisory Report is the latest accepted
assessment, using the SCA model through 2023. DFO identifies the 2019 run as the
last full assessment and says the same population model was used again in 2024.
The 2024 database record is current but partial: it contains the reported 2023
SSB estimate and does not borrow observation series or estimated outputs from
the 2019 run. Structural settings are cross-referenced to the 2019 detailed
report only where the 2024 summary confirms the same population model was used.
An official DFO dataset was found for annual SSB medians and the 2.5th, 25th,
75th, and 97.5th percentiles from the assessment to 2023, in thousands of
tonnes. Its [SSB CSV](https://api-proxy.edh-cde.dfo-mpo.gc.ca/catalogue/records/fe51e3da-0e0b-11ef-90aa-8b219c568296/attachments/NAFO-4T4VN-Atlantic-Cod-spawning-stock-biomass-estimates-1950-2023.csv)
is named `NAFO-4T4VN-Atlantic-Cod-spawning-stock-biomass-estimates-1950-2023.csv`.
The dataset landing page reports an update date of 2026-04-17. The attachment
could not be downloaded in this environment because outbound TLS connections
fail, so no annual SSB records have been added yet. The reported 2023 median
and 95% interval already stored from the 2024 report remain the only current-run
output values in the database.

The detailed 2019 SCA assessment remains as a historical record under the 2012
framework. It contains 54 annual stock-catch values for 1965-2018
(tonnes), plus 480 source-reported landed numbers-at-age values for 1971-2018
(ages 3-12+, in thousand fish). These age counts do not recover the full fitted
catch-composition series, which the model describes as ages 2-12+. It also
retains 759 annual maturity rows carried forward between source-listed change
years. Its survey records keep the model's aggregate index separate from its
age-composition input. Both are explicitly marked reconstructions from
published age-specific source tables, and no derived age-specific abundance
series is entered. The reconstructed RV aggregate biomass series has 46 years because
weights-at-age are not tabulated for 1980 and 1985; the mobile-sentinel series
covers 2003-2018; and the longline index covers 1995-2017. These inputs remain
partial because the report does not publish all original composition samples,
all RV biomass-index years, or every input over the model's 1950-2018 span. The
recorded population outputs are maximum-likelihood estimates from Tables
21-23; the report's other population summaries are generally posterior
medians. The report also states terminal estimates for M at ages 5-8 and 9+
and fully recruited q for the RV and mobile-sentinel surveys; these are stored
as grouped or time-invariant outputs with their source identified as report
text.
