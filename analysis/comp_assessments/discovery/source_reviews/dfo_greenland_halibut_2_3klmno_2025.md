# DFO Greenland halibut, NAFO Subarea 2 + Divisions 3KLMNO

## Unit and current advice

The current unit is Greenland halibut (*Reinhardtius hippoglossoides*) in NAFO
Subarea 2 and Divisions 3KLMNO. It is a transboundary DFO/NAFO unit, distinct
from DFO's Gulf of St. Lawrence 4RST stock and from Greenland turbot units in
Alaska. No corresponding canonical stock record was found. The DFO
Newfoundland and Labrador inventory entry is provisionally treated as a new
unit pending completion of the cross-organization inventory crosswalk.

NAFO adopted a new Management Procedure (MP) for this stock in 2024. The 2025
Scientific Council report applies that procedure to calculate a 2026 TAC of
14,729 t. This is the current advice pathway: it uses aggregate survey mean
weight-per-tow indices, normally from five series. The EU-Spain 3L series was
excluded for the 2025 calculation because of missing recent values; the other
four indices were used. The latest survey year available for the procedure is
2024. The procedure also uses a relative age-4 recruitment indicator to assess
exceptional circumstances. It does not fit a new age-structured population
model to produce the annual TAC.

## Recent age-structured model evidence

Two detailed age-structured model updates were published in 2023 for the
management-strategy evaluation (MSE): an updated SCAA base case and an updated
state-space model (SSM). Their fitted data run through 2021. The SCAA follows
the methodology used in 2017, with data updated for 2017–2021. The SSM update
was explicitly prepared to supply updated operating-model outputs for the 2023
MSE review. It evaluated a base model using the EU-Spain 3NO survey, a
sensitivity that combined EU-Spain 3LNO, and a sensitivity that added EU-Spain
3L separately. These are model variants for one stock unit, not additional
stock cases. They should be represented as the last detailed age-structured
model events, separately from the current annual MP advice.

The 2023 SSM report documents:

- Commercial catch numbers-at-age in thousands of fish for ages 1–10+, by year
  from 1975 to 2021. The report also provides landings in kilotonnes; the SCAA
  uses landings from 1960 onward.
- Survey mean numbers per tow at ages 0–10+, with year coverage varying by
  survey. The base model uses Canadian Autumn 2J3K, Canadian Spring and Autumn
  3LNO, EU-Spain 3NO, and EU 3M surveys. EU-Spain 3L is only included in
  sensitivities. Canadian age-8+ indices are omitted for 2018, 2019 and 2021
  when deep strata were not covered.
- Catch weights-at-age in kilograms and survey biomass indicators as mean
  weight per tow. The SCAA report also gives proportions mature-at-age and
  catch weights-at-age; pre-1975 weights use the 1975–1979 average.
- Fixed natural mortality of 0.12 per year in the SSM base case. Fishing
  mortality is fixed at zero for ages 1–2 and estimated for ages 3–10+.

The SSM observation model corrects abundance for survey timing using the
midpoint of survey dates as a fraction of year. The reviewed report describes
this convention but does not provide the numerical year-by-survey midpoint
series. Recover or derive those values from the source survey dates or native
model files before making an exact tinyAM translation.

The 2023 SSM report includes age-specific abundance and biomass proportions,
fishing mortality-at-age, recruitment, exploitable biomass ages 5–9 and mean F
ages 5–9, with 95% intervals for key outputs. The detailed SCAA report provides
biomass, fishing-mortality and recruitment trajectories and tabulated estimates
with Hessian-based CVs. Numeric age-specific fitted arrays are shown in report
figures rather than supplied as a machine-readable output file in the reviewed
sources. The 2025 MP report's TAC and survey indicators must not be described as
fitted age-specific population outputs.

## Screening decision

This unit is eligible for source curation of the recent 2023 age-structured
SCAA/SSM assessment updates: catch-at-age, age-specific survey observations,
weights, maturity, a stated M assumption and age-structured outputs are
available in detailed official reports. Eligibility is for those detailed
age-based model events, not for a direct age-structured replication of 2025
annual advice. The canonical source record should preserve the current 2025 MP
advice and its survey inputs separately from the 2023 age-structured MSE models.
Do not combine 2025 MP observations with the 2023 fitted age-model outputs.

Exact within-year survey timings and machine-readable fitted N/F-at-age values
remain source-recovery tasks. No numerical inputs or outputs have been imported
to the canonical database during this screen.

## Sources

- [2025 NAFO Scientific Council Report, SCS Doc. 25/13](https://www.nafo.int/Portals/0/PDFs/sc/2025/scs25-13.pdf)
- [2024 NAFO Conservation and Enforcement Measures, revised, including the Greenland halibut MP](https://www.nafo.int/Portals/0/PDFs/com/2024/comdoc24-01REV.pdf)
- [2023 updated SCAA base-case assessment, NAFO SCR Doc. 23/008REV](https://www.nafo.int/Portals/0/PDFs/sc/2023/scr23-008REV.pdf)
- [2023 updated state-space model, NAFO SCR Doc. 23/014](https://www.nafo.int/Portals/0/PDFs/sc/2023/scr23-014.pdf)

Reviewed PDFs, retrieval dates, source URLs, and SHA-256 checksums are in the
gitignored local cache `source_cache/discovery_dfo_greenland_halibut_2_3klmno_2025/`.