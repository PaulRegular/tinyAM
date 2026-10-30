# Assessment discovery inventory

`candidates.csv` is a discovery register, not a canonical assessment table. Its
rows are source-list entries, some reviewed and others not yet reviewed: they are not
automatically verified distinct new stock units, eligible age-structured assessments,
or a selected import batch. Screening, identity matching, source verification and
exclusion decisions remain outstanding. `candidate_unit_id` links inventory rows only
when reviewed evidence confirms a shared stock unit; blank values are unresolved or
not yet crosswalked. Do not overwrite the register after review has begun.

The initial snapshot was retrieved on 2026-10-28 from:

| Source | Entries | Local cached source and provenance | Coverage limit |
|---|---:|---|---|
| NOAA Stock SMART via NOAA-EDAB | 432 stock IDs | `source_cache/discovery_noaa_stocksmart/stock_assessment_summary.rda`, revision `664c4237e26520b260a090e1821a2cee5d1cdf07`, SHA-256 `AD47F9EEB240C2E77B2B952809A12C4FA739E3F0D39BEAB283B413FDBAEF76A8` | Latest metadata row per stock ID, not a verified official production run or proof of age observations. Regional official sources require separate review. |
| ICES SID | 274 published stock entries | `source_cache/discovery_ices_sid/stock_list_2026.json`, SHA-256 `913D6A3EEAB9449AE3CD0DD15661A150D36FF79B6E6A9BA78D86D46BA1035F42` | Published 2026 stock list includes assessments and advice that may not use age observations. |
| DFO Sustainability Survey | 198 selected stocks | `source_cache/discovery_dfo_sustainability/survey_2024.html`, SHA-256 `A10E65EEF028C2A524973FFCFA155915D673AE3AB9A9F1312581449E29644A28` | The 2024 survey covers only selected managed stocks, including non-finfish. Supplement with regional CSAS inventories. Official IDs and scientific names require review. |

## Current inventory status

The inventory has 905 source-list rows. These are not unique stock counts.
`inventory_progress.csv` breaks row-level identity, screening, and eligibility
status down by organization and source region; `inventory_units.csv` summarizes
the rows already linked to a stock or candidate unit.

| Organization | Rows | Matched to database | New candidate rows | Ambiguous identity | Identity unreviewed | Eligible rows | Ineligible rows | Eligibility unresolved | Screened | Not screened |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DFO | 199 | 20 | 48 | 10 | 121 | 21 | 133 | 45 | 176 | 23 |
| ICES | 274 | 59 | 84 | 2 | 129 | 51 | 76 | 147 | 145 | 129 |
| NOAA | 432 | 28 | 153 | 4 | 247 | 54 | 101 | 277 | 187 | 245 |
| **Total** | **905** | **107** | **285** | **16** | **497** | **126** | **310** | **469** | **508** | **397** |

The crosswalk currently links 312 distinct unit keys across 331 inventory rows:
100 link to canonical stocks and 212 are candidate units: 210 are new, one is
ambiguous and one has mixed identity status. Of the 210 new
candidate units, 30 are eligible for source curation, 149 are ineligible for
the age-disaggregated set, and 31 remain unresolved. A further 574 inventory
rows have no unit crosswalk yet. Row-level status and unit-level deduplication are intentionally
reported separately.

The 2024 DFO WCVI natural-origin Chinook assessment is now recorded as a new
unit, with age-input eligibility unresolved. Its detailed report describes
run-reconstruction samples and supporting data, but a complete, aligned age-
specific input set for the natural-origin unit has not been verified. See
`source_reviews/dfo_wcvi_chinook_2024.md`.

The 2026 ICES boarfish unit `boc.27.6-8` is newly screened ineligible for the age-disaggregated set: its benchmarked SS3 assessment fits total catch, biomass indices and length compositions, while its 2026 outputs include age-derived quantities. See `source_reviews/ices_boarfish_6_8_2026.md`.

The 2025 ICES red gurnard advice (`gur.27.3-8`) describes a survey trends-based assessment with an aggregate biomass index, no reliable catch data, and no age-specific inputs or population outputs. It is screened ineligible. See `source_reviews/ices_red_gurnard_3_8_2025.md`.

NOAA Stock SMART Pacific bluefin tuna (ID 11687) is a new Pacific-wide ISC/PBFWG assessment unit. The 2024 fully accepted SS3.30 run uses catch totals, size compositions, and annual fishery CPUE indices; its reported catch-at-age is estimated internally, not an observed age input. The age-input candidate is screened ineligible. See `source_reviews/noaa_pacific_bluefin_2024.md`.

NOAA Stock SMART ID 11827 is the Main Hawaiian Islands Deep 7 multispecies complex. Its accepted 2024 JABBA assessment uses total catch and aggregate fishery/survey indices; a supplementary Opakapaka age-structured model uses length and weight compositions but no documented age-composition likelihood. The complex is screened ineligible for age-disaggregated inputs. Its seven component species remain one assessment unit. See `source_reviews/noaa_mhi_deep7_2024.md`.

NOAA Stock SMART IDs 11997, 11998, 12003, and 12005 have newer 2026 American Samoa species assessments than the 2023 rows in the pinned snapshot. The WPSAR-reviewed update fits four species using aggregate catch, aggregate boat-based CPUE, and length compositions; it reports SSB, age-0 recruitment, and annual F, but not age-specific observed inputs or numeric N/F-at-age outputs. They are screened ineligible for the age-disaggregated set. Five other 2023 Samoa species entries are outside the revised 2026 assessment list, while the older 2019 multispecies complex has an unresolved one-to-many crosswalk. *Etelis boweni* and *E. carbunculus* are current management species but were not modeled separately because their data are confounded; they have no separate Stock SMART inventory rows. See `source_reviews/noaa_american_samoa_bottomfish_2026.md`.

NOAA Stock SMART IDs 11696 (Guam, 2024) and 11731 (CNMI, 2025) are current 13-species bottomfish-complex JABBA assessments. Both use aggregate catch and standardized creel CPUE, report biomass and harvest-rate outputs, and provide no age-specific likelihood inputs or numeric N/F-at-age outputs; each complex is one ineligible assessment unit. The 2026 CNMI life-history survey report is not an input source for the 2025 model. See `source_reviews/noaa_marianas_bottomfish_2024_2025.md`.

NOAA Stock SMART IDs 12452 and 14617-14623 are eight species-specific Main Hawaiian Islands reef-fish assessments in the fully accepted 2016 PIFSC report. The report uses aggregate catch, diver abundance and size composition; it derives status metrics with an age-structured calculation but does not provide the age-specific observations and annual N/F-at-age outputs required by this candidate set. See `source_reviews/noaa_mhi_reef_fish_2016.md`.

NOAA Stock SMART ID 11686 links to the fully accepted 2022 ISC benchmark for North Pacific blue shark. The SS3 model is age- and sex-structured but uses aggregate catch, CPUE, and length compositions; the public assessment does not provide the age-disaggregated observations and annual N/F-at-age arrays required for this set. See `source_reviews/noaa_north_pacific_blue_shark_2022.md`.

NOAA Stock SMART 10434, Gulf of Alaska sharpchin rockfish, is the Tier 4 component of the GOA Other Rockfish complex. Stock SMART lists a fully accepted 2023 REMA Model 15.2 run; the 2024 SAFE update carried forward Tier 4 biomass because no new survey data were available to refit, and the 2025 SAFE catch report says no assessment was conducted and prior values were rolled forward. The assessment uses aggregate survey biomass and catch, not age-specific population likelihoods or outputs, so the unit is screened ineligible. See `source_reviews/noaa_goa_sharpchin_2024.md` and its local source cache manifest.


Regional supplementation: the [2026 DFO Pacific groundfish plan](https://www.pac.dfo-mpo.gc.ca/fm-gp/mplans/ground-fond-ifmp-pgip-sm-eng.html), retrieved 2026-10-28, led to a British Columbia yellowtail rockfish candidate absent from the initial survey list. The cached regional page SHA-256 is `2CB8792B9AEBE6FE5B27F3887BEFF57F02E2E7CE15D9344353AD705AED18475B`. This partial regional review is not an exhaustive CSAS inventory; further regional sources must still be audited.

The [Pacific Fishery Management Council groundfish assessment listing](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/) was retrieved 2026-10-28; cached HTML SHA-256 `77FEB82D8A674A018FD3FC7E02CE16579B36F3101B8603FF61286927B1D8B32D`. It provides reports and model files for recent West Coast assessments. The chilipepper identity has been reviewed, but the remainder of this regional listing still needs reconciliation and screening.

The [2026 WGBIE working-group report](https://doi.org/10.17895/ices.pub.32790321.v1) is a regional assessment listing and source. A partial review of its megrim and sole sections is recorded in `source_reviews/ices_iberian_wgbie_2026.md`; the remaining sections still need reconciliation. The locally cached sections and working-document annex have individual URLs and checksums in `source_cache/discovery_ices_wgbie_2026/sources.csv`.

The full [2026 WGCSE report](https://doi.org/10.17895/ices.pub.32625759.v1) has been cached and used to screen its method categories. Eleven candidate entries are ineligible under the current non-age-structured methods; 14 category-1 age-based stocks remain for stock-specific source screening. See `source_reviews/ices_wgcse_2026.md` and `source_reviews/ices_sol_27_7a_2026.md`; the remaining age-based sections still require reconciliation and input/output review.

The [2026 WGNSSK report](https://doi.org/10.17895/ices.pub.32676345.v2) is a North Sea and Skagerrak regional assessment listing. Its witch section has been screened and cached; other sections still require reconciliation. See `source_reviews/ices_witch_northsea_2026.md`.

ICES plaice in Division 5.a (ple.27.5a) is a distinct current Icelandic unit screened against the 2026 SAM production assessment. Catch-at-age, age-specific spring-survey indices and age-specific biology are documented, along with annual recruitment, SSB and F outputs with uncertainty. Numerical input matrices are mainly figure-based, so this is eligible for source curation but not yet a ready numerical import. The 2026 fitted series ends in 2025; the 2026 assessment-graph row is a projection. See source_reviews/ices_plaice_iceland_2026.md.
ICES Atlantic wolffish in Division 5.a (caa.27.5a) is a distinct Icelandic unit assessed with SAM in 2026. It uses commercial catch-at-age, age-specific spring and autumn survey indices, age-specific biology and fixed M = 0.15; ICES reports annual age-4 recruitment, SSB and F-at-ages 10-15 with uncertainty. The current report's catch-age coverage differs from the 2022 Stock Annex, and current input matrices are primarily figure-based. See source_reviews/ices_atlantic_wolffish_iceland_2026.md.
The 2026 Celtic Sea herring assessment (`her.27.irls`) is a distinct provisional new eligible ICES case. Its report provides numerical catch-at-age, age-specific acoustic survey inputs, biology, timing and annual SSB/F/recruitment outputs; raw bootstrap files and generated numeric N/F-at-age outputs are absent from the public TAF snapshot. The report also documents unquantified mixing with Irish Sea herring. See `source_reviews/ices_celtic_herring_2026.md`. The northern Irish Sea herring assessment (`her.27.nirs`) is a separate eligible unit with catch-at-age, ISAS age-specific abundance, an ISSS SSB index, annual biology and annual SSB/recruitment/F outputs. Its report flags mixed-stock observations and survey catchability scaling (2.7–4.9); numerical fitted N/F-at-age arrays were not recovered. See `source_reviews/ices_irish_sea_herring_2026.md`.
Cached files and their sources.csv manifests are gitignored. The register retains
inventory source URLs and leaves unknown identities, run details and data evidence
blank. `scripts/database/186_build_discovery_inventory.R` builds the initial
snapshot and refuses to overwrite an existing register, protecting later reviews.
Candidate discovery does not modify the canonical database.

DFO Atlantic halibut in 3NOPs4VWX5Zc is a distinct DFO management unit from the separately defined 4RST unit. The latest full assessment is the 2022 SISCAL model, with detailed inputs published in 2025. It is ineligible for the age-indexed candidate set because observed likelihoods are based on catch-at-length and aggregate survey indices, despite the model being age-structured internally. The relationship to NOAA Stock SMART 10056 is not one-to-one, so the two units are not coalesced. See source_reviews/dfo_atlantic_halibut_3nops4vwx5zc_2022.md.

DFO Atlantic halibut 4RST is a separate provisional new unit. Its 2024 assessment (reviewed in 2025) is classified as a delay-difference model and reports aggregate exploitable biomass, fishing mortality and recruitment; it is ineligible for the age-indexed candidate set. The separate 2025/028 research document reviews 2022 data and describes otoliths aged for growth, not an age-indexed likelihood for the 2024 assessment. The 2024 SAR's assessment-date field conflicts with its title and method-change description; see `source_reviews/dfo_atlantic_halibut_4rst_2024.md`.

DFO Greenland halibut in NAFO Subarea 2 + Divisions 3KLMNO is a provisional new unit eligible for source curation of its detailed 2023 age-structured SCAA/SSM updates. Current 2025 advice uses the 2024-adopted Management Procedure, so the age-model updates from the MSE work must be kept separate from the current MP event. Exact survey timing and machine-readable age-specific output arrays remain to be recovered; see `source_reviews/dfo_greenland_halibut_2_3klmno_2025.md`. The distinct Gulf of St. Lawrence 4RST unit is ineligible: its latest detailed assessment is index-based and the 2025 update reports aggregate biomass above 40 cm; age and maturity analyses are supporting biology, not age-specific assessment inputs. Assessment-date fields differ across the 2025 and 2026 reports; see `source_reviews/dfo_greenland_halibut_4rst_2025_update.md`.

The DFO Cambridge Bay Arctic Char survey entry covers five river-specific
management units. Its latest detailed reports include data-limited DB-SRA work
for Halokvik and Jayko, with biological age samples but no catch-at-age model
likelihood or age-specific population outputs. The aggregate entry remains
identity-ambiguous and is excluded from the age-disaggregated set; see
`source_reviews/dfo_arctic_char_cambridge_bay.md`.

The DFO Cumberland Sound Arctic Char entry also spans separate management
units. The latest detailed assessment located is the 2017 review published as
Science Advisory Report 2018/021 for Ijaruvung Lake, Iqalujjuaq Fiord and Irvine
Inlet. It uses aggregate harvest and data-limited methods; biological age
samples and aggregate abundance/biomass outputs do not provide age-structured
inputs or population-at-age outputs. The row remains identity-ambiguous and is
ineligible for the age-disaggregated set. Ikaluit Lake (Robert Peel Inlet) is a
separate 2023 assessment. See
`source_reviews/dfo_arctic_char_cumberland_sound.md`.


The DFO North Slope Dolly Varden entry combines several river populations with
different assessment histories. The register row remains identity-ambiguous and
eligibility-unresolved; Rat River is a distinct follow-up candidate because
DFO's 2025 research document describes an integrated statistical catch-at-age
model. Its numerical inputs, link to the formal advice run, and the newer
2015–2022 monitoring report still need verification. See
`source_reviews/dfo_north_slope_dolly_varden.md`.

The DFO Gulf and Newfoundland and Labrador Atlantic salmon rows each combine
many river-level populations and are not single stock units. Their current
advice is based on river counts, fishery or juvenile indices, and egg-deposition
calculations; age and size biology in some programs does not establish
catch-at-age and population-at-age inputs and outputs for an age-structured
assessment. Both regional rows are screened out as composite candidates, while
river-specific units remain possible subjects for separate review. See
`source_reviews/dfo_atlantic_salmon_gulf_2025.md` and
`source_reviews/dfo_atlantic_salmon_nl_2026.md`.

The DFO Greenland halibut Cumberland Sound Turbot Management Area is a separate
management unit from NAFO Subarea 0, though the 2019 review documents movement
and emigration between Cumberland Sound and offshore waters. Its latest review
was held in 2019 and published as Science Advisory Report 2026/007. The report
provides aggregate survey abundance and biomass indices, but the generic
population model had poor convergence and was not used for advice; matched
age-specific inputs and accepted population-at-age outputs were not located.
The entry is ineligible for the age-disaggregated set. The report PDF redirected
to an archived-information notice during local retrieval; available supporting
documents and checksums are in the ignored source cache. See
`source_reviews/dfo_greenland_halibut_cumberland_sound_2019.md`.

NOAA Stock SMART 11082 and DFO's Pacific halibut entry identify one shared
IPHC Convention Area stock, new to the canonical database. It is eligible for
source curation: public commercial age frequencies and area-specific landings,
FISS age samples and NPUE, biological data, and annual spawning-biomass outputs
support documented age-input reconstruction. Exact Stock Synthesis input files
require IPHC permission, and the public commercial-age table was revised after
the 2025 assessment cutoff; see
`source_reviews/noaa_dfo_pacific_halibut_iphc_2025.md`.

ICES Bay of Biscay sea bass (`bss.27.8ab`) is screened eligible for source
curation from the 2026 update (data through 2025) and the 2025 benchmark input
annex. Annual conditional age-at-length data support catch-at-age
reconstruction, while two juvenile indices, adult LPUE, age-based mortality,
growth and maturity assumptions, and recruitment/SSB/Fbar outputs are
documented. Exact weight-at-age and index timing still need validation; no
canonical data were added. See
`source_reviews/ices_seabass_8ab_2026.md`.

ICES Bay of Biscay whiting (`whg.27.89a`) is screened ineligible from its 2026
Category 3.2 assessment: the current method reports an aggregate biomass index
and length-based fishing-pressure proxy, not age-disaggregated inputs or
outputs. No canonical data were added; see
`source_reviews/ices_whiting_89a_2026.md`.

ICES eastern Baltic cod (`cod.27.24-32`) is reconciled to the existing eastern
Baltic cod stock record. The canonical extent now reflects subdivisions 24–32
and retains the historical Charbonneau–Keith 25–32 label as an alias. Its
current 2026 index assessment remains distinct from the historical 2024
age-structured assessment; see
`../source_reviews/ck_ices_wgbfas_morhua_bs_25_32.md`.

## Inventory progress snapshot

Snapshot dated 2026-10-30: the register contains 905 source-list entries. These are not deduplicated stock units.

| Organization | Entries | Matched | Provisional new | Ambiguous | Identity unreviewed | Eligible | Excluded | Eligibility unresolved | Screened | Not started |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DFO | 199 | 20 | 48 | 10 | 121 | 21 | 133 | 45 | 176 | 23 |
| ICES | 274 | 59 | 84 | 2 | 129 | 51 | 76 | 147 | 145 | 129 |
| NOAA | 432 | 28 | 153 | 4 | 247 | 54 | 101 | 277 | 187 | 245 |
| Total | 905 | 107 | 285 | 16 | 497 | 126 | 310 | 469 | 508 | 397 |

These are inventory-entry counts. The linked-unit crosswalk in
inventory_units.csv represents 312 distinct units across 331 inventory rows:
100 canonical stock IDs and 212 candidate-unit IDs. Of the 210 new candidate
units, 30 are screened eligible, 149 are screened ineligible, and 31 remain
eligibility-unresolved. Two candidate-unit IDs have ambiguous or mixed identity
status. The other 574 entries do not yet have a crosswalk unit. Regional
counts by organization are in inventory_progress.csv. These are provisional
audit counts, not the selected expansion batch.

The distinct ICES English Channel sprat unit (`spr.27.7de`) is screened
ineligible for age-disaggregated coverage. The latest 2026 Category 3.3 method
uses aggregate PELTIC biomass and a constant harvest-rate rule. Survey age and
length samples are collected but are not used in the assessment. The 2025
biomass index and 2026 catch advice are documented in
`source_reviews/ices_sprat_english_channel_2026.md`; no canonical data were
added.

The separate Celtic Seas and West of Scotland sprat unit (`spr.27.67a-cf-k`)
is also screened ineligible. Its 2025 Category 5.2 advice uses catch/landings
data only and the published graph contains catch history only. See
`source_reviews/ices_sprat_celtic_seas_2025.md`.

ICES English Channel small-eyed ray (`rje.27.7de`) is a distinct current unit,
separate from the Bristol Channel/Celtic Sea North unit. Its 2024 Category 5
advice uses landings only, with no age-specific population inputs or outputs,
so it is screened ineligible for the age-disaggregated set. See
`source_reviews/ices_small_eyed_ray_7de_2024.md`.

The ICES alfonsino unit (`alf.27.nea`) is screened ineligible: the 2026
Category 3.2 advice reports aggregate landings, a biomass index in kg/day and a
length-based fishing-pressure proxy, without age-specific population inputs or
outputs. The graph's “SSB” heading is identified as `B_index` in its metadata
and is not treated as spawning biomass. See
`source_reviews/ices_alfonsinos_nea_2026.md`.

ICES black-bellied anglerfish (`ank.27.78abd`) is also ineligible for the
age-disaggregated set. The 2026 production run is a length-based age-structured
Stock Synthesis model, with catch and survey observations at length and public
outputs limited to age-0 recruitment, SSB and F ages 3–10. No canonical data
were added. See `source_reviews/ices_black_bellied_anglerfish_78abd_2026.md`.

Both ICES white anglerfish units added from WGBIE 2026 are screened ineligible
for age-disaggregated coverage. Subarea 7/8.a/8.b/8.d (`mon.27.78abd`) is a
length-based age-structured SS3 assessment with length observations; divisions
8.c/9.a (`mon.27.8c9a`) explicitly has no age information. Their public graphs
provide recruitment, SSB and age-group F summaries rather than annual
age-specific population surfaces. See
`source_reviews/ices_white_anglerfish_78abd_2026.md` and
`source_reviews/ices_white_anglerfish_8c9a_2026.md`.

The 2025 NOAA Acadian redfish Gulf of Maine / Georges Bank Run 5 source screen
found native catch and age-composition inputs, spring and fall survey indices
and compositions, age-specific biology, and comparable fitted N-at-age and
F-at-age outputs. It is eligible for source curation; catch-age gaps and the
native survey-index unit code still need to be resolved during import. See
[`noaa_acadian_redfish_2025.md`](source_reviews/noaa_acadian_redfish_2025.md).

The 2025 NOAA Atlantic Bluefish WHAM update (Stock SMART 10388) is eligible
for source curation. It provides digitizable commercial and recreational age
compositions, relative survey-index and age-composition plots, age-specific
biology, and numeric fitted N-at-age and F-at-age tables through 2024. The
SASINF bundle lacks numeric input matrices and model code; exact sample sizes,
the production M vector, and fractional survey timing still require review.
See [`noaa_bluefish_2025.md`](source_reviews/noaa_bluefish_2025.md).

The Gulf of Mexico Gray Snapper candidate maps to the accepted 2022 SEDAR 75
operational Stock Synthesis assessment through 2020. Its sources include fishery
age compositions, age-0 and age-1 FWRI abundance indices, an age-specific fixed M
vector, and annual aggregate biomass, SSB, abundance, and recruitment outputs.
The age-composition inputs and biological surfaces need source-by-source
transcription, and annual survey timing will require documented approximations.
See [`noaa_gray_snapper_2022.md`](source_reviews/noaa_gray_snapper_2022.md).

The Gulf of Mexico vermilion snapper candidate is eligible for source curation
using the 2020 SEDAR 67 operational Stock Synthesis assessment through 2017. Its
report prints the final model inputs and numeric annual aggregate biomass, egg
production, abundance, and recruitment. The 2023 Stock SMART row is a monitoring
update, and the separate Gulf interim analysis used G-FISHER to consider catch
advice; it is not a new accepted age-structured run. The report does not provide
age-specific survey abundance indices or numeric N-at-age/F-at-age series. See
[`noaa_vermilion_snapper_gulf_2020.md`](source_reviews/noaa_vermilion_snapper_gulf_2020.md).
The Gulf Council's February 2025 scheduling summary recommends an update in
2027; it is a future recommendation, not a completed assessment. A separate
August 2025 vermilion-snapper schedule entry concerns the South Atlantic unit.

NOAA Stock SMART Gulf king mackerel (10171) maps to the accepted SEDAR 38 Update
through FY 2017-2018. The report documents conditional age-at-length inputs and
reports recruitment, egg-based spawning output, exploitation rate, and
landings, but numeric run-matched age-composition arrays were not recovered.
Growth and length-weight parameters appear to support deriving weights-at-age;
eligibility remains unresolved pending the source arrays needed to verify
catch-at-age reconstruction. SEDAR 99 Gulf assessment is still underway; its
October 2026 schedule lists a
final report for September 2027. Do not mix its preliminary source papers with
SEDAR 38 Update outputs. See
[`noaa_gulf_king_mackerel_sedar38_2020.md`](source_reviews/noaa_gulf_king_mackerel_sedar38_2020.md).

NOAA Stock SMART Southern Atlantic king mackerel (10172) is a separate unit
from the Gulf stock, despite their shared winter mixing zone. Its latest
completed run is the accepted 2020 SEDAR 38 Update through FY 2017. The report
provides annual egg-based spawning output, recruitment and exploitation, and
documents conditional age-at-length data, but the run-matched numeric
age-at-length and length-composition inputs needed to verify catch-at-age
reconstruction were not recovered. Growth and length-weight parameters appear
to support deriving weights-at-age. Eligibility remains unresolved. A draft February 2026 terms of
reference plans a future update through FY 2026–2027. See
[`noaa_south_atlantic_king_mackerel_sedar38_2020.md`](source_reviews/noaa_south_atlantic_king_mackerel_sedar38_2020.md).

NOAA Stock SMART Scamp rows 10324 (Gulf) and 10326 (South Atlantic) are
separate regional units, although both SEDAR 68OA assessments model Scamp and
Yellowmouth Grouper jointly. The Gulf's 2025 update contains projections based on
the 2022 fitted assessment; the South Atlantic fit ends in 2021 with 2023 SSC
advice. Both reports provide age-specific population outputs, but their numerical
age-composition inputs and run-matched biology remain unresolved. See
`source_reviews/noaa_scamp_sedar68oa_2022.md`; the source PDFs are cached locally
with checksums.

NOAA Stock SMART Yellowedge Grouper row 10323 maps to one Gulf of Mexico
management unit. The 2023 SEDAR 85 assessment fits commercial and bottom-longline
survey age compositions, aggregate abundance indices, and age-based biology; its
2024 SSC outputs and projections derive from that fit. Age compositions are plotted
and sample counts are tabulated, but complete numerical composition arrays and
run-matched model files were not recovered. Eligibility remains unresolved. See
`source_reviews/noaa_yellowedge_sedar85_2023.md`; the reviewed PDFs are cached locally
with checksums.

NOAA Stock SMART Gulf Gag Grouper row 10348 is a separate unit from South Atlantic
Gag (row 10349). The latest accepted production run located is the SEDAR 72
September 2022 amendment using calibrated Florida SRFS private-mode catch inputs;
its inputs and outputs cover 1963–2019. Annual catch and index tables, age-composition
figures and sample counts, fixed age-specific biology, and SSB/recruitment outputs
support source curation. Exact composition arrays and native run files were not
recovered, so any age data import will need a documented figure-based reconstruction.
SEDAR 105 is underway in 2026 and has not superseded the completed production run.
See `source_reviews/noaa_gag_sedar72_2022.md`; five reviewed PDFs are cached locally
with checksums.

NOAA Stock SMART Atlantic yellowfin and bigeye entries 10166 and 10169 are
crosswalked to their single Atlantic-wide ICCAT assessment units. Yellowfin's
conditional age-at-length observations inform growth, while bigeye's reported
age-specific exploitation rates are derived from fitted states; neither report
provides a complete age-disaggregated observation matrix. Their source data may
still support a limited reconstruction, so eligibility remains unresolved. See
`source_reviews/noaa_iccat_atlantic_tunas_2024_2025.md`; the four ICCAT reports
are cached locally with checksums.

NOAA Stock SMART 10165 and the DFO Western Atlantic bluefin listing are
crosswalked as one ICCAT unit. The 2021 Stock Synthesis assessment was used for
advice, but the report's numeric catch-at-age and N/F-at-age tables belong to a
separate VPA rejected for quantitative status advice. The advice model's
conditional age-at-length inputs could support a reconstruction if its numeric
run-matched data are recovered; current eligibility remains unresolved. See
`source_reviews/iccat_west_atlantic_bluefin_2021.md`.
NOAA Stock SMART 10158 is the Western Atlantic skipjack unit assessed by ICCAT
in 2022. Although the Stock Synthesis model is age-structured internally, its
fitted observations are catch in mass, catch-at-size/size compositions and
aggregate CPUE; the official data-preparatory report says catch-at-age was not
planned. It is screened ineligible for the current age-disaggregated source
set. The assessment reports age-specific fitted quantities graphically, but
those states do not replace age-based observations. See
`source_reviews/iccat_west_atlantic_skipjack_2022.md`.

NOAA Windowpane Stock SMART units 10192 (Gulf of Maine/Georges Bank) and 10193
(Southern New England/Mid-Atlantic) are separate current units. The 2023 northern
update uses an empirical swept-area biomass and catch/biomass exploitation method;
the southern update uses aggregate AIM observations. Neither uses age-disaggregated
assessment inputs or reports age-specific population outputs. See
`source_reviews/noaa_windowpane_2023.md`.

NOAA Stock SMART 10160 maps to the IATTC Eastern Pacific skipjack benchmark
assessment conducted in 2024. The assessment uses quarterly fleet catches,
aggregate abundance/biomass indices and length compositions, while its internal
age classes are derived from assumed growth because reliable aging is
unavailable. No observed catch-at-age, age compositions or age-specific survey
series were found, so the unit is ineligible for the age-disaggregated input
set. The report and NOAA summary provide official population outputs; the
generated IATTC run files carry a restriction against redistribution and were
not cached. See
`source_reviews/noaa_eastern_pacific_skipjack_2024.md`.

NOAA Gulf red snapper (Stock SMART 10230) is a new, eligible source-curation
candidate based on the 2026 SEDAR 98 final report, whose terminal data year is
2023. The assessment uses total catches and age-at-length compositions and
reports annual biomass, SSB, age-1+ abundance, recruitment, and exploitation;
numeric N-at-age and F-at-age outputs were not found in the published report.
The report says Council SSC review remains pending, and Stock SMART's 2021
assessment metadata are stale for this run. This is one Gulf unit with three
assessment areas, distinct from Southern Atlantic red snapper. See
`source_reviews/noaa_gulf_red_snapper_2026.md`.

NOAA Gulf gray triggerfish (Stock SMART 10062) is a new, eligible source-
curation candidate based on SEDAR 100, with data through 2024. The assessment
uses fishery conditional age-at-length data and aggregate indices; age-specific
catch and survey inputs for tinyAM would require documented reconstruction.
The Review Panel revised-base report includes annual biomass, SSB, age-1+
abundance, recruitment and a biomass harvest-rate output, but not annual
N-at-age or F-at-age arrays. NOAA's peer-review listing has no actual
dissemination date and formal management adoption was not confirmed in the
sources reviewed. See
`source_reviews/noaa_gulf_gray_triggerfish_2026.md`.

`inventory_progress.csv` also reports distinct matched stock IDs,
crosswalked candidate units and entries without an identity link by source
region. Region names follow each source inventory; ICES entries with multiple
ecoregions remain grouped under their original combined label. Cross-agency
units can appear once under each organization's region, so use
`inventory_units.csv` for global linked-unit counts.
The DFO Maritimes herring component rows are crosswalked to the existing
4VWX canonical management-unit record, not counted as three additional stocks.
Current component indicators and historical age-based material do not establish
a run-matched age-structured production input/output set. The Southwest New
Brunswick migrant-juvenile fishery is a mixed aggregation, not a separate
assessed stock. See `source_reviews/dfo_4vwx_herring_component_crosswalk_2025.md`.

The DFO Atlantic mackerel inventory row matches the existing northern-contingent
record. Its pinned 2025 run has age-based catch inputs and fitted N/F-at-age,
SSB and recruitment outputs; the aggregate egg-production index is not directly
supported by tinyAM, and output uncertainty remains partial. The NOAA mackerel
record covers a broader assessment unit and remains separate; see
`source_reviews/dfo_northern_mackerel_2025.md`.
The DFO Western Component Pollock entry is a distinct unit (4Xopqrs + 5Yb + 5Zc)
from the overlapping NOAA 5YZ assessment. Its accepted 2022 ASAP model and 2026
research reports document catch-at-age, two age-specific survey series, biology,
time-varying M and population outputs through 2020. It is eligible for source
curation, but numeric series are largely figure-based and require source recovery
or documented digitization; nothing has been imported. The 2024 HCR's later SSB
extension is not a continuation of the 2022 population fit. See
`source_reviews/dfo_western_component_pollock_2022.md`.
The DFO 2J+3KL Capelin complex is a distinct unit with fishery catch-at-age
reconstructed from biological samples and age-disaggregated spring acoustic
indices. The current 2026 assessment is index-based and provides no fitted
population N/F/M-at-age or SSB surfaces, so the case is outside the current
age-structured comparison set. Its latest advice reports a 2025 acoustic
biomass index of 766 kt (90% CI 560–1,070 kt); detailed 2025 assessment inputs
are not yet available. See `source_reviews/dfo_capelin_2j3kl_2025.md`.
The separate DFO 4RST Capelin unit is also index-based. Its 2025 assessment
reports mature-biomass and fishing-mortality indices but confirms that age data
are unavailable; it is ineligible for the age-structured comparison set. See
`source_reviews/dfo_capelin_4rst_2025.md`.
The Quebec North Shore herring inventory row is one management-area listing for
spring and fall stocks that DFO assesses separately. The latest update is
index-based, reports no mortality data, and does not provide fitted population
states. It remains identity-ambiguous and ineligible for the age-structured set;
see `source_reviews/dfo_quebec_herring_4s_2025.md`.

NOAA Stock SMART 10573 and the DFO Pacific Sardine management entry are
crosswalked to one northern-subpopulation unit. The 2025 NOAA update confirms
age-composition likelihoods and reports catch and aggregate biomass/recruitment,
but its run-matched numeric age compositions and age-specific outputs were not
recovered. The public Stock Synthesis archive is for the 2024 benchmark and is
not substituted for the latest run. DFO's reviewed advice uses the NOAA biomass
assessment rather than a separate DFO age-structured model. Keep eligibility
unresolved pending the 2025 electronic files; see
`source_reviews/pacific_sardine_northern_subpopulation_2025.md`.

NOAA Stock SMART 10574 (“Southern Pacific Coast”) maps to the 2021 PFMC central-subpopulation northern anchovy assessment adopted in 2022. Its native model archive supplies catches and age compositions, with catch-at-age and age-specific survey abundance reconstructed using report-documented weights and compositions. The report also provides numbers-at-age, biomass, SSB, recruitment and an estimated M. It is eligible for source curation; nothing has been imported. See `source_reviews/noaa_northern_anchovy_central_subpopulation_2021.md`.

NOAA South Atlantic red porgy (Stock SMART 10100) is a new eligible source-curation candidate. The latest completed full assessment found is SEDAR 60 (2020; data through 2017), which reports numeric abundance-at-age, SSB, recruitment and fishing mortality. Catch-at-age and age-specific survey inputs can be reconstructed from the reported fleet totals, age/length compositions, survey indices and fixed biology, but digitizing figures and deriving trawl ages are explicit approximations. No successor assessment was completed in the reviewed 2026 planning documents; see `source_reviews/noaa_south_atlantic_red_porgy_2020.md`.

On 2026-10-29, twelve already-matched ICES records were screened against
their latest reviewed assessments: six are eligible for source curation,
four are ineligible because the current official methods are not
age-structured, and two remain unresolved because current-run age inputs
could not be verified. Their source reviews are linked in the candidate
register. `imported` means some assessment data are already represented in
the canonical tables; it does not mean the assessment record is complete.

Northern Shelf cod (`ICES:cod.27.46a7d20`) now has a completed source screen. Official 2026 advice and SAG summaries are public, but the 2026 run-specific inputs and fitted object are absent from its pinned public TAF repository. The canonical database holds the distinct 2025 detailed event; it is not treated as a substitute for 2026. Current-run eligibility is unresolved. See `source_reviews/ices_cod_north_sea.md` and the gitignored 2026 source cache.

DFO Unit 3 redfish is a provisional new current management unit, but is excluded from the age-disaggregated candidate set: its latest status update uses catch-at-size and survey biomass/abundance-at-length, and the source reports that Acadian and Deepwater redfish are not separated in catches. See source_reviews/dfo_unit3_redfish_2022.md.

Four additional DFO Pacific herring regions (Central Coast, Haida Gwaii,
Prince Rupert District and West Coast of Vancouver Island) are provisionally
new units, screened against their latest regional 2025 assessments and the
SISCAH framework. Public catch and biological-sample inputs are available, but
the abundance likelihood is an aggregate spawn-biomass index and no survey-age
composition was recovered. They are therefore ineligible for the current
age-indexed tinyAM translation; source curation remains possible. See
`source_reviews/dfo_pacific_herring_areas_discovery.md` and the local source
cache manifest.

DFO NAFO 4T yellowtail flounder is a distinct new candidate screened against
the 2025 state-space age assessment. Its catch and survey compositions,
age-specific biology, and official SSB, recruitment and mortality outputs make
it eligible for source curation. Age-specific catch and survey series require
documented reconstruction, and the detailed 2026/054 PDF could not be cached
because its download URL redirected to an archive notice. See
`source_reviews/dfo_yellowtail_flounder_4t_2025.md`.

DFO NAFO 4T white hake is a distinct current unit, separate from the NOAA Gulf
of Maine/Georges Bank stock. Its 2024 assessment uses fishery catch-at-age,
survey age information, age-specific biology, estimated mortality and reported
population trajectories. The detailed report includes age-bin and non-RV
weight-vector gaps that must be resolved before import. See
source_reviews/dfo_white_hake_4t_2022.md.

The DFO Sustainability Survey's broad Rougheye Rockfish row maps to two units in the 2020 assessment: REBS north (5DE) and REBS south (3CD5AB). The published report documents age proportions and model outputs, but complete numeric series were not recovered in machine-readable form. Identity is ambiguous and import eligibility remains unresolved. See `source_reviews/dfo_rougheye_blackspotted_bc_2020.md`.

Fifty-six entries are both provisionally new and eligible (seven DFO, nine ICES and forty NOAA rows); 70 eligible entries are matched to existing records. The NOAA review has screened 176 of 432 inventory rows. NOAA:10231 South Atlantic red snapper is eligible based on the detailed 2024 SEFSC update to SEDAR 73, with reported N-at-age and F-at-age through 2023; its egg-based spawning output and unavailable fractional survey timing are noted in `source_reviews/noaa_south_atlantic_red_snapper_2024.md`. NOAA:10239 Gulf lane snapper is ineligible for the age-disaggregated set because its 2026 SEDAR 49 interim analysis uses an aggregate G-FISHER index and does not directly use ages or lengths; see `source_reviews/noaa_gulf_lane_snapper_2026.md`. The three NOAA Stock SMART rows for Gulf of Alaska rex sole are linked to one candidate stock unit, not counted as three cases. NOAA Pacific Coast rex sole (Stock SMART 10059) is a separate provisional eligible unit: its 2023 locked assessment includes WCGBTS conditional age-at-length data and numeric age-specific outputs, while age-specific catch and survey inputs require documented length-based reconstructions. See `source_reviews/noaa_pacific_rex_sole_2023.md`. The latest 2026 full assessment recommends Model 25.1 and reports numbers-at-age by year, area and sex. The report documents fishery age compositions and survey conditional age-at-length data, but the final-run input archive and the number scale of the age-output table need verification before import; see `source_reviews/noaa_goa_rex_sole_2026.md`. The latest Alaska flatfish batch screened five entries: one new eligible BSAI flathead sole-complex case, three new ineligible cases lacking fishery age inputs, and one existing BSAI arrowtooth match. Northern-rockfish screening added two distinct eligible units: BSAI, whose 2026 report is a harvest projection based on its 2023 assessment, and GOA, whose latest detailed fit located is the 2024 m24 assessment. The GOA Pacific ocean perch candidate is marked eligible: the latest detailed fitted model located is the 2023 Model 20.1 operational update, with input/output files in a pinned NOAA repository revision. Its catch numbers and age-specific survey index can be reconstructed from source totals, age compositions and type-specific weights, subject to explicit likelihood, age-grid and timing approximations. The US West Coast arrowtooth candidate is screened as a distinct stock unit, but its input and age-specific output files remain unresolved; the report documents age and age-at-length data but its named native files were not located. US West Coast NOAA Dover sole (Stock SMART 10022) is a new eligible case: the 2021 full assessment bundle includes native age-composition data, survey age information, aggregate catches and indices, and fitted numbers-, catch-, F- and biomass-at-age outputs. Its age-based tinyAM inputs still require documented reconstruction from compositions and derived weights, and its conditional age-at-length survey likelihood is not a direct age-specific index. In contrast, GOA Dover sole (Stock SMART 10021) is screened ineligible for the current age-input set: the 2023 Model 19.3.1 assessment reports no fishery otoliths, so fishery catch-at-age is unavailable; survey age-at-length samples are conditional and do not supply that missing input. NOAA widow rockfish is another new eligible case: the 2025 full update bundle through 2024 contains fishery age compositions, WCGBTS conditional age-at-length data and age-specific fitted outputs. The assessment report, supplemental review and bundled output give inconsistent Stock Synthesis version labels; the source review records this and identifies the bundled run by its own metadata. NOAA West Coast Yelloweye Rockfish (Stock SMART 10426) is a new eligible source-curation candidate: the adopted 2025 single-sex, two-area assessment includes fishery and survey conditional age-at-length data and native age-specific outputs. It does not provide direct annual catch-at-age or age-specific abundance observations, so tinyAM inputs would need explicit reconstructions. NOAA Cape Cod/Gulf of Maine yellowtail flounder (Stock SMART 10035) is screened as a provisional new unit. The final 2025 peer-review panel accepted its WHAM assessment, but the recovered support archive contains reports and age-surface plots rather than numeric age inputs and outputs, so eligibility remains unresolved; see `source_reviews/noaa_yellowtail_ccgom_2025.md`. NOAA Atlantic halibut is identity-ambiguous with the two DFO halibut inventory units and ineligible for the age-disaggregated set because the current 2024 FSD assessment is an aggregate index-trend method. See `source_reviews/noaa_atlantic_halibut_2024.md`. NOAA BSAI Greenland turbot is provisionally new and unresolved for the age-input screen: the current Model 16.4c uses fishery/survey length compositions, and its available shelf age compositions are excluded from the likelihood; numerical age inputs for a run-matched translation were not recovered. September 2026 alternatives remain preliminary. See `source_reviews/noaa_greenland_turbot_bsai_2024.md`. DFO B.C. coastwide Sablefish is a distinct new unit, screened but unresolved for age-based import: the latest advice uses the 2025-26 management procedure, while the latest detailed age-structured operating model is the 2022 update. Public StRS age records are available, but the commercial trap age-composition arrays and numeric age-specific fitted outputs were not recovered. DFO B.C. coastwide Bocaccio is also a distinct new unit. The latest 2024 update supplies a comparable terminal SSB estimate and uses the detailed 2019 age-structured assessment, but its numeric age-composition arrays are not tabulated and the non-trawl catch stream has no age compositions; age-based import eligibility remains unresolved. NOAA Pacific Coast Canary Rockfish (Stock SMART 10415) is a provisionally new U.S. unit, separate from the Canadian BC stock. Its official 2023 Stock Synthesis archive through 2022 includes catch and biomass-index inputs, age compositions, biological settings and fitted age-specific outputs; age-specific catch and survey series require documented reconstruction. The assessment PDF could not be cached, but the run-matched native archive is available locally. See source_reviews/noaa_canary_rockfish_pacific_coast.md. NOAA Pacific Coast English Sole (Stock SMART 10046) is another new coastwide unit, but its latest full assessment is the 2013 catch-and-index-only data-moderate model. The 2019 update supplied catch-only projections; older age-structured assessments are not substituted for the current run. It is ineligible for this age-specific input set. See source_reviews/noaa_english_sole_pacific_coast.md. NOAA Pacific Coast Darkblotched Rockfish (Stock SMART 10562) is a new eligible source-curation candidate. The 2025 Stock SMART year refers to a catch-only projection; the latest full age-structured assessment is the 2017 update using data through 2016. The cached Stock Synthesis archive has fishery and survey age compositions, aggregate catch and active biomass-index series, biological settings and age-specific model outputs. Age-specific catch/index series require documented reconstruction, and the archive has a fourth survey fleet whose NWCBO data rows are commented out. See source_reviews/noaa_darkblotched_rockfish_pacific_coast.md. NOAA Gulf of Alaska thornyhead (Stock SMART 10457) is provisionally new but ineligible for the age-disaggregated set. The 2024 SAFE assesses the GOA thornyhead complex with a Tier 5 biomass-based REMA model; the age-structured method was suspended in 2003, current fishery age samples are few and none were aged, and the model reports biomass rather than age-specific population outputs. See source_reviews/noaa_goa_thornyhead_complex_2024.md. NOAA Pacific Coast shortspine thornyhead (Stock SMART 10458) is a distinct coastwide unit, provisionally new and ineligible for the age-disaggregated input set. The latest full detailed assessment is the 2023 length-based age-structured Stock Synthesis model; the later 2025 update is catch-only. The 2023 report and native data file explicitly confirm that no age-composition observations were used, despite the model's age-specific derived outputs. See `source_reviews/noaa_shortspine_thornyhead_pacific_coast.md`. NOAA Pacific Coast longspine thornyhead (Stock SMART 10460) is provisionally new but ineligible: the latest full assessment is the 2013 length-based model, which used no age-composition data; the Stock SMART 2019 year is a later catch-only update, and the 2025 SAFE identifies no newer full assessment. A DFO inventory entry for Canadian longspine thornyhead remains separate pending its own screen. See `source_reviews/noaa_longspine_thornyhead_pacific_coast.md`. NOAA Southern California California scorpionfish (Stock SMART 10461) is a provisional new eligible source-curation candidate: the latest full 2017 assessment has fishery length compositions, NWFSC conditional age-at-length observations for 2005–2016, biomass indices and age-specific fitted outputs. Reconstructing catch-at-age and an age-specific index requires an explicit age-length-key approximation; no fishery age-composition likelihood is supplied. See `source_reviews/noaa_california_scorpionfish_southern_california_2017.md`. NOAA Oregon kelp greenling (Stock SMART 10462) is a provisional new eligible source-curation case. Its 2015 assessment report prints the run-specific SS input files, including catch and age-composition data, biology and fixed sex-specific M; aggregate SSB and recruitment outputs are reported. The three CPUE indices are fishery-dependent, no fitted numeric age surfaces were located, and the prose and data-file header disagree on the age grid. See `source_reviews/noaa_oregon_kelp_greenling_2015.md`. Several latest assessments remain unresolved where age-specific inputs for the selected production run were not recovered. These counts do not represent a selected batch, and new entries may still overlap across organizations. Region-level counts and the audit trail are in `inventory_progress.csv` and `candidates.csv`. Screening and identity reconciliation remain incomplete.

The 2024 NOAA Atlantic cod screen reconciles Stock SMART 12805 Eastern Georges Bank to the existing DFO transboundary EGB record and identifies the 2024 WGOM, EGOM and SNE assessments as distinct current units. The three accepted WHAM runs include age-composition inputs and age-specific figures, but the public support archives contain plots rather than numeric age arrays; eligibility remains unresolved. The reviews and cached reports are in `source_reviews/noaa_atlantic_cod_2024_units.md` and `source_reviews/noaa_eastern_georges_bank_2023_identity.md`.

The NOAA Atlantic halibut 2024 screen is ineligible for this age-disaggregated set because its FSD method is aggregate and empirical; its relationship to DFO halibut units is unresolved. See `source_reviews/noaa_atlantic_halibut_2024.md`.
ICES `cod.27.5a` is now matched and screened against the existing Icelandic cod assessment record. Its 2026 MFRI source tables include commercial catch-at-age, two age-disaggregated survey series, age-specific biology, and N/F/SSB/recruitment outputs. Exact 2026 survey timing and production-run settings remain unresolved; the source review also records conflicting maturity-history and terminal SSB values. See `source_reviews/ices_cod_iceland_2026.md`.
Five more matched ICES North Sea rows are now screened against existing 2026 assessments: autumn-spawning herring, haddock, plaice, saithe and whiting. The review confirms direct catch-at-age, age-specific survey information, biological inputs and official fitted outputs for each. Remaining differences are explicit: herring includes partial larval likelihoods and approximate timing; haddock has uncalibrated survey scales and conditional state uncertainty; plaice has approximate timing and unresolved survey units; saithe lacks exact current sampling fractions and its run object; whiting's survey and biology plus groups require expansion in the canonical data. See the corresponding ICES North Sea source reviews.
Four matched DFO inventory entries are now screened against their curated source records: Atlantic cod 3Ps, Atlantic cod 4TVn, 4X5Y haddock and joint coastal Pacific hake. The 3Ps 2025 update and 4TVn/4X5Y newer advice are kept separate from the latest detailed age-structured runs represented by earlier assessment IDs. The Pacific hake 2026 catch-only update is likewise separate from the detailed 2025 fit. The 3Ps, 4TVn and 4X5Y age-specific fitted output coverage is incomplete; Pacific hake has age compositions and biomass outputs but no recovered instantaneous F-at-age series.
Five further matched Iceland/Faroe ICES assessments are screened: Icelandic saithe and haddock (2026), Faroe Plateau cod, haddock and saithe (2025). The current source reviews verify age-based catch and survey inputs and official fitted outputs for those runs. Icelandic data retain unresolved scales/configuration details; the Faroe SAM runs have native objects and stronger run-level checks. See the corresponding species source reviews.

The WGCSE 2026 screening also reconciled five non-age-structured ICES units as provisional distinct stock codes (plaice 7.h-k and 7.b-c, haddock 7.a, cod 6.b, and sole 7.b-c) and retained two megrim entries as identity-ambiguous because their inventory names specify only *Lepidorhombus* spp. Their screening decisions remain ineligible under current methods; none is an import candidate. Counts and row-level rationale are recorded in `candidates.csv` and `inventory_progress.csv`.

The DFO Gulf American plaice inventory entry is matched and screened against the existing 2016 4T VPA record. Catch-at-age and survey indices, source biology, and age-specific N/F/M plus SSB outputs are already represented in the canonical database, although the output coverage remains partial. The 2025 indicator update is kept separate from the detailed population assessment.
ICES Division 7.d plaice is now screened against the 2026 SAM assessment. Native age-specific catch, survey and biology inputs, plus official SAG recruitment, SSB and Fbar, make it eligible for translation review; the native fitted-object endpoint and official N/F-at-age outputs remain unavailable, and the local refit is not treated as an official output.
Northeast Arctic haddock is matched and screened against the 2026 JRN-AFWG SAM run. The final native fit and report verify catch, four survey series, biology, fixed/reconstructed M, and official population/mortality outputs; grouped plus-age observations and the native index scales remain documented limitations.
North Sea Norway pout is matched to the 2026 WKBWNG benchmark and screened eligible. Its native SESAM fit supplies quarterly age-based catch, surveys, fixed age/quarter M, population outputs and fitted-observation predictions; the source's quarterly process and grouped age-2–3+ survey row remain explicit translation limits.
Norwegian spring-spawning herring is matched and screened against the official 2026 WGWIDE run. The source provides numeric catch-at-age, four survey index series, age-specific biological inputs, and age-specific plus aggregate outputs; annual mortality deviations, RFID numeric values and exact fitted settings remain incomplete.
Northeast Arctic golden redfish is matched but screened ineligible for age-structured work because its current 2026 assessment uses aggregate SPiCT, which replaced the historical age-length model. The official relative B/BMSY and F/FMSY series are retained without converting them to absolute quantities or inventing age-specific data.
Gulf of Riga herring is matched and screened against the current 2026 SAM run. Its native fit and report verify catch-at-age, an age-specific acoustic index, annual biology, fixed M, and age-specific plus aggregate outputs; the index scale remains native because no physical multiplier is stated.
ICES Sandeel Areas 1r, 2r and 3r are matched and screened against their 2026 HAWG assessments. Each has seasonal catch-at-age, age-specific survey input, seasonal biology, supplied M and official age/aggregate outputs; source-specific survey intervals, plus-group semantics and unavailable age-state uncertainty remain recorded in the individual reviews.
Sandeel Area 4 is screened eligible using its latest detailed age-based 2025 assessment. The later 2026 category-5 summary remains a separate event and is not mixed into the age-based run's inputs or outputs.
North Sea sole is matched to the current 2026 SS3 assessment and screened unresolved: current-run catch, biology, fixed mortality and fitted age-specific outputs are available, but a run-matched age-specific survey input or age composition was not recovered. Its 2025 age-specific survey series is kept separate rather than substituted into the 2026 run.
Common sole in ICES subdivisions 20-24 is matched and screened eligible using the official 2026 SAM assessment. It provides catch-at-age, three age-specific tuning series, source biology and fitted N/F-at-age. The exact survey timing and precision settings and the named native run object remain unavailable, so any translation must expose those approximations.
Northeast Atlantic blue whiting is matched and screened eligible against the official 2026 BW-2026 SAM fit. The native object supplies catch and IBWSS numbers-at-age, model biology and fitted population outputs; catch and survey units differ, and advice-adjusted recruitment is kept separate from the fitted value.
Icelandic summer-spawning herring is screened unresolved. The canonical 2026 record contains age-specific MFRI tables through 2023 and age-specific population estimates through 2024, but their correspondence to the exact production SAM input mapping, infection-related mortality and later report catch/acoustic series is not verified.
Baltic sprat is matched and screened eligible against the 2026 SAM assessment, which reports age-specific catch, four survey fleets, annual biology and mortality, and population outputs. Sampling times are approximated; age-specific uncertainty and a fitted q surface are unavailable.
Southern horse mackerel remains unresolved against its latest detailed 2024 assessment: age-specific catch and IBTS inputs plus aggregate outputs are available, but exact survey timing, numeric commercial CPUE and age-specific N/F outputs were not recovered. The later 2026 SAG series is kept summary-only.
Northeast Arctic cod is matched and screened eligible using the current official 2026 JRN-AFWG fit, although the ICES SID inventory row still labels its assessment year as 2021. The native run verifies age-based catch, five age-specific indices, biology, timing, fixed M and fitted population surfaces; the distinct advice RCT3 recruitment forecast is not substituted for the fitted SAM recruitment.
Eastern English Channel sole is screened eligible based on the current detailed 2026 assessment, which reports catch-at-age, multiple age-specific survey indices, biology and age-specific and aggregate outputs. Three aggregate commercial biomass indices and several SAM variance details remain translation limitations; a source-based tinyAM fit has not converged, which is tracked separately from source eligibility.
Western Baltic spring-spawning herring is screened eligible against the official 2026 SAM run. The source fit supplies catch, four survey series, biology, spawning fractions and fitted outputs; its survey plus groups differ from the population ages, so the translation's 0–5+ aggregation remains a documented approximation.
Northeast Arctic beaked redfish is matched to the current 2026 JRN-AFWG SCAA run rather than the older ICES SID assessment year. Its age data and population outputs are imported, but exact catch fleet assembly, some survey likelihood selections and index scaling remain unresolved, so it is not yet classified eligible.

NOAA starry flounder required a stock-boundary reconciliation: legacy ID 10025 is a broad 2005 assessment entry, while current IDs 14804 and 14805 represent Northern and Southern assessment units. The current 2017 DB-SRA remains the basis for 2027-28 harvest specifications but uses aggregate catch and life-history inputs, without age-disaggregated observations or outputs. The legacy ID is retained as an ambiguous historical row, not counted as a third current stock; all three entries are screened ineligible for the age-structured candidate set. See `source_reviews/noaa_starry_flounder.md`.

The NOAA winter flounder screen reviewed separate current units. Gulf of Maine (ID 10028) uses a 2025 area-swept 30+ cm biomass model and is ineligible for age-structured comparison. Southern New England/Mid-Atlantic (ID 10030) uses a 2025 ASAP update, but the available report omits the run's age-specific catch and survey arrays, numerical biology and timing files; eligibility remains unresolved. Georges Bank is matched to an existing canonical record. See `source_reviews/noaa_winter_flounder_gom_sne.md`.

NOAA West Coast Sablefish (ID 10466) is screened as a distinct eligible unit. Its 2025 native model bundle contains fishery and survey age compositions, empirical weight-at-age, monthly survey timing, and fitted N/F-at-age outputs. Catch-at-age and age-specific survey indices require explicit transformations from the source observations. The report and native output identify different Stock Synthesis version strings; see `source_reviews/noaa_west_coast_sablefish_2025.md`.

NOAA West Coast Yelloweye Rockfish (ID 10426) is screened as a distinct eligible unit. Its 2025 report and native bundle use the same Stock Synthesis version and provide age-composition inputs, source biology, timing and age-specific fitted output. The assessment is single-sex with two areas; catch-at-age and age-specific survey indices are not direct observed inputs and would require documented reconstruction. See `source_reviews/noaa_west_coast_yelloweye_2025.md`.

NOAA Northern Yellowtail Rockfish (ID 10393) is screened as a distinct eligible unit. Its 2025 single-area, two-sex benchmark includes fishery and survey age data, age-based maturity, estimated mortality, and native fitted age-specific outputs. Catch-at-age and age-indexed observations would require reconstruction. See `source_reviews/noaa_yellowtail_north_2025.md`.

NOAA Rougheye/Blackspotted Rockfish (ID 14516) is screened as one eligible U.S. West Coast complex assessment, not two species-level cases. Its 2025 model includes fishery conditional age-at-length data, survey age information, and age-specific fitted outputs. Catch-at-age and age-specific survey indices would require reconstruction. See `source_reviews/noaa_rougheye_blackspotted_2025.md`.

NOAA California Quillback Rockfish (ID 15342) is a new eligible California assessment unit. Its 2025 assessment has fishery age compositions, fishery-independent survey age data, biological inputs, and fitted age-specific outputs. Catch-at-age and age-specific indices would require reconstruction. See `source_reviews/noaa_quillback_california_2025.md`.

Recent DFO screens also reconciled Strait of Georgia Pacific herring to its existing `dfo_herring_sog` record. Although its database inputs include fishery age compositions and age-specific model outputs, the official survey input is aggregate spawning biomass and native survey age samples are absent; it is therefore ineligible for the current age-indexed tinyAM translation screen.

The three DFO Pacific ocean perch inventory entries map to the existing coastwide `dfo_pop_bc_2023` assessment and represent its 5ABC, 3CD and 5DE subareas, not three new stocks. The detailed assessment has fishery and survey age compositions, aggregate survey indices, fixed biological inputs, and official sex-specific mortality and terminal-biomass outputs. A source-supported age-indexed translation is possible, but the complete matching coastwide historical output series and biomass definition remain unresolved; see `source_reviews/dfo_pop_bc.md`.

The DFO Canary Rockfish inventory row matches the curated coastwide `dfo_canary_bc_2022` Run24/24.01 assessment. Catch totals, fishery and survey age compositions, fixed biology, survey timing, a female SSB posterior subset and sex-specific M outputs are documented; see `source_reviews/dfo_canary_bc.md` for translation limits.

Northwest Atlantic spiny dogfish is a transboundary population represented by separate NOAA and DFO inventory rows. The latest NOAA management-track report states age and growth data are lacking and reports aggregate spawning output, F and recruitment; DFO's accepted model is stage-based, and its later assessment model was not accepted. The DFO 4VWX5 and NOAA Atlantic Coast definitions are not one-to-one, so the cross-agency identity remains unresolved while both are screened ineligible for the age-structured candidate pool. See `source_reviews/noaa_dfo_spiny_dogfish_northwest_atlantic.md`.

Norway pout in ICES Division 6.a (`r`nop.27.6a`) is distinct from the North Sea/Skagerrak stock already in the database. Official ICES scheduling identifies its 2024 advice output as stock status only and schedules the next assessment for 2027; the row is screened out of the current age-structured candidate pool without asserting that age observations are absent from all source datasets. See `source_reviews/ices_norway_pout_6a.md`.

ICES North Sea horse mackerel (`hom.27.4bc7d`) is a distinct current unit with a 2026 SAM assessment and public SSB, recruitment and F summaries. The detailed 2024 benchmark gives age-specific catch and biological inputs for its prior run, but those values are not substituted for the latest run; its stock key also changed from `hom.27.3a4bc7d`, so Division 3.a requires boundary verification. Eligibility remains unresolved pending current-run inputs and survey age information. See `source_reviews/ices_north_sea_horse_mackerel_2026.md`.

ICES whiting in Division 3.a (`whg.27.3a`) is a separate stock unit from the North Sea, Irish Sea and West of Scotland whiting entries. Its 2026 official assessment is Category 3 and the SAG source reports an aggregate biomass index, not age-specific population surfaces; it is screened out of the age-structured candidate pool without claiming that raw age samples do not exist elsewhere. See `source_reviews/ices_whiting_3a_2026.md`.

DFO/NAFO redfish in Divisions 3LN and 3O are separate managed units. The latest detailed 3LN assessment located is the 2024 index-based status assessment; its 2025 operating-model update is MSE work, not a new production fit. The latest detailed 3O assessment is the 2025 survey-index assessment. Both are screened out of the age-structured candidate pool because their current official methods do not report age-specific population surfaces. See `source_reviews/dfo_redfish_3ln_2024.md` and `source_reviews/dfo_redfish_3o_2025.md`.

DFO's two redfish inventory rows represent distinct stock species (*S. mentella* and *S. fasciatus*) across management Units 1 and 2. The latest detailed report is Research Document 2026/037 for the 2024 assessment. The approach is index-based and uses management-unit landings and survey indices; it reports SSB and small-fish recruitment biomass, but no age-structured population inputs or outputs. Science Response 2026/009 updates 2025 indicators, but no Unit 2 survey was available for that year. See `source_reviews/dfo_redfish_units_1_2_2024.md` and the ignored local source cache.

The Southern Gulf spring- and fall-spawning herring entries also match existing, separate stock records. Both have detailed 2024 age-specific inputs and outputs, while 2026 summary-only events are preserved separately. The fall mobile-gear catch-at-age table duplicates the fixed-gear table and is flagged rather than treated as verified data. See `source_reviews/dfo_herring_4tvn_spring.md` and `source_reviews/dfo_herring_4tvn_fall.md`.

NOAA Stock SMART rows for Gulf of Alaska northern and southern rock sole are crosswalked by species: legacy Gulf-wide IDs 10836 and 10040 and current central/western components 15412/15413 and 15414/15415 represent two species-specific assessment units, not six stocks. The latest full models are from 2021; later monitoring records continue through 2024, and the 2025 SAFE catch report says no assessment was conducted that year. The documented models use total catch and length compositions rather than fishery catch-at-age, so both units are screened ineligible for the current tinyAM age-structured translation criteria. See `source_reviews/noaa_goa_northern_rock_sole.md` and `source_reviews/noaa_goa_southern_rock_sole.md`.

The DFO Maritimes Haddock 5Zjm entry and NOAA Stock SMART 12806 are crosswalked
to candidate unit `haddock_eastern_georges_bank`. NOAA Stock SMART 10520 remains
matched to the broader U.S. Georges Bank assessment in the canonical database.
The DFO unit's latest assessment outputs are public, but numerical age-input
arrays and age-specific outputs for the assessment through 2024 remain unresolved
pending access to the detailed 2026/048 report. See
`source_reviews/dfo_noaa_eastern_georges_haddock.md`.

`inventory_progress.csv` summarizes the candidate register by organization and source region. NOAA regions are the science-centre labels in the pinned Stock SMART snapshot; DFO regions are the original Sustainability Survey regions; ICES regions are the exact SID ecoregion labels, including combinations. A row with multiple ICES regions remains one inventory row. These groupings are source reporting categories, not standardized ocean basins.

Counts are **inventory entries**, not deduplicated stock units. `matched`, `r`new`, `ambiguous` and `identity_unreviewed` describe current identity-review decisions; a `r`new` label is provisional until all potential cross-organization matches and boundary changes are resolved. Eligibility and screening columns describe separate decisions and must not be added to identity counts. An unresolved or unstarted entry is not an eligible stock or an exclusion. Because regional listings and identity and source reviews remain incomplete, the snapshot cannot establish final organization coverage or the number of eligible new stocks.

`eligible_new` counts entries provisionally marked both `r`new` and `eligible`; it is not a distinct-unit count. Use `inventory_units.csv` to count only units with an explicit crosswalk ID. `eligible` alone also includes matched existing assessments. Identity, screening and eligibility are separate decisions: an entry can have an evidence screen while its cross-inventory identity remains unresolved.

Regenerate the regional snapshot with `scripts/database/195_summarize_discovery.R` and the linked-unit crosswalk with `scripts/database/196_summarize_inventory_units.R` after updating the register. Both use the cached pinned NOAA metadata described above; the regional summary checks that each status set sums to the source-list count, and the unit summary checks that each linked entry appears exactly once. The candidate register remains the audit trail for individual decisions.

The Newfoundland and Labrador DFO inventory's 4R spring- and fall-spawner rows now reflect the current 4RSw assessment boundary. The 2024 assessment is a WHAM state-space catch-at-age model, but the detailed framework/model report is still cited as in preparation; current-run numerical inputs and age-specific fitted outputs remain unresolved. Spring and fall components remain separate candidate units. See `source_reviews/dfo_herring_4rsw_2024.md`; the six supporting PDFs and checksum manifest are in the ignored local source cache.

The DFO Sustainability Survey's NAFO 3LNO yellowtail flounder entry is a distinct current stock unit, but the latest detailed 2025 NAFO assessment is a Bayesian surplus-production model. It uses annual catch and aggregate survey biomass indices, and the report states age distributions have been unavailable since aging stopped in 2001. It is screened ineligible for the age-disaggregated set; see `source_reviews/dfo_yellowtail_flounder_3lno_2025.md` and the ignored cache `source_cache/discovery_dfo_yellowtail_3lno_2025/`.
The DFO Sustainability Survey's Witch Flounder - NAFO 3NO entry refers to a distinct NAFO stock unit. The latest detailed full assessment located is the 2024 Bayesian surplus-production model; 2026 advice updates status and says the next full assessment is due in 2028. It does not use age-based likelihoods, and Canadian survey ageing has been unavailable since 1994. See `source_reviews/dfo_witch_flounder_3no_2024.md` and the ignored cache `source_cache/discovery_dfo_witch_flounder_3no_2026/`.
The DFO Maritimes NAFO 5Z yellowtail inventory entry is crosswalked to the existing combined U.S.-Canada Georges Bank record and NOAA Stock SMART 10033. It is one stock unit, not a second assessment case. The current approved 2025 WHAM run remains only partially represented in the database; its run-matched age inputs and age-specific output surfaces are unresolved. See source_reviews/dfo_yellowtail_5z_crosswalk.md.

The DFO Newfoundland and Labrador 2J3KLPs herring inventory row is screened
as identity-ambiguous and ineligible for the age-structured population-model
candidate set. The 2025 detailed assessment describes commercial numbers-at-age
and age-4 research-gillnet indices, but the current official assessment is
index-based and reports no fitted population-at-age or mortality outputs. The
broad management area spans five separate 3KLPs complexes and includes 2J,
which is outside those assessment complexes. See
`source_reviews/dfo_nl_herring_3klps_2025.md`; the detailed 2025 report, current
update and framework PDF are cached locally with checksums.


NOAA Stock SMART 10063 is the South Atlantic gray triggerfish unit, separate from Gulf unit 10062. SEDAR 82 provides a detailed 2024 research-track BAM with age inputs and reported egg-production/F outputs, but the review panel described it as work to carry forward to an operational assessment; current production adoption is unconfirmed. Its egg-based measure called SSB is not biomass SSB, and numeric N-at-age/F-at-age tables were not verified. It remains unresolved for the production-assessment set. See `source_reviews/noaa_south_atlantic_gray_triggerfish_2024.md`; source PDFs and checksums are in the ignored local cache.


NOAA Stock SMART 10207 is the Gulf greater amberjack unit, provisionally new relative to the canonical database. Its latest detailed operational fit located is SEDAR 70 (2020; terminal year 2018), accepted by the Gulf Council SSC and still cited by the 2023 framework action. SEDAR reports annual aggregate abundance, biomass and mature-female SSB with uncertainty, but the public age-composition report gives figures rather than numeric input matrices; no species-specific survey age composition was verified. Eligibility remains unresolved pending those age inputs. See `source_reviews/noaa_gulf_greater_amberjack_2021.md` and the ignored local source cache.

NOAA Stock SMART 11682 and the DFO Pacific albacore listing are crosswalked to one North Pacific biological stock and the 2026 ISC benchmark (data through 2024), despite the older inventory metadata. The current assessment fits catch totals, aggregate fishery CPUE and length compositions, not observed age-specific inputs; its modeled age structure and age-specific plotted outputs do not supply those missing observations. Both entries are screened ineligible for the current age-input tinyAM set. See `source_reviews/noaa_dfo_north_pacific_albacore_2026.md` and the ignored cached report under `source_cache/discovery_north_pacific_albacore/`.

NOAA Stock SMART 10163 is the U.S.-managed Pacific Coast Pacific mackerel
unit. The latest Council-adopted 2025 catch-only update retains age-composition
and survey likelihoods from the 2023 age-structured benchmark, alongside
updated catches through model year 2023. Fishery and survey age compositions,
weights, maturity, mortality assumptions, survey timing and biomass outputs
are documented; age compositions end in 2021 and annual numeric N-at-age and
F-at-age outputs were not located. This is a new eligible source-curation
candidate, not a canonical import. The cached Stock SMART metadata conflicts
with the detailed current report, so the run-specific report and adoption
record determine eligibility. See
`source_reviews/noaa_pacific_mackerel_2025.md`.


NOAA Stock SMART 10168 is the Eastern Pacific yellowfin tuna unit assessed
by IATTC in the 2025 benchmark using data through 2023. The integrated model
is age-structured internally but fits catches, aggregate indices and length
compositions; the age-at-length information informs external growth, not
fishery-specific catch-age observations. It is ineligible for the current
age-disaggregated input set. The report gives probabilistic status results,
while run-matched numeric N/F-at-age outputs were not recovered. IATTC
restricts redistribution of generated run files; only the public report is
cached. See `source_reviews/noaa_eastern_pacific_yellowfin_2025.md`.

NOAA Stock SMART 10413 is the Southern Pacific Coast bocaccio unit. Its latest
fitted assessment is the 2017 update, carrying forward age-at-length data from
the 2015 assessment; the 2025 document is a catch-only projection, not a new
fit. The source likelihood uses conditional age-at-length observations rather
than direct marginal catch-at-age, alongside aggregate abundance indices and
fixed biological relationships. It is eligible for source curation, but a
catch-at-age reconstruction requires verification of the native age and
expansion records. Numerical annual N-at-age and F-at-age tables were not
located. See `source_reviews/noaa_bocaccio_southern_pacific_2017.md`.

NOAA Stock SMART 10512 is the separate Bogoslof pollock unit. Its accepted
2024 production assessment uses Tier 5 and aggregate acoustic-trawl biomass;
the report's historical age table supports separate research on M and is not
an age-structured production likelihood. A September 2026 SAFE chapter is
later draft evidence for 2027-2028 specifications. NOAA Stock SMART 10514 is
the separate Southeast Outside pollock unit east of 140° W; its 2024
assessment uses a Tier-5 random-effects model of aggregate bottom-trawl
biomass, while the age-structured C/W/WYK model applies west of 140° W. Both
are screened ineligible for the age-disaggregated candidate pool, not as
claims that all supporting age data are absent. See
`source_reviews/noaa_bogoslof_pollock_2024.md` and
`source_reviews/noaa_pollock_southeast_outside_2024.md`.

NOAA's northern and southern silver hake entries (Stock SMART 10521 and 10522) are distinct U.S. management units. Their latest detailed update is the 2023 empirical catch-and-survey-biomass assessment through 2022; survey numbers-at-age appear in older documentation but are not current model likelihood inputs, and fishery catch-at-age and population-at-age outputs were not verified. The DFO 4VWX Scotian Shelf silver hake record is kept separate; its latest 2025 update (published in 2026) uses a logistic biomass-dynamics method. All three are screened ineligible for the age-disaggregated candidate set. Historical 2012 DFO framework and assessment PDFs are cached with hashes; direct local retrieval of current NOAA/DFO PDFs was blocked, so only complete source files are retained. See `source_reviews/noaa_silver_hake_2023.md` and `source_reviews/dfo_silver_hake_4vwx_2025.md`.

The DFO thorny-skate entry is crosswalked to the single NAFO 3LNOPs biological
stock, while retaining its 3LNO source-list name and the separate Canadian
management treatment of 3Ps. The latest detailed assessment located is the
qualitative 2024 NAFO review; 2026 advice could not be updated because the
designated-expert position was vacant. It does not use age-disaggregated inputs
or report population-at-age outputs. The DFO white-hake entry is crosswalked to
NAFO 3NOPs; its latest detailed full assessment is the 2023 qualitative
catch-and-survey review, and no full assessment or monitoring report was
completed in 2026. Both are screened ineligible for the current age-disaggregated
candidate set. Their advice sheets and detailed supporting documents are cached locally.
The thorny-skate report is marked DRAFT, and the white-hake report says not to cite it without prior reference to the authors; see
`source_reviews/dfo_nafo_thorny_skate_2026.md` and
`source_reviews/dfo_nafo_white_hake_2026.md`.

The DFO Inside and Outside Yelloweye Rockfish inventory entries are now screened as separate provisional units. Inside has survey age compositions in the detailed 2021 rebuilding analysis, but its commercial age sample is from one trip; a representative annual catch-at-age series and several run-matched inputs and outputs remain unverified. Outside received a 2024 status update published in 2026, but the linked official PDF returned an archive placeholder during review. Both remain eligibility-unresolved and no canonical data were imported. See `source_reviews/dfo_yelloweye_inside_2020.md` and `source_reviews/dfo_yelloweye_outside_2024.md`.

DFO Yellowmouth Rockfish has been screened as a distinct coastwide provisional unit. The latest detailed assessment found is the 2021 two-sex Stock Synthesis assessment. It documents annual catch, commercial and survey age compositions, four aggregate survey indices, age-specific biology and official SSB/recruitment outputs. The source curation can proceed with documented reconstructions, but exact age matrices and fractional survey timing remain to be recovered; numeric N-at-age and instantaneous F-at-age outputs were not found. See `source_reviews/dfo_yellowmouth_2021.md`.
The DFO Outside Pacific Spiny Dogfish inventory entry is screened ineligible: the latest 2024 assessment is age-structured but fits aggregate catch and indices with length compositions, while available age samples are not linked to fishery or survey events. No canonical values were imported. See source_reviews/dfo_dogfish_outside_2024.md.

NOAA South Atlantic Spanish mackerel (Stock SMART 10174) maps to the fully accepted SEDAR 78 operational BAM assessment, with data through 2020. The report supplies numeric landings, age compositions, abundance indices, age-specific biology, and N-at-age, F-at-age, SSB and recruitment outputs. Handline and cast-net age compositions are pooled across years; commercial catch totals are weights while other removals are numbers. Survey timing is not one fixed fraction, and SSB timing descriptions conflict within the report. It is eligible for source curation, with these limitations to resolve during extraction. See `source_reviews/noaa_south_atlantic_spanish_mackerel_2022.md`. No canonical values have been imported.

DFO NAFO 4T winter flounder is eligible for source curation from its latest full assessment, the 2017 statistical catch-at-age model (published as DFO Research Document 2019/026). It reports 1973-2016 catch and RV survey age proportions, weights-at-age, maturity-at-age, estimated M, and fitted population summaries. Catch numbers-at-age require a documented reconstruction from tonnes, age proportions, and annual fishery weights-at-age. The 2025 interim update is indicators-only and says a new full assessment was planned for 2026; verify whether that assessment has since been published before importing. See `source_reviews/dfo_sgsl_winter_flounder_4t_2017.md`.

The DFO Pacific Dogfish inside-water unit is screened ineligible for the age-disaggregated candidate set. The latest detailed inside assessment located is the 2010 surplus-production analysis; its authors explicitly state that catch-age data were too sparse for an age-structured model. The 2024 assessment is for outside waters. See `source_reviews/dfo_dogfish_inside_2010.md` and the ignored local cache `source_cache/discovery_dfo_dogfish_inside/`.

NOAA Stock SMART 10190 is the PFMC coastwide Pacific sanddab unit. Its detailed 2013 age-structured Stock Synthesis assessment includes retained-fish age compositions and age-at-length survey data, but the Council did not use its biomass scale for harvest specifications; status determination only was accepted. The latest detailed public evidence does not support total catch-at-age because discards were not age-sampled, and the current harvest specification remains Category 3. Eligibility is unresolved pending source-supported catch-age reconstruction and run-matched inputs. See `source_reviews/noaa_pacific_sanddab_2013.md`; source PDFs and checksums are cached under the ignored local source cache.

ICES megrim `meg.27.7b-k8abd` is a distinct Bay of Biscay/Celtic Seas unit with a 2025 a4a catch-at-age assessment and advice for 2026. WGBIE tabulates two age-specific survey-index series and annual aggregate outputs, but the numerical run-matched catch-at-age matrix is only shown in figures and remains unrecovered. The 2025 report links the revised assessment to advice key 21222; its WD12 says the initial base run was rejected before corrections, so exact run linkage remains to be verified. The downloadable stock annex still says it was last updated in 2018 despite the report describing later updates. Eligibility remains unresolved; see `source_reviews/ices_megrim_7b_k8abd_2025.md` and the ignored source cache.

NOAA Butterfish (Stock SMART 10187) is a provisional new eligible unit. The accepted 2024 Management Track archive provides numerical total-catch proportions by weight, weights-at-age, four seasonal age-index series, and annual SSB, fully selected F, recruitment, and CVs through 2023. The 2021 Research Track report documents maturity ogives for its final ASAP3 run, while the 2024 production update uses WHAM; do not assume those exact ogives are the updated WHAM inputs. The 2024 report states M = 1.278 was carried forward. Current WHAM maturity linkage, complete survey timing, and the 2024 fitted object remain source follow-ups. The panel accepted the assessment for use despite the archived report retaining a draft cover label. See `source_reviews/noaa_butterfish_2024.md` and the ignored local source cache.


DFO Outside Lingcod is a distinct current unit covering PMFC Areas 3CD + 5ABCDE, excluding inside Area 4B. The latest full assessment (2025) used catch, four abundance indices, and age compositions from three fleets and three surveys. The public Science Advisory Report provides fleet catch totals and posterior summary outputs, but the detailed Research Document is still cited as in preparation; numerical age-composition inputs, survey series, and fitted age-specific output arrays were not located. Eligibility remains unresolved pending those source details. See source_reviews/dfo_lingcod_outside_2025.md and the ignored local source cache.

NOAA Atlantic herring (Stock SMART 10572) is screened against the 2026 WHAM Management Track assessment for the U.S. Gulf of Maine/Georges Bank stock complex. Its ignored run archive contains fixed- and mobile-fleet catch-age proportions for ages 1-8 in 1987-2025, plotted survey and biology inputs, and current assessment outputs. Numeric fleet catch totals, survey indices and compositions, full biological arrays, and exact timing remain unresolved; the candidate is not matched to DFO's 4VWX herring components. See `source_reviews/noaa_atlantic_herring_2026.md`.

NOAA black sea bass Stock SMART units 10310 (South Atlantic) and 10311 (Mid-Atlantic) are distinct stocks with age-structured assessment models. The South Atlantic SEDAR 76 Update reports N-at-age and F-at-age through 2023, but its annual age-composition input matrices remain unverified and a revision is underway. The accepted Mid-Atlantic 2025 multi-WHAM update reports aggregate indicators through 2024; its final input package and age-specific output arrays were not recovered. Both remain eligibility-unresolved. See `source_reviews/noaa_black_sea_bass_south_atlantic_2025.md` and `source_reviews/noaa_black_sea_bass_mid_atlantic_2025.md`.

Recent discovery source reviews: [ICES Iberian sardine](source_reviews/ices_iberian_sardine_2025.md), [ICES northern hake](source_reviews/ices_northern_hake_2026.md), [ICES turbot in Division 3.a](source_reviews/ices_turbot_3a_2026.md), [NOAA northern and southern Windowpane](source_reviews/noaa_windowpane_2023.md), [NOAA Atlantic yellowfin and bigeye](source_reviews/noaa_iccat_atlantic_tunas_2024_2025.md), [NOAA Gulf of Mexico red grouper (Stock SMART 10331)](source_reviews/noaa_red_grouper_gom_2025.md), [NOAA Gulf Spanish mackerel (Stock SMART 10173)](source_reviews/noaa_gulf_spanish_mackerel_2023.md), [DFO Gulf witch flounder, NAFO 4RST](source_reviews/dfo_witch_flounder_4rst_2022.md), [DFO Inside Quillback](source_reviews/dfo_quillback_inside_2021.md), [DFO Outside Quillback](source_reviews/dfo_quillback_outside_2021.md), [DFO British Columbia Longspine Thornyhead](source_reviews/dfo_longspine_thornyhead_bc.md), and [DFO Greenland halibut in NAFO Subarea 0](source_reviews/dfo_greenland_halibut_0a0b_2025.md).


ICES turbot in Division 3.a (`tur.27.3a`) is a distinct current advice unit with a 2026 SPiCT assessment. The model uses total catches and three quarterly aggregate exploitable-biomass indices; its use of length data to construct those indices does not provide age-disaggregated likelihood inputs. It is screened ineligible for the age-specific source set, while its relative biomass and fishing-pressure outputs remain available for aggregate comparisons. The WGNSSK report also notes that biological connectivity with North Sea and Baltic populations remains under study. See `source_reviews/ices_turbot_3a_2026.md`.


ICES northern hake (`hke.27.3a46-8abd`) is assessed with a sex-separated Stock Synthesis model. Its current 2026 observations are catches, aggregate abundance indices, and length-frequency data, without catch-at-age or age-specific survey series. Public outputs include historical recruitment, SSB, aggregate F1-7 and biomass; the report shows numbers-at-age only for a short-term forecast. The detailed TAF repository cited in the report was not publicly accessible during review. The current assessment is ineligible for the age-disaggregated set; see `source_reviews/ices_northern_hake_2026.md`.

ICES southern hake (`hke.27.8c9a`) is assessed with a 2026 length-based production model. Standard Graphs reports age-0 recruitment, SSB, and F averaged over ages 1-7 with 90% intervals, but no historical N-at-age, F-at-age, or M-at-age surfaces. It is ineligible for the age-disaggregated set; see `source_reviews/ices_southern_hake_8c9a_2026.md`.

ICES southern black-bellied anglerfish (`ank.27.8c9a`) is assessed with a 2026 surplus-production model. Its graph provides aggregate recruitment, SSB, catch, and F with 90% intervals, not age-specific inputs or population outputs. It is ineligible for the age-disaggregated set; see `source_reviews/ices_black_bellied_anglerfish_8c9a_2026.md`.


ICES Iberian sardine (pil.27.8c9a) is assessed using a 2025 Stock Synthesis run. Its March 2026 Stock Annex documents catch-at-age, spring acoustic age indices and compositions, age-specific biology, and fixed M; Standard Graphs provides recruitment, B1+, and aggregate F outputs. The actual annual input matrices and numeric fractional-year survey timing were not recovered because the assessment TAF repository was not publicly accessible. Eligibility remains unresolved pending those data. See source_reviews/ices_iberian_sardine_2025.md.

NOAA Stock SMART row 10349 is the separate South Atlantic Gag Grouper unit. Its
latest completed assessment located is SEDAR 71 (2021; data through 2019);
SEDAR 106, planned with terminal year 2024, was delayed and has no new schedule.
The report provides numeric aggregate catch and index series, age-composition
figures, biological inputs and age-specific abundance/biomass outputs. Exact
fishery composition arrays and native fit files were not recovered, so
catch-at-age requires documented figure digitization. Its numeric chevron-trap
age compositions informed the SERFS video selectivity but were not themselves
fitted as an age-specific abundance index. This unit must remain separate from
Gulf Gag Grouper (Stock SMART 10348). See
source_reviews/noaa_gag_sedar71_2021.md.

NOAA Stock SMART row 10339 is a distinct South Atlantic Snowy Grouper unit. The latest completed assessment located is SEDAR 36 Update, revised January 2021, using data through 2018. Its report provides removal totals, aggregate indices, age-composition figures and sample sizes, age-specific biology, and numeric abundance-at-age, biomass-at-age, recruitment, SSB, F-at-age, and estimated-removals-at-age outputs. Reconstructing observed catch-at-age requires digitizing composition figures. The 2013 SEDAR 36 WP15 model files are for the superseded benchmark and are not a substitute for update inputs. See source_reviews/noaa_snowy_grouper_sedar36u_2021.md.

NOAA Stock SMART blueline tilefish entries 10381 (South Atlantic) and 14922 (Mid-Atlantic) are distinct management-unit candidates with nonidentical boundaries relative to SEDAR 92. The 2025 southern ASPIC assessment was rejected for management; the replacement southern DLM and northern DLM use aggregate removals, length/life-history data, and TAC methods, not age-disaggregated likelihood observations or population-at-age outputs. Both are screened ineligible for the age-input set. See source_reviews/noaa_blueline_tilefish_sedar92_2025.md; four official assessment and SSC PDFs are cached locally with checksums.

NOAA Stock SMART identifies South Atlantic (10384; SEDAR 89, 2025 run through 2022) and Mid-Atlantic (10385; 2024 WHAM run through 2023) golden tilefish as separate, fully accepted units. Both assessments use age compositions and age-specific biological inputs and report population-at-age or related outputs, so each is eligible for source curation. The South Atlantic input/output matrices are partly tabulated and partly figure-based; the Mid-Atlantic series are largely chart-based in the public SASINF reports. Exact numeric age matrices still require extraction, and the Mid-Atlantic report is publicly labelled draft despite the peer-review panel's acceptance for management. See `source_reviews/noaa_golden_tilefish_2024_2025.md`; the reviewed PDFs and SASINF archive are cached locally with checksums.

ICES sandeel Area 5r (san.sa.5r) is screened ineligible for the age-disaggregated set. Its 2025 advice is Category 5; the HAWG report lists the Area 5r Category 5.4 analysis as last updated in 2019 and describes it as catch-only. The published graph provides annual catches but no age-specific population inputs or outputs. See source_reviews/ices_sandeel_area_5r_2025.md.

ICES sandeel Area 6 (san.sa.6) and Area 7r (san.sa.7r) are screened ineligible for the age-disaggregated set. The latest advice records located are from 2023 and cover 2023–2024; the 2026 HAWG report lists their catch-only Category 5.2 and 5.3 analyses as last updated in 2019. The associated stock graphs provide catch entries only. See `source_reviews/ices_sandeel_area_6_2023.md` and `source_reviews/ices_sandeel_area_7r_2023.md`.


ICES capelin cap.27.2a514 is a new distinct unit with unresolved age-data eligibility. The 2026 MFRI/GINR update is survey-index based; current public age tables cover selected recent surveys, but a run-matched annual age-specific acoustic input series and mortality assumptions remain unverified. See source_reviews/ices_capelin_2a514_2026.md.

ICES cod in Kattegat (`cod.27.21`) is a distinct 2026 Category 3, age-based SAM
trend assessment. Current public outputs provide relative recruitment, SSB and
mortality summaries, but the run-matched age-disaggregated inputs remain
unverified, so eligibility is unresolved. See
`source_reviews/ices_cod_kattegat_2026.md`.

ICES greater silver smelt in Subarea 14 and Division 5.a (`aru.27.5a14`) is a distinct eligible candidate for source curation. The 2026 MFRI report documents the age-structured inputs, but the annual matrices are mostly figure-based; official outputs include recruitment, SSB and F summaries with uncertainty. See `source_reviews/ices_greater_silver_smelt_5a14_2026.md`.

ICES ling in Division 5.a (`lin.27.5a`) is a distinct, eligible candidate for source curation. The 2026 MFRI report documents a SAM assessment with catch-at-age, three age-specific survey series, biology, and fixed M = 0.15; catch-age inputs and much of the age-specific data are figure-based, with a documented sampling gap. See `source_reviews/ices_ling_5a_2026.md`.

ICES blue ling in Subarea 14 and Division 5.a (`bli.27.5a14`) is a distinct current stock screened ineligible for the age-disaggregated set. The latest assessment is length-based GADGET; the report says reliable age composition and weights-at-age are unavailable, while maturity is modeled at length. See `source_reviews/ices_blue_ling_5a14_2026.md`.

ICES lemon sole in Subarea 4 and divisions 3.a and 7.d (lem.27.3a47d) is a distinct current stock screened ineligible under the 2026 advised method. The WGNSSK report describes a SURBAR relative-SSB stock-size indicator and chr advice rule, not an age-structured population model. See source_reviews/ices_lemon_sole_3a47d_2026.md.
NOAA Stock SMART IDs 10079, 10081, and 10083 are now linked to their distinct Gulf of Mexico brown, pink, and white shrimp units. They remain outside the initial marine-finfish scope and are not assessment-input candidates. See source_reviews/noaa_gulf_shrimp_inventory_2026.md.
NOAA Stock SMART IDs 10080–14842 covering 30 shrimp, crab, lobster, scallop, clam, and squid management units are now linked to distinct candidate units. They are screened out only because they fall outside the initial marine-finfish scope; their assessment inputs and outputs have not been audited. See source_reviews/noaa_invertebrate_inventory_2026.md.
ICES grey gurnard (`gug.27.3a47d`), brill (`bll.27.3a47de`), and North Sea flounder (`fle.27.420`) are current distinct stocks screened ineligible under their 2026 advised methods. The report describes a category-3 trends method for grey gurnard and SPiCT assessments for brill and flounder, without age-disaggregated population inputs or outputs. See `source_reviews/ices_wgnssk_n_age_methods_2026.md`.
ICES SID stock codes for 32 Norway lobster and northern shrimp entries are now linked to distinct candidate units. They remain screened outside the marine-finfish scope; assessment inputs and outputs have not been audited. See `source_reviews/ices_invertebrate_inventory_2026.md`.
