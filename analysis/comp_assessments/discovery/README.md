# Assessment discovery inventory

`candidates.csv` is a discovery register, not a canonical assessment table. Its
rows are source-list entries, some reviewed and others not yet reviewed: they are not verified distinct new stock
units, eligible age-structured assessments, or a selected import batch. Screening,
identity matching, source verification and exclusion decisions remain outstanding.
Do not overwrite it with a fresh inventory after review has begun.

The initial snapshot was retrieved on 2026-10-28 from:

| Source | Entries | Local cached source and provenance | Coverage limit |
|---|---:|---|---|
| NOAA Stock SMART via NOAA-EDAB | 432 stock IDs | `source_cache/discovery_noaa_stocksmart/stock_assessment_summary.rda`, revision `664c4237e26520b260a090e1821a2cee5d1cdf07`, SHA-256 `AD47F9EEB240C2E77B2B952809A12C4FA739E3F0D39BEAB283B413FDBAEF76A8` | Latest metadata row per stock ID, not a verified official production run or proof of age observations. Regional official sources require separate review. |
| ICES SID | 274 published stock entries | `source_cache/discovery_ices_sid/stock_list_2026.json`, SHA-256 `913D6A3EEAB9449AE3CD0DD15661A150D36FF79B6E6A9BA78D86D46BA1035F42` | Published 2026 stock list includes assessments and advice that may not use age observations. |
| DFO Sustainability Survey | 198 selected stocks | `source_cache/discovery_dfo_sustainability/survey_2024.html`, SHA-256 `A10E65EEF028C2A524973FFCFA155915D673AE3AB9A9F1312581449E29644A28` | The 2024 survey covers only selected managed stocks, including non-finfish. Supplement with regional CSAS inventories. Official IDs and scientific names require review. |


Regional supplementation: the [2026 DFO Pacific groundfish plan](https://www.pac.dfo-mpo.gc.ca/fm-gp/mplans/ground-fond-ifmp-pgip-sm-eng.html), retrieved 2026-10-28, led to a British Columbia yellowtail rockfish candidate absent from the initial survey list. The cached regional page SHA-256 is `2CB8792B9AEBE6FE5B27F3887BEFF57F02E2E7CE15D9344353AD705AED18475B`. This partial regional review is not an exhaustive CSAS inventory; further regional sources must still be audited.

The [Pacific Fishery Management Council groundfish assessment listing](https://www.pcouncil.org/stock-assessments-star-reports-stat-reports-rebuilding-analyses-terms-of-reference/groundfish-stock-assessment-documents/) was retrieved 2026-10-28; cached HTML SHA-256 `77FEB82D8A674A018FD3FC7E02CE16579B36F3101B8603FF61286927B1D8B32D`. It provides reports and model files for recent West Coast assessments. The chilipepper identity has been reviewed, but the remainder of this regional listing still needs reconciliation and screening.

The [2026 WGBIE working-group report](https://doi.org/10.17895/ices.pub.32790321.v1) is a regional assessment listing and source. A partial review of its megrim and sole sections is recorded in `source_reviews/ices_iberian_wgbie_2026.md`; the remaining sections still need reconciliation. The locally cached sections and working-document annex have individual URLs and checksums in `source_cache/discovery_ices_wgbie_2026/sources.csv`.

The full [2026 WGCSE report](https://doi.org/10.17895/ices.pub.32625759.v1) has been cached and used to screen its method categories. Eleven candidate entries are ineligible under the current non-age-structured methods; 14 category-1 age-based stocks remain for stock-specific source screening. See `source_reviews/ices_wgcse_2026.md` and `source_reviews/ices_sol_27_7a_2026.md`; the remaining age-based sections still require reconciliation and input/output review.

The [2026 WGNSSK report](https://doi.org/10.17895/ices.pub.32676345.v2) is a North Sea and Skagerrak regional assessment listing. Its witch section has been screened and cached; other sections still require reconciliation. See `source_reviews/ices_witch_northsea_2026.md`.
Cached files and their sources.csv manifests are gitignored. The register retains
inventory source URLs and leaves unknown identities, run details and data evidence
blank. `scripts/database/186_build_discovery_inventory.R` builds the initial
snapshot and refuses to overwrite an existing register, protecting later reviews.
Candidate discovery does not modify the canonical database.

## Inventory progress snapshot

Snapshot dated 2026-10-29: the register contains 905 source-list entries. These are not deduplicated stock units.

| Organization | Entries | Matched | Provisional new | Ambiguous | Identity unreviewed | Eligible | Excluded | Eligibility unresolved | Screened | Not started |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DFO | 199 | 15 | 2 | 0 | 182 | 2 | 98 | 99 | 103 | 96 |
| ICES | 274 | 57 | 6 | 0 | 211 | 4 | 44 | 226 | 53 | 221 |
| NOAA | 432 | 22 | 2 | 1 | 407 | 2 | 37 | 393 | 39 | 393 |
| Total | 905 | 94 | 10 | 1 | 800 | 8 | 179 | 718 | 195 | 710 |

Two entries are both provisionally new and eligible (NOAA); the six other eligible entries are already matched to existing DFO or ICES stocks. These counts do not represent a selected batch, and new entries may still overlap across organizations. Region-level counts and the audit trail are in `inventory_progress.csv` and `candidates.csv`. Screening and identity reconciliation remain incomplete.

`inventory_progress.csv` summarizes the candidate register by organization and source region. NOAA regions are the science-centre labels in the pinned Stock SMART snapshot; DFO regions are the original Sustainability Survey regions; ICES regions are the exact SID ecoregion labels, including combinations. A row with multiple ICES regions remains one inventory row. These groupings are source reporting categories, not standardized ocean basins.

Counts are **inventory entries**, not deduplicated stock units. `matched`, `new`, `ambiguous` and `identity_unreviewed` describe current identity-review decisions; a `new` label is provisional until all potential cross-organization matches and boundary changes are resolved. Eligibility and screening columns describe separate decisions and must not be added to identity counts. An unresolved or unstarted entry is not an eligible stock or an exclusion. Because regional listings and identity and source reviews remain incomplete, the snapshot cannot establish final organization coverage or the number of eligible new stocks.

`eligible_new` counts entries provisionally marked both `new` and `eligible`; it is the closest current register measure of potentially eligible new units, but can still include duplicate stock units across organizations until cross-inventory reconciliation is complete. `eligible` alone also includes matched existing assessments.

Regenerate the snapshot with `scripts/database/195_summarize_discovery.R` after updating the register. It uses the cached pinned NOAA metadata described above and checks that each set of categories sums to the source-list count. The candidate register is the audit trail for individual decisions.
