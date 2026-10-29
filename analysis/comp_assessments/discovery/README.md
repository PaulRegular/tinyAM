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

Cached files and their `sources.csv` manifests are gitignored. The register retains
inventory source URLs and leaves unknown identities, run details and data evidence
blank. `scripts/database/186_build_discovery_inventory.R` builds the initial
snapshot and refuses to overwrite an existing register, protecting later reviews.
Candidate discovery does not modify the canonical database.

## Inventory progress snapshot

`inventory_progress.csv` summarizes the candidate register by organization and source region. NOAA regions are the science-centre labels in the pinned Stock SMART snapshot; DFO regions are the original Sustainability Survey regions; ICES regions are the exact SID ecoregion labels, including combinations. A row with multiple ICES regions remains one inventory row. These groupings are source reporting categories, not standardized ocean basins.

Counts are **inventory entries**, not deduplicated stock units. `matched`, `new`, `ambiguous` and `identity_unreviewed` describe current identity-review decisions; a `new` label is provisional until all potential cross-organization matches and boundary changes are resolved. Eligibility and screening columns describe separate decisions and must not be added to identity counts. An unresolved or unstarted entry is not an eligible stock or an exclusion. Because regional listings and identity and source reviews remain incomplete, the snapshot cannot establish final organization coverage or the number of eligible new stocks.

Regenerate the snapshot with `scripts/database/195_summarize_discovery.R` after updating the register. It uses the cached pinned NOAA metadata described above and checks that each set of categories sums to the source-list count. The candidate register is the audit trail for individual decisions.
