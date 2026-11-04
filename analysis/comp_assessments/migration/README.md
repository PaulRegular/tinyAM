# Assessment-source migration history

This history belongs to tinyAM. SAIO is an independent source database: its checkout, instructions, importers, tests and canonical records do not require tinyAM or a translation workflow.

The initial manifest is preserved byte for byte in [saio_initial_migration.json](saio_initial_migration.json). The former SAIO snapshot note is preserved in [SAIO_SNAPSHOTS.md](SAIO_SNAPSHOTS.md). They describe frozen source revisions, not current validation results. Historical revisions remain recorded in these migration files and Git history. [saio_revisions.csv](../saio_revisions.csv) now selects one manually refreshed local CSV snapshot shared by both tiers; see the current [workflow](../README.md).

## Former SAIO migration note

## Migration and provenance

The initial snapshot is an exact copy of the canonical CSVs at tinyAM comp-assessments commit 34470212e1102f3f31bf3f08c26a8a72a985c948. migration.json records table sizes and SHA-256 hashes. Historical commits remain available in the [original repository](https://github.com/PaulRegular/tinyAM/tree/34470212e1102f3f31bf3f08c26a8a72a985c948/analysis/comp_assessments); history has not been rewritten. This repository begins with the preserved source snapshot.

The reviewed and preliminary tinyAM translations remain in tinyAM. Their [translation protocol](https://github.com/PaulRegular/tinyAM/blob/main/analysis/comp_assessments/TINYAM_TRANSLATION.md) and [stock model notes](https://github.com/PaulRegular/tinyAM/tree/main/analysis/comp_assessments/translation_notes) are maintained there. Source curation does not require tinyAM.

## Attribution and reuse

Retain the source attribution and any source-specific reuse terms. This migration does not assign a new blanket licence to third-party assessment material. Code originating in tinyAM retains its original attribution and GPL (>= 3) licence, preserved in CODE_LICENSE.md.


## Source-review reconciliation

# Assessment source reviews

Keep official run identity, source definitions, units, transformations, verification and unresolved source gaps here. Downstream tinyAM settings and fit diagnostics live in [tinyAM translation notes](https://github.com/PaulRegular/tinyAM/tree/main/analysis/comp_assessments/translation_notes). Source curation and completeness do not depend on a successful translation.

The 30 pre-migration `legacy_main` review copies were reconciled with their primary reviews during cleanup. The newer recoveries take precedence over superseded row counts and missing-data claims. Distinct 2025 Norwegian herring source details are retained in their own historical review; the 2025 Northern Shelf cod SAG identifiers remain with that run. Earlier review versions remain available in [the initial SAIO snapshot](https://github.com/PaulRegular/saio/tree/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/legacy_main).


## Independent-source cleanup

The source repository was reviewed from SAIO commit `8eeb542`. Project-specific references remained in documentation, generated annotations, discovery decisions and importer controls. The subsequent cleanup removes those references and makes the dry-run controls source-owned (`SAIO_DATABASE_DRY_RUN` and `SAIO_IMPORT_DRY_RUN`). The former `tinyam_surveys` assumption records recovered survey coverage; it is renamed `recovered_survey_coverage`, preserving its value, missing-survey caveat and source citation. Only explanatory text and that recovery-inventory label change in the canonical tables; source observations, estimates, units, age/year keys, identifiers and citations are preserved.

## Source-test follow-up

The five reported SAIO test failures were resolved by validating current discovery accounting and using no-write source-import previews. The inventory test no longer replays historical screening decisions; stock-specific source tests remain in their own files. Inventory membership uses stable candidate identifiers, allowing either candidate or canonical-stock keys after curation.

All five affected tests, the inventory-lookup regression checks and canonical validation passed. During this follow-up, all 900 protected source files—including the five canonical tables, candidate register, selected batch, source reviews and retrieval provenance—remained byte-for-byte unchanged. Only the two derived inventory summaries were refreshed. This does not advance tinyAM's data pins, change its model specifications or rerun any assessment fits.

## Source-review content separation, 2026-11-03

The follow-up content audit covers all 476 Markdown files in SAIO's source-review and discovery-review folders. It separates secondary analysis choices from official assessment facts, including passages whose subjects were lost during the earlier name-only cleanup. Additional model recipes, pooling choices, initializations, fit diagnostics and agreement results are preserved in [translation notes](../translation_notes/README.md), linked to their original SAIO revision `e723863b83e12b2f9ca7142eb147d1ae082bdb98`. These are inherited records, not new fitting results or confirmation of current recipes.

SAIO retains official model methods, source definitions, units, run checks, biological inputs, canonical derived-output calculations and unresolved recovery gaps. Its README describes source inputs and outputs without fitting instructions. Checksums confirm that all 1,330 pre-existing tracked files outside the source reviews and README are unchanged, including every canonical table, discovery decision, importer, source test and provenance file. No data pin, model specification or fitted result changes.

The new source-documentation scope check, six existing source/inventory checks and canonical database validation pass. Documentation links, source citations, code fences and both repository whitespace checks pass. The package suite and `R CMD check` were not repeated for this documentation-only separation.
