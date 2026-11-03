# Assessment-source migration history

This history belongs to tinyAM. SAIO is an independent source database: its checkout, instructions, importers, tests and canonical records do not require tinyAM or a translation workflow.

The initial manifest is preserved byte for byte in [saio_initial_migration.json](saio_initial_migration.json). The former SAIO snapshot note is preserved in [SAIO_SNAPSHOTS.md](SAIO_SNAPSHOTS.md). They describe frozen source revisions, not current validation results. The downstream data pins remain in [saio_revisions.csv](../saio_revisions.csv); this cleanup does not advance them.

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
