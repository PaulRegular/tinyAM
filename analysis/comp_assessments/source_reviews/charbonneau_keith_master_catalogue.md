# Charbonneau–Keith Master catalogue import

## Source and scope

The metadata snapshot is the `Master.csv` file from the [Age-structured marine fish database](https://github.com/JulieCharbonneau/Age-structured-marine-fish-database), pinned to revision `7e4daca3cbfc14ce7fa0b28ed824eebe7113249b`. It was retrieved on 2026-10-17 from the [revision-specific raw file](https://raw.githubusercontent.com/JulieCharbonneau/Age-structured-marine-fish-database/7e4daca3cbfc14ce7fa0b28ed824eebe7113249b/Master.csv). The local gitignored copy is `source_cache/charbonneau_keith/Master.csv`; its MD5 is `430d85ba1f4144fc87bca098c8090699` and SHA-256 is `52FA158B197574B45B92D485D738E1DC7D62A4904315AEEF90C0EDB526E449ED`. The manifest records the URL, retrieval date, checksum and revision.

The pinned file contains 102 rows and 102 unique `unique_ID` values. This differs from the 101-row count anticipated in the task plan; the import follows the source file and represents all 102 identifiers once. This is catalogue metadata only. The catalogue’s historical model flags and labels do not establish that any record is a current official production assessment.

## Import and crosswalk

The importer is `scripts/database/071_import_charbonneau_master.R`. It uses the source scientific name, assessor code, `unique_ID` area suffix and PDF label. It does not import the catalogue’s model flags, estimated quantities or status fields. All 102 source identifiers are now present exactly once in `stocks.csv`: 27 were linked to existing stock rows and 75 were added as metadata-only rows, increasing the stock table from 37 to 112 rows. The other four canonical database tables were not changed.

For the 27 existing records, the Master `unique_ID` is now the `charbonneau_id`; any previous identifier is retained in the stock note as an alias. The reviewed crosswalk is explicit in the importer and requires both scientific-name and regional-basin agreement. The importer rejects stale or duplicate matches rather than dropping them silently. Its focused test checks the crosswalk, area preservation, provenance notes, duplicate protection, rerun stability and that imports leave assessment, assumption, input and output tables unchanged.

New `stock_id` values are stable slugs of the Master `unique_ID`, prefixed with `ck_`. Common names and authority-specific stock codes are left blank when the catalogue does not provide them. Area is taken from the `unique_ID` suffix, retaining the source text. Regional basin is assigned from assessor family: AFSC to Northeast Pacific; ICES to Northeast Atlantic; and DFO/NEFSC families to Northwest Atlantic. AFSC rows use `Alaska` as region; other rows use the assigned basin as their broad region. These regional labels are a reproducible catalogue classification, not a claim that the Master file supplies a normalized basin field. No Northwest Pacific entries occur in this snapshot.

## Identity limitations retained

- Nine catalogue entries are marked as possible overlaps in their stock notes and remain separate because exact identity was not established: three non-sex-specific cases (NEFSC summer flounder area, NEFSC 5Y haddock variant and ICES western Baltic herring subdivisions 22–24) and six sex-specific entries for AFSC arrowtooth flounder, Alaska plaice and northern rock sole. The sex-specific rows remain separate catalogue metadata; they are not counted as additional assessment cases.
- `ICES-WGCSE_merlangus_CS7a` has the source assessor code `ICES-WGSCE`. Both values are retained in its note rather than silently correcting the identifier.
- `AFSC_ macrocephalus_BSAI` is linked to the existing Eastern Bering Sea cod record based on its legacy identifier and PDF label, although the Master area suffix says `BSAI`. This discrepancy is explicitly recorded in the stock note; confirm it against an assessment source before using the catalogue row as assessment identity.

These exceptions are visible in `stocks.csv` and the importer crosswalk. They do not add assessments, assumptions, numerical inputs or outputs. Stock-specific production assessment curation remains governed by the corresponding stock source review.
