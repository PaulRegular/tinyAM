# Comparing stock assessments with tinyAM

Source data are maintained in the independent **SAIO — Stock Assessment Inputs and Outputs Database**. This folder contains tinyAM translations and comparisons.

## Curated and preliminary translations

- scripts/translation/curated/: the 35 scripts retained from main at 78f3d347331b0e770d55a855de4c5b90307fd400. This preserves the existing paper working set, including its documented blocked cases. Membership does not imply convergence or scientific equivalence.
- scripts/translation/preliminary/: additional translations and changed alternatives from comp-assessments. They have not been promoted into the paper working set.
- translations.csv: the explicit registry, including each script's originating commit. An assessment can have a curated version and a preliminary alternative.

Curated runs are the default. Preliminary runs require an explicit tier. Extraction completeness, translation review and numerical convergence are separate. Move a script and update the registry only after reviewing its sources, assumptions and comparisons. If both tiers contain the same assessment, selecting its ID without a tier uses the curated version.

## SAIO data

Clone [SAIO](https://github.com/PaulRegular/saio) beside tinyAM. The default lookup uses the primary tinyAM checkout's parent directory, so it also works from Git worktrees. On this laptop the locations are `C:/Users/paulr/Documents/tinyAM` and `C:/Users/paulr/Documents/saio`.

For another location, set SAIO_PATH to your local SAIO checkout, or use:

```r
options(tinyAM.saio_path = "/path/to/saio")
source("analysis/comp_assessments/R/run_assessment.R")
database <- read_database()
```

The default reader uses the exact commit for its tier in saio_revisions.csv, including locally available Git LFS objects. It does not silently follow SAIO's current branch or uncommitted changes. Fetch each pinned revision and its LFS objects when setting up another machine, using git -C /path/to/saio lfs fetch origin <revision>. The reader does not require checking out the historical revision. To deliberately inspect a working-tree database, use read_database(file.path(saio_path(), "database")).

Change saio_revisions.csv deliberately and recheck affected curated examples before refreshing paper figures. Diagnostics record the SAIO revision and tinyAM revision independently. DATABASE_STRUCTURE.md and PROTOCOL.md now point to their maintained SAIO documents.

## Run or preview a translation

Open scripts/translation/run_stock.R, select translation_tier and assessment_id, then source it. For a batch, use scripts/translation/run_all.R. Its translation_tier defaults to curated; use preliminary deliberately for the wider library. Batch membership comes from translations.csv, rather than SAIO's is_current flag.

```r
result <- run_assessment("afsc_cod_ebs_2024", tier = "curated",
  database = database, fit = FALSE)
preliminary <- assessment_translations("preliminary")
```

The full scientific translation instructions remain in [TINYAM_TRANSLATION.md](TINYAM_TRANSLATION.md). Stock-specific assumptions and inherited fit reviews are in [translation_notes/](translation_notes/README.md); official source evidence remains in SAIO. Shared converters and dashboard helpers serve both tiers. These are analysis-local functions; package model mathematics and public interfaces are unchanged.

## Local outputs

All generated results are ignored by Git and separated under results/curated/ and results/preliminary/, including cache/, dashboards, initialization attempts, summaries and figures. Plot scripts use the selected tier and filter assessment IDs against its registry. Preliminary results never enter curated summaries by default.

Open results/curated/index.html or results/preliminary/index.html for migrated dashboards. Older comp-assessments outputs for curated stocks are under curated/legacy_comp_assessments/ and require revalidation before paper use. Older mixed batch summaries remain in preliminary/ and are historical outputs, not a newly completed preliminary batch. results/output_migration.json records the source and destination of relocated files.

## Verification and purpose

Source recovery and translations used substantial agentic-AI assistance. Not all records and translations have been exhaustively human-verified. Check authoritative sources and diagnostic warnings before research use. tinyAM translations are illustrative abstractions, not official assessments or management advice.

The broad library supports future standardized comparative modelling. A stock-specific approximation of an official assessment does not automatically define the common formulation for a depensation study.

## Preserved paper inputs

The curated pin preserves the five source tables from the original main working set byte for byte. The preliminary pin uses the expanded SAIO source database. Both are immutable commits in the same SAIO history, not separately maintained databases. Select the same tier for the database and the translation; switching to newer source data is a deliberate scientific update.

## Migration history

The source-database split and its initial snapshot manifest are documented in [migration/](migration/README.md). SAIO is maintained independently; tinyAM owns its data pins, translations and migration history.
