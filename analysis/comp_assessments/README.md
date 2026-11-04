# Comparing stock assessments with tinyAM

Source data are maintained in the independent **SAIO — Stock Assessment Inputs and Outputs Database**. This folder contains tinyAM translations and comparisons.

## Curated and preliminary translations

- scripts/translation/curated/: the 35 scripts retained from main at 78f3d347331b0e770d55a855de4c5b90307fd400. This preserves the existing paper working set, including its documented blocked cases. Membership does not imply convergence or scientific equivalence.
- scripts/translation/preliminary/: additional translations and changed alternatives from comp-assessments. They have not been promoted into the paper working set.
- translations.csv: the explicit registry, including each script's originating commit. An assessment can have a curated version and a preliminary alternative.

Curated runs are the default. Preliminary runs require an explicit tier. Extraction completeness, translation review and numerical convergence are separate. Move a script and update the registry only after reviewing its sources, assumptions and comparisons. If both tiers contain the same assessment, selecting its ID without a tier uses the curated version.

The first 20-stock expansion beyond Charbonneau–Keith is tracked in
[expansion_batch_01.csv](expansion_batch_01.csv), using the selected stock units
in SAIO. Its new recipes are preliminary. The registry's `source_commit`
records the origin of migrated scripts; new recipes leave it blank and use Git
history and run diagnostics for their model revision. The SAIO data revision
is recorded separately.

The batch register also records translation limitations. `source_readiness_blocked`
means the reviewed source record does not yet support a defensible fit; it does
not establish that data can never be recovered. `reconstruction_pending` means
native observations are present but their age/length, sex, area or removal
definitions still need a reviewed conversion. Neither status counts as a fit.

## SAIO data

Clone [SAIO](https://github.com/PaulRegular/saio) beside tinyAM. The default lookup uses the primary tinyAM checkout's parent directory, so it also works from Git worktrees. On this laptop the locations are `C:/Users/paulr/Documents/tinyAM` and `C:/Users/paulr/Documents/saio`.

For another location, set SAIO_PATH to your local SAIO checkout, or use:

```r
options(tinyAM.saio_path = "/path/to/saio")
source("analysis/comp_assessments/R/run_assessment.R")
database <- read_database()
```

Both translation tiers read the latest **committed local SAIO main** by default. The reader resolves main once and records that exact commit; it ignores uncommitted table edits and another checked-out branch. Curated membership remains the same 35 scripts, while their inputs now follow current SAIO data.

Update SAIO before an analysis with `git -C /path/to/saio switch main`, `git -C /path/to/saio pull --ff-only` and `git -C /path/to/saio lfs pull`. Reading or fitting never fetches data automatically. If the fetched origin/main is ahead of or diverges from local main, reading current data stops with an update instruction. Locally committed changes ahead of origin/main are included and identified by their actual commit.

For a reproducible historical run, supply `read_database(revision = "<SAIO commit>")` or `read_committed_database(revision = "<SAIO commit>")`. The commit and its LFS objects must be available locally; fetch missing objects explicitly with `git -C /path/to/saio lfs fetch origin <revision>`. To deliberately inspect edited tables, use `read_database(file.path(saio_path(), "database"))`; those reads have a separate content fingerprint.

Diagnostics and cached fits record the actual SAIO revision, selection and tinyAM revision. Recheck affected translations before refreshing paper figures when source data advance. DATABASE_STRUCTURE.md and PROTOCOL.md point to their maintained SAIO documents.

## Run or preview a translation

Open scripts/translation/run_stock.R, select translation_tier and assessment_id, then source it. For a batch, use scripts/translation/run_all.R. Its translation_tier defaults to curated; use preliminary deliberately for the wider library. Batch membership comes from translations.csv, rather than SAIO's is_current flag.

```r
result <- run_assessment("afsc_cod_ebs_2024", tier = "curated",
  database = database, fit = FALSE)
preliminary <- assessment_translations("preliminary")
```

The full scientific translation instructions remain in [TINYAM_TRANSLATION.md](TINYAM_TRANSLATION.md). Stock-specific assumptions and inherited fit reviews are in [translation_notes/](translation_notes/README.md); official source evidence remains in SAIO. Shared converters and dashboard helpers serve both tiers. These are analysis-local functions; package model mathematics and public interfaces are unchanged.

## Local outputs

All generated results are ignored by Git and separated under results/curated/ and results/preliminary/. Fits, dashboards and initialization attempts use `cache/<SAIO commit>/<assessment_id>/`; batch checkpoints share that revision folder. Explicit edited-table reads use a `working-<fingerprint>` key. Previous unversioned caches and historical revision folders are preserved and are never substituted for current fits. Plot scripts reject result tables or caches from another or unknown data revision. Preliminary results never enter curated summaries by default.

Open results/curated/index.html or results/preliminary/index.html for dashboards. Regenerated indexes label the selected source revision separately from historical or unknown-revision outputs, without rewriting old HTML. Fresh previews show source provenance and reject a known mismatch between the supplied fit and assessment context. Older comp-assessments outputs for curated stocks are under curated/legacy_comp_assessments/ and require revalidation before paper use. Older mixed batch summaries remain in preliminary/ and are historical outputs, not a newly completed preliminary batch. results/output_migration.json records the source and destination of relocated files.

## Verification and purpose

Source recovery and translations used substantial agentic-AI assistance. Not all records and translations have been exhaustively human-verified. Check authoritative sources and diagnostic warnings before research use. tinyAM translations are illustrative abstractions, not official assessments or management advice.

The broad library supports future standardized comparative modelling. A stock-specific approximation of an official assessment does not automatically define the common formulation for a depensation study.

## Historical source snapshots

`saio_revisions.csv` retains historical snapshot references, not active defaults. The curated snapshot preserves the original paper-input tables byte for byte; the preliminary snapshot preserves the expanded database at migration. Request one explicitly with `read_database(revision = .saio_snapshot_revision("curated"))`. Old cached estimates keep their original provenance and are not relabeled as results from current data.

## Migration history

The source-database split and its initial snapshot manifest are documented in [migration/](migration/README.md). SAIO is maintained independently; tinyAM owns its translation provenance and migration history.
