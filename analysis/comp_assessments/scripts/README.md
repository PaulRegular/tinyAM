# Assessment scripts

Scripts are grouped by purpose.

## `database/`

Seed, import, review, validate, and export assessment records here. The fit-readiness check also belongs here because it checks whether the recorded inputs are complete enough to translate. Run these scripts from the repository root.

For example:

```sh
Rscript analysis/comp_assessments/scripts/database/002_validate_database.R
Rscript analysis/comp_assessments/scripts/database/003_fit_readiness.R
```

## `translation/`

Converter checks and focused, stock-specific translation scripts belong here. Run the converter check from the repository root:

```sh
Rscript analysis/comp_assessments/scripts/translation/004_test_translation.R
```