# Assessment scripts

Scripts are grouped by purpose. Run them from the repository root.

## `database/`

Use these scripts to seed, import, review, validate, and export assessment records. The fit-readiness check reports whether each record is ready for translation.

## `translation/`

`run_translations.R` loads the reviewed database and runs the matching stock specifications in `translation/stocks/`. Each stock file defines one `translate_stock()` function with its observation conversion, fit settings, comparison scales, and background table.

Run the translation workflow and its focused converter checks with:

```sh
Rscript analysis/comp_assessments/scripts/translation/run_translations.R
Rscript analysis/comp_assessments/tests/test_translation.R
```
