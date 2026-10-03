# Assessment scripts

Scripts are grouped by purpose. Run them from the repository root.

## `database/`

Use these scripts to seed, import, review, validate, and export assessment records. The fit-readiness check reports whether each record is ready for translation.

## `translation/`

`translation/run_translations.R` loads the reviewed database and runs current
accepted assessments with matching specifications in `translation/stocks/`.
Each stock file is named for its `assessment_id` and defines one
`translate_stock()` function with its data conversion, fit settings, and
background text. These files are intentionally unnumbered: the assessment ID
identifies the stock, and the driver selects the current records.

Run the translation workflow and its focused converter checks with:

```sh
Rscript analysis/comp_assessments/scripts/translation/run_translations.R
Rscript analysis/comp_assessments/tests/test_translation.R
```
