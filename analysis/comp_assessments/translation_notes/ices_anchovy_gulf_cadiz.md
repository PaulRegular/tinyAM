# Gulf of Cadiz and southern Portugal anchovy — translation notes

## Additional inherited translation context

Preserved from the [SAIO review at `e723863`](https://github.com/PaulRegular/saio/blob/e723863b83e12b2f9ca7142eb147d1ae082bdb98/source_reviews/ices_anchovy_gulf_cadiz.md) during the 2026-11-03 content separation. These passages describe earlier analysis choices or checks; they were not refitted or revalidated as current recipes. Use the registered stock script and `translations.csv` for current settings and status.

Missing current
empirical stock weights and native quarterly compositions prevent a
defensible replication-oriented fit. A future annual approximation must
also address July settlement, strong seasonal mortality and seasonal
survey/composition timing. It must not pretend July recruitment occurred
on 1 January, or use the older MSE states to manufacture observations. No optimization or starting-value retry is justified before those input
gaps are resolved. All three canonical source-completeness statuses remain
partial.

The stock recipe stops at source readiness for both preparation and fitting;
it runs no optimizer and constructs no starting values. A standalone local
dashboard at
`results/cache/ices_anchovy_gulf_cadiz_2025/dashboard_official_only.html`
shows the available official SSB, July recruitment and age-3 F with source
uncertainty. Missing observations, population patterns and other states
remain unavailable rather than being reconstructed from the older MSE run.
