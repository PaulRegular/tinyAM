# Curated comparison: SAIO November 16, 2026

This snapshot records the November 26, 2026 rerun of all curated tinyAM recipes after the recruitment argument changed to `R_settings`. See the [comparative assessment workflow](../../README.md) for instructions to reproduce the fits and figures.

- SAIO release: [`saio-2026-11-16-curation`](https://github.com/PaulRegular/saio/releases/tag/saio-2026-11-16-curation).
- SAIO [manifest](https://github.com/PaulRegular/saio/blob/14e86b223f04c376486b610a625efd2f35a48130/database/releases/saio-2026-11-16-curation.csv) commit: `14e86b223f04c376486b610a625efd2f35a48130`; manifest blob and result key: `00dd3b63d9934f01425ad1fa678f751a539c827b`.
- tinyAM revision used for fitting and plotting: [`a05957ba`](https://github.com/PaulRegular/tinyAM/commit/a05957ba6b6883d3b8cf9007dbfe9990d74a2e57).
- Prior comparison: [November 16 results on the same SAIO release](https://github.com/PaulRegular/tinyAM/tree/ff674e12/analysis/comp_assessments/snapshots/saio-2026-11-16-curation).

Of 41 cases, 28 converged, seven failed numerical checks, one reached the 20-minute limit, and five were blocked by source inputs or reference data. Every case had the same outcome as in the previous batch. Across 101 comparisons with matching years, units, and definitions, the largest change in median absolute percentage difference was 0.001 percentage points for British Columbia Yellowtail Rockfish recruitment. Changes stayed below the review thresholds of 10 percentage points in median absolute percentage difference, 0.10 in trend correlation, and 20% in median or terminal fitted SSB.

- Figures: [Agreement heatmap](agreement_heatmap.png), [SSB](ssb.png), [recruitment](recruitment.png), [abundance](abundance.png), and [average F](F_bar.png).
- Tables: [Fit diagnostics](fit_diagnostics.csv), [agreement summaries](comparison_summary.csv), and [comparison to the previous run](comparison_to_previous.csv).

The figures compare tinyAM fits with official assessment results assembled by [SAIO](https://github.com/PaulRegular/saio). See [original assessment sources](assessment_sources.csv) for attribution. Fitted objects, annual comparison records, and local archives remain under the local results workflow.
