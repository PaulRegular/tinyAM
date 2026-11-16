# Curated comparison: SAIO November 16, 2026

This snapshot records a complete rerun of the curated tinyAM recipes. See the [comparative assessment workflow](../../README.md) for instructions to reproduce the fits and figures.

- SAIO release: [`saio-2026-11-16-curation`](https://github.com/PaulRegular/saio/releases/tag/saio-2026-11-16-curation).
- SAIO [manifest](https://github.com/PaulRegular/saio/blob/14e86b223f04c376486b610a625efd2f35a48130/database/releases/saio-2026-11-16-curation.csv) commit: `14e86b223f04c376486b610a625efd2f35a48130`; manifest blob and result key: `00dd3b63d9934f01425ad1fa678f751a539c827b`.
- tinyAM revision used for fitting and plotting: [`f204832f`](https://github.com/PaulRegular/tinyAM/commit/f204832fe18429eac8b081243edc23c77da286ec).
- Prior comparison: November 12 results on the historical SAIO revision `c36cc22df29dc0efb9a96b7e2500caf4bec36c86`.

Of 41 cases, 28 converged, seven failed numerical checks, one reached the 20-minute limit, and five were blocked by source inputs or reference data. Every case had the same outcome as in the previous batch. No result crossed the material-change thresholds on common years, units, and definitions: 10 percentage points in median absolute percentage difference, 0.10 in trend correlation, or 20% in median or terminal fitted SSB change. The largest change in median absolute percentage difference across 101 like-for-like comparisons was 0.001 percentage points.

- Figures: [Agreement heatmap](agreement_heatmap.png), [SSB](ssb.png), [recruitment](recruitment.png), [abundance](abundance.png), and [average F](F_bar.png).
- Tables: [Fit diagnostics](fit_diagnostics.csv), [agreement summaries](comparison_summary.csv), and [comparison to the previous run](comparison_to_previous.csv).

The figures compare tinyAM fits with official assessment results assembled by [SAIO](https://github.com/PaulRegular/saio). See [original assessment sources](assessment_sources.csv) for attribution. Fitted objects, annual comparison records, and local archives remain under the local results workflow.
