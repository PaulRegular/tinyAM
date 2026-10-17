# British Columbia petrale sole — candidate preflight

## Official assessment

The coastwide 2024 assessment is documented in
[Science Advisory Report 2025/002](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/4127961x.pdf).
It uses a two-sex Stock Synthesis catch-at-age model: catches for 1938–2023,
three synoptic survey indices for 2003–2023, and fishery/survey age compositions
for 2004–2019. The official result is an ensemble with fixed M of 0.154, 0.124
and 0.184, weighted 50%, 25% and 25%, respectively. A single base-M fit would
be a component approximation, not the official ensemble.

No final matching detailed research document was located in this search.
The 2024 State of the Pacific Ocean technical report cites the Mazur, Holt,
Fisch and English research document as in press. That citation does not
establish its present publication status. An older public working template
is not used as a production assessment source.

## Recovered numerical output

The [gftrends repository](https://github.com/pbs-assess/gftrends/tree/7f1c59b8cc3895a3477204ac16da6cf15e2466c6)
provides `data-raw/petrale-sole-bc-mcmc-2024.rds` and its processing script.
The cached object has 174,000 rows: **2,000 complete posterior SSB trajectories
for 1938–2024**, with corresponding BMSY and reference points. There are no
duplicate draw/year cells, and all biomass values are finite and positive.
The script extracts `SSB_1938` through `SSB_2024` from medium-, low- and high-M
nuisance-posterior CSVs, joins BMSY by draw and downsamples to 2,000 draws.
The underlying CSVs are not included in the inspected public tree.

| Quantity | Cached 2.5% | Cached median | Cached 97.5% | Official report 2.5% / median / 97.5% |
|---|---:|---:|---:|---|
| Female SSB in 2024 (t) | 10,512.47 | 17,274.80 | 28,107.54 | 10,830 / 17,330 / 27,977 |
| BMSY (t) | 2,798.42 | 5,729.49 | 8,342.80 | 2,793 / 5,704 / 8,323 |

These are close, but not identical, summaries. Downsampling is documented;
the original component sample sizes and the intended ensemble weights still
need confirmation. Do not relabel these as the exact published quantiles.
The 2024 state is the assessment boundary after catches through 2023, not an
additional observed catch year.

## Remaining readiness requirements

The [pbs-assess/petrale repository](https://github.com/pbs-assess/petrale/tree/fc214bd26288dd98633d3e6b2c0c7f1bddca5db0)
contains environmental-analysis scripts for assessment appendices. They read
external SS3 outputs and `mcmcresults.RData`; those objects and the native
production input arrays are not in the inspected tree. The recovered gftrends
object contains no N-at-age, F-at-age or recruitment. Original catch/survey
composition arrays, biological inputs, timing, likelihood weights and final
configuration must be recovered before a faithful tinyAM translation.

BC petrale is a **promising Canadian replacement candidate**, with a numerical
SSB comparison already recoverable. It is not yet an added canonical case,
runnable translation or verified replication. No fit or source-based
observation reconstruction was attempted. Any addition remains a catalogue
exception unless an authoritative catalogue identifier is verified.

## Cache and checks — 2026-10-17

Sources, pinned repository inventories and checksums are stored locally under
`source_cache/dfo_petrale_bc_2024_preflight/`. The validated advisory PDF is the
Waves copy; the publications.gc.ca download returned HTML and is excluded.
Offline checks verify complete year/draw coverage, positive finite SSB/BMSY,
and the 0.4/0.8 BMSY reference-point calculations. No package code or canonical
database table changed.
