# Gulf of Bothnia herring, 2026

## Assessment selected

The database record represents ICES assessment key 22377 for the combined Gulf
of Bothnia stock `her.27.3031` (subdivisions 30 and 31). WGBFAS 2026 reports
the production Stock Synthesis assessment, version 3.30.22.1. It fits total
catch through 2025, uses catch age compositions from 1980, and reports a 2026
forecast. The database keeps the fitted history through 2025 and excludes the
forecast. Older catalogue rows for subdivisions 30 and 31 remain separate;
they are not treated as this current combined assessment.

## Sources and local cache

- [WGBFAS 2026 report, section 4.3](https://literatur.thuenen.de/digbib_extern/dn071152.pdf), especially Tables 4.3.3–4.3.7.
- [ICES SAG source data, key 22377](https://sg.ices.dk/ViewSourceData.aspx?key=22377) and [the annual XML export](https://sag.ices.dk/SAG_API/StockXmlDownload/22377).
- [WKCLuB 2021 benchmark report](https://doi.org/10.17895/ices.pub.5989), used for the fixed natural-mortality age vector.
- [IBPCluB 2019 report](https://jukuri.luke.fi/bitstreams/adfbdc03-5e3b-45f2-8492-ec5bf7085724/download) provides the preceding inter-benchmark context.

The 2026 report PDF, extracted text, SAG XML, acoustic table transcription,
source-page images, and tidy transcriptions of catch-at-age, catch weights,
and maturity are in the gitignored
`analysis/comp_assessments/source_cache/ices_herring_gulf_bothnia_2026/`.
`manifest.csv` records each cached file's URL, checksum, and verification
date. The WKCLuB PDF could not be cached from this shell because the remote
socket was unreachable; an earlier response was HTML, not a PDF, and was
removed. The failure and alternate source URL are recorded in the manifest.
The reported M breakpoints were verified against the benchmark's searchable
published report, but the local cache still lacks that PDF.

## Recovered inputs and assumptions

- ICES SAG total catch is retained in tonnes for 1963–2025. The published
  catch-at-age table is retained in thousand fish for ages 1–14 and 15+,
  1980–2025. Within-year number proportions are reconstructed from those
  counts for the SS3 age-composition input; no composition sample sizes or
  effective likelihood weights were recovered.
- Table 4.3.4 gives annual catch weights in kg per fish for ages 0–16 and
  17+, 1963–2025. A separate stock-weight surface was not recovered. The
  tinyAM recipe uses catch weights as a documented stock-weight proxy and
  applies the reported age-15 weight to its 15+ group.
- Table 4.3.5 gives annual maturity proportions for ages 1–16 and 17+,
  1980–2025. The published ogives are retained without smoothing or
  extrapolation.
- The 2021 benchmark uses a time-invariant, age-varying M curve with values
  0.563, 0.472, 0.332, 0.290, 0.267, and 0.257 at age midpoints 0.5, 1.5,
  3.5, 5.5, 8.5, and 15.5. The database reconstructs ages 0–20 by linear
  interpolation and holds the terminal value above age 15.5. This is a fixed
  benchmark input, not an estimated M output from the 2026 assessment.
- Table 4.3.6 provides the BIAS acoustic abundance index for 2007–2025,
  ages 0–15+, in millions of fish. Published values are already area-corrected;
  the table's `AREA_CORR_FACTOR` is preserved in the local transcription and
  is not applied a second time. The survey samples late September to early
  October, covers most of subdivision 30, and does not survey subdivision
  31. Coverage was reduced in 2021–2025. The report also describes a trapnet
  index for 1990–2006, ages 3–9, but its numeric input series was not
  recoverable from the report figures.
- The report describes hybrid F and an unweighted Fbar over ages 3–7. The
  ICES SAG provides annual recruitment at age 0, SSB, and Fbar for 1963–2025
  with a 90% central interval. It does not provide fitted N-at-age, F-at-age,
  or estimated M-at-age surfaces in the recovered output tables. Therefore
  total abundance and age-specific mortality comparisons cannot be derived.

## Translation limits

The initial tinyAM model uses years 1980–2025 and ages 1–15+, matching the
reported catch-composition age range. Recruitment is therefore represented
at age 1, not the source assessment's age 0, and recruitment is not treated
as a like-for-like comparison. The BIAS index is translated to thousand fish
by the shared unit converter; the script assigns within-year time 0.75 as an
approximation for the reported late-September/early-October survey.

tinyAM reconstructs catches by combining annual total biomass, age-number
proportions, and catch weights, then fits the resulting age-specific catches
with its lognormal likelihood. This does not reproduce SS3's total-catch and
age-composition likelihoods or their unknown precision settings. tinyAM uses
one shared estimated error for the acoustic index and age-specific
catchability; the assessment's exact index weighting is unavailable. SSB is
compared only as an approximate aggregate because tinyAM starts at age 1,
uses catch weight as its stock-weight proxy, and the report does not expose
the exact spawning-time settings. Source Fbar is unweighted, whereas tinyAM's
Fbar is abundance-weighted; it is not a direct comparison.

## Initial tinyAM fit

The package-default initialization converged (optimizer code 0, objective
561.3017, maximum absolute gradient 0.000285, positive-definite Hessian;
26 fixed and 1,365 random effects). The only numeric aggregate comparison is
approximate SSB: mean absolute percent difference 28.6%, terminal difference
+34.6%, and annual trend correlation 0.20 over 1980–2025. This is a sizable
disagreement under the documented differences above; no source observations
or starting values were changed to improve it. Recruitment and Fbar are marked
non-equivalent, and the M match is only the reused fixed input. The cached
dashboard and fitted object are available locally under
`analysis/comp_assessments/results/cache/ices_herring_gulf_bothnia_2026/`.
