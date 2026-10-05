# Northern Gulf cod (3Pn4RS) source review

The canonical record represents [DFO Research Document 2026/010](https://publications.gc.ca/site/eng/9.958344/publication.html),
the detailed assessment peer-reviewed in February 2025 and published in
February 2026. Its data and reported estimates run through 2024, although the
English title calls it the assessment in 2025. The [2026/012 Science Advisory
Report](https://publications.gc.ca/collections/collection_2026/mpo-dfo/fs70-7/Fs70-7-2026-012-eng.pdf)
is a newer summary update and does not replace this detailed assessment record.

The [2025/074 framework report](https://publications.gc.ca/collections/collection_2025/mpo-dfo/fs70-5/Fs70-5-2025-074-eng.pdf),
[2022/033 catch and removals report](https://publications.gc.ca/collections/collection_2022/mpo-dfo/fs70-5/Fs70-5-2022-033-eng.pdf),
[2022/049 survey-input report](https://publications.gc.ca/collections/collection_2022/mpo-dfo/fs70-5/Fs70-5-2022-049-eng.pdf),
and [2024/045 weight and model report](https://publications.gc.ca/collections/collection_2024/mpo-dfo/fs70-5/Fs70-5-2024-045-eng.pdf)
are cached with checksums in the local source manifest.

The report provides catch-at-age in thousands of fish for ages 2–11+ in
1974–2024 and model outputs for abundance, biomass, fishing mortality, natural
mortality, recruitment, and two Fbar age ranges. Fixed M values for ages 2–11+
through 1983 are also reported. From 1984 onward, M is estimated for ages 4–11+;
ages 2 and 3 remain fixed. These values and the reported output series are
stored with their original units and age groups. The six survey series and
their year/age coverage are recorded as assumptions.

The assessment uses six age-specific abundance indices: DFO August (1985–2024,
ages 2–11+), Sentinel mobile (1995–2024, ages 2–11+), Sentinel gillnet
(1995–2024, ages 4–11+), Sentinel summer longline (1995–2024, ages 3–11+),
Sentinel fall longline (1995–2020, ages 3–11+), and the Minet bottom-trawl
series (1973–1976, ages 3–11+). The main report describes these observation
series but does not tabulate their numerical values. Public DFO survey files
contain source survey data, but they are not the assessment's processed index
series; reproducing those indices would require the source survey selection,
standardization, age-reading, and calibration steps. They cannot be substituted
with the report's fitted population surfaces. DFO's public [Teleost
survey](https://open.canada.ca/data/en/dataset/40381c35-4849-4f17-a8f3-707aa6a53a9d),
[Cabot survey](https://open.canada.ca/data/en/dataset/7001783a-4dc0-41bc-8932-1dae1e699d91),
[Alfred Needler survey](https://open.canada.ca/data/en/dataset/4eaac443-24a8-4b37-9178-d7cce4eb7c7b),
and [mobile Sentinel survey](https://open.canada.ca/data/en/dataset/929fe07f-ab8e-4b3c-8ee3-1aa7a9ea0b1a)
files were checked as possible sources. These contain underlying survey data,
not the model-ready indices, and do not by themselves recover the assessment's
preprocessing and calibration.

Annual beginning-of-year stock weights and female maturity ogives are also
required for the assessment's biomass calculations. The report describes the
weight and maturity models but does not publish their full numeric age-year
matrices. The separate [2025 maturity technical report](https://publications.gc.ca/collections/collection_2025/mpo-dfo/fs97-6/Fs97-6-3671-eng.pdf)
provides figures and methods rather than a complete numeric matrix. It could
not be retrieved into the local cache because the source download failed in
this environment; no maturity values were transcribed from figures. A full accepted-model input
record and a defensible tinyAM fit therefore remain blocked pending numerical
survey indices, stock weights, and maturity values or the native model data.
No output estimates have been reused as observations or biological inputs.

The cached source files and checksums are listed in
`source_cache/dfo_cod_3pn4rs_2025/manifest.csv`. Table extraction is checked by
`scripts/database/039_validate_dfo_cod_3pn4rs_2025.R`; the general database
validator also passes. This stock is recorded as partial and is not counted as
a completed database-to-tinyAM analysis.
