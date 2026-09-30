# SAM input bridge and compatibility audit

Translating SAM data does not reproduce a SAM assessment. This branch separates
input conversion, mathematical compatibility, and extraction of public reference
results. It does not alter tinyAM's model or defaults.

The first case is North Sea cod from the public
[fishfollower/SAM repository](https://github.com/fishfollower/SAM), pinned to
`c6cfd035c7de59f7b3421dde31901efbca4cb0e8`. Inputs and explicit configuration
changes are in `testmore/nscod/`. A fitted regression-test object is also public
at `stockassessment/tests/nscod/fit.expected.Rdata`; its provenance and agreement
with the selected case must be checked before comparison. It is a repository
example, not necessarily a current ICES advice assessment.

Source review covers `stockassessment/R/reading.R` (`read.ices`, `read.surveys`,
`read.data.files`, `setup.sam.data`), `conf.R` (`defcon`, `loadConf`), `run.R`
(`sam.fit`), `tables.R`, and the C++ implementations under
`stockassessment/inst/include/SAM/`. The
[current package manual](https://fishfollower.r-universe.dev/stockassessment/doc/manual.html)
is a secondary reference; source resolves discrepancies with older documentation.

The bridge must preserve these distinctions:

- Survey values are divided by the row's effort denominator. Sampling time is
  the mean of the two timing endpoints, as in `setup.sam.data`.
- Negative survey entries represent missing values; explicit zeros remain zeros
  in the translated observations. tinyAM later treats zero catch/index values
  as missing, like SAM's setup step.
- q and observation-SD grouping can use factor columns and formulas. SD sharing
  across catch and index tables cannot be enforced by their separate parameter
  vectors.
- Sharing a mean through a formula does not share a latent F state. SAM's
  correlated random-walk innovations are also different from tinyAM's stationary
  AR1 process.
- SAM's Fbar is an arithmetic average over its configured ages; tinyAM's
  reported F_bar is abundance weighted. Compare arithmetic Fbar separately.
- Spawning-time mortality and catch weight must be retained as source metadata;
  they must not silently replace tinyAM's beginning-year SSB or stock-weight
  yield conventions.

No installed `stockassessment` package is required. Offline fixtures will test
the supported file formats and mappings. Public inputs will be downloaded by a
pinned, reproducible workflow; tests must not access the internet.
