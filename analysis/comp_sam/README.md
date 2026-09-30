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

## Run the example

From the repository root, after installing tinyAM (or loading the development
package with `pkgload::load_all()`):

```r
source("analysis/comp_sam/001_download.R")
source("analysis/comp_sam/002_nscod.R")
# Optional: validate and refit with an installed stockassessment package.
source("analysis/comp_sam/003_validate_sam.R")
```

The first script verifies 14 source files against `source_manifest.csv`, including
the public fitted object. Raw downloads and serialized working objects are
ignored by Git; readable input, audit, reference and diagnostic CSVs are retained
under `results/`. `nscod.cfg` records the pinned `defcon()` defaults with the
explicit changes in `testmore/nscod/script.R`. It excludes the alternate
`scriptcc.R` prediction-dependent observation-variance model.

The core bridge reads ICES table codes 1/2/3/5, standard survey blocks, and SAM's
`$field` configuration files. Unknown configuration fields are retained and
marked `not_checked`; absent settings are never filled with inferred defaults.
Executable file attributes, unsupported observation fleets and silent catch-fleet
aggregation are deliberately excluded. This is a first-case reader, not a
replacement for every SAM import format.

## Exact mappings for this case

The modeled period is 1963–2015, ages 1–6, with age 6 a population plus group.
Catch observations end in 2014; the missing 2015 catch is represented by `NA`.
Both surveys retain their original age ranges, fleet identity and sampling time.
Biological input grids include 2015. Source values and units are unchanged;
these files do not independently establish physical units, so no unit conversion
is assumed.

| SAM input/assumption | tinyAM mapping |
|---|---|
| `cn.dat` | `catch$obs`, with year/age/fleet metadata |
| `survey.dat` | `index$obs`, survey identity, effort and mean timing |
| `sw.dat`, `mo.dat` | `weight$obs`, `maturity$obs` |
| `nm.dat` | `weight$M_assumption`; `M_settings = list(process = "off", mu_form = NULL, mu_supplied = ~ M_assumption)` |
| `keyLogFpar` | `q_block` from global keys; `q_form = ~ 0 + q_block` (nine coefficients) |
| `keyVarObs` | `sd_block` from global keys; `sd_form = ~ 0 + sd_block` (three catch and four index coefficients) |
| `keyVarLogN = c(0,1,1,1,1,1)` | IID survival residuals with one SD; separate recruitment SD |
| recruitment code 0 | Random walk on log recruitment |
| independent LN observation errors | Existing independent lognormal likelihood |
| zero spawning-mortality timing | Beginning-year SSB |

The reusable converter also handles q fixed at 1: active `-1` q cells receive
zero rows in numeric `q_key_*` indicator columns. A formula without an intercept
then estimates only the other blocks. With every q fixed, `q_form = ~ 0` gives
zero log q. Unit tests exercise these mappings through `make_dat()`.

For this case there is no shared observation-SD key across catch and surveys.
Such sharing would be only partially supported because tinyAM estimates separate
catch and index coefficient vectors. Shared q across surveys is exact. Neither
`mu_form` nor `mono()` can reproduce equality of two fitted latent F states or
correlation between process innovations.

## Barriers to exact assessment replication

Three active process/observation assumptions in the selected case are unsupported:

1. `corFlag = 2`: AR1 correlation across the **F random-walk increments**.
   tinyAM's `rw` has independent increments; its stationary `ar1` process is a
   different model.
2. `keyVarF = c(0,1,1,1,1,1)`: two F-process SDs. tinyAM currently has one.
3. Estimated catch scaling in 1993–2005: thirteen shared scale parameters.
   Adjusting raw observations by unknown fitted scales would not replicate this
   likelihood.

Initial abundance is **partially supported**, rather than fully equivalent:
`initState = 0` and `N_settings$init = "free"` both omit an initial-state density.
However, SAM integrates all first-year log N states, whereas tinyAM estimates
`log_r0` and free `log_n0` as fixed parameters. The conditional population
equations agree, but the fitted marginal likelihood differs. This additional
inference distinction must be addressed before claiming exact replication.

`002_nscod.R` constructs only an explicitly simplified structural tinyAM model:
independent F increments, one F SD and no catch scaling. It retains the raw data,
q/observation-SD blocks, supplied M and IID N process. Its initial joint objective
and gradient are finite. It is **not optimized** and its initial objective must
not be compared with SAM's fitted marginal objective.

Two reported-quantity differences remain even where observation/population
equations agree: SAM reports arithmetic Fbar over ages 2–4, whereas tinyAM reports
an abundance-weighted F_bar; SAM uses catch mean weight for catch biomass, whereas
tinyAM yield uses stock mean weight. Catch weight and spawning timing remain
metadata for explicit calculations. They do not silently change tinyAM summaries.

## Public and current SAM reference results

`sam_reference()` extracts fitted N/F/q, stored reported SSB/recruitment/Fbar, and
stored observation predictions. It never calls an old serialized optimizer and
rejects unfitted starting objects. `SAM_*.csv` contains the public saved reference;
`reference_availability.csv` records which outputs were stored. N/F table `obs`
columns contain fitted state estimates, not measurements. No uncertainty is
invented during extraction.

The installed package was **stockassessment 0.12.0**, built from the same pinned
source commit. `003_validate_sam.R` checks inputs against `read.ices()`, constructs
data with `setup.sam.data()`, checks settings with `loadConf()`, and compares both
references with `ntable()`, `faytable()`, `qtable()`, `ssbtable()`, `rectable()` and
`fbartable()`. SAM's `qtable()` returns log q; the bridge reports its exponential.
All checks passed.

The fresh fit used SAM's default optimizer and three Newton steps:

| Diagnostic | Result |
|---|---|
| Optimizer | code 0, relative convergence (4) |
| Objective | 145.516678075642 |
| Maximum absolute gradient | 4.07 × 10⁻¹¹ |
| SD estimates finite | yes |
| Fixed-effect Hessian positive definite | yes |
| Public saved objective | 145.516678075519 |
| Largest absolute SSB difference from saved fit | 4.88 × 10⁻⁶, in source units |

Current outputs are named `SAM_current_*.csv`; diagnostics include package
version/source, elapsed fit time and comparison with the saved fit. This numerical
agreement supports the input interpretation; it is not evidence of equivalence
between SAM and the simplified tinyAM model.

Scientific decisions for a subsequent comparison are whether to retain the F
innovation correlation, two F-process SDs and catch scaling, how to handle the
initial-state integration, and which biomass and Fbar definitions to compare.
Changing the three process/observation assumptions together would obscure the
cause of differences. This branch documents those decisions without changing
tinyAM's core model.
