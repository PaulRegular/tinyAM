# Remaining assessment capabilities

## Pacific expansion: current coverage and findings

The database now contains **35 current detailed stock records**, including the
twelve Pacific additions documented in [PACIFIC_REVIEW.md](PACIFIC_REVIEW.md).
Completeness remains partial where source information is unresolved. Nine new
cases received optimization attempts: five passed numerical checks, three did
not converge and one was interrupted without usable reporting. Three additional
cases are blocked before fitting. Eight new dashboards are cached locally.
This is an explicit shortfall of three actually attempted translations against
the requested twelve, not twelve successful replications.

- **Original likelihoods:** the new assessments generally fit total biomass and
  multinomial, Dirichlet-multinomial or correlated age/length compositions.
  tinyAM uses separate lognormal age observations. Conversion changes their
  covariance and precision. Source zero categories remain canonical, but zeros
  are treated as missing by tinyAM's lognormal likelihood; no epsilon was added.
- **Catch-length data:** EBS cod permits a documented pooled empirical fishery
  age–length allocation covering more than 99.4% of observed length mass; the
  uncovered fraction is left unallocated. Arrowtooth lacks verified numerical
  compositions/age–length keys. Fitted N/F cannot manufacture observations.
- **Aggregate surveys:** SoG herring's observed spawning-biomass index has no
  survey age compositions. The age-specific index equation cannot represent it
  directly. This is a confirmed observation-model limitation, distinct from
  the source's density-dependent M and discrete fishery timings.
- **Sex, area and fleet structure:** female subsets or explicit sex/area pooling
  support some approximations, but do not reproduce sex-specific survival,
  spatial recruitment or normalized fleet selectivity. BC POP is one three-area
  case; all source distinctions remain canonical. Complete-area catch pooling
  leaves missing-area years missing, rather than creating partial totals.
- **Priors and identification:** unconstrained M approaches zero in the numerically
  converged hake and BC POP translations, with log-M SEs about 298 and 1925.
  These are practical identification failures despite finite numerical
  covariance. The informative source M priors/bounds are omitted. Starting
  values are not priors, and agreement alone did not trigger another fit.
- **Biological definitions:** original fixed SS3 growth/CV/allometry can recover
  BC POP weights without fitting source states. EBS cod weights are estimated
  outputs used as labelled plug-in biology. Its worksheet combined mature
  biomass differs by a factor of two from reported female SSB. BC POP's public
  reproductive-output configuration and reported biomass need further source
  confirmation; AI cod also has unresolved single-sex SSB scaling. These are
  definition cautions, not proof of an error in tinyAM mathematics.
- **Reference uncertainty:** native MCMC medians/means and 90% or 95% credible
  limits remain labelled. Marginal age/sex medians and component modes do not
  determine posterior means/intervals of a pooled state. SoG component MPD N is
  not ensemble posterior N; no such series or uncertainty was invented.
- **Source gaps:** AI cod's public September native run differs from its final
  November fit. Arrowtooth lacks matching likelihood arrays/keys and timing.
  SoG's public DDM ensemble ends in 2023, distinct from the interim 2025 DIM;
  the two are not combined. Unverified activation/sharing/point conventions
  remain unresolved, rather than being described as absent features.

Full-library review flags use >20% mean absolute or terminal differences and
correlation <0.90. Fixed M is checked as an input, unavailable comparisons are
listed separately, and failed numerical/uncertainty checks are not hidden by
good individual trends. No new core model feature was implemented for this
expansion. The earlier prevalence inventory below is historical, not a recount
over the enlarged library.

## Historical 23-assessment capability inventory

**Process-expansion follow-up:** F now supports AR1 age-correlated RW increments
(`cor_rw`), and N/F/M support ordinary age-based log-SD sharing formulas.
Stock-wide spawning-time SSB is also available. The prevalence counts below
remain those of the source review; these three rows now describe supported
options requiring source-definition checks, rather than absent features. See
[the process study](../formula_expansion/PROCESS_EXPANSION.md) for recovery and application evidence.

The original review covered **23 then-current detailed assessment records: 21 stock recipes
and two blocked translations**. Counts below are confirmed minimums, once per
assessment, not fleet or database-row counts. The explicit membership sets are
in the [historical review script](https://github.com/PaulRegular/tinyAM/blob/782034b/analysis/comp_assessments/scripts/translation/review_translations.R);
source reviews link the reports,
frameworks and pinned code used to verify them. Unrecorded or unresolved
specifications are not counted as absence. Norwegian herring's exact process
configuration, Icelandic haddock's native sharing settings, and some published
survey-key assignments remain unresolved. The two blocked records lack enough
source inputs to fit; their limited known specifications remain in the inventory.

| Feature family | Confirmed / 23 | Examples and source evidence | Current workaround / support | Potential benefit | Complexity |
|---|---:|---|---|---|---|
| Correlated F random-walk innovations across ages | 10 | [North Sea haddock](source_reviews/ices_haddock_north_sea.md), [blue whiting](source_reviews/ices_bluewhiting_northeast_atlantic_2026.md) | `cor_rw` now supports AR1 age-correlated increments. Verify source correlation and first-state conventions; arbitrary covariance remains unsupported. | Retain common annual mortality changes without replacing an RW with a stationary process. | Implemented for F |
| Correlated observation errors, including full biomass covariance | 9 | [NE Arctic cod](source_reviews/ices_cod_northeast_arctic.md), [EBS pollock](source_reviews/afsc_pollock_ebs.md) | Independent lognormal observations; mean/SD random effects are not equivalent. | Represent joint information and uncertainty more faithfully. | Moderate–high |
| Multiple age-sharing groups for latent F/N process SDs | 9 | [Plaice](source_reviews/ices_plaice_north_sea.md), [Northern Gulf cod](source_reviews/dfo_cod_3pn4rs.md) | Ordinary age-based `sd_form` is now supported for N/F/M, with M designs constant within fitted state blocks. Recruitment remains separate. | Different process variability for young ages or plus groups without extra observation noise. | Implemented |
| Totals plus age/length composition likelihoods | 9 | [GOA cod](source_reviews/afsc_cod_goa.md), [Northern Shelf cod](source_reviews/ices_cod_north_sea.md) | Convert recoverable totals/compositions to age observations, or omit unsupported components. Multinomial, Dirichlet-multinomial and logistic-normal source likelihoods differ. | Use original observation units and sample-size information. | High |
| Fleet, substock and normalized selectivity structures | 8 | [GOA pollock](source_reviews/afsc_pollock_goa.md), [summer flounder](source_reviews/nefsc_summer_flounder.md) | Aggregate removals and flexible F-at-age. Rising mono/logistic F means are supported with IID/AR1 residuals; they do not reproduce fleet-specific, normalized or double-logistic selectivity processes. | Preserve fishery distinctions, selectivity normalization and associated observations. | High |
| Survival to spawning time in SSB | 6 | [Baltic sprat](source_reviews/ices_sprat_baltic.md), [EBS pollock](source_reviews/afsc_pollock_ebs.md) | `data$maturity$spawn_time` supports age/year timing, shared by F/M. Separate F/M fractions and seasonal biology remain unsupported. | Match parent SSB and reported biomass timing. | Implemented for supplied age/year fractions |
| Priors on starting mortality, q or selectivity | 3 | [Southern Gulf cod](source_reviews/dfo_cod_4t4vn.md), [GOA pollock](source_reviews/afsc_pollock_goa.md) | Starting values or q logit restriction are not priors. Fixed Gaussian effect SDs are supported, including a standalone M mean RW; the source starting-level priors remain different. | Reproduce external information and stabilize genuinely weak boundaries. | Moderate |
| Recruitment curves beyond BH/Ricker | 2 | [Western Baltic herring](source_reviews/ices_herring_western_baltic_2026.md), [Southern Gulf cod](source_reviews/dfo_cod_4t4vn.md) | RW or stationary recruitment; hockey-stick and SSB times an AR1 recruitment rate remain distinct. | Match those recruitment equations without substituting another curve. | Moderate |

These are feature families, not eight independent package defects. Some have
partial mathematical support, and their implementation would require decisions
about identifiability, missing data and reporting. Other documented limitations
include quarterly Norway-pout dynamics, fitted biological GMRF observations,
aggregate-only indices, density-dependent q powers, ageing-error likelihoods
and spawning-component indices. Their presence does not justify inventing
missing observation allocations. Broader source gaps are separate from package
capabilities: Northern Gulf cod lacks annual stock weights and fitted maturity;
Georges Bank haddock lacks recoverable current model-grid observations/native
files. New formulas do not solve those blockers.

## What the expanded formulas already cover

IID/RW/AR1 recruitment and BH/Ricker median curves, fixed recruitment SDs,
Gaussian catchability and mean effects, monotone/logistic survey curves,
categorical/numeric grouping and observation-SD effects are available. Many
SAM baselines already use their supported recruitment and sharing structures.
GOA cod fixes BH steepness at one: a free BH curve is not a faithful upgrade.
The Northern cod curve describes age-zero recruitment, while the present
translation starts at age two. Its BH candidates are survivor approximations,
not exact source relationships. See the review table for actual attempts and
retention decisions; fitting an available feature is not itself an improvement.

## Priorities from the source review

Items 1 and 3 below are now implemented within the limits described above.
The remaining items are proposals, not authorized additions.

1. **Age-correlated F RW innovations**, followed by simple **latent SD sharing
   groups**: common needs that fit the existing state-space architecture.
2. **Observation correlation**: begin with age AR1 before considering arbitrary
   supplied covariance matrices. Keep this separate from random mean effects.
3. **Spawning-time SSB**: a comparatively contained biological option, with
   consistent recruitment-parent and reporting definitions.
4. **Narrowly specified priors**, if justified by source use and uncertainty
   testing. Known M increment SDs can already be supplied through standalone
   mean effects; a dedicated residual-SD interface would mainly ease that
   specification. Avoid a general prior language initially.

Composition likelihoods and explicit fleets are useful but would substantially
expand tinyAM. Additional recruitment curves have lower confirmed prevalence
in this inventory. None of these recommendations assumes that a new feature
would necessarily improve agreement with an official assessment.
