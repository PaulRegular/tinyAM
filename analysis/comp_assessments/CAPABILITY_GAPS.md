# Remaining assessment capabilities

This review covers **23 current detailed assessment records: 21 stock recipes
and two blocked translations**. Counts below are confirmed minimums, once per
assessment, not fleet or database-row counts. The explicit membership sets are
in `scripts/translation/review_translations.R`; source reviews link the reports,
frameworks and pinned code used to verify them. Unrecorded or unresolved
specifications are not counted as absence. Norwegian herring's exact process
configuration, Icelandic haddock's native sharing settings, and some published
survey-key assignments remain unresolved. The two blocked records lack enough
source inputs to fit; their limited known specifications remain in the inventory.

| Feature family | Confirmed / 23 | Examples and source evidence | Current workaround / support | Potential benefit | Complexity |
|---|---:|---|---|---|---|
| Correlated F random-walk innovations across ages | 10 | [North Sea haddock](source_reviews/ices_haddock_north_sea.md), [blue whiting](source_reviews/ices_bluewhiting_northeast_atlantic_2026.md) | Independent age RWs or stationary 2D AR1. Shared mean effects do not give the source innovation covariance. | Retain common annual mortality changes without replacing an RW with a stationary process. | Moderate |
| Correlated observation errors, including full biomass covariance | 9 | [NE Arctic cod](source_reviews/ices_cod_northeast_arctic.md), [EBS pollock](source_reviews/afsc_pollock_ebs.md) | Independent lognormal observations; mean/SD random effects are not equivalent. | Represent joint information and uncertainty more faithfully. | Moderate–high |
| Multiple age-sharing groups for latent F/N innovation SDs | 8 | [Plaice](source_reviews/ices_plaice_north_sea.md), [North Sea haddock](source_reviews/ices_haddock_north_sea.md) | One residual-process SD per F/N/M. Recruitment already has its own SD; ordinary recruitment-versus-survival sharing is **not** counted as a gap. Observation-SD formulas are supported. | Different process variability for young ages or plus groups without extra observation noise. | Moderate |
| Totals plus age/length composition likelihoods | 9 | [GOA cod](source_reviews/afsc_cod_goa.md), [Northern Shelf cod](source_reviews/ices_cod_north_sea.md) | Convert recoverable totals/compositions to age observations, or omit unsupported components. Multinomial, Dirichlet-multinomial and logistic-normal source likelihoods differ. | Use original observation units and sample-size information. | High |
| Fleet, substock and normalized selectivity structures | 8 | [GOA pollock](source_reviews/afsc_pollock_goa.md), [summer flounder](source_reviews/nefsc_summer_flounder.md) | Aggregate removals and flexible F-at-age. Logistic q is supported, but is not logistic fishery F, double-logistic selectivity or fleet-specific F. | Preserve fishery distinctions, selectivity normalization and associated observations. | High |
| Survival to spawning time in SSB | 6 | [Baltic sprat](source_reviews/ices_sprat_baltic.md), [EBS pollock](source_reviews/afsc_pollock_ebs.md) | Start-year mature biomass; recompute common-definition references when N and biology permit, otherwise label approximate. | Match parent SSB and reported biomass timing. | Low–moderate |
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

## Suggested priorities

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
would necessarily improve agreement with an accepted assessment.
