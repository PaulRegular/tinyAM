# Assessment capabilities and remaining gaps

## Historical 23-assessment capability inventory

**Process-expansion follow-up:** F now supports AR1 age-correlated RW increments
(`cor_rw`), and N/F/M support ordinary age-based log-SD sharing formulas.
Age/year-specific spawning-time SSB is also available. The prevalence counts below
remain those of the source review; these three rows now describe supported
options requiring source-definition checks, rather than absent features. See
[the validation studies](../formula_expansion/README.md) for recovery and application evidence.

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

The [final standalone West of Scotland cod assessment](source_reviews/ices_cod_west_scotland_2022.md)
also confirms an unsupported observation-mean feature: one fitted catch-scaling
factor per year in 1995–2006. Changing observation SDs is not equivalent to that
mean correction, and fitted source factors must not be applied to original
observations. Its translation retains the observations and omits the correction.
An actual fit with latent initial abundance passes numerical checks only after
a source-start retry, but gives extreme abundance, near-zero F and extensive
parameter correlation. This is a flagged identifiability/approximation outcome,
not evidence of successful replication or proof that any one omitted feature
caused it. The historical prevalence counts above are unchanged.

[Georges Bank yellowtail flounder](source_reviews/ck_dfo_nefsc_ferruginea_gb.md)
has verified 2025 production summary outputs, but the recovered numerical
2024 archive belongs to the superseded constant-M candidate. The selected
age-based-M run lacks corresponding numerical input/output tables in that
archive, and the 2025 update revised maturity. This is a source-data blocker:
candidate biology or age states must not fill the production record. The
partial summary import is not a completed tinyAM translation.

[Northeast Arctic beaked redfish](source_reviews/ck_ices_afwg_mentella_deep_1_2.md)
now has published N, reconstructed rounded F and source biological vectors in
the database. Its stock-specific SCAA uses two fleet selectivity processes,
explicit older cohorts and a separate deep-pelagic composition likelihood.
Missing young-age biology, numerical pelagic compositions and exact overlapping
catch-series assembly still prevent translation. Those source gaps must not be
treated as model capabilities or filled by borrowing another run's inputs.

For [4VWX herring](source_reviews/ck_dfo_mar_harengus_4vwx.md), the recovered
native example inputs end in 2018, while the production reference operating
models were conditioned through 2020. The 2025 component indicator update is
recorded separately. Its acoustic spawning-biomass index is not a stock-wide
fitted SSB series; a single tinyAM comparison would first need the production
inputs, reference outputs and a justified component/ensemble definition.

[Atlantic striped bass](source_reviews/ck_nefsc_garmiii_saxatilis_gom_chatt.md)
has a separate 2018 historical SCA translation using published inputs and
outputs, but its 2024 production observation arrays and annual weights remain
unrecovered. The 2018 approximation converges after source-based starts, but
F agreement is poor; it does not stand in for the current 2024 update. The
source uses totals plus compositions, fleet selectivity blocks and constrained
recruitment deviations. Missing 2024 input arrays are a source-data blocker;
published fitted states and unrelated simulated data cannot replace them.
Female SSB requires the separate female fraction and maturity inputs.

[Eastern Georges Bank cod](source_reviews/ck_nefsc_trac_cert_morhua_gb.md)
uses a supplied Lorenzen M surface and decoupled older-cohort process error;
its inferred apparent mortality is not a separately estimated M output.
Recruitment is a random walk with correlated innovations, which is distinct
from both independent RW increments and AR1 levels. The selected NTS catch-age
data are not the CATCH age estimates printed in the input-review table.
Missing selected observations and complete biology currently block fitting;
separate fleets and aggregate/composition likelihoods would still require
documented approximations after source recovery. Historical review counts
above are unchanged.

[Southern New England/Mid-Atlantic yellowtail flounder](source_reviews/ck_nefsc_garmiii_ferruginea_sne_ma.md)
now has the published catch-age and three trawl-survey source tables, plus a
separate 2025 official summary. Its selected framework uses age-specific M,
not the original constant-M candidate; the update removes fleet-selectivity
random effects and changes biology and GSI inputs. Exact production biology,
GSI configuration and fitted age surfaces remain missing, so no translation
is claimed. The selected 2024 age-specific M vector is now recorded from its
supporting model slide. The 2025 update does not restate whether that vector
was carried over unchanged. Aggregate/composition likelihoods and correlated
older-cohort N would still require documented approximations. Historical
counts are unchanged.

The [2026 Sandeel Area 1r assessment](source_reviews/ices_sandeel_area_1r.md)
adds a verified seasonal example: July recruitment, half-year mortality and
biology, effort-proportional F, and a survey prediction averaged over an
interval. The retained annual ages-1–4+ approximation sums seasonal mortality
hazards and catch counts, uses January biology and a survey midpoint, and
substitutes IID recruitment/F around documented mean structures. Annual F
means can share ages and periods, but shared mean coefficients do not force
latent ages 3–4+ to be identical. Its default-start fit converges with usable
uncertainty; recruitment and abundance remain flagged for disagreement.
Exact catch flags and two seasonal effort rows remain source gaps. This
limited translation does not establish support for native seasonal dynamics,
integrated survey predictions or the downweighted hockey-stick likelihood.
The historical inventory counts above are unchanged.

The corresponding [Area 2r translation](source_reviews/ices_sandeel_area_2r.md)
also converges from defaults with usable uncertainty, retaining four F mean
periods and one age-1 Dredge series. Its terminal SSB, recruitment and abundance
differences remain above review thresholds. These two numerical successes
support documented annual approximations; they do not establish support for
the native seasonal likelihood. Source code also explains the different
recruitment averaging windows in summary and forecast outputs; those future
assumptions are not historical observations.
The public configurations also use SD argument names absent from their
declared pinned builder. That is an unresolved production-software provenance
gap, distinct from a tinyAM capability gap or a demonstrated error in the
official numerical outputs.

The [Area 3r source record](source_reviews/ices_sandeel_area_3r.md) adds an
exact zero-effort/catch/F terminal year (2025). Finite log-F states cannot
represent this structural closure; zero lognormal catch observations also
do not enforce it. Source zeros must remain unchanged, and terminal F
percent differences are unavailable rather than passes. An annual model
with correlated F increments failed from both default and aligned source
starts as the estimated age correlation approached its boundary. This is a
documented numerical failure of that approximation, not evidence that the
official seasonal assessment failed or that all annual alternatives fail.
An independent-increment trial passes numerical checks from defaults but
drifts toward essentially zero F, enormous abundance and unbounded SSB/q
intervals. It is retained as exploratory, with all comparable core metrics
flagged. Finite curvature alone is insufficient evidence of an identifiable
abundance scale. This case is an unsuccessful annual replication; no
pseudo-effort, altered catches, terminal-year deletion or agreement-driven
starting-value retry is used to conceal that outcome.

The [detailed Area 4 case](source_reviews/ices_sandeel_area_4.md) recovers the
same seasonal framework, with two Dredge eras and estimated survey log-SDs
subject to a 0.3 lower bound. The annual recipe retains juvenile indices
without changing their values: constant supplied age-0 M and zero age-0 F
allow their prediction to be expressed against next January N1 with a
redefined free q. Original observation dates remain in provenance; reference
recruitment stays in the same calendar year. This conditional mapping is not
a general solution for estimated/time-varying juvenile mortality.
Both default and aligned-start fits fail with non-finite gradients, so no
tinyAM comparison statistics are available. The precise numerical cause is
unresolved; source differences are not proof of a particular cause. The
Official-only dashboard is source context, not a completed replication.

Area 4 also verifies that native `Catchtot` includes only predicted catch
cells with positive observed catch. The source's conditional yield can be
zero despite some positive F. It must not be treated as unrestricted
predicted yield or observed total landings. Its second-half effort table has
reversed Area 4/5r headers, reconciled against both independent annual totals
and explicitly labelled as a source reconstruction. These are source-data
interpretation issues, distinct from model capability or convergence.

## What the expanded formulas already cover

The [Celtic Sea plaice review](source_reviews/ices_plaice_celtic_sea_7fg.md)
illustrates a different limit: the current production assessment is a
category-3 biomass/length advice rule, despite its historical age-structured
catalogue label. Verified survey-at-age tables are ancillary; current
catch-at-age and biology are not recovered, and the official rule does not
report fitted N/F/SSB states. This is a source and assessment-definition gap,
not evidence that an additional tinyAM process is required. Its biomass
index and length-based fishing-pressure ratio are kept distinct from SSB
and instantaneous F. No fit or agreement statistics are claimed.

[Northeast Arctic golden redfish](source_reviews/ices_golden_redfish_northeast_arctic.md)
likewise moved away from an age-length assessment in 2026, to SPiCT using
total landings and a compiled adult-biomass index. Its B/BMSY and F/FMSY
states are recovered with reported bounds, but current age-based inputs
and numerical survey data are not. The benchmark also documents an
alternative high-biomass solution under jittering. This source and
assessment-definition limitation must not be represented as successful
tinyAM agreement or addressed by substituting older Gadget observations.

[Western English Channel plaice](source_reviews/ices_plaice_western_channel_7e.md)
uses an MSE-tuned chr rule, not annual population estimates. Its harvest
indicator uses landings plus half of discards under a 50% discard-survival
assumption; reported total catch includes surviving discards. This distinction
must precede any future translation of removals. Stock-level catches also
include migration into 7.d. Current catch-at-age/biology are not recovered,
so aggregate indices and MSE operating models do not establish a tinyAM
comparison. The 2026 monitoring update issues no new advice and is kept
separate from the 2025 production inputs.

[Gulf of Cadiz anchovy](source_reviews/ices_anchovy_gulf_cadiz.md) has an
age-based SS3 assessment despite category-3 advice. It settles age-0 recruits
on 1 July and uses quarterly fishing fleets, high age-specific M, biomass
indices and age compositions. These seasonal population and observation
definitions require an explicit annual approximation before tinyAM can be
compared. The immediate blocker is source data: the public 2025 TAF tree
omits its empirical stock-weight and quarterly composition inputs. A public
1989–2024 MSE reference differs from the current official estimates and is
not a replacement. Published survey ages are expanded abundances, not sample
counts; they must not be added as an independent likelihood alongside the
same biomass/composition observations. No fit or numerical agreement is
claimed for this source record.

[US Georges Bank cod](source_reviews/nefsc_cod_georges_bank.md) illustrates
why model assumptions must be checked against the selected update: the
2024 run replaces historical IID survival errors with decoupled 2DAR1
errors, estimates initial N freely, and turns off lognormal adjustments.
It also fixes fall-survey selectivity to one at ages 1–3, an unusual pattern
explicitly flagged by the official review. Current original age arrays and
biology remain unrecovered, so no tinyAM agreement is claimed from its
published aggregate catch/status table or fixed M alone.

The [White Hake Run_005 audit](source_reviews/ck_nefsc_garmiii_tenuis_gom_gb.md)
likewise fixes BH steepness at one, but gives the recruitment penalty zero
weight. A constant IID tinyAM mean preserves the flat expectation, not the
source's unpenalized sum-to-zero fixed-effect recruitment vector. The source
also penalizes initial abundance around first-year equilibrium survival and
uses deterministic older-cohort survival. These distinctions require explicit
translation notes; a supplied CV alone does not establish an active process.

IID/RW/AR1 recruitment and BH/Ricker median curves, fixed recruitment SDs,
Gaussian catchability and mean effects, monotone/logistic survey curves,
categorical/numeric grouping and observation-SD effects are available. Many
SAM baselines already use their supported recruitment and sharing structures.
GOA cod fixes BH steepness at one: a free BH curve is not a faithful upgrade.
The Northern cod curve describes age-zero recruitment, while the present
translation starts at age two. A BH curve at age two would be a survivor approximation,
not the source relationship. Fitting an available feature is not itself an improvement.

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
