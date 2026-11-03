# Eastern Bering Sea Pacific cod: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/afsc_cod_ebs.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/afsc_cod_ebs.md).

## Expanded-library review

### Translation

The original likelihood uses catch-at-length. A tinyAM approximation needs an
explicit allocation from length to age. The native inactive fishery age samples
provide a potential pooled empirical key, without using fitted N to adjust an
observation. Its covered length bins contain 99.46–99.99% of annual fishery
number proportions in 2000–2024. Uncovered lengths must remain unallocated,
not be assigned invented ages or silently redistributed. Any such translation
must disclose pooling across years, omitted mass, age-error differences,
weight–length bin approximation and plug-in fitted population biology.

Original maturity-at-age and empirical stock weights were not supplied to SS.
Do not label a weight-derived maturity proxy as the official maturity curve.
Canonical completeness remains partial; source preservation and numerical
model readiness are separate claims.

The retained recipe fits 2000–2023, ages 1–12+, using the pooled fishery key,
original VAST numbers and active survey ages. Catch total numbers use the
fixed weight–length equation at five-cm bin midpoints. Population weights are
plug-in official outputs; the maturity proxy applies the original length
ogive to weight-equivalent length. Pool ages 12–20+ with official N weights
for biological reporting. This is not the source's length-integrated maturity.
SSB comparisons therefore use a labelled common proxy; native female SSB and
its reported SD remain separate. Neither fitted q nor fitted N changes an
observation. Zero source proportions remain zero; tinyAM's lognormal
observations treat these zero entries as missing.

#### Proposed paired-model candidate (not yet fitted)

The updated translation continues to use deterministic older-age survival (`N process = "off"`),
free first-year older-age abundance (`init = "free"`), and IID log recruitment around
an estimated constant mean. **The recruitment log SD is now estimated** (the
`sd` setting was removed), whereas the accepted assessment fixes it at 0.6908.
This deliberately adds flexibility and is not an exact match to the source
recruitment penalty. Natural mortality remains fixed at 0.3866; January spawning
timing and the existing proxy biology are unchanged.

Relative to the earlier paired candidate, the new change is the estimated
recruitment SD. The other two non-operational process approximations remain:
F-at-age follows **correlated temporal random walks** (`F_settings = list(process = "cor_rw", mu_form = NULL)`), in which
annual log-F innovations are correlated among ages. This is a simple proxy for
collective changes in fishing intensity and selectivity, not an exact representation
of length-based fishing or the official selectivity penalties. Second, the survey
index uses **rising logistic age-specific catchability with an annual common log-q
random walk** (`q_form = ~ logistic(age) + rw(year)`). The logistic component
approximates the ascending survey selectivity curve. The RW permits the whole survey
q level to change among years; it does **not** reproduce the official annual changes
in length-selectivity width or the official scalar-q specification.

The translation now ends in a single `if (do_fit)` block calling `tinyAM::fit_tam()`.
The previous `fit_default_first()` helper, fitted-official-N/F retry starts and
initialization-attempt cache have been removed. Original catch/survey reconstruction,
source data, biological proxy and comparison definitions are retained.
The estimated recruitment SD should be reported from the new fit rather than
assumed to equal 0.6908. Its value and the new fit diagnostics have not been
provided, so no numerical conclusions are attributed to this candidate.

#### Historical baseline results (not results for this candidate)

The **superseded** independent age-specific F random walk / unrestricted age-q
configuration converged with optimizer code 0, objective 372.2591, maximum gradient
0.0013577, positive-definite Hessian and finite covariance (approximately 6.45 s).
Its source comparisons reported 42.1% mean absolute N-at-age difference, -43.7%
terminal N difference, within-age N correlation 0.681; 70.1% mean absolute F-at-age
difference, +107.0% terminal F difference, within-age F correlation -0.200; and
21.7% mean absolute proxy-SSB difference. None of these figures should be quoted
as performance of the newly specified `cor_rw` F / logistic-plus-year-RW q candidate.

Substantial limitations remain: fishery catch numbers are reconstructed with a
pooled empirical age–length key; the source fits catch-at-length and age/length
compositions; the source selectivity varies in length space; age error and survey
precision are simplified; and age-0 recruitment/SSB definitions are not fully
matched. Official fixed M is not an independent validation of M estimation.
