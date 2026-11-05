# Simulated-cod example exploration

The primary selection criterion is absolute mean annual relative SSB error below
10% over all displayed years, with successful numerical checks. This is a
single-realization teaching criterion, not an estimate of repeated-sampling bias.
Mean absolute relative error and terminal error are reported alongside the signed
mean. Abundance and recruitment errors are secondary checks. No annual years are
removed to improve the criterion.

The original 80 x 80 simulation, seed 438, packaged inputs and local raw caches
are retained. Trials use separate ignored files under
`data-raw/cache/sim_cod_exploration/`; each attempt records settings, diagnostics,
warnings/failures and numerical results. Run one fitting/simulation job at a time.

1. On the original data, replace IID N deviations and constant M with exact
   cohort survival and age/year AR1 M, retaining the supplied M baseline of 0.3.
   Test the existing separable and AR1 F specifications separately with default
   starts and fitting controls.

   Both numerical fits passed. SSB mean relative errors were -63.41% (separable
   F) and -56.59% (AR1 F); terminal errors were -83.38% and -80.79%. Thus shifting
   process variation to M alone did not solve the recovery problem.
2. Retain the closed N and AR1 M models, but replace the broad q age groups with
   RV logistic and Sentinel double-logistic curves. The previous q group 1--3
   pooled ages with very different RV length-selectivity; it could distort
   inferred recruitment and apparent survival. This is a separate q hypothesis,
   evaluated from default starts, not a retry of the earlier model.
   These curved-q fits also failed numerical checks. Mean SSB errors were
   -49.53% and -30.28%; neither qualifies regardless of those point estimates.
3. Give each survey a separate q at each modeled age (with 14+ sharing the last
   q), eliminating the arbitrary young-age grouping without forcing nonlinear
   curves into flat solutions. Compare fixed M log-mean baseline 0.3 with an
   estimated intercept; otherwise retain separable F, AR1 M and closed N.
   The fixed-mean model converged but retained -55.94% mean SSB error. Estimating
   the M intercept produced an implausibly large population (SSB roughly
   8.2e7 times truth), despite numerical checks passing, and 15 extreme parameter
   correlations. This weak-scale solution is rejected.
4. Retain age-specific q, adding independent AR1 year trajectories for each
   survey's ages 1, 2, 3, 4 and 5+ groups. The grouping follows the declared
   spatial operating-model structure, not fitted abundance or hidden numerical
   availability values. Compare supplied M baseline 0.3 and 0.4 as explicit
   alternatives, keeping separable F and N off. The latter is a design-based
   diagnostic motivated by the remaining-M summary, not an estimated M mean.
   Both converged, but mean SSB errors were -57.03% and -44.49%. The temporal
   terms also generated sparse-support/persistence warnings. Neither qualifies.
5. Diagnostic only: supply realized M, use all biological ages 1--20, and
   compare the same separable and AR1 F structures with temporal q. These fits
   deliberately use hidden truth and are ineligible as the final teaching
   example. They test whether remaining-M assumptions or another observation/
   population mismatch drives the error. Empty true population cells use Z as
   supplied M; their true catches are zero. Finite terminal-age loss remains.
   The separable-F fit failed numerical checks. The AR1-F fit passed, but its
   mean SSB error was +65.51% (mean absolute error 65.51%), with abundance +34.3%
   and recruitment +33.95%. It also had weak gillnet-curve parameters and extreme
   correlations. Supplying realized M therefore did not remove all discrepancies.

An abundance-scale ambiguity is important: for admissible positive k, scaling
N and recruitment by k, q and F by 1/k, and replacing M by Z-F/k preserves the
Baranov catches, index predictions and cohort survival. Process assumptions or
external mortality information select a point on this observational ridge;
allowing more random effects does not independently anchor abundance.

## Prospective simpler operating model

`stable_scenario.R` declares a new scenario before fitting: seed 438, 40 x 40
cells, one shared spatial field across ages and years, spatial SD 2.8 with no
unused age/year persistence, and the original recruitment, Z and growth
parameters. This is a static-availability teaching benchmark, not an unchanged
realization of the original spatial experiment. A single positive fleet-share
multiplier sets the geometric mean of the reference M baseline to 0.3 using
declared Z log-mean and integrated selectivities. No fitted estimates or SSB
errors enter this calibration. Realized M is reported separately.
The multiplier is 2.035524; its largest local base share is 0.6107 and largest
combined reference fishing fraction is 0.5683, leaving positive remaining M.
Physical access allocations use the simulated spatial population; the target
calculation uses the declared Z centre, not a target for realized remaining M.

Known growth-derived mean length is an explanatory column, not hidden population
abundance. Logistic/double-logistic relationships are evaluated on this length
scale. Both candidate assessments have closed N, AR1 F and AR1 M; compare scalar
baseline 0.3 with the declared reference age-specific M0, using all 20 biological
ages. M0 is a process-mean assumption, not realized M states. Reference F/M
coupling, finite tail, capture and sampling still differ from tinyAM's model.

The original remaining M is not centered at 0.3: its geometric mean over ages
5--10 is 0.359 (arithmetic mean 0.390), and abundance-weighted across all ages the
geometric mean is 0.455. A reference-removal preflight suggests multiplying
fishery allocations by about 1.48 would bring the ages 5--10 geometric mean
toward 0.3; this is only a design approximation, not a new simulated dataset.

Further trials and their outcomes will be recorded here before selecting the
final example. No exploratory trial changes package likelihoods or settled
assessment specifications.

### Completed static 40 x 40 trials

All 60 fleet-year targets were met, monitoring was complete, and generation
retained no warnings. Fisheries took 3,586 seconds; population and surveys each
took about five seconds. Peak resident memory was approximately 9.8 GB. Reducing
the grid alone does not remove the cost of simulating operations and samples.

| F / q means | M specification | Numerical checks | Mean SSB error |
|---|---|---|---:|
| Curves / curves | Default shared block, scalar 0.3 | Pass | -43.77% |
| Curves / curves | Default shared block, M0 | Preparation error: incompatible means | — |
| Curves / curves | Every age, scalar 0.3 | Fail | -29.90% |
| Curves / curves | Every age, M0 | Fail | -25.97% |
| Curves / free q | Four blocks, scalar 0.3 | Pass | -34.03% |
| Curves / free q | Four blocks, block-mean M0 | Pass | -34.97% |
| Free age / free q | Four blocks, block-mean M0 | Pass | -28.36% |
| Correlated RW / free q | Four blocks, block-mean M0 | Optimization error | — |
| Curves / free q | Deterministic M0 | Pass | -26.03% |
| Curves / curves | Deterministic M0 | Pass | -22.66% |

The supplied M0 age curve requires compatible M blocks; using the default
shared-actual-M block was an incorrect initial specification, retained above.
Four-block models use ages 1–3, 4–5, 6–10 and 11–20 and geometric block means.
Free q uses separate ages through 13, with older observations sharing the tail.
None qualifies. Failed-fit point estimates do not establish accuracy.

Realized geometric M was 0.366, but SSB-weighted geometric M was 0.278 versus
reference M0 of 0.264. Thus the unweighted mean discrepancy alone does not
explain the SSB error. Capture and depletion also alter the realized age profile:
at ages 6–10, geometric M exceeded the Z-adjusted reference by about 19%.
Meeting annual fleet targets does not preserve reference F at age. Separating
F, M and q remains difficult even with stable regional availability.

`teaching_scenario.R` predeclares a 20 x 20 scenario with the same seed and
recruitment settings, static common availability, Z innovation SD 0.1 and
length-at-age log-SD 0.03. These are deliberate operating-model simplifications,
not a controlled one-factor comparison with the previous simulation. The two
queued assessments share all assumptions except deterministic versus AR1 M.
