# Historical cheatsheet experiments

This is an archival snapshot, not a supported workflow. It is intentionally
removed in the following cleanup commit; retrieve it from Git when needed.

The three R scripts record trials with the original cod inputs, catch-scaled
initial abundance, explicit recruitment starts, and alternative F processes.
Their absolute paths and calls are preserved as run; use the package revision
from their provenance when available rather than assuming current compatibility.
`attempts.csv` records stored outcomes and warnings without rerunning fits.
Unknown historical revisions remain missing. Bulk fitted objects are not tracked.

Earlier five-model examples (including BH, F curves and monotone q) are already
recoverable from cheatsheet history before e5e35fa. Their failed fits were not
replaced by undocumented tuning. N_iid remains numerically unsuccessful in the
current example and is unranked.

Catch-conditioned, survey-only scoring used 211 observations in 18 of the 20
requested forecast years; two years had no scored survey observations. The
shared M block scored 0.9470344695. Adding breaks c(2,4,6,8,14) yielded
0.7654028667 versus Base 0.7738511352. All 20 age-blocked M forecast fits passed
numerical checks; boundary/correlation advisories and an optimizer warning
were retained. This is an example comparison, not evidence of universal ranking.

AR1 M forecasts were already conditional on the terminal state. The earlier
claim that they reset immediately to their mean was incorrect; f0fa93a added
regression tests without changing those mathematics. The status-quo scores in
the CSVs distinguish the old examples from the catch-conditioned study.

Disposable inspection, plotting, cache-migration and test-launch scripts are
not retained here. Catch-solver, zero-catch, AR1 M, parallel hindcast and
recruitment-boundary checks were promoted to package regression tests in
01c3d16, f0fa93a, 909cfbe and d91852e.
