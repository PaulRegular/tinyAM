# Compare SAM with tinyAM

Load/install this branch and run `source("analysis/comp_sam/001_compare.R")`
from the repository root. The script loads the cached `WKCOD_combined_99` SAM
fit (or retrieves it with `fitfromweb`), translates its inputs/settings, fits
tinyAM, summarizes percent differences, and renders the dashboard.

The comparison uses **1983–2022**, because original maturity is missing before
1983. Original inputs are retained; neither observations nor biology are filled
using fitted SAM values. tinyAM uses ordinary starting values and the existing
translated settings. No model equations or assumptions were changed by cleanup.

The outputs in `results/` are:

- `percent_differences.csv`: annual SSB, recruitment, N/F/M at age, and arithmetic
  Fbar estimates and `100 * (tinyAM / SAM - 1)`. Positive means tinyAM is higher.
  A zero SAM estimate gives an undefined percent difference (`NA`).
- `summary.csv`: mean signed difference, mean absolute difference (avoiding
  cancellation), and terminal-year difference, by metric and age. All differences
  are percentages; missing annual differences are excluded from the averages.
- `audit.csv`: exact mappings, approximations, unsupported features and unknowns.
- `diagnostics.csv`: optimizer code/message, gradient and Hessian status.
- `SAM_tinyAM_dashboard.html`: detailed states, q, predictions, available
  uncertainty, original input tables and model calls; tables can be downloaded.

Fbar is arithmetic over the same ages (2–4 here) for both models. SSB retains each
assessment's native definition: SAM estimates maturity, while tinyAM uses the
original maturity inputs. This difference must be stated when reporting SSB.
M is supplied identically here, so differences are zero to numerical precision.
No extra transformed summaries or confidence intervals are manufactured.

The reference RDS and provenance/checksum files retain the source. The tinyAM
fit is saved for later use. RDS files and the generated HTML are ignored by Git;
the scripts and four review CSVs are tracked.

The SAM reference retains its pre-1983 fitted history, uses correlated F increments,
and integrates initial abundance. tinyAM uses independent F increments and free
initial abundance. Legacy `initState`/`logNMeanAssumption` remain unresolved in
the audit. SAM's cached optimizer code 1 ("false convergence") is retained despite
its tiny gradient and positive-definite Hessian. The tinyAM fit converges; these
results establish similar historical trends, not interchangeability for advice.
