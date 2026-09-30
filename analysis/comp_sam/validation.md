# Validation of the fitted-SAM comparison workflow

Validated on 2026-09-30 with R 4.4.1, RTMB and stockassessment 0.12.0
(SAM source SHA recorded in reference_provenance.csv).

- Offline conversion/settings/audit/comparison tests: 62 assertions passed.
- Focused mixed SAM/tinyAM and missing-report dashboard test: passed.
- Full suite, run by the final installed-package check: 1,128 assertions passed,
  no failures or test warnings, one interactive-browser test skipped.
- `R CMD check --no-manual`: **Status: OK** (zero errors, warnings or notes).
- Documentation regenerated from roxygen. Repository search finds no references
  to the removed parser or old sam_tam_assumptions/sam_reference interfaces.
- Final diff against the branch's main base confirms no changes to make_dat,
  make_par, nll_fun, process densities, fit_tam or mono model mathematics.

SAM and RTMB together expose an upstream ambiguity in TMB's parallel one-step
residual helper. tinyAM now runs that requested calculation sequentially when
stockassessment is loaded. The equations are unchanged; a regression test loads
SAM before the existing residual test.

Windows text-mode file reading stops at Ctrl-Z embedded in dashboard JavaScript.
The dashboard regression test reads binary-mode HTML, verifying actual reference
notes rather than mistakenly treating the file as truncated. Missing native
reports render as unavailable and do not acquire invented SEs or optimizer data.
Both real dashboards rendered successfully. An automated visual browser preview
was unavailable because the browser policy disallows local file URLs; rendering,
HTML-content and mixed-model tests supplied the verification.

## Assessment comparison

Reference: WKCOD_combined_99, retrieved through fitfromweb without refitting.
The downloaded reference was verified equal to the earlier cached fit and was
cached with retrieval provenance and checksum. Original maturity is missing in
1963–1982. Per the explicit period choice, tinyAM is fitted on original inputs
from 1983–2022 and comparisons use that common period. The SAM reference retains
its earlier fitted history; this is an explicit source of non-equivalence.

The tinyAM fit converged (code 0, max gradient 0.000403, positive-definite Hessian,
finite reported SEs). Its 27 fixed parameters and 553 random parameters were
freely estimated after initialization from SAM estimates. No rescue changes to
settings or constraints were applied.

The native SAM reference records code 1, message "false convergence (8)", stored
maximum fixed-effect gradient 8.4e-9, positive-definite Hessian and finite reported
SEs. It is not silently classified as converged or refitted. Native and tinyAM
objectives are exported solely as optimizer diagnostics, not comparable scores.

Shared-definition SSB: correlation 0.9974, mean absolute relative difference 3.7%,
terminal difference -13.1%. Recruitment: correlation 0.9978, mean absolute relative
difference 4.5%, terminal difference -2.3%. Arithmetic Fbar: correlation 0.9701,
mean absolute relative difference 9.4%, terminal difference +51.7%. Predicted
catch biomass: correlation 0.9966 over 39 years, last available (2021) difference
+34.8%; original catch weights and SAM catch predictions are unavailable for 2022.
See exported tables for native definitions, age-specific results and uncertainty.

These results demonstrate similar historical trends but meaningful recent F/catch
differences. They do not establish interchangeability for advice or projections.
Independent versus correlated F increments, initial-state integration, the longer
SAM history, modeled maturity, and native summary definitions remain documented
assumption differences. Their causal contributions have not been isolated by
controlled sensitivities. No arbitrary agreement threshold was imposed.
