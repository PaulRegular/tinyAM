# Validation

Validated on 2026-09-30 with R 4.4.1 and stockassessment 0.12.0.

- The simplified WKCOD workflow fitted 1983–2022 and rendered the dashboard.
  Ordinary tinyAM starting values reached the same objective as the previous
  SAM-based starts (442.7456, difference below 1e-5). Optimizer code: 0;
  maximum absolute gradient: 0.000857; positive-definite Hessian.
- Verified all 960 annual comparisons and 24 summary rows against their source
  estimates. M differences are zero. Additional checks cover cancellation of
  signed differences, terminal-year selection and undefined zero-reference ratios.
- The focused mixed SAM/tinyAM dashboard test passed, including metric-specific
  confidence-interval and SE headers and absence of generic `.x`/`.y` SE headers.
- Full installed-package suite: 1,142 assertions passed, no failures or test
  warnings, one interactive-browser test skipped.
- `R CMD check --no-manual`: **Status: OK**, no errors, warnings or notes.

No package model equations, settings mappings or public interfaces changed.
SSB retains each model's native definition; arithmetic Fbar uses ages 2–4 for
both. The cached SAM optimizer code 1 remains visible in diagnostics. These
comparisons do not isolate the causes of disagreement or certify equivalence.
