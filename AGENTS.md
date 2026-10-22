
# tinyAM development conventions

tinyAM is an R package built with RTMB. Keep changes small,
transparent, and consistent with existing package structure.

## General package conventions

- Treat this as a standard R package, not a collection of scripts.
- Put exported/user-facing functions in `R/`.
- Put internal helpers in `R/` with a leading `.` where appropriate.
- Do not source R files manually from package code.
- Do not use `library()` or `require()` inside package functions.
- Use namespace-qualified calls for non-imported functions where practical,
  e.g. `stats::model.matrix()`, `RTMB::dnorm()`.
- Add package dependencies only when necessary.
- Update `DESCRIPTION` if a new dependency is genuinely required.
- Do not hand-edit `NAMESPACE`; generate it from roxygen.
- Do not hand-edit generated `.Rd` files unless explicitly requested;
  update roxygen comments and regenerate documentation.
- Preserve compatibility with R >= 4.1.0 unless explicitly changing the
  package requirement.

## Documentation

- Every exported function must have roxygen2 documentation.
- Use roxygen markdown syntax.
- Document new arguments in the function where they are introduced.
- Include mathematical definitions for model parameters and processes when useful.
- Keep terminology consistent across `prepare_tam()`, `make_par()`, `nll_fun()`,
  simulation, tidying, and documentation.
- Regenerate documentation after modifying exported interfaces.
- Write user-facing documentation primarily for fisheries biologists and ecologists
  with basic statistical training.
- Use plain biological and modelling language first.
- Explain what an option does, what assumption it represents, and when a user might
  choose it before describing implementation details.
- Keep equations and statistical precision where useful, but place more technical
  material such as parameterizations, likelihood details, transformations,
  covariance structures, and implementation terminology in `Details` or clearly
  identified technical sections.
- Define statistical jargon when it cannot reasonably be avoided.
- Write examples to emphasize specifying, fitting, and interpreting assessment
  models rather than package internals.
  
## Errors and messages

- Use `cli::cli_abort()` for user-facing errors.
- Use `cli::cli_warn()` for warnings.
- Use `cli::cli_inform()` for informational messages.
- Prefer concise, actionable messages.
- Do not use base `stop()`, `warning()`, or `message()` for new user-facing
  conditions unless there is a specific technical reason.

## Testing

- Add or update `testthat` tests for every behavioural change.
- Test both intended behaviour and important invalid inputs.
- Prefer small deterministic unit tests for model mathematics.
- Add regression tests when fixing bugs.
- Run targeted tests first, then the full test suite.
- Run `R CMD check` for substantial changes when practical.
- Do not weaken existing tests merely to make new code pass.

## RTMB / model conventions

- Preserve the latent-state convention:
  `log_r`, `log_n`, `log_f`, and `log_m` are absolute latent log states.
- Process deviations should be calculated explicitly as `eta_*`.
- Expected/mean log surfaces should use `log_mu_*`.
- Estimation and simulation must use the same mathematical model.
- A new process option must update:
  `prepare_tam()`, `make_par()`, `nll_fun()`, simulation, tidying,
  documentation, and tests as applicable.
- Avoid hidden special cases that cause model semantics to differ between
  fitting and simulation.
- Do not add unused parameters or parameters with flat likelihood directions.

## API design

- Prefer compact, explicit settings over adding many top-level arguments.
- Follow existing structures such as:
  `N_settings`, `F_settings`, `M_settings`, `catch_settings`,
  and `index_settings`.
- Defaults should favour stable, identifiable models.
- Scientific assumptions should be visible in the interface rather than hidden.
- Preserve backward compatibility unless the task explicitly authorizes an
  API change.
- If retiring an experimental argument, update all package code, tests,
  examples, analyses, and documentation rather than leaving compatibility
  hacks unless explicitly requested.

## Code style

- Match the surrounding code style.
- Prefer clear vectorized/base-R code over clever abstractions.
- Keep helpers small and single-purpose.
- Dashboard chunks should primarily call exported plotting functions, so users
  can reproduce the displayed plots outside the dashboard. Keep small
  dashboard-specific overlays in the template, with styling visible there.
- Avoid deep nesting.
- Use descriptive object names consistent with existing model notation.
- Do not introduce a new dependency just for code style or convenience.

## Repository-wide changes

When changing parameter dimensions, names, settings, or semantics, search the
entire repository for downstream assumptions, including:

- `R/`
- `tests/`
- `analysis/`
- `inst/examples/`
- `inst/rmd/`
- `README.md`
- roxygen/man documentation

Warm starts, retrospectives, simulation, plotting, tidying, and saved analysis
objects are especially likely to depend on parameter structure.

## Analysis code

- Package functionality belongs in `R/`; analysis-specific code belongs in
  `analysis/`.
- Do not move one-off analysis logic into package code unless it is genuinely
  reusable.
- Keep competing models identical except for the assumption being tested.
- Do not silently change multiple scientific assumptions in the same model
  comparison.

## Generated files and experimental work

- Prefer stable output names and existing entrypoints to new scripts for each run.
- Keep generated fits, dashboards, figures and build intermediates in ignored
  output/cache directories. Track only designated deliverables, such as the
  canonical cheatsheet PDF; generate other assets from their existing sources.
- Use temporary directories for disposable probes and render intermediates;
  clean them on completion. Retain one useful failure log rather than numbered
  copies of every build.
- Record substantive model trials, failed approaches and their outcomes in
  commits before removing their development files. Git history is the archive;
  do not leave retired experiments in the current checkout solely for history.
- Keep current-run warnings, failures and provenance with cached results.
  Preserve useful history before pruning obsolete caches.
- Promote a scratch check to a permanent test only when it protects meaningful
  package behavior that existing tests do not cover.
- Retain the current scientific workflow, canonical imports/validators, final
  stock recipes, concise source provenance and meaningful regression tests.
  Keep compact recovery studies and their current report sources; retire
  stock-specific experiment runners and duplicate delivery reports.
- Judge files by behavior rather than names: a review script that populates
  canonical assumptions is provenance and must remain reproducible.
- Before deleting a report, move lasting definitions or current limitations
  into maintained guidance/source notes and repair its consumers and links.
  Remove tests of retired review artifacts without weakening scientific checks.
- Ask about genuinely ambiguous material. Preserve local source/review caches
  and explicitly designated summary figures/tables unless authorized otherwise.
- Finish each milestone by removing stale scratch files. Preserve user edits
  and report open/locked files rather than closing applications forcibly.

## Versioning and releases

- Increment versions for releases, not individual commits or pushes.
- Use patch releases for fixes without intentional API or model changes;
  use minor releases for new capabilities or substantial pre-1.0 API changes.
- Use a development suffix (e.g. `0.1.1.9000`) for unreleased work.
- Maintain concise `NEWS.md` entries and identify breaking changes.
- The evaluation and fleet expansion targets a future `0.2.0`; do not tag,
  release, or merge a branch without the user's request.

## Before finishing a task

Report:
- files changed;
- API or mathematical changes;
- tests added/updated;
- targeted test results;
- full test-suite result;
- `R CMD check` result if run;
- any remaining warnings or limitations.

Do not merge branches unless explicitly asked.
