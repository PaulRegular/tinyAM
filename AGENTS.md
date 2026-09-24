
# tinyAM development conventions

tinyAM is an experimental R package built with RTMB. Keep changes small,
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
- Include mathematical definitions for model parameters/processes when useful.
- Keep terminology consistent across `make_dat()`, `make_par()`, `nll_fun()`,
  simulation, tidying, and documentation.
- Regenerate documentation after modifying exported interfaces.

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
  `make_dat()`, `make_par()`, `nll_fun()`, simulation, tidying,
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
