Row
----------------------------------------------------------------------------------------

### Dashboard from tinyAM cheatsheet example

This dashboard accompanies the [tinyAM cheatsheet](https://paulregular.github.io/tinyAM/cheatsheet/tinyAM_cheatsheet.pdf). It shows how three model structures fitted to the same Northern cod example data can produce different population estimates and diagnostics.

The fits use `cod_obs` for 1983-2025 and model ages 2-14, with age 14 as the plus group. The examples resemble statistical catch-at-age (SCAA), state-space assessment model (SAM), and bespoke Northern cod assessment model (NCAM) structures. They illustrate alternative assumptions rather than reproduce those assessments.

| Model | General structure |
|---|---|
| **SCAA** | Deterministic cohort survival, annual recruitment parameters, separable fishing mortality, and fixed natural mortality. |
| **SAM** | Abundance deviations, random-walk recruitment, correlated random walks in fishing mortality, and fixed natural mortality. |
| **NCAM** | Deterministic cohort survival, random-walk recruitment, and AR1 variation in fishing mortality and age-blocked natural mortality. |

### Navigating the dashboard

Use **Fishery** and **Surveys** to compare observations with fitted values and inspect residual plots. **Trends** compares population estimates such as spawning stock biomass. The **Parameters** pages show fixed parameter estimates, including coefficients and process or observation standard deviations, and random parameter estimates, including latent states and process deviations across years and ages.

On **Random**, **Process deviations** shows fitted log-scale departures from mean or cohort predictions, or annual increments for random walks. **Joint-sample residuals** standardizes these innovations using the fitted process correlations and variation, then propagates one shared draw from Gaussian parameter uncertainty. This local Gaussian approximation provides reproducible diagnostic values (default seed `123456`). Missing or invalid joint uncertainty is identified explicitly; fitted deviations remain available. Each viewer separates recruitment line plots from N/F/M matrices or bubbles.

The **Tables** pages contain numerical and statistical diagnostics, input data, model settings and calls, and fitted outputs such as population estimates and parameter tables. Use these checks to assess model fit before interpreting differences. The separate hindcast call in the cheatsheet assesses survey forecast skill conditional on observed catches.

Hover over Plotly figures to inspect values and units. Drag to zoom, hold Shift while dragging to pan, and double-click to reset the view. Use the legend to show or hide series and use any animation controls to move between models or years. Shaded bands show 95% confidence intervals where available.
