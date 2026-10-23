## Companion to the tinyAM cheatsheet

This dashboard accompanies the [tinyAM cheatsheet](https://paulregular.github.io/tinyAM/cheatsheet/tinyAM_cheatsheet.pdf). It shows how three model structures fitted to the same Northern cod example data can produce different population estimates and diagnostics.

The fits use `cod_obs` for 1983-2025 and model ages 2-14, with age 14 as the plus group. Weights and maturity retain the older biological ages. The examples resemble statistical catch-at-age (SCAA), state-space assessment model (SAM), and bespoke Northern cod assessment model (NCAM) structures. They illustrate alternative assumptions rather than reproduce those assessments.

| Model | General structure |
|---|---|
| **SCAA** | Deterministic cohort survival, annual recruitment parameters, separable fishing mortality, and fixed natural mortality. |
| **SAM** | Abundance deviations, random-walk recruitment, correlated random walks in fishing mortality, and fixed natural mortality. |
| **NCAM** | Deterministic cohort survival, random-walk recruitment, and AR1 variation in fishing mortality and age-blocked natural mortality. |

## Exploring the results

Use **Fishery** and the survey pages to compare observations with fitted values, **Trends** to compare population estimates such as spawning stock biomass, and **Parameters** to inspect model components. The **Tables** menu includes numerical and statistical diagnostics. Hover over plots for values and units; shaded bands show 95% confidence intervals where available.

This dashboard shows the main fits. Use the numerical and statistical diagnostics to assess model fit before interpreting differences. The separate hindcast call in the cheatsheet assesses survey forecast skill conditional on observed catches.

The complete specifications are printed in the cheatsheet and maintained in [the case-study models](https://github.com/PaulRegular/tinyAM/blob/main/analysis/case_study/001_models.R). To reproduce this example dashboard with its background text, run `Rscript cheatsheet/examples.R dashboard` from the repository root.
