# Eastern English Channel common sole (`sol.27.7d`), 2026: tinyAM translation notes

Transferred from the pre-cleanup SAIO review without refitting. Check the registered stock script and review tier for current settings. [Official source review](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_sole_eastern_channel_7d_2026.md).

## Translation boundary

The fitted tinyAM approximation can use catch-at-age, age-specific surveys, and source M, but cannot reproduce the SAM commercial aggregate FSB indices or all SAM variance-sharing and cross-age model details. Those are documented omissions, not reasons to alter source observations. No fitted SAM object, native bootstrap data, or age-specific uncertainty was found in the public sources reviewed.

The translation follows the age-specific observation-variance keys in Table 17.13 and the pinned TAF configuration: UK BTS age 1, age 2, and ages 3–6 share three variance parameters. Table 17.15 reports different BTS SD estimates for ages 3, 4, and 5–6 (0.0966, 0.0780, and 0.0730), which conflicts with those keys. The recipe does not use those SD estimates as starting values; the report and public TAF snapshot do not resolve the discrepancy. The index scale also remains unresolved, so source q estimates are not used as starts.

The tinyAM recipe fits the full 1982–2025 age 1–11+ grid, uses an IID N process, random-walk recruitment and F, fixed M=0.1, and the source q and observation-variance groupings where representable. BTS timing is set to 0.625, the midpoint of the reported 0.5–0.75 window. The three aggregate commercial FSB indices are retained in the database but excluded from the fit.

The package-default fit did not converge (optimizer code 1, false convergence, objective 1958.75, maximum gradient 31.61, and a non-positive-definite Hessian). An `sdreport` object was returned, but its fixed-effect standard errors were non-finite. One retry used source N-at-age and F-at-age as starting states while retaining package defaults for other parameters; its gradient evaluation was non-finite. An exploratory dashboard was rendered from the failed fit, but its estimates and comparison statistics are not suitable for inference. The comparison derives SSB-at-age from source N × weight × maturity because the published 2021–2025 age-specific SSB rows do not sum to reported annual SSB. Canonical outputs are unchanged.
