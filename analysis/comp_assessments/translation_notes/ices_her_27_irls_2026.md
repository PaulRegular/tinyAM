# Celtic Sea herring: preliminary 2026 translation

The recipe retains April–March model-year labels and fits 2002–2024, the
acoustic-era overlap with complete source biology. Published counts are
already in thousands of fish. The acoustic index is at model-year fraction
0.5, with the official 2004, 2017 and 2024 exclusions. The unresolved 2019
age-2 value of 98,229 is preserved even though it exceeds the published total.
Missing 2025 age-9+ weights are not filled. Source records remain unchanged.

The default attempt converged: optimizer code 0, objective 395.6306, maximum
gradient 0.0014380, positive-definite Hessian and available standard errors.
No source-start retry was used; fitted N/F arrays are unavailable anyway.
The local dashboard is under
`results/preliminary/cache/ices_her_27_irls_2026/`.

The SSB comparison is approximate because annual spawning weights are used
without recovered spawning-survival fractions. Its mean absolute percent
difference is 15.9%, terminal difference −24.6% and trend correlation 0.980.
The terminal difference triggers review. Native total biomass and Fbar are
not treated as matched common-definition metrics. Recruitment remains
non-equivalent in automated reporting because the canonical recruitment-age
metadata is unresolved there, despite the review describing age-1 recruitment.
No missing metric is counted as a pass.

Other approximations include deterministic survival, IID recruitment around
a mean, free initial abundance, static fishery selectivity and independent
age errors in place of aggregate/composition likelihoods. Source facts and
conflicts remain in the
[SAIO review](https://github.com/PaulRegular/saio/blob/757b59c/discovery/source_reviews/ices_celtic_herring_2026.md).
