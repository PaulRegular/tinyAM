# Step 5 source-readiness review

Reviewed from the local SAIO restore identified by manifest `saio-local-2026-11-14-step5-sole-herring`. The manifest is tracked in SAIO commit `e6f2d7e0ba07f228ab3aa9ab34370777f1275573`; its checksums identify the restored CSV tables. The three registered preliminary recipes below stop before constructing observations or fitting. The three assessment records are included in the published SAIO expansion and curation releases. The historical 28-stock snapshot predates these records. The recipes record where the official source inputs do not support tinyAM's age-based observation structure. No source values are substituted, and no fitted source states are used as observations.

| Assessment | Translation blocker |
|---|---|
| `noaa_goa_rex_sole_2026` | The final-run numerical age-composition inputs are not available; the public appendix is marked draft. The report does not provide an age-specific survey index, and numeric run-matched age biology is incomplete. Reconstructing marginal age observations from conditional age-at-length samples would require unsupported source mappings. |
| `iphc_pacific_halibut_2025` | Public tables provide aggregate index and output series, but not exact catch-age likelihood vectors, age-specific survey observations, or run-matched numeric weight and maturity surfaces. The native input files are permission-restricted. |
| `noaa_chilipepper_2025` | The official likelihood uses catch biomass, aggregate indices, conditional age-at-length and length compositions. It does not provide observed catch-at-age or an age-specific abundance index. Model-estimated age surfaces cannot be reused as observations, and the reported biology is length-based. |

These records are not ready for a tinyAM fit. Revisit only if a source adds the missing age-resolved observations and compatible biological inputs for the same official assessment run. The 2026 Irish Sea sole and Atlantic herring imports were not repeated in this review.
