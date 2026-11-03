# Baltic plaice in ICES Subdivisions 21-32: WGBFAS 2026: tinyAM translation notes

Transferred from the [SAIO source review](https://github.com/PaulRegular/saio/blob/4b347a10b95e8d1f9b3d34156bf58ed63c82f775/source_reviews/ices_plaice_baltic_21_32_2026.md) during the repository cleanup. These are previously recorded settings, caveats and diagnostics, not new fitting results. The registered stock script defines the current model; check its curated or preliminary tier before using these notes.

Official assessment evidence is maintained in [SAIO](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_plaice_baltic_21_32_2026.md).

## Expanded-library review

### tinyAM translation

The translation retains the 2002-2025, ages 1-7+ observation grid; Q1 and Q4 timing; source q and observation-variance sharing; fixed time-varying M; and source maturity and weight surfaces. It uses random-walk recruitment, IID N-process errors, an AR(1) F process, exponential initial abundance, and the official Fbar ages 3-5.

The translation cannot reproduce the full SAM state-space structure, age-correlated F random walk, estimated stock-weight process, or all source likelihood details. These are stated in the stock recipe and comparison dashboard. No observations are adjusted to improve agreement, and no explicit starting values are used.
