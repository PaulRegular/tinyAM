# First source-curation batch

Selected 2026-11-02 from the screened NOAA, ICES and DFO register. The
selection uses the 90 distinct candidate units classified `new` and `eligible`
in `inventory_units.csv`; mixed-status units are held out until their
eligibility is resolved. The selection balances the three authorities while
prioritizing detailed assessment inputs and official outputs. This is a
source-curation batch, not a translation-readiness ranking.

| Priority | Organization | Candidate unit | Stock | Inventory assessment year |
|---:|---|---|---|---:|
| 1 | DFO | `dfo_greenland_halibut_2_3klmno` | Greenland halibut, NAFO 2 + 3KLMNO | 2026 |
| 2 | DFO | `dfo_sff_2024_gulf_white_hake_nafo_4t` | White hake, NAFO 4T | 2022 |
| 3 | DFO | `dfo_sff_2024_gulf_yellowtail_flounder_nafo_4t` | Yellowtail flounder, NAFO 4T | 2025 |
| 4 | DFO | `dfo_western_component_pollock_4xopqrs5` | Pollock, western component | 2024 |
| 5 | DFO | `dfo_yellowmouth_bc` | Yellowmouth rockfish, coastwide BC | 2021 |
| 6 | DFO | `dfo_sgsl_winter_flounder_4t` | Winter flounder, NAFO 4T | Unspecified in inventory |
| 7 | ICES | `ices_atlantic_wolffish_iceland_5a` | Atlantic wolffish, Iceland | 2026 |
| 8 | ICES | `ices_cod_west_greenland_inshore` | West Greenland inshore cod | 2026 |
| 9 | ICES | `ices_cod_west_greenland_offshore` | West Greenland offshore cod | 2026 |
| 10 | ICES | `ices_golden_redfish_561214` | Golden redfish, subareas 5, 6, 12 and 14 | 2026 |
| 11 | ICES | `ices_greater_silver_smelt_5b6a` | Greater silver smelt, divisions 5.b and 6.a | 2026 |
| 12 | ICES | `ices_her_27_irls` | Herring, Irish and Celtic Seas | 2026 |
| 13 | ICES | `ices_seabass_8ab` | Seabass, divisions 8.a and 8.b | 2026 |
| 14 | ICES | `ices_whg_27_7b-ce-k` | Whiting, divisions 7.b, 7.c and 7.e-k | 2026 |
| 15 | NOAA | `noaa_acadian_redfish_gom_gb` | Acadian redfish, Gulf of Maine / Georges Bank | 2025 |
| 16 | NOAA | `noaa_stocksmart_10466` | Sablefish, Pacific Coast | 2025 |
| 17 | NOAA | `noaa_stocksmart_10568` | Widow rockfish, Pacific Coast | 2026 |
| 18 | NOAA | `noaa_stocksmart_10426` | Yelloweye rockfish, Pacific Coast | 2025 |
| 19 | NOAA | `noaa_stocksmart_10415` | Canary rockfish, Pacific Coast | 2025 |
| 20 | NOAA | `noaa_black_rockfish_oregon_2023` | Black rockfish, Oregon | 2023 |

The two West Greenland cod units are separate current ICES assessment units,
not separate model variants. Candidate screening has identified public
age-composition or age-specific inputs, biology, mortality/timing evidence,
and official outputs for each selected unit. Before each import, verify the
latest assessment used for advice, the detailed run and terminal year, the
native observation basis, units, source coverage, and whether the source files
are recoverable. The year shown above is the inventory's assessment year; it
does not replace that run-level verification. Preserve partial latest-unit
records when a series cannot be recovered. Do not substitute an older stock
definition to fill a gap.

The selected batch is recorded in `selected_batch_01_2026.csv`. Candidate
years are being refreshed against current official advice before each source
import. The Greenland halibut source review and partial database import are
recorded separately; its 2023 age-structured source is not represented as the
model used for the 2026 management-procedure advice.
