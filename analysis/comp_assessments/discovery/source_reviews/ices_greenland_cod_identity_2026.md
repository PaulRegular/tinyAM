# Greenland cod identity and source screen

**Review date:** 2026-11-01  
**Inventory:** ICES Stock Information Database snapshot retrieved 2026-10-28  
**Scope:** Identity and first-pass age-input/output screen; no canonical records added.

## Identity and assessment boundaries

ICES lists three Greenland cod advice units in 2026: West Greenland offshore spawning cod (`cod.21.1.osc`), West Greenland inshore spawning cod (`cod.21.1.isc`), and East Greenland-Iceland offshore spawning cod (`cod.21.27.1.14`). Their shared fishery does not make them one assessment unit.

The West Greenland inshore advice code is one advice unit, but the 2026 report describes two distinct assessment components: northern inshore cod (`cod.21.1a-c`, N-WISC) and southern inshore cod (`cod.21.1d-f`, S-WISC). Do not count those components as two additional inventory stocks. Preserve their separate assessment evidence if this unit is curated.

These three units are distinct from the ICES cod stocks already in the canonical database, which cover the North Sea/West of Scotland and Northeast Arctic; the DFO and NOAA cod records cover other areas.

## First-pass screening

- **West Greenland offshore (`cod.21.1.osc`) — eligible for source curation.** The current 2026 report describes an analytical SAM assessment over 2000-2025, commercial catch-at-age, age-specific bottom-trawl survey indices, catch and stock weights-at-age, maturity fixed between years, and natural mortality fixed at 0.2 for all ages. The report identifies the fit as `NWWG26_WOSC_v1` and provides recruitment at age 2, SSB, catch, and F ages 5-8 with uncertainty intervals. Exact input matrices, age coverage, source units, survey timing fractions, and fitted age surfaces still need to be checked during curation.
- **West Greenland inshore (`cod.21.1.isc`) — eligible for partial source curation, with a component limitation.** The current report retains a Category 1 SAM assessment for N-WISC (`cod.21.1a-c`) and reports comparable age-2 recruitment, SSB, catch, and F ages 5-8. It identifies commercial catch-at-age, a gillnet survey targeting ages 2-3, age-based weights, fixed maturity, and M=0.2 for all ages. The same report says the S-WISC (`cod.21.1d-f`) assessment was changed to Category 5: its 2026 model was not accepted, the gillnet index was excluded after a timing shift, and calculated catch-at-age and weights-at-age were not used. Advice is issued under the combined `cod.21.1.isc` code. Keep N-WISC and S-WISC separate in any source record and label the current S-WISC model results unavailable/non-comparable. The N-WISC report names run `NWWG25_N-WISC` while reporting 2026 results; verify the exact run/data revision during source curation.
- **East Greenland-Iceland offshore (`cod.21.27.1.14`) — ineligible for the age-based comparison batch.** Age-disaggregated commercial and survey data are described, but the 2026 report says the SAM model failed to fit observed catches for some ages (up to 92%) and is not a robust basis for advice. ICES uses Category 5 catch-based advice with a 20% precautionary buffer; no comparable accepted population or F outputs were identified. Keep this as a distinct current stock unit. This screening decision does not prevent recording verified partial current-unit inputs if later source curation finds them.

## Sources and cached files

The official [ICES 2026 Northwestern Working Group report](https://doi.org/10.17895/ices.pub.32414577) was retrieved as the stock-specific chapter PDFs below. The three official advice pages are also cached. All files are in the gitignored `source_cache/discovery_ices_greenland_cod_2026/`; the local `sources.csv` manifest records URLs, retrieval date, bytes and SHA-256 values. The script `491_screen_ices_greenland_cod_2026.R` checks these before applying the screening decisions.

| Report chapter | Source file | Figshare file URL | SHA-256 |
|---|---|---|---|
| 13, West Greenland offshore spawning cod | `nwwg_2026_west_greenland_cod.pdf` | [65495862](https://ndownloader.figshare.com/files/65495862) | `AA905234DF5792401A83CBB6AEECE710143FC4F8B3A8A16DBF0DC8959642EC46` |
| 14, Southwest Greenland inshore spawning cod | `nwwg_2026_southwest_greenland_cod.pdf` | [65495868](https://ndownloader.figshare.com/files/65495868) | `B1486382F7456DF8DDF71BF29409BA004E1B180BB8CEC916F58FDE955CFC7137` |
| 15, Northwest Greenland inshore spawning cod | `nwwg_2026_northwest_greenland_inshore_cod.pdf` | [65495871](https://ndownloader.figshare.com/files/65495871) | `FF4F837590C8115E496F9C002BF296363DE3E458C76361FBBC9BFBB9692952F5` |
| 16, East Greenland-Iceland offshore spawning cod | `nwwg_2026_east_greenland_iceland_offshore_cod.pdf` | [65495874](https://ndownloader.figshare.com/files/65495874) | `D32EBED8DB3F91B0E2F13C1F2FCEE286AC3F93B6C039A51ACFEE1ADC00C16EB1` |

The 2026 advice pages are [5088](https://asd.ices.dk/viewAdvice/5088), [5089](https://asd.ices.dk/viewAdvice/5089), and [5091](https://asd.ices.dk/viewAdvice/5091). The official SID snapshot and advice establish identity and advice status; neither alone proves age-specific inputs. No canonical assessment, input, or output rows have been added.
