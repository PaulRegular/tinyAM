# DFO Atlantic herring, NAFO 4VWX

## Assessment and source selection

The Charbonneau-Keith metadata row `ck_dfo_mar_harengus_4vwx` identifies DFO Atlantic herring in NAFO areas 4VWX. Its common name is blank in the canonical stock table. DFO treats 4VWX as a management unit containing several distinct spawning and fishery components.

The latest published advice located is DFO Science Response 2025/026, *4VWX Herring Stock Status Update to 2024*. It reports landings, component-level acoustic indicators, and a 2022-2024 three-year mean acoustic SSB index of 296,430 t against the 317,846 t limit reference point. It gives no catch advice. The current CSAS schedule lists a March 31, 2027 assessment for the SWNS/BoF component with an update of indices for 4VWX; no results from that process were published as of October 23, 2026.

The 2022 assessment used a management strategy evaluation (MSE) framework for the Southwest Nova Scotia/Bay of Fundy (SWNS/BoF) component. The 2022 Science Advisory Report says other spawning components continue to be reported using indicators from the 2018 assessment. The detailed 2018 evaluation (Research Document 2020/062) describes the stock-wide advice as indicator-based rather than a single analytical age-structured assessment. The 2019 data-input report (Research Document 2020/028) documents commercial catch and age data, acoustic surveys, biology, and other candidate inputs.

The 2022 MSE is component-specific. Its public conditioning report (Research Document 2023/022) links to a GitHub repository with input RDS/CSV files and scripts to build reference and sensitivity operating models. The pinned local snapshot is commit `515857664de2719f7fd1b227b41862b82673ac56`. The repository does not contain saved fitted operating-model state surfaces. DFO's 2024 update says the reference operating models were conditioned through 2020; the source code/data snapshot and any run outputs need to be reconciled to that production reference set before deriving official states.

## Translation readiness

The available materials do not yet identify one production fitted N/F/M-at-age surface for the complete 4VWX management unit. The 2022 framework is an ensemble of alternative operating models for only SWNS/BoF, while other components are represented by separate indicators. Treating one reconstructed operating model as a stock-wide accepted assessment would change the assessment scope and select one scenario from the reference set without an official basis.

No assessment rows or numerical database values are added in this source review. The published 2025 update can support a separate summary-only record of component indicators, but that record would not provide a direct tinyAM fit or the full outputs needed for a comparison dashboard. A component-specific tinyAM approximation may be considered after the official production reference OM and its state outputs are reproduced and selected without ambiguity. No fitting, source-based starting values, or dashboard has been attempted.

The source PDFs and pinned code/data archive are cached under `source_cache/ck_dfo_mar_harengus_4vwx/`; checksums and URLs are recorded in `accessible_sources_manifest.csv`. The cache remains local and gitignored.

## Sources

- [DFO Science Response 2025/026: 4VWX Herring Stock Status Update to 2024](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/4130293x.pdf)
- [DFO Science Response 2024/030: Stock Status Update of 4VWX Herring for the 2024 Fishing Season](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41251349.pdf)\n- [DFO Science Advisory Report 2022/050: 2022 Assessment of 4VWX Herring](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41087719.pdf)
- [DFO Research Document 2020/062: 2018 Evaluation of 4VWX Herring](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/40951686.pdf)
- [DFO Research Document 2020/028: 2019 SWNS/BoF Herring Framework Data Inputs](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/40946265.pdf)
- [DFO Research Document 2023/022: Framework for Conditioning SWNS/BoF Herring Operating Models](https://waves-vagues.dfo-mpo.gc.ca/library-bibliotheque/41206903.pdf)
- [Pinned public operating-model repository](https://github.com/z5a1n/herring_OM_conditioning/tree/515857664de2719f7fd1b227b41862b82673ac56)
- [Current CSAS assessment schedule](https://csas-scas.dfo-mpo.gc.ca/)

