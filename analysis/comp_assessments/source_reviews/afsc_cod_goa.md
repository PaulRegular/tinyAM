# Gulf of Alaska Pacific cod: accepted-assessment source review

Status: accepted-run/source inventory in progress; no canonical record imported yet.
Charbonneau identifier: AFSC_GOA_Gadus_macrocephalus.

## January 2026 assessment

The January 2026 update (presented February 2026) uses Model 24.0 accepted in 2024, updated through 2025. The 54-page detailed update and 2024 full model-defining chapter are cached under source_cache/afsc_cod_goa_2026/.
Current detailed source: https://meetings.npfmc.org/CommentReview/DownloadFile?fileName=C1+GOA+Pcod+Assessment.pdf&p=f00593eb-12f5-458c-842e-a5cdd45306bb.pdf
Framework/full 2024 source: https://files.npfmc.org/SAFE/2024/GOApcod.pdf

The Council's February 6, 2026 announcement confirms that the updated assessment was reviewed by the Plan Team and SSC and superseded December 2025 specifications. Thus the 2025 shutdown catch-only rollover is not the current accepted assessment for this stock.
Acceptance source: https://www.npfmc.org/council-recommends-substantial-increase-to-gulf-of-alaska-pacific-cod-catch-limits-for-2026-27/

The report states 2026 ABC 41,520 t and OFL 49,782 t, projected total age-0+ biomass 182,156 t, and female SSB 52,772 t. These are projection values, not historical 2025 estimates. Historical values must be extracted from their own tables without mixing columns or assessments.

## Native model files

The report page 2 directly links to https://afsc-assessments.github.io/goapcod/2025_Assessment/January_Model/.
Corresponding repository: https://github.com/afsc-assessments/goapcod
Inspected revision: e632807e4947686c16caf99b864bfe8466f8dbca.
The native archive docs/2025_Assessment/January_Model/model_files/M24.0_SS3_files.zip is cached and unpacked locally. It contains GOAPcod2025Dec08.dat, Model24_0.ctl, starter.ss and forecast.ss. No fitted output file is included in this archive; linked diagnostic/figure pages may provide additional outputs.

The starter names those exact data/control files and has retro_yr=0. The data specify 1977–2025, one annual season, two subseasons, spawning month 1, one sex and area, maximum modeled age 10, and seven fleet definitions. The three removal fleets are FshTrawl, FshLL and FshPot. The report explicitly says the longline category includes jig catches. Survey definitions also include Srv, LLSrv, ADFG and Seine; determine which actually contribute positive-weight observations before assuming four fitted surveys. Report inventory identifies bottom trawl and longline surveys as fitted abundance/composition streams.

Native inputs include length compositions and conditional age-at-length observations. Do not collapse these into age-only catches or infer an age composition from model predictions. A faithful representation will need to retain conditioning length bins and native weighting, with a narrow schema extension if required.

The control uses estimated natural mortality with a temporal block including 2014. Its parameter starting values are not estimated M outputs and must not be exported as fixed M inputs. Maturity is length-logistic in the control; do not substitute a guessed age-maturity vector. Detailed parameter, temperature-covariate, biology, survey timing and uncertainty interpretation remains pending.

An older author repository, https://github.com/pete-hulson/goa_pcod at facf41573f9a0b609d0096611bc9302aaab43abe, was inspected initially; the current report's directly linked repository is preferred. Earlier 2024 numerical files are context only, not substitutes for the current accepted run.

## Native observation inventory

The source-specific inventory reader cached in inventory_native_observations.py checks column counts and section sentinels, and stages every source cell in native_observation_sections_raw.csv. It found 150 catch rows: three initial-equilibrium rows and 147 historical fleet/year catches (49 years for each of trawl, longline/jig and pot). The index section has 110 source rows, of which 52 have positive year and fleet codes: 17 bottom-trawl indices and 35 longline indices. ADFG/Seine definitions do not imply active fitted indices.

Length compositions have 182 positive year/fleet rows: trawl 48, longline/jig 46, pot 36, bottom trawl 17, longline survey 35. The age section has 923 source rows, of which 857 have positive year/fleet codes: trawl 205, longline/jig 190, pot 168 and bottom trawl 294. These are composition records, not counts of unique observation years. Their conditioning length intervals, sample sizes, age-error codes, partition and sex fields must be retained. Remaining age/index records include negative codes and must not be counted as fitted streams.

The official Stock Synthesis manual explains that a negative fleet code excludes a composition observation's likelihood contribution, even though predictions and diagnostics may still be calculated. It also confirms that the observation's month determines survey timing; the fleet-definition timing field is not sufficient.
Manual: https://nmfs-ost.github.io/ss3-doc/SS330_User_Manual_release.html

The repository worktree had concurrent user edits to database_to_tiny_obs.R and 004_test_translation.R during this review; these are outside the source extraction work and must remain untouched.

## Historical summary staging

stage_historical_summaries.py cached 147 current-run values (49 years each for female SSB, total age-0+ biomass and age-0 recruitment) from Tables 2.7 and 2.8. Row counts passed. Earlier-assessment columns and the separate 2026 forecast were excluded. Reported SDs are retained separately without inventing intervals.

Recruitment units need resolution before canonical import: Table 2.8 labels values such as 0.48 as millions, whereas Table 2.6 reports log(mean recruitment)=13.09 and Stock Synthesis normally represents numbers in thousands. exp(13.09) is approximately 484,000 model units, or 484 million fish; this suggests the printed table may actually be in billions, but it is not yet direct evidence of each reported recruitment's scaling. Obtain native time-series outputs or verify the table-generation scale before changing the unit or importing a guessed correction. Preserve the original printed values in the source staging table.

## Recruitment scale resolved and accepted output verified

The author repository DOES contain the current management run under 2025/mgmt/24.0. Its data and control numeric contents match the report-linked January archive exactly after removing comments/whitespace; byte hashes differ because formatting/comments differ. The cached ss3.rep is a compact fitted report, not a full Report.sso, and contains historical spawning output/recruitment and selectivity information.

Direct reporting code 2025/R/safe_tbls.R at facf41573f9a0b609d0096611bc9302aaab43abe divides recruitment Value and StdDev by 1,000,000. Since native SS numbers are thousands, Table 2.8 is in BILLIONS of fish despite its printed millions heading. It also divides native spawning output Value and StdDev by 2 for reported female SSB. Thus native unadjusted spawning output is not directly the published female SSB.

verify_historical_summary.py checked all 49 historical years: native recruitment / 1e6 matches Table 2.8 to 0.0051 printed units, and native spawning output / 2 matches Table 2.7 female SSB to one tonne (compact-report rounding). Both checks passed. The verified staging file records billion fish units; the original raw extraction retains the misleading printed label and remains cached. No input observations were adjusted using fitted parameters.

Source code: https://github.com/pete-hulson/goa_pcod/blob/facf41573f9a0b609d0096611bc9302aaab43abe/2025/R/safe_tbls.R
Native output: https://github.com/pete-hulson/goa_pcod/blob/facf41573f9a0b609d0096611bc9302aaab43abe/2025/mgmt/24.0/ss3.rep
## Composition and timing checks

All 857 positive-year, positive-fleet age records are conditional age compositions. Each has identical lower/upper conditioning labels, ranging from 4.5 to 104.5, an age-error code of 1 or 2, and partition 0. Their ten age proportions sum to between 0.99999 and 1.00001, consistent with source rounding. These records must retain the conditioning labels rather than being collapsed into annual unconditional age compositions.

The data file sets Lbin_method=1 (population-length bins), while the observation labels are half-integers matching the length-composition labels. The interpretation of these labels needs verification against the fitted SS version/implementation before treating them as physical centimetre bounds. Preserve native values; do not silently reinterpret them as integer bin indices.

Every positive index observation uses month 7. Length and age composition records use months 1 or 7. Survey timing must therefore be derived from the observation month, not the fleet-definition default. Fishery composition timing and continuous catch timing remain distinct.

The composition schema validator now checks finite numeric dimensions, non-negative integer age-error/partition codes, positive supplied sample sizes, observation IDs and ordered conditioning labels before announcing validation success. The existing canonical database passes structural validation; the Pacific cod composition rows are not yet imported.
## Canonical coverage and validation

Assessment afsc_cod_goa_2026 now contains 12,643 input cells: 147 catch totals, 52 index observations, 52 index log SDs, 3,822 length-proportion cells and 8,570 conditional-age-proportion cells. The 1,039 composition observations retain supplied sample sizes, native length labels/conditioning bounds, age-error codes and partitions. The canonical table also contains 147 historical outputs (49 each for SSB, total biomass and recruitment) and 16 initial assumptions.

The source-specific validator checks each composition cell against the cached native data, verifies observation counts, sample sizes and timing, and confirms that adding optional composition dimensions leaves all earlier input values unchanged. These checks and full database structural validation pass. All completeness statuses remain partial: equilibrium inputs, numerical ageing errors, mean size-at-age, temperature covariates, complete selectivity/q assumptions and N/F/M-at-age outputs still require work. This is an intermediate import, not a completed assessment.
## Additional material input review

The control links longline catchability to environmental variable 1 (env_var&link=101). All 46 supplied annual temperature-covariate values, 1979–2024, are now represented verbatim as environmental covariates, including negative values. Native scaling units remain unresolved; no restandardization was applied. Source-specific validation checks every year/value against the data section.

All 16 mean-size-at-age observations have negative fleet code -4 and therefore no fitted likelihood contribution. Their absence from canonical fitted observations is intentional; source records remain cached. Three initial-equilibrium catch values and their supplied SEs are retained in fleet-specific assumptions because native year -999 is not a calendar year. The two ageing-error mean/SD definitions (age 0 through age 10) are retained once as numerical assumption vectors. Definition 2's -1 mean flags must not be interpreted as measured negative ages or an empirical misclassification matrix.

Pacific cod coverage is now 12,689 input cells, 147 historical outputs and 25 assumptions. Structural and source-specific validation pass. Biological parameter semantics, full selectivity/q structure and fitted N/F/M surfaces remain incomplete; all statuses remain partial. The earlier gap list is superseded for temperature, equilibrium catches and source ageing-error vectors.
## Estimated mortality and fixed biological parameters

Table 2.6 reports estimated M=0.50 (SD 0.023) outside the 2014–2016 block and M=0.84 (SD 0.053) inside it. The native control specifies age-constant M (natM_type=0) with replacement block 4 covering exactly 2014–2016. The database now expands these rounded published estimates to 1977–2025, ages 0–10+, as 539 mortality output cells. These are repeated representations of two shared parameters, not independently estimated age/year values; notes state the rounding and shared uncertainty. No confidence intervals are manufactured.

Fixed weight-length coefficient/exponent, maturity length/slope, stock-recruit steepness and recruitment sigma are now retained as numerical assumptions from negative-phase control parameters. They are input parameters, not estimated age-specific biological surfaces. The compact ss3.rep does not contain N-at-age or M-at-age tables; its exploitation and length-selectivity sections cannot be relabeled as full F-at-age.

Coverage is now 12,689 input cells, 686 outputs and 31 assumptions. Source validation verifies the M block/year/age mapping and absence of invented intervals, alongside native composition/covariate checks; full structural validation passes. N/F-at-age and remaining process/selectivity interpretation still require work. Status remains partial.
