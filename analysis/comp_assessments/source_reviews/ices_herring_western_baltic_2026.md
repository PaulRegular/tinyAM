# Western Baltic spring-spawning herring (her.27.20-24), 2026

## Formula review (2026-10-10)

An IID older-cohort survival candidate was tested separately from the retained
recruitment RW, q sharing and F process. It reported false convergence
(optimizer code 1, gradient 0.042), a non-positive-definite Hessian and a
near-zero variance component. Scale disagreement also increased. Retain the
baseline; the source hockey-stick recruitment and correlated observation/F
processes remain different. No BH/Ricker curve was substituted for hockey-stick.

## Assessment represented

The canonical record represents the accepted 2026 ICES HAWG assessment, fitted for 1991–2025. The report identifies the production run as WBSS_HAWG_2026, a single-fleet SAM assessment available from stockassessment.org. The official fit object is cached locally with the full HAWG 2026 report and the ICES Stock Assessment Graphs source page. The model object reports SAM 0.12.0, optimizer convergence code 0, and objective 414.9946.

## Inputs and assumptions

The cached fit contains all 315 total-catch-at-age observations and 372 finite survey observations (five HERAS ages 2-6 are missing in 1999), annual stock and catch weights, maturity, fixed natural mortality, spawning fractions, and the finite observation precision weights supplied to SAM and their equivalent relative SD factors. Survey names, age ranges, plus groups, and sampling fractions are taken from the accepted object. The report and native configuration specify the source age sharing and likelihood structures.

The report states that M is time-invariant and age-specific, derived from North Sea autumn-spawning herring M, and profiled during the 2025 benchmark. The production model uses mortalityModel=0; the M values in the database are fixed inputs, not estimated by this assessment. Maturity is time-invariant. Fractions of annual fishing and natural mortality before spawning are 0.168 and 0.25.

SAM treats recruitment with a hockey-stick relationship (code 61), gives recruitment a separate process variance, shares the survival-process variance across ages 1–8, and correlates F states across ages with AR(1). The accepted F state for ages 7 and 8 is shared. Observation residuals use AR(1) across ages for HERAS, GERAS and IBTS/BITS Q1, while N20 residuals are independent. The model object's keyLogFpar specifies catchability sharing by survey and age group.

## Output checks

Annual recruitment, SSB, Fbar ages 2–5, and total stock biomass extracted from report Table 3.6.11 agree with the native fit after applying the keyLogFsta mapping to expand fitted F to ages. Differences are less than one unit for counts and biomass and less than 0.00051 for Fbar (report rounding). The native stock-number and F-at-age surfaces agree with report Tables 3.6.12–3.6.13 to their displayed precision. The 2025 reported values are recruitment 20,038,259 thousand fish, SSB 110,526 tonnes (95% CI 85,439–142,978), Fbar 0.011 (95% CI 0.008–0.015), and total stock biomass 266,183 tonnes.

## Source references

- ICES HAWG 2026, Report of the Herring Assessment Working Group, Section 3.6, especially pp. 173–174 and Tables 3.6.4–3.6.13, pp. 247–256: [cached report](../source_cache/ices_herring_western_baltic_2026/HAWG_2026_Full_Report.pdf).
- ICES Benchmark Workshop on Selected Herring and Sprat Stocks (WKBHERSPRAT), 2025, DOI [10.17895/ices.pub.31538287](https://doi.org/10.17895/ices.pub.31538287).
- Accepted native SAM object WBSS_HAWG_2026/run/model.RData, cached locally at [model.RData](../source_cache/ices_herring_western_baltic_2026/model.RData).
- ICES Stock Assessment Graphs, source-data key 21294: [official source page](https://sg.ices.dk/ViewSourceData.aspx?key=21294), cached locally at [HTML copy](../source_cache/ices_herring_western_baltic_2026/ICES_Stock_Assessment_Graphs_key_21294.html).

## Remaining limitations

The report does not tabulate age-specific uncertainty for the population or mortality surfaces. The native object contains fitted states but tinyAM will use the source estimates only as starting values. tinyAM cannot reproduce the hockey-stick stock-recruitment relationship, the correlated cross-age F process, or AR(1) age residuals in three survey likelihoods. The translation will retain source M, survey timing, catchability-sharing keys, and observation variance groups where tinyAM permits.
