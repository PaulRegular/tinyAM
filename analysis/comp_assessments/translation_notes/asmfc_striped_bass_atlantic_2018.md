# Atlantic striped bass 2018 translation

The registered recipe fits pooled Chesapeake Bay and Coast removals from the SCA in Appendix B10 of the 2018 benchmark assessment. It shares Maryland spawning-stock survey (MD SSN) catchability across ages 4–15+. The original fit estimated separate catchability at each supported survey and age. On SAIO `cecb45ad`, that fit passed numerical checks over 1982–2017. The median absolute difference for female SSB was 28.9%, with trend correlation 0.986 and a terminal difference of −27.4%. The translation approximates female SSB using the available stock biology.

The trial used separate Bay and Coast catches with the source's shared catch-weight table. It used fleet-specific AR1 F around time-blocked means. Bay selectivity used a double logistic in the source periods starting in 1982, 1985, 1990 and 1996. Coast selectivity used an ascending logistic in the source periods starting in 1982, 1985, 1990 and 1997. The source uses a Thompson dome for the Bay and a Gompertz curve for the Coast. These tinyAM curves represent the corresponding broad shapes.

The trial prepared successfully. The default fit exceeded the intended 20-minute limit and failed the gradient check (0.0356 versus 0.01). The earliest Bay block's estimated curve midpoint was outside the observed age range. The source-start retry timed out before returning. The hybrid has no validated comparison series. The pooled recipe remains preliminary.

## Time-shared selectivity sensitivity

The follow-up kept the same fleet catches, time-block F levels and AR1 deviations. It fitted one normal dome for Chesapeake Bay and one ascending logistic for the Coast, with each shape shared across periods. The optimizer returned code zero and a maximum gradient of 0.00199, but the Hessian was not positive definite. The Bay curve was effectively flat, with a fitted width above 11,000 ages and poorly determined peak. The weakest curvature direction was almost entirely the Bay width. The Coast curve and the shared F innovation SD and correlations had finite intervals.

In the pooled baseline, terminal total abundance was 4.4% above the reference, while abundance at ages 14 and 15+ was only 55% and 12% of the reference. Applying the translation's biology to the official N changed terminal reference SSB by only 3.8%. The main discrepancy is therefore mature age composition. The time-shared curves did not yield a numerically acceptable replacement for the pooled recipe.

## Source-informed start on SAIO `312bc721`

The source record and prepared model matched the previous time-shared trial. A single restart initialized each fleet's F from catches and official N and Z using the Baranov equation. The Chesapeake Bay normal curve started at a peak of 5.5 years and a width of 2.01 years. Other starts used the same source states and survey data. All states and selectivity parameters were re-estimated.

The restart passed numerical checks and estimated a finite Bay dome, with peak 8.66 and width 2.24 years. The objective declined from 4850.45 to 4833.79 for the same specification. Across 36 years the median SSB difference was 82.5% with trend correlation 0.909. Terminal SSB was 145.8% above the official series, with age-15+ abundance 6.26 times the reference. The startup change resolved the numerical problem and left the age-composition and SSB discrepancies. The pooled model remains the registered preliminary recipe.

## Survey precision sensitivity

Published index CVs matched all 412 fitted survey-years. Annual composition effective sample sizes were unavailable for the 202 fitted composition survey-years, so this comparison used annual relative precision only. The pooled fit supplied `sqrt(log1p(CV^2))` and estimated one multiplicative log-SD scale for each survey. It converged from defaults. Median SSB difference was 28.9% and trend correlation was 0.990, with a terminal difference of −34.0% versus −27.4% in the baseline. This weighting left the scale mismatch.

## Maryland spawning-stock survey catchability

The 2018 source fixes MD SSN selectivity at one above age 3 after correcting for mesh-size selectivity. The translation estimates one catchability for ages 4–15+ and retains separate estimates at ages 2 and 3. The index likelihood keeps a constant log-SD within each survey. MD SSN contributes 28 positive observations at age 15+, including 2017. [SAW-66 assessment, section B4.21.1.2 and Appendix B10](https://repository.library.noaa.gov/view/noaa/23031/noaa_23031_DS1.pdf)

The change converged from defaults on SAIO `312bc721`. Median absolute differences over the same 36 years fell from 28.9% to 8.8% for SSB, 7.9% to 3.8% for abundance, 7.6% to 5.5% for recruitment and 9.2% to 5.1% for average F. Terminal SSB moved from 49,682 to 74,046 t, versus 68,477 t in the comparison reference. SSB trend correlation declined from 0.986 to 0.968, and the terminal total-abundance difference increased from +4.4% to +12.7%. Terminal abundance at age 15+ rose from 69,312 to 751,526 fish, versus 590,437 in the official assessment. The fitted age-15+ catchabilities fell in all six surveys with positive observations at that age, consistent with the higher fitted abundance shared across surveys.
