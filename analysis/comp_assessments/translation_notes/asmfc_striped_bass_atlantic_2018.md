# Atlantic striped bass 2018 translation

The registered recipe fits pooled Chesapeake Bay and Coast removals from the SCA in Appendix B10 of the 2018 benchmark assessment. On SAIO `cecb45ad`, it passed numerical checks over 1982–2017. The median absolute difference for female SSB was 28.9%, with trend correlation 0.986 and a terminal difference of −27.5%. The translation approximates female SSB using the available stock biology.

The trial used separate Bay and Coast catches with the source's shared catch-weight table. It used fleet-specific AR1 F around time-blocked means. Bay selectivity used a double logistic in the source periods starting in 1982, 1985, 1990 and 1996. Coast selectivity used an ascending logistic in the source periods starting in 1982, 1985, 1990 and 1997. The source uses a Thompson dome for the Bay and a Gompertz curve for the Coast. These tinyAM curves represent the corresponding broad shapes.

The trial prepared successfully. The default fit exceeded the intended 20-minute limit and failed the gradient check (0.0356 versus 0.01). The earliest Bay block's estimated curve midpoint was outside the observed age range. The source-start retry timed out before returning. The hybrid has no validated comparison series. The pooled recipe remains preliminary.
