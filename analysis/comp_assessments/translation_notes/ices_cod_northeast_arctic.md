# ices cod northeast arctic: tinyAM translation notes

Inherited settings and diagnostics; no models were rerun. Check the registered recipe and review tier for current settings.

## Initial abundance

The cached official SAM fit has `initState = 0` and no `initN` parameters.
The translation retains `init = "free"`, estimating older initial abundances
as fixed effects. This differs from integrating the source initial states
without an initial density; first-year recruitment remains an estimated fixed
effect in tinyAM.

[Official source review](https://github.com/PaulRegular/saio/blob/main/source_reviews/ices_cod_northeast_arctic.md).
