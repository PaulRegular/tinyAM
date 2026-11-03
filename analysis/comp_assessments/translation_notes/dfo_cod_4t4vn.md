# Southern Gulf cod: tinyAM translation notes

Transferred from SAIO without refitting. The registered stock script is authoritative for the current model.

## Mortality and recruitment approximation

The official assessment estimates log-M random walks for ages 2–4, 5–8, and
9+. Each group's initial M level is held constant through 1971 and has a
normal prior on the natural M scale (means 0.65, 0.15, and 0.15; SD 0.05);
log-M innovations begin in 1972 and have a fixed SD of 0.075. The tinyAM
translation now includes 1971 as the free first process state, so the first
penalized change is from 1971 to 1972. The prior means are used only to
initialize the tinyAM states: tinyAM cannot apply the source's initial-M prior
or fix the innovation SD, so it estimates that SD and leaves the first M state
unpenalized. Because the direct random-walk fit was numerically unstable from
default starting values, the stock script first fits an IID-M approximation and uses
its estimates only to initialize the random-walk fit. The official assessment
reports terminal M of 0.81 for ages 5-8 and 0.85 for ages 9+.

The source models age-2 recruitment as lagged SSB times an AR1 recruitment
rate. The retained tinyAM age-3 entry follows a RW; it is not the same
stock–recruit relationship.

[Official source review](https://github.com/PaulRegular/saio/blob/main/source_reviews/dfo_cod_4t4vn.md).
