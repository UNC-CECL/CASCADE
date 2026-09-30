# 2026-09-28 — option A ends re-solved against the LOWESS-7 target

Hannah, 2026-09-28: "switch the runner to 7 and re-solve the ends". The runner's CoastSat scoring target went from a 10-domain to a 7-domain LOWESS (the group's smoothing range), and the end rates solved against the LOWESS-10 target no longer matched it.

| | |
|---|---|
| waves | option A: Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5 |
| scenario | full management, edgeBE, no groin, relocations off (as the 09-27 solve) |
| target | each window's CoastSat LRR: GIS 1 against the raw domain mean (unaffected by the window), GIS 90 against the LOWESS-7 value |
| solve | `HAT_resolve_ends_metres.py --tag … --seed "1996=4.8394,17.545;2010=18.8,24.535"`: step 1 probes the LOWESS-10 ends, then the safeguarded secant; converged at \|residual\| ≤ 0.02 m/yr |

## Result — adopted in `HATTERAS_BE_EDGE_ONLY` the same day

| window | GIS 1 | GIS 90 | residuals | steps | LOWESS-10 values |
|---|---|---|---|---|---|
| 1996–2010 | +4.8394 | +18.2545 | +0.006 / +0.016 | 5 | +4.8394 / +17.545 |
| 2010–2024 | +18.8657 | +24.2358 | −0.005 / +0.001 | 3 | +18.8 / +24.535 |

GIS 90's target moved +0.125 m/yr (1996) and −0.066 m/yr (2010); with a response of about 0.2 m/yr per m/yr imposed, GIS 90 moved +0.71 and −0.30 m/yr. The edgeBE matrix was re-run on these values the same day (the LOWESS-10 runs are in `raw_runs/archive/2026-09-28-loess10-ends/`).

`tables/ends.json` is the record; `tables/solve_log_*.csv` every probe; `runs/` on disk only.
