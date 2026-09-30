# 2026-09-26 — wave shortlist with the end domains solved per setting

> **Record.** The settings to use: [`../2026-09-27-wave-recommendation/README.md`](../2026-09-27-wave-recommendation/README.md).

Hannah, 2026-09-26: "when you solve for the ends, are there different wave
parameters that perform the best?" The 2026-09-25 grid
(`../2026-09-25-wave-grid-smoothed-score/`) ran zeroBE: nothing imposed at the
two end domains. The stored edgeBE values (`HATTERAS_BE_EDGE_ONLY`) were solved
at Hs 2.5 on the old ÷10 offset and do not carry over; what the ends must carry
depends on the waves, so each setting gets its own solve.

| | |
|---|---|
| shortlist | the top 10 zeroBE settings (smoothed score) per window × scenario: 40 chains (`tables/shortlist.csv`) |
| ends | solved against each window's CoastSat LRR, as the matrix end values were: GIS 1 against the raw domain mean, GIS 90 against the LOWESS-10 value |
| solve | step 0 = the setting's zeroBE grid run; each end stepped on its own, first by the 2026-09-11 response (≈0.09 m/yr residual per m/yr imposed at GIS 1, 0.13 at GIS 90), then by the secant; all chains in lockstep; converged at \|residual\| ≤ 0.02 m/yr at both ends, at most 5 probes |
| score | as the grid: share of the alongshore variation explained by the model smoothed like CoastSat, interior GIS 2–89; the zeroBE rank and score beside it (`tables/all_runs.csv`) |
| fixed | metres offset (dune line v1), edgeBE through `HAT_BE_OVERRIDE`, no groin, no relocations, Barrier3D route_overwash fix |

```
tables/shortlist.csv     the 40 settings and their zeroBE scores
tables/solve_log.csv     every probe: imposed ends and residuals
tables/chains.json       each chain's final state
tables/all_runs.csv      converged runs scored, zeroBE rank vs ends-solved rank
logs/<scenario>/step<k>/<period>_<settings>.log
runs/<scenario>_step<k>/<window>/edgeBE/<run_name>/     on disk only
```

Reproduce: `python scripts/hatteras_ms/experiments/HAT_wave_shortlist_ends_solved.py run`

## Results

Run 2026-09-26 17:20 – 2026-09-27 01:10: a first pass (plain secant, 5 probes)
converged 6 of 40 chains -- GIS 1 in 2010-2024 does not respond smoothly
(0-25 m/yr gives residuals near zero, above ~25 gives +5 to +15 whatever the
value), the secant jumped to -466 / +335 m/yr and four probes drowned the
barrier -- then `resume` (safeguarded, steps 6-9). Final: 1996-2010 16 of 20
converged, the rest within 0.08 m/yr at GIS 1 and 0.26 at GIS 90; 2010-2024 2 of
20 converged, median residual ~0.1 m/yr, a few off by 6-7 at GIS 1.

| | best with ends solved | smoothed | ends GIS 1 / 90 (m/yr) | zeroBE best |
|---|---|---|---|---|
| natural 1996–2010 | Hs 1.25, Tp 9, asym 0.7, high-angle 0.45 | +27% | +8.6 / +16.2 | Hs 1.25, Tp 10, asym 0.8: +24% |
| managed 1996–2010 | Hs 1.5, Tp 8–10, asym 0.7, high-angle 0.45 | +21–22% | +9 / +20–22 | the same waves, +23% |
| managed 2010–2024 | Hs 2, Tp 7, asym 0.5, high-angle 0.55 | −118% | +39 / +7 | −124% |
| natural 2010–2024 | Hs 2, Tp 10, asym 0.5, high-angle 0.5 | −511% | +11 / +27 | −530% |

- The best waves barely move: Hs 1.25–1.5, Tp 8–10, **asymmetry 0.7,
  high-angle 0.45** in both 1996–2010 scenarios; the asymmetry-0.8 settings
  lose 3–4 points once their ends are solved.
- Solved ends lift natural (+24% → +27%), not managed.
- The ends need +6 to +9 m/yr at GIS 1 and +13 to +30 at GIS 90 in 1996–2010,
  against the stored +32 / +10 (solved at Hs 2.5 on the ÷10 offset).
- 2010–2024 improves by 5–20 points at most and stays far below a flat line.
- A preliminary reading (unconverged chains) had put natural Hs 2 / Tp 10 /
  asym 0.6 / high-angle 0.5 first at +27%; converged it is +24%, rank 5.
