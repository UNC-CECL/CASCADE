# 2026-09-28 — Hs 2.5 against Hs 2.0 on the adopted setup, each on its own ends

Hannah, 2026-09-28: "set up the Hs 2.5 check". The wave sensitivity sweep on the
adopted setup (`raw_runs/sensitivity/`) left option A best on high-angle,
asymmetry and Tp, but Hs kept improving above 2.0 (interior RMSE 1996
1.168 → 1.139 at 2.5, 2010 2.123 → 2.084). Those cells ran on ends solved at
Hs 2.0. On the pre-ceiling setup the same trend was the ends mismatch: with
ends re-solved at Hs 2.5, 1996 managed dropped to +18.1% from +20.2%
(`../2026-09-27-wave-recommendation/`). The adopted setup changed the dunes and
storms, which is what Hs acts on, so the check is repeated here.

| | |
|---|---|
| question | does Hs 2.5 beat Hs 2.0 once each has its own end rates, on the adopted setup? |
| waves | Hs 2.5, Tp 7.5, asymmetry 0.6, high-angle 0.5 (option A with Hs 2.5) |
| setup | Barrier3D `hatteras/adopted`, storms `v3_trim24`, the dune-cap fix (`bdm_dune_cap_applies_to` = added sand only; every run here has it), full management, edgeBE, no groin |
| ends | solved by `HAT_resolve_ends_metres.py` exactly as `../../end-domain-boundaries/2026-09-28-ends-resolved-adopted/` (GIS 1 against the raw domain mean, GIS 90 against LOWESS 7; \|residual\| ≤ 0.02 m/yr) |
| seed | 1996 +4.0 / +33.0, 2010 +8.0 / +38.5 m/yr: the adopted Hs-2 ends moved by what Hs 2.5 needed on 09-27 |
| compared against | the matrix full-management runs, Hs 2.0 on its own ends (1996 +4.3509 / +19.0935, 2010 +8.0 / +21.2582 since the dune-cap fix), re-run with the fix at 22:18 and 22:56. The sweep's Hs-2.5 cells were meant as a third row but predate the fix, so they are left out |
| scored on | interior (GIS 2–89) RMSE and mean bias, and the raw share of alongshore variation explained, as in the matrix |

The config (`HATTERAS_BE_EDGE_ONLY`) is not changed. Runs under
`runs/full_management_step<k>/`, the solve in `tables/`, the log
`logs/solve.log`.

Reproduce:

```
python scripts/hatteras_ms/experiments/HAT_resolve_ends_metres.py \
  --tag wave-climate/2026-09-28-hs2p5-check-adopted-setup \
  --hs 2.5 --tp 7.5 --asym 0.6 --ahf 0.5 \
  --seed "1996=4.0,33.0;2010=8.0,38.5"
```

## Result (2026-09-28, 23:00): Hs 2.0 stays

Ends solved at Hs 2.5 (`tables/ends.json`): 1996 +4.0 / +33.1847 (closest probe,
residuals −0.019 / −0.021; the GIS 90 response is noisy at ±0.05 m/yr and the
last four probes repeated), 2010 +6.2249 / +34.65 (converged).

Interior GIS 2–89 against the LOWESS-7 CoastSat target (`tables/comparison.csv`):

| window | run | ends GIS 1 / 90 | RMSE | bias | share explained, raw | smoothed |
|---|---|---|---|---|---|---|
| 1996–2010 | Hs 2.0 (matrix) | +4.351 / +19.093 | 1.172 | +0.06 | +22.5% | +21.1% |
| 1996–2010 | Hs 2.5 | +4.000 / +33.185 | 1.174 | +0.09 | +22.3% | +23.3% |
| 2010–2024 | Hs 2.0 (matrix) | +8.000 / +21.258 | 2.068 | −1.36 | −76.9% | −69.3% |
| 2010–2024 | Hs 2.5 | +6.225 / +34.650 | 2.053 | −1.32 | −74.4% | −60.1% |

- On the raw score (the one the matrix and the 09-27 picks use) the two are
  tied: 1996 −0.3 points for Hs 2.5, 2010 +2.5 points, RMSE within 0.015 m/yr.
- The smoothed score leans to Hs 2.5 (1996 +2 points, 2010 +9), but 2010 is far
  below a flat line either way, and the bias does not move.
- Hs 2.5 needs about +14 m/yr more imposed at GIS 90 in both windows for that
  tie. The sweep's gain above Hs 2.0 was the ends mismatch, as on 09-27.
- Both 2010 runs drown NC-12 once (GIS 14, as in the whole adopted sweep).

**Option A (Hs 2.0) stays.** Nothing here justifies the change; Hs 2.0–2.5 is
a flat optimum on the adopted setup too.
