# 2026-09-29 — end rates re-solved on the split12 storms

Hannah, 2026-09-29: "use split12 as the storm series going forward", then "re-run the matrix and re-solve the ends". The storms changed from `v3_trim24` to `v3_split12_trim24`: grouped events are split at ≥12 h below the berm, so Fran 1996 and Jose 2017 are back. Evidence: `experiments/storms-and-overwash/2026-09-29-event-splitting/`.

The protocol is as `2026-09-28-ends-resolved-dunecap/`: option A waves, full management, edgeBE, no groin, CoastSat LRR target (GIS 1 raw, GIS 90 LOWESS-7). The secant was seeded at the config's ends. The driver was:

    HAT_resolve_ends_metres.py --periods 1996 2010 --hs 2.0 --tp 7.5 --asym 0.6 --ahf 0.5
        --accept 0.03 --tag end-domain-boundaries/2026-09-29-ends-resolved-split12
        --seed "1996=4.3509,19.0935;2010=8.0,21.2582"

| window | step | GIS 1 | GIS 90 | residuals GIS 1 / GIS 90 |
|---|---|---|---|---|
| 1996–2010 | 1 (the trim24 ends) | +4.3509 | +19.0935 | −0.033 / −0.008 |
| | **2** | **+4.3888** | +19.0935 | **+0.005 / −0.008** |
| 2010–2024 | 1 (the trim24 ends) | +8.0 | +21.2582 | −0.040 / −0.018 |
| | **2** | **+8.0405** | +21.2582 | **+0.003 / −0.018** |

**Adopted in `HATTERAS_BE_EDGE_ONLY`: 1996 (+4.3888, +19.0935), 2010 (+8.0405, +21.2582).**

- Only GIS 1 moved, by +0.04 m/yr in each window.
- The extra storm hours added a little erosion at the southern end, which a slightly larger imposed rate offsets.
- GIS 90 was already inside tolerance.
- The 2010 GIS 1 response is monotonic around +8 (the 09-28 direct probes: +8.0 → +0.002, +8.8 → +0.27), so the secant step is on safe ground. The flat region is above about +10.

`tables/ends.json` is the record; `tables/solve_log_*.csv` holds the probes and `runs/` holds the runs.
