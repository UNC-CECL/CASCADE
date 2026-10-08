# 2017 Buxton fill shifted two domains north (2026-10-08)

**Question:** does placing the 2017 Buxton fill on GIS 8–18 instead of the reported GIS 6–16 improve the 2009–2025 test run?

**Setup:** `scripts/hatteras_ms/experiments/HAT_buxton_2017_shift_north.py` (run / score) and `_plot.py`. Test period, domainBE (set 1), blocking groin, full management, relocations off. The unchanged runner is used; only the 2017 project's `gis_domains` is swapped in the child process. The volume (2.6 M cy) and the width (11 domains) are unchanged, so the fill is the same 361.4 m³/m (checked in the shifted run's `nourishment_log.csv`: GIS 8–18). The 2022 Buxton fill stays on GIS 6–16. The reported arm matches the matrix run exactly (max difference 0 m).

## Result

| score | reported | shifted |
|---|---|---|
| interior GIS 2–89, LOWESS-7: bias / RMSE / r | −18.2 / 32.1 / 0.36 | −18.2 / 31.2 / 0.43 |
| Buxton reach GIS 1–20 raw RMSE | 52.7 | 51.0 |
| GIS 1–5 raw bias | −69.7 | −71.3 |
| GIS 6–7 raw bias (sand removed) | +47.7 | +25.2 |
| GIS 17–18 raw bias (sand added) | −14.0 | +15.1 |

- The model keeps the sand where it is placed: the change made by the shift is almost entirely in GIS 6–7 (−21 to −24 m) and GIS 17–18 (+26 to +33 m), with only a few metres leaking to GIS 4–5, 8, 16 and 19.
- GIS 6–7 halves its overshoot but is still +25 m (the observations show −40 m there). GIS 17–18 goes from 14 m too low to 15 m too high, so that end swaps one misfit for the other.
- GIS 1–5, where the observed sand went (south past the groin), gets slightly worse.
- The island-wide gain (RMSE −0.9 m, r +0.07) comes from removing sand at GIS 6–7.

**Reading:** the shift helps by taking sand off GIS 6–7, not by putting it where it was observed. It mimics part of the south bypass without the sand reaching GIS 1–5, and it has no support in the placement record (GIS 6–16, from the Haulover point and the groin). Experiment only: the config still has GIS 6–16 and nothing in the main code changed; adoption is Hannah's call. Related: `../2026-10-06-fill-footprint-coastsat/` (observed footprints, worse), and the known mismatch in `data/hatteras_init/4-mgmt-forcing/README.md`.

Figure: `figures/buxton_2017_reported_vs_shifted_north_net_change_2009_2025.png`, in the run net-change layout (shifted run solid, reported run dashed, both footprints as bars; caption in CAPTIONS.md). Tables: `tables/scores.csv`, `tables/per_domain_net_change.csv`.
