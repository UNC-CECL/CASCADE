# 2026-09-28 — end rates re-solved on the adopted model

Hannah, 2026-09-28: "include the overwash fixes, keep option A, go ahead". The model changed under the ends:

- **Barrier3D** `hatteras/adopted`: the three overwash fixes and per-cell dune ceilings.
- **Storms** `v3_trim24`.

So GIS 1 and GIS 90 were re-solved on it. Protocol as `2026-09-28-ends-resolved-lowess7/`:

| | |
|---|---|
| waves | option A |
| scenario | full management, edgeBE, no groin |
| target | CoastSat LRR (GIS 1 the raw domain mean, GIS 90 the LOWESS-7 value) |
| driver | `HAT_resolve_ends_metres.py --seed` from the LOWESS-7 ends |

## Adopted in `HATTERAS_BE_EDGE_ONLY`

| window | GIS 1 | GIS 90 | residuals | how |
|---|---|---|---|---|
| 1996–2010 | +4.3509 | +19.0935 | −0.006 / +0.007 | secant, 4 steps |
| 2010–2024 | +8.0 | +22.4937 | +0.002 / −0.004 | GIS 90 by the secant; GIS 1 by direct probes |

**2010 GIS 1.** The secant stalled at about +12 to +13.5: the response above about +10 m/yr is flat and even reverses (+13.5 to +19 all leave about +0.8 m/yr). Direct probes at GIS 90 = +22.4937 mapped it (runs `runs/direct_gis1_<rate>/`):

| GIS 1 imposed | −6 | −2 | +2 | +6 | **+8.0** | +8.8 | +9.5 | +10 |
|---|---|---|---|---|---|---|---|---|
| residual | −13.02 | −9.81 | −5.82 | −1.27 | **+0.002** | +0.27 | +0.50 | +0.52 |

The last LOWESS-7 value was +18.8657. A drop to +8.0 means the adopted model needs far less sand imposed at the southern end.

`tables/ends.json` is the record (the 2010 GIS 1 entry is the direct-probe value); `tables/solve_log_*.csv` holds the secant probes.
