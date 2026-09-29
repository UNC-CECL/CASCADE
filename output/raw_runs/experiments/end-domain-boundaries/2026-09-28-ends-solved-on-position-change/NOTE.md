# Ends solved on position change (2026-09-28)

Hannah, 2026-09-28: solve the edge source/sink on the observed end-minus-start
position change instead of the LRR (option 1: an experiment; solver and config
unchanged). GIS 90 against the LOESS-7 value (Hannah), GIS 1 raw.

Scripts: `scripts/hatteras_ms/experiments/HAT_resolve_ends_on_position_change.py`
(solve), `HAT_position_change_ends_figure.py` (natural runs at the ends + figure).

| window | GIS 1 | GIS 90 | LRR-solved (LOESS 7) | residual (m/yr; x14 = m) |
|---|---|---|---|---|
| 1996-2010 | +1.231 | +16.708 | +4.839 / +18.255 | -0.001 / -0.001, converged in 4 |
| 2010-2024 | +22.332 | +18.572 | +18.866 / +24.236 | +0.050 / +0.006, not converged |

2010-2024 GIS 1 is not monotonic between +21.5 and +22.3 (residual +-0.1 m/yr);
steps 6-8 re-probed +21.933 (deterministic, -0.088). Closest probe kept.

Interior (GIS 2-89) barely moves: managed LRR RMSE 1.19 -> 1.20 (1996),
2.31 -> 2.31 (2010); position RMSE 20.8 -> 20.9 m, 27.9 -> 27.8 m. The change
is at the ends: 1996 GIS 1 model LRR drops below the target LRR (+3.2) while
its net change now matches the observed +11.5 m.

Figure: `figures/recommended_vs_coastsat_position_change_ends.png`; the LRR-ends
version is `../../wave-climate/2026-09-27-wave-recommendation/figures/3_recommended_vs_coastsat.png`.
