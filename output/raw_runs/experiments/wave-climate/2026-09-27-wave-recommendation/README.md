# 2026-09-27 — recommended hindcast wave climate

> **Adopted 2026-09-27: option A** is now the model default (waves in `HAT_hindcast_config` / `hat_run.yaml`, ends in `hatteras_site_config.HATTERAS_BE_EDGE_ONLY`). Option B is recorded as `HATTERAS_WAVE_OPTION_B` in the site config, not wired.

Synthesis of every wave test of 24–27 September (no new runs except the
recommended setting under the natural scenario). Write-up with the figures:
https://claude.ai/artifact/L3ezkgxG7aQNDDPMYazDjG (private).

**Recommended, both windows, both scenarios: Hs 2.0 m, Tp 7.5 s, asymmetry 0.6,
high-angle 0.5**, with the end rates solved for it (full management, against
CoastSat): 1996–2010 GIS 1 +4.84 / GIS 90 +17.55, 2010–2024 +18.8 / +24.535 m/yr
(`../../end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/tables/ends.json`).
Raw share explained 1996–2010: managed +20% (bias +0.09), natural +24% (−0.18);
2010–2024: managed −135%, natural −699% (no setting fits that window).

Figures (`figures/`, drawn by `scripts/hatteras_ms/experiments/HAT_wave_recommendation_figures.py`):
1 score by parameter (zeroBE coarse grid), 2 agreement across the five searches,
3 recommended vs CoastSat, 4 what the high-angle fraction does, 5 the 2010 windows.

## Hs above the grid edge (27 Sep, Hannah: "run the Hs 2.25 and 2.5 check")

`hs_check_2p25_2p5.csv`, `hs_check_2_to_3.csv`: Hs 2.25–3.0 at Tp 7.5 / asym 0.6 /
high-angle 0.5, ends fixed at the Hs-2 values: managed 1996–2010 raw 20.2 → 25.1%
at 2.5–2.75, 24.1% at 3.0; natural ±3 points of noise. The ends mismatch grows with
Hs (GIS 90 −2.5 m/yr at 2.5). Ends re-solved at Hs 2.5 (1996 +4.379 / +31.404,
2010 +8.0 / +40.399, `ends.json`), then Hs 2.25/2.5/2.75 on them
(`hs_check_on_hs2p5_ends.csv`): managed 1996 **+18.1%** at 2.5 (below Hs 2.0's
+20.2% on its own ends), natural +20.4% (vs +24.5%). The apparent gain came from
under-supplying GIS 90. **Hs 2.0 kept**; 2.0–2.5 is a flat optimum. `ends.json` was
set back to the Hs-2 values (1996 +4.8394 / +17.545, 2010 +18.8 / +24.535); the
Hs-2.5 values are under `history`.

## Same vs period-specific (27 Sep; report version 3 at the link above)

| option | 1996–2010 | 2010–2024 | ends (GIS 1 / 90) | 1996 M · N | 2010 M · N |
|---|---|---|---|---|---|
| A same (default) | Hs 2.0, Tp 7.5, asym 0.6, ha 0.5 | same | 1996 +4.8/+17.5, 2010 +18.8/+24.5 | +20 · +24% | −135 · −699% |
| B one change | same | **Hs 2.5** | 2010 +8.0/+40.4 | unchanged | −122 · −538% |
| C best per period | M Tp 9 / N as A | M 2.0,7,0.5,0.55 / N 2.0,10,0.5,0.5 | Hs-2 ends | +21 · +24% | −116 · −516% |

Every one-parameter change for 2010–2024 was run (`oneparam_*` in the fixed-ends
study): only Hs helps both scenarios; asymmetry and high-angle make it worse.
Figure 6 (`figures/6_hs_change_2010.png`) compares A and B.
