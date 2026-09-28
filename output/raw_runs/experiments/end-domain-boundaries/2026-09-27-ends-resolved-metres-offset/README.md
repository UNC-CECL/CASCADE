# 2026-09-27 — the end domains re-solved under the metres offset

> **Current values (2026-09-27):** 1996–2010 GIS 1 **+4.8394** / GIS 90 **+17.545**; 2010–2024 **+18.8** / **+24.535** m/yr, solved at Hs 2 / Tp 7.5 / asym 0.6 / high-angle 0.5 (`tables/ends.json`, key `ends_m_yr`). The first Results section below is the Hs-1 solve, since replaced; the Hs-2 solve and the Hs-2.5 check follow it. For Hs 2.5 in 2010–2024 (option B) use +8.0 / +40.399.

Hannah, 2026-09-26/27: the wave sweep should hold the end source/sink terms
**fixed**, and the ends need re-solving now that the island offset is in
metres. The stored values (`HATTERAS_BE_EDGE_ONLY`: 1996 +32.2 / +10.0,
2010 +72.6 / +31.3 m/yr at GIS 1 / GIS 90) were solved at Hs 2.5 on the ÷10
offset.

| | |
|---|---|
| reference waves | the step-2 baseline, Hs 1.0 m, Tp 8 s, asymmetry 0.8, high-angle 0.45: neutral, not the winner of either search, so the fixed ends do not pre-favour the sweep |
| scenario | full management (as the matrix end values always were); the pair is used for natural and managed runs alike |
| target | each window's CoastSat LRR: GIS 1 against the raw domain mean, GIS 90 against the LOESS-10 value |
| solve | step 0 a fresh zeroBE run, then safeguarded Newton steps (first from the metres response measured 2026-09-26, then secant capped at ±30 m/yr until bracketed, then interpolation inside the bracket); converged at \|residual\| ≤ 0.02 m/yr |
| output | `tables/ends.json` (read by `HAT_wave_grid_fixed_ends.py`), `tables/ends.csv`, `tables/solve_log.csv` |

The config was not changed by this study. On 2026-09-27 the Hs-2 values (option A) were written into `HATTERAS_BE_EDGE_ONLY`, and the option B pair into `HATTERAS_WAVE_OPTION_B`.

Driver: `scripts/hatteras_ms/experiments/HAT_resolve_ends_metres.py`

## Results

Run 2026-09-27 01:12–01:32 (`tables/ends.csv`, `tables/solve_log.csv`).

| window | GIS 1 (m/yr) | GIS 90 (m/yr) | residuals GIS 1 / 90 | state | stored (÷10 era) |
|---|---|---|---|---|---|
| 1996–2010 | **+7.02** | **+9.83** | −0.017 / +0.017 | converged, 4 steps | +32.2 / +10.0 |
| 2010–2024 | +132.9 to +137.6 | **+14.73** | +0.07 to +0.10 / +0.013 | GIS 90 converged; GIS 1 stalled at ~0.1 m/yr | +72.6 / +31.3 |

- 1996–2010: GIS 1 needs a fifth of the stored value; GIS 90 about the same.
- 2010–2024, GIS 1: the residual goes −0.27 (112.5), −0.17 (121.1), +0.07
  (137.6), +0.10 (132.9): at this level the end's response is not monotonic
  to better than ~0.1 m/yr, so the bracket stopped shrinking. It is solved
  to ~0.1 m/yr, not to the 0.02 criterion, and the value is very large: the
  2010–2024 CoastSat rate at GIS 1 is +6.9 m/yr (strong accretion at the
  Cape Point end), which at Hs 1 takes ~+135 m/yr of imposed supply.
- **Accepted (Hannah, 2026-09-27): 2010–2024 GIS 1 = +137.581 m/yr**, the
  closest probe (step 6, residual +0.069), with GIS 90 = +14.731. Written to
  `tables/ends.json` with the note; the fixed-ends sweep uses it.

## 2010–2024 re-solved at Hs 2 (2026-09-27, Hannah: "re-solve the 2010 ends at Hs 2 and rerun")

At the adopted managed waves, **Hs 2 / Tp 7.5 / asym 0.6 / high-angle 0.5**,
full management: step 0 residuals GIS 1 −9.61 / GIS 90 −4.43 m/yr. GIS 90
converged at **+24.535** (residual −0.008). The safeguarded solve stalled at
GIS 1 +26.26 (residual +0.151: probes from +23.6 to +40 all leave +0.15 to
+0.5, a noisy flat region), so GIS 1 was probed directly below it:

| GIS 1 | +8 | +14 | +18.8 | +19.3 | +20 | +23.6 | +26.3 |
|---|---|---|---|---|---|---|---|
| residual | −1.192 | −0.313 | **−0.045** | +0.099 | +0.062 | +0.212 | +0.151 |

**Taken: GIS 1 +18.8 m/yr** (residual −0.045; the response is not monotonic
within ~0.1 m/yr). The Hs-1 value it replaces (+137.581) overshot by about
7× at Hs 2: larger waves carry more sand alongshore on their own. `ends.json`
keeps the replaced value under `history`; the 2010–2024 fixed-ends runs made
on it are in `archive/2026-09-27-fixed-ends-2010-ends-solved-at-hs1/`.
