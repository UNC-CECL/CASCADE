# 2026-09-27 — the wave grid with the end domains fixed

Hannah, 2026-09-26/27: sweep the four wave parameters with the end
source/sink terms held **fixed**, rather than zero
(`../2026-09-25-wave-grid-smoothed-score/`) or re-solved per setting
(`../2026-09-26-wave-shortlist-ends-solved/`). The fixed pair per window comes
from `../../end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/`
(solved at the step-2 baseline waves, full management, against CoastSat, on
the metres offset) and is imposed in every run, natural and managed alike.

Everything else is the zeroBE grid's, unchanged: coarse grid Hs 0.75, 1, 1.5, 2 ×
Tp 7, 8, 10 × asymmetry 0.5, 0.7, 0.9 × high-angle 0.3, 0.45, 0.55, both
windows, both scenarios; refine around each best; cross-runs; the smoothed
score. The two studies compare one to one.

Layout as the zeroBE grid (`tables/`, `figures/`, `logs/`, `runs/<phase>_<scenario>/`).

Driver: `scripts/hatteras_ms/experiments/HAT_wave_grid_fixed_ends.py run all --jobs 8`

## Results

**1996–2010 coarse grid** (run 2026-09-27 01:35–03:58, ends fixed at
GIS 1 +7.02 / GIS 90 +9.83 m/yr): 216 cells, 198 scored, 18 drowned (every
Hs 0.75 / Tp 10 cell). 2010–2024 not run yet: its GIS 1 solve stalled at
~+133–138 m/yr (residual ~0.1), awaiting Hannah's decision. Refine and
cross-runs not run yet.

| | best with ends fixed | smoothed | bias (m/yr) | zeroBE rank / score |
|---|---|---|---|---|
| natural | Hs 1.5, Tp 10, asym 0.7, high-angle 0.45 | +28.5% | −0.46 | #1 / +23.4% |
| natural #2 | Hs 1.5, Tp 7, asym 0.7, high-angle 0.45 | +25.2% | −0.58 | #22 / +13.7% |
| managed | Hs 2.0, Tp 10, asym 0.7, high-angle 0.45 | +28.4% | −0.34 | #27 / +11.2% |
| managed #2 | Hs 1.5, Tp 10, asym 0.9, high-angle 0.45 | +28.0% | −0.14 | #7 / +17.1% |
| managed #5 | Hs 1.5, Tp 8, asym 0.7, high-angle 0.45 | +26.6% | −0.18 | #3 / +21.6% |

- Fixed ends raise every setting: +15 points on average (natural), +9
  (managed). Ranking correlation with zeroBE 0.93 natural, 0.83 managed.
- Natural keeps its zeroBE winner. Managed moves up the ridge to larger
  waves (Hs 1.5–2, Tp 8–10); the top five are within 2 points.
- High-angle 0.45 in every top setting; asymmetry 0.7 (or 0.9).
- The gain is at the ends (GIS 1–5 now near CoastSat's +2 to +3) and at
  Tri-Village; the accreting peaks at GIS 18, 29 and 42 are still missing.

## Final (sweep finished 2026-09-27 ~11:40)

2010–2024 ends accepted by Hannah at GIS 1 +137.581 / GIS 90 +14.731. Coarse
432 (36 drowned: every Hs 0.75 / Tp 10 cell, year 2 in 2010–2024), refine 118,
cross 14; no crashes. Figures `figures/best/grid/`; pairs in
`tables/one_parameter_between_periods.csv`.

| | best with ends fixed | smoothed | bias (m/yr) | zeroBE best |
|---|---|---|---|---|
| managed 1996–2010 | Hs 1.75, Tp 9, asym 0.7, high-angle 0.45 | **+29.1%** | −0.26 | Hs 1.5, Tp 10, asym 0.7, ha 0.45: +22.7% |
| natural 1996–2010 | Hs 1.5, Tp 10, asym 0.7, high-angle 0.45 | **+28.5%** | −0.46 | Hs 1.25, Tp 10, asym 0.8, ha 0.45: +23.9% |
| managed 2010–2024 | Hs 2, Tp 7.5, asym 0.6, high-angle 0.5 | −114.7% | −1.63 | −123.5% |
| natural 2010–2024 | Hs 2, Tp 10, asym 0.5, high-angle 0.5 | −514.7% | −3.24 | −530.3% |

**One setting for both windows** (lowest mean RMSE ÷ flat line):
- managed: **Hs 2, Tp 7.5, asym 0.6, high-angle 0.5** — 1996–2010 +24.4%
  (bias +0.07), 2010–2024 −114.7% (that window's best). With the ends fixed
  the shared pick no longer gives up 1996–2010 (it cost it down to +8% with
  zeroBE).
- natural: Hs 2, Tp 10, asym 0.5, high-angle 0.5 — 1996–2010 +16.9%.

**One parameter changing between windows** (1996–2010 kept within 2 points of its best):
- natural: **Hs 2, Tp 10, high-angle 0.5; asymmetry 0.6 → 0.5** — +26.6% and
  −514.7% (2010's best). The same pair the zeroBE grid picked.
- managed: Tp 10, asym 0.7, high-angle 0.45; **Hs 2 → 1.5** — +28.4% and
  −147.1%; or Hs 1.5, Tp 10, high-angle 0.45, asym 0.9 → 0.7 (+28.0 / −147.1).

- Fixed ends lift 1996–2010 by 5–6 points at the top and move managed up the
  ridge (Hs 1.75–2); high-angle 0.45–0.5 and asymmetry 0.6–0.7 throughout.
- 2010–2024 still cannot be fit; the ends change it by 10–15 points.
- GIS 1 in 2010–2024 was solved at the reference waves (Hs 1); at Hs 2 the
  +137.6 m/yr overshoots (model ≈ +10 against +6.9 observed, ~+300 m of
  position change at GIS 1).
- Still no setting makes the accreting peaks at GIS 18, 29, 42.

## Picked on the raw (unsmoothed) score — 2026-09-27, Hannah

The runner's RMSE and every matrix run score the per-domain model, so the
picks were redone on the raw score (`figures/best/grid_raw/`, dark model line,
no smoothed curve; `tables/best_shared_grid_raw.csv`). The raw score penalises
domain-scale roughness, which high-angle 0.5 damps, so the picks move from
high-angle 0.45 to 0.5 and the per-window and shared answers nearly coincide:

| | raw best | raw | smoothed |
|---|---|---|---|
| managed 1996–2010 | Hs 2, Tp 7.5, asym 0.6, high-angle 0.5 | +24.6% | +24.4% |
| natural 1996–2010 | Hs 1.75, Tp 10, asym 0.6, high-angle 0.5 | +25.8% | +25.3% |
| managed 2010–2024 | Hs 2, Tp 7, asym 0.5, high-angle 0.55 | −115% | |
| natural 2010–2024 | Hs 2, Tp 10, asym 0.5, high-angle 0.5 | −522% | |

- **Shared, managed: Hs 2, Tp 7.5, asym 0.6, high-angle 0.5** — it is also
  the managed 1996–2010 best; 2010–2024 −126%.
- **Shared, natural: Hs 2, Tp 10, asym 0.6, high-angle 0.5** — 1996–2010
  +24.5%; with one change, asymmetry 0.6 → 0.5, 2010–2024 takes its best.
- The smoothed-score winners (high-angle 0.45) score only +15.7% raw.

## Final — ends re-solved at the adopted waves, targeted reruns (2026-09-27 afternoon)

The ends were re-solved at **Hs 2 / Tp 7.5 / asym 0.6 / high-angle 0.5**
(full management): 1996–2010 **+4.84 / +17.55** (converged), 2010–2024
**+18.8 / +24.535** (GIS 1 probed, residual −0.045). All runs made on the Hs-1
ends went to `archive/2026-09-27-fixed-ends-{1996,2010}-ends-solved-at-hs1/`
(with this study's superseded figures and pick tables under
`study_outputs_before_resolve/`). Rather than the full grid, a **targeted** set
was rerun per window (`tables/targeted_<period>_cells.csv`: top 15 on the raw
score from the superseded runs and from the zeroBE grid, the other window's top
5, the adopted candidates): 1996–2010 51 settings, 2010–2024 40 (+8 coarse).
Figures `figures/best/grid_raw/` (raw score); pairs
`tables/one_parameter_between_periods_raw.csv`.

| | best (raw score) | raw | bias (m/yr) |
|---|---|---|---|
| managed 1996–2010 | Hs 2, Tp 9, asym 0.6, high-angle 0.5 | +20.6% | +0.06 |
| natural 1996–2010 | Hs 1.75, Tp 10, asym 0.6, high-angle 0.5 | +23.6% | −0.20 |
| managed 2010–2024 | Hs 2, Tp 7, asym 0.5, high-angle 0.55 | −115.8% | −1.52 |
| natural 2010–2024 | Hs 2, Tp 10, asym 0.5, high-angle 0.5 | −515.6% | −3.33 |

- **Managed, one setting: Hs 2, Tp 10, asym 0.6, high-angle 0.5** — +20.5% /
  −132.8% (Tp 7.5 is equivalent: +20.2% / −135.4%).
- **Natural, one setting: Hs 2, Tp 10, asym 0.6, high-angle 0.5** — +22.1% /
  −550.8%; **one change, asymmetry 0.6 → 0.5** in 2010–2024: −515.6%
  (that window's best).
- 1996–2010 scores are 3–4 points lower than on the Hs-1 ends (+7.0 / +9.8,
  against +4.8 / +17.5 now); the ranking is the same ridge. 2010–2024 unchanged within 1–2 points: the ends
  do not reach its interior.
