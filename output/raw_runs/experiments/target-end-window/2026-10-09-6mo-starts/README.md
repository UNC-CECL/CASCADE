# 2026-10-09-6mo-starts

**Question.** Does a ±6-month start shoreline (the model's island offset and the target's start) change the runs and their scores, compared with the ±1-yr one in use?

**Why.** Hannah's definition of the start is ±6 months around the DEM date, but the config has used ±1 yr since 2026-09-29 (offset v2).

**Method.** `six_month_starts.py`, the unchanged runner, 4 at a time. Two arms × three presets (zeroBE and edgeBE no groin, domainBE + blocking groin) × both periods:

- **matrix**: the production runs in `matrix/1996_2009/` and `matrix/2009_2025/`, on offset v2 (±1 yr: 1995-10-12 → 1997-10-12, 2008-08-17 → 2010-08-17).
- **6mo**: the same runs on offset v3 (±6 months: 1996-04-12 → 1997-04-12, 2009-02-17 → 2010-02-17), pinned with `HAT_OFFSET_VERSION_<year>_SHORELINE=v3`, filed here under `1996_2009/` and `2009_2025/`.

Each arm is scored against the target that starts from its own line; the end windows are the config's (calibration 2007-07-01 → 2009-06-30, test 2025-02-17 → 2026-02-17). Interior GIS 2–89, both sides LOWESS-7. Runs use the config as of 2026-10-09: the test runs 17 years (2009–2025).

**The ±6-month lines are thin.** At the usual minimum of 10 positions per transect they do not cover the island (1996: 81 domains, GIS 8–14, 76 and 81 missing; 2009: 69 domains, GIS 1–22 mostly missing). They were built with a minimum of 5 (905 of 906 transects each, median 10–11 positions against 20–28 at ±1 yr). See `data/hatteras_init/2-brie-offset/<year>/shoreline/v3/PROVENANCE.md`.

## Results

**The start lines.** In the model frame v3 sits +1.3 m (1996; −3.7 to +7.6) and +3.0 m (2009; −3.3 to +9.9) seaward of v2 on average.

**What it does to the model** (net change, 6mo minus matrix, interior LOWESS-7): mean within ±0.4 m, RMS 0.8–1.6 m for every run. The model forgets most of its start-line difference.

**What it does to the target** (6mo minus matrix): calibration +1.0 m (RMS 2.1), test −0.4 m (RMS 3.2).

**Scores** (`scores.csv`, bias / RMSE / r)

| run | calibration, ±1 yr | calibration, ±6 mo | test, ±1 yr | test, ±6 mo |
|---|---|---|---|---|
| zeroBE, no groin | +1.5 / 16.6 / 0.56 | +0.6 / 16.7 / 0.55 | −13.4 / 28.3 / 0.27 | −12.7 / 28.4 / 0.30 |
| edgeBE, no groin | +2.4 / 17.0 / 0.53 | +1.3 / 16.9 / 0.54 | −11.6 / 25.1 / 0.42 | −11.0 / 25.5 / 0.44 |
| domainBE + groin | −3.3 / 5.6 / 0.98 | −4.4 / 6.8 / 0.97 | −20.2 / 34.0 / 0.36 | −19.8 / 33.8 / 0.39 |

domainBE (BE set 1) was solved on the old calibration target (±1-yr start, DEM-centred end), so its calibration scores on either line are not evidence.

**Reading.** The start window moves bias by about 1 m and RMSE by less than 0.5 m (1.2 m for the fitted domainBE). That is small against the misfits (RMSE 17 m calibration, 25–34 m test), so the choice between ±6 months and ±1 yr is a definitional one, not one the fit can decide. The ±1-yr line averages twice the positions, so it is the less noisy of the two.

**Status.** Diagnostic. The config still uses ±1 yr (v2 CURRENT). Adopting ±6 months means setting `CURRENT` to v3 for both starts, changing `SHORELINE_WINDOW_FOR_YEAR` and the start windows in `NET_CHANGE_WINDOWS`, and re-running the matrix.

## Also checked

- **No drift.** The calibration matrix re-run (inputs unchanged) is bit-identical to the runs archived in `archive/2026-10-09-pre-test-2025/`. The 17-year test runs match the archived 16-year ones exactly through 1 Jan 2025.
- **What 2025 adds.** Over the extra year the test runs lose about 1.3 m island-wide (edgeBE model mean +1.32 → +0.06 m), so the test bias grows by about 1.2 m (edgeBE −10.4 → −11.6 m).

**Note.** `score` was re-run after the edgeBE matrix runs moved to the 2026-10-09 end rates; the 6mo edgeBE runs keep the old ends. The matrix edgeBE scores moved by 0.01 m, so the table stands.

Outputs: `scores.csv`, `domain_net_change.csv` (per GIS: each arm's target and each run, LOWESS-7), `logs/` (one runner log per run).

Back to [the topic](../README.md).
