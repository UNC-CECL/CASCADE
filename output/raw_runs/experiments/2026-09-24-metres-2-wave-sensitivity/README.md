# Wave-climate sensitivity in metres, natural and full management (2026-09-24)

Step 2 of 3 in the 2026-09-24 chain: `../2026-09-24-metres-INDEX.md`. Step 1
chose the offset unit (`../2026-09-24-metres-1-offset-units/`); this study
found the Barrier3D bug that step 3 fixed (`../2026-09-24-metres-3-barrier3d-overwash-fix/`).
Filed as `sensitivity/2026-09-24-natural-waves/` until the 2026-09-24
reorganisation (runs retagged; the name predates the full-management sweep).

## The question

With the island offset in metres (the default since 2026-09-24), how does the
model's alongshore shoreline change respond to each wave parameter when
nothing human acts on the island, in both canonical windows?

The design was set by interview with Hannah on 2026-09-24; every choice below
is hers.

## Design

| setting | value |
|---|---|
| scenario | `natural`: no road management, no beach/dune management, no fills, no relocations |
| groin | none |
| background erosion | `zeroBE` (no imposed rates; the edgeBE end rates were solved under ÷10) |
| island offset | metres, dune line, `2-brie-offset/<start>/duneline/v1` |
| periods | 1996–2010 and 2010–2024, each against its own CoastSat LRR target |
| baseline | Hs 1.0 m, Tp 8 s, asymmetry 0.8, high-angle fraction 0.45 |

The baseline is the best metres point of
`experiments/2026-09-24-metres-1-offset-units/` (found there under
full management, 1996–2010).

**Stage 1, one parameter at a time around the baseline** (baseline value in bold):

| parameter | folder | values |
|---|---|---|
| significant wave height, Hs (m) | `wave_height/` | 0.65, 0.75, **1.0**, 1.25, 1.5, 2.0, 2.5, 3.0 |
| high-angle fraction | `high_angle/` | 0.1, 0.2, 0.3, 0.4, **0.45**, 0.5, 0.55 |
| asymmetry | `asymmetry/` | 0.3, 0.4, 0.5, 0.6, 0.7, **0.8**, 0.9 |
| wave period, Tp (s) | `wave_period/` | 6, 7, **8**, 10, 12 |

The baseline itself is run once per period (`baseline/`), and once more under
full management (`baseline_full_management/`) for the natural-vs-managed
comparison: 25 runs per period, 50 in all. Asymmetry below 0.5 reverses the
direction of net drift; high-angle 0.55 steps past the ~0.5 limit where most
of the coast turns anti-diffusive, to map where the model breaks. Hs below
0.65 drowns the barrier (the offset-scale study).

**Stage 2, a grid**, chosen automatically from stage 1: the two parameters
whose range moves the share of alongshore variation explained the most,
averaged over both periods, with a drowned run counted as the worst score of
its period. Each gets the 5 stage-1 values centred on its best value (shifted
inward at the ends of its range); the other two stay at the baseline. 5 × 5
per period; cells stage 1 already ran are reused, not re-run.
`tables/stage2_selection.csv` records every parameter's effect, the pair
chosen and the grid values.

## Labels

| where | example |
|---|---|
| folder | `high_angle/2010_2024/zeroBE/<run_name>/` (parameter, period, preset) |
| log | `logs/high_angle/2010_2024/Hs1.0_period8.0_asymmetry0.8_highangle0.3.log`: every setting |
| tables | `tables/all_runs.csv`: every setting and score as a column, drowned runs with the reason |
| run name | `HAT_2010_2024_zeroBE_offsetmetres_noroad_nobdm_nogroin_waveHs1asym0p8ahf0p3` |

The run name is the runner's and leaves out whatever is at the model's
default: Tp 8 has no `Tp` token, and the baseline's Hs 1.0, asymmetry 0.8 and
high-angle 0.45 always appear because they are off the model's defaults. Use
the log name or the tables for a complete label. `noroad_nobdm` is the natural
scenario; a full-management run reads `road_bdm` there instead.

## Layout

```
README.md
tables/        all_runs.csv, observed_targets.csv, stage2_selection.csv
figures/
  stage1/      scores_by_<parameter>_both_periods.png
  alongshore/<period>/rate_and_position_change_by_<parameter>_<period>.png
  management/  natural_vs_full_management_baseline_both_periods.png
  stage2/      grid_<p1>_x_<p2>_both_periods.png
  (each with supporting/: PDF, CAPTIONS.md, the data CSV)
logs/<group>/<period>/<settings>.log, logs/drivers/
baseline/ baseline_full_management/ wave_height/ high_angle/ asymmetry/
wave_period/ grid_<p1>_x_<p2>/          <period>/zeroBE/<run_name>/
```

**What is in git.** The run folders are on disk only: this scenario's run name
is about 85 characters and appears twice in each file path, so every run path
is past Windows' 260 and git cannot index it (the same choice as the
offset-scale study). Committed: this README, `tables/` and each figure's PDF,
caption and data CSV. PNGs are ignored under `output/raw_runs/` repo-wide.

## Scores

The same as the offset-scale study: interior (GIS 2–89) mean bias and RMSE of
the modelled LRR rate against the window's CoastSat LRR target (LOESS 10
domains), and the share of the observed alongshore variation explained,
1 − Σ(model − obs)² / Σ(obs − obs mean)², with 0% meaning no better than a
flat line at the observed mean (`tables/observed_targets.csv` gives each
window's flat-line RMSE). The target is rebuilt the way the runner builds it
and every run's RMSE is checked against the runner's. Position change is the
modelled end-minus-start position against the observed CoastSat change (the
mean position over the last calendar year minus that over the first,
`5-scr/3-rates/coastsat/total_change/<window>/`, smoothed at 10 domains).

## What happened while it ran

- **Paused and resumed.** Stopped at Hannah's request at 15:27 after 26 of 50
  stage-1 cells; the five in flight were removed (each held only its
  parameters file) and re-run on resume at 16:14.
- **Stage 2 re-chosen.** The first selection rule took each parameter's range
  of variation explained over both periods and counted every unscored run as
  its period's worst score. It chose Hs × Tp, for the wrong reasons: every
  2010–2024 natural run explains −800% to −2800%, so that window's ranges
  measure how badly a setting fails; the drowned Tp 12 became each period's
  worst score; and two crashes were counted as worst although the rule named
  drownings only. Hannah stopped that grid (12 runs finished, kept in
  `grid_wave_height_x_wave_period/`, selection in
  `tables/stage2_selection_first_rule.csv`) and chose the rule now in the
  driver: the range over runs that survived in 1996–2010, best value from the
  same runs, grid in both periods. It selects Hs (effect 11.0) and high-angle
  fraction (3.5) over Tp (1.9) and asymmetry (0.6): `tables/stage2_selection.csv`.
- **Two crashes, no groin.** 2010–2024 high-angle 0.2 (year 8) and asymmetry
  0.7 (year 13) died with no Python error: the silent access violation in
  Barrier3D's jitted `route_overwash`, which reads out of bounds once a domain
  has prograded (memory note of 2026-09-11, where it was found with the
  groin). Not confirmed with `NUMBA_BOUNDSCHECK=1` here. They are tables rows
  with status `process crashed …`, and open diamonds in the figures; a crash
  is not a model result, unlike a drowning.

## Results

Observed targets (`tables/observed_targets.csv`): 1996–2010 mean −0.42 m/yr,
flat-line RMSE 1.17; 2010–2024 mean **+1.06** m/yr, flat-line RMSE 1.46.

**1996–2010: a ridge, not a point.** Share of the alongshore variation
explained, Hs × high-angle fraction (asymmetry 0.8, Tp 8):

| Hs \ high-angle | 0.3 | 0.4 | 0.45 | 0.5 | 0.55 |
|---|---|---|---|---|---|
| 0.75 | −44% | +0% | +8% | +3% | +2% |
| 1.0 | −91% | −3% | +6% | +2% | +1% |
| 1.25 | −158% | −18% | +9% | **+12%** | +9% |
| 1.5 | −307% | −51% | +1% | +10% | +9% |
| 2.0 | −998% | −227% | −54% | +5% | +8% |

- A larger Hs needs a larger high-angle fraction: the good region runs from
  (0.75, 0.45) to (2.0, 0.55). Off it, the error rises steeply, fastest at low
  high-angle fraction and high Hs.
- On the ridge the bias is −0.2 to −0.5 m/yr and the best share explained in
  this grid is 12% (Hs 1.25, high-angle 0.5; RMSE 1.10 against the flat line's
  1.17): the model beats a flat line, but not by much. **The best natural
  1996–2010 run of the study is in the Hs × Tp grid below (+18%).**
- One-at-a-time, Tp 6–7 s scored best of any single change (+12%, +13%); Tp 10
  lost 181% and Tp 12 drowned the barrier in year 4, in both windows.
  Asymmetry mattered least (−50% at 0.3 to +6% at 0.8); reversing the drift
  direction (asymmetry < 0.5) made things worse, not better.
- Full management at the baseline did better than natural in 1996–2010 (17%
  vs 6%, bias −0.04 vs −0.44 m/yr).

**1996–2010: Hs × wave period** (asymmetry 0.8, high-angle 0.45). The grid the
first stage-2 rule chose; stopped when that rule was replaced, then finished
for 1996–2010 only (Hannah, 2026-09-24; 4 cells added, the rest reused):

| Hs \ Tp | 6 | 7 | 8 | 10 | 12 |
|---|---|---|---|---|---|
| 1.0 | +12% | +13% | +6% | −181% | drowned |
| 1.25 | +7% | +5% | +9% | **+18%** | drowned |
| 1.5 | −8% | −7% | +1% | +8% | +14% |
| 2.0 | −89% | −78% | −54% | −30% | −16% |
| 2.5 | −244% | −233% | −212% | −142% | −90% |

- **A diagonal ridge: a larger Hs needs a longer period** — (1.0, 7 s) +13%,
  (1.25, 10 s) +18%, (1.5, 12 s) +14%; bias −0.36 to −0.49 m/yr along it. Past
  Hs 1.5 no period reaches the flat line, though longer periods always help.
- **Long periods at low Hs fail**: Tp 10 at Hs 1.0 loses 181%, and Tp 12 drowns
  the barrier at Hs 1.0 and 1.25 but not at 1.5. BRIE's shoreface response
  grows with Tp^2.5 while its depth is set by Hs (8.9·Hs), so a long period on
  a shallow shoreface pulls the barrier down.
- **Best combination found, natural 1996–2010: Hs 1.25, Tp 10, asymmetry 0.8,
  high-angle 0.45** — +18%, bias −0.39 m/yr, RMSE 1.06 (flat line 1.17). The
  +12% to +18% cells along both ridges differ by less than the fit on one
  window can separate: read it as a region, not a value. Full management's
  best (its baseline, +17%, bias −0.04) was only swept one parameter at a time.

**2010–2024: no wave setting works.** Every scored natural run explains
between −739% and −2824%, with bias −3.9 to −7.3 m/yr (−3.9 to −5.2 across
the grid): at the baseline the model erodes about 3.4 m/yr where CoastSat
shows +1 m/yr. The best grid cell (Hs 2.0,
high-angle 0.55) still explains −739%. Full management halves the bias (−1.87
at the baseline) but the pattern is still wrong (pattern-only −49%). The
natural runs in this window also zig-zag from domain to domain, which the
observations do not. This is a problem of the window or its setup, not of the
wave climate: fills, the 2021 CoastSat step (memory note: a real +17 m jump),
and the 2010–2024 storm record are the places to look.

### Full management, stage 1 (added 2026-09-24)

After stage 1 showed management halving the 2010–2024 bias, Hannah asked for
the same stage-1 sweep under `full_management` (road, beach and dune
management, the historical fills; no relocations, no groin; zeroBE, metres),
filed as `full_management_<parameter>/` beside the natural folders and sharing
`baseline_full_management/`. 46 new runs: 44 scored, 2 drowned (Tp 12, year 4,
both windows, as in natural). 8 of them first crashed and were re-run on the
fixed Barrier3D (below).

Share explained / mean bias (m/yr), one parameter at a time:

| | natural 1996–2010 | managed 1996–2010 | natural 2010–2024 | managed 2010–2024 |
|---|---|---|---|---|
| baseline | +6% / −0.44 | **+17% / −0.04** | −1008% / −4.47 | −212% / −1.87 |
| best Hs | +9% (1.25) | +17% (1.0) | −821% (1.5) | −197% (1.25) |
| best high-angle | +6% (0.45) | +17% (0.45) | −980% (0.5) | −201% (0.5) |
| best asymmetry | +6% (0.8) | +17% (0.8) | −957% (0.6) | −193% (0.4) |
| best Tp | +13% (7) | +17% (8) | −1008% (8) | −212% (8) |

- **1996–2010: full management is better at every setting**, by about 5–10
  points of share explained, and its bias sits near zero (−0.2 to +0.0 on the
  good settings). The baseline is the best managed point of all four sweeps:
  every one-at-a-time change leaves it the same or worse. The parameters rank
  the same way as in natural: Hs and high-angle fraction matter, asymmetry
  least, Tp 10–12 fails.
- **2010–2024: management cuts the error by a factor of about five, but no
  setting reaches the flat line**: the best managed run explains −197%, bias
  −1.9 m/yr. Wave tuning barely moves it (−197% to −263% over Hs 0.75–2.0).
  So even with the fills, the model erodes where CoastSat shows gain.
- **Eight managed 2010–2024 runs first crashed, all in year 13**: every
  high-angle value below the baseline (0.1–0.4) and every asymmetry below 0.7
  (0.3–0.6). The cause turned out to be an indexing bug in Barrier3D, not an
  event in the record (next section); re-run on the fix, all eight scored.
  With them the managed 2010–2024 picture is unchanged: the best asymmetry run
  (0.4) explains −193%, lower high-angle fractions are worse (−408% at 0.1).

### The crashes: an index bug in Barrier3D (confirmed 2026-09-24)

**Resolved.** All 10 crashed cells (natural 2010–2024 high-angle 0.2 and
asymmetry 0.7; managed 2010–2024 high-angle 0.1–0.4 and asymmetry 0.3–0.6)
were re-run on the fixed Barrier3D (branch `fix/route-overwash-axis-swap`,
now the one in use) and all are scored. `tables/all_runs.csv` says which
Barrier3D each run used (`barrier3d_route_overwash_fix`: True for these 10,
False for the rest, which the fix changes negligibly). The crash logs are kept
under `logs/<group>/2010_2024/crashed_before_fix/`. One re-run
(managed high-angle 0.1) finished its 14 years and then failed to replace the
shared `run_index.csv` while a parallel run held it (a Windows lock); its
outputs were complete and are scored, and the index write now retries.

Re-running crashed and scored cells with `NUMBA_BOUNDSCHECK=1` (logs under
`../2026-09-24-metres-3-barrier3d-overwash-fix/logs/diagnostics/`, moved there in the reorganisation):

| run | bounds-checked result |
|---|---|
| managed 2010–2024, high-angle 0.4 (crashed in year 13 unchecked) | `IndexError` in year 2 |
| managed 2010–2024 baseline (**scored**) | `IndexError` in year 2 |
| natural 1996–2010 baseline (**scored**) | `IndexError` around year 10 |
| managed 2010–2024 baseline, offset ÷10 on the pre-metres build | `IndexError` in year 2 |

With the JIT off the error is exact: `IndexError: index 45 is out of bounds
for axis 1 with size 45` at `Barrier3D/barrier3d/barrier3d.py:1092`, in
`route_overwash`:

```python
if Elevation[TS, d, i] > SL or np.sum(np.greater(Elevation[TS, i, d + 1: d + 10], SL)) > 0:
```

Everywhere else in the routine `d` is the cross-shore row (axis 1) and `i` the
alongshore column (axis 2), e.g. `SedFluxIn[TS, d + 1, i]` two lines down. The
second clause has them swapped: it reads row `i`, columns `d+1…d+10`, where the
cells nine rows landward of the flow, `Elevation[TS, d+1:d+10, i]`, are
evidently meant. Upstream code (`git blame`: I. Reeves, 2024-03-01,
UNC-CECL/Barrier3D), unmodified in this copy.

- **Every run with overwash reads the wrong cells** wherever `i` is less than the
  domain's row count: the subaerial test that switches how overwash sediment is
  routed looks across the island instead of landward. This is in every run in
  the project, ÷10 and metres alike.
- **Where `i` is not less than the row count** (a domain has 50 columns; at
  initialisation 2 of 90 real domains have fewer than 50 rows, 32 and 29), the
  read is outside the array. Compiled without bounds checking it returns
  whatever is in memory, and crashes only when it hits protected memory: the 10
  crashes in this study.
- The 2026-09-11 finding (the groin crash) was this same line; prograding a
  domain was one way to reach it, not the cause.

**Measured 2026-09-24** (`experiments/2026-09-24-metres-3-barrier3d-overwash-fix/NOTE.md`): the
one-line fix, on the local Barrier3D branch `fix/route-overwash-axis-swap`,
removes every out-of-bounds read and the crashes, and moves scores only in the
third decimal (per-domain rates by at most 0.11 m/yr except in natural 2010–2024,
where 14 domains move by more than 0.1 and one by 2.2 m/yr). No result here
changes. Barrier3D `master` is still unpatched and checked out.

## Figures

| figure | shows |
|---|---|
| `stage1/scores_by_<parameter>_both_periods.png` | share explained, bias, RMSE against each parameter, both windows, natural (solid) and full management (dashed) |
| `alongshore/<window>/rate_and_position_change_by_<parameter>_<window>.png` | rate and position change along the island for each value, against CoastSat, natural |
| `alongshore/<window>/rate_and_position_change_by_<parameter>_full_management_<window>.png` | the same under full management |
| `management/natural_vs_full_management_baseline_both_periods.png` | the baseline, natural against full management, both windows |
| `stage2/grid_wave_height_x_high_angle_both_periods.png` | the Hs × high-angle grid as lines: share explained and bias against Hs, a line per high-angle fraction |
| `stage2/grid_wave_height_x_wave_period_1996_2010.png` | the Hs × Tp grid, 1996–2010, a line per period |

## Not yet done

- Why 2010–2024 fails: a natural run with no fills cannot match an observed
  gain, but full management fails too; check the fills, the 2021 CoastSat step
  and the storm series for this window.
- Confirm the two crashes with `NUMBA_BOUNDSCHECK=1`.
- No search covered all four parameters together, and full management was swept
  one parameter at a time only.

## Reproduce

```
python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py run stage1
python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py score
python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py run stage2
python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py score
python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity_plot.py
```
