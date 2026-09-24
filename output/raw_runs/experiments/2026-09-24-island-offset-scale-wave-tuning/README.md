# Island-offset scale against wave-climate tuning (2026-09-24)

## The question

BRIE's shoreline position (`brie.x_s`) is in metres, and so is the island-offset
file. But every calibrated run hands BRIE the offset **divided by ten**
(`offset_mode: asrun`), a units error: `cascade/brie_coupler.py:369` adds the
offset with no conversion, and the only decametre conversion is on the
Barrier3D side (line 344). The Ocracoke scripts in
`scripts/input_prep/2-brie-offset/1-produce/colleague_old_version/` hand
CASCADE metres (a decametre file × 10).

The 2026-08-21 test put the offset back at full scale and skill collapsed, but
that test kept the wave climate tuned for the shrunken planform. **With the
offset at full scale, can the model be re-tuned through its wave climate to
match the ÷10 runs?** Hs is treated as a tuning parameter, not a physical one.

## Two sweeps

| sweep | varied | fixed | offset scales | offset sources | runs |
|---|---|---|---|---|---|
| `wave_height` | Hs 0.5, 0.6, 0.65, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0 m | asymmetry 0.7, high-angle 0.1 | all three | dune line, shoreline | 54 |
| `wave_angle` | asymmetry 0.5–0.8 × high-angle fraction 0.1–0.4 at Hs 1.0; 0.45 and 0.5 for metres only; metres at Hs 0.75, asymmetry 0.8, high-angle 0.3–0.5 | — | all three | dune line | 60 |

Both sweeps: 1996–2010, zeroBE, full management, no groin, relocations off,
base geometry, Tp 8 s. zeroBE means no end-domain rates are imposed, so no rate
solved under one setting is carried into another.

## How to read the labels

**Offset scale.** The folder label, and the value the runner reads:

| folder label | `HAT_OFFSET_MODE` | what BRIE receives |
|---|---|---|
| `div10` | `asrun` | offset ÷ 10: every calibrated run (0.6 km span) |
| `metres` | `metres` | the measurement as is, including the island's lean (6.2 km span) |
| `metres-detrended` | `detrended` | metres with the linear trend removed; curvature kept (2.1 km span) |

**Offset source.** `duneline` = `2-brie-offset/1996/duneline/v1` (CURRENT);
`shoreline` = `2-brie-offset/1996/shoreline/v1` (CoastSat 1995–1997 mean).

**Where each setting is written:**

| place | example | complete? |
|---|---|---|
| folder | `runs_wave_angle/metres_duneline/1996_2010/zeroBE/` | sweep, scale, source, period, preset |
| log file | `logs/wave_angle/metres_duneline/Hs1.0_asymmetry0.6_highangle0.3.log` | yes: every tuned value |
| tables | `tables/all_runs.csv` | yes: every setting is a column |
| run name | `HAT_1996_2010_zeroBE_offsetmetres_road_bdm_nogroin_waveHs1asym0p6ahf0p3` | **no**, see below |

**Run names leave out anything at its default.** The name is the runner's
own. A `div10` run has no `offset` token, asymmetry 0.7 has no `asym` token and
high-angle fraction 0.1 has no `ahf` token. So `…_road_bdm_nogroin_waveHs1` is
÷10, asymmetry 0.7, high-angle 0.1. Use the folder, the log name or the
tables for a complete label. `score` reads every setting from each run's
metadata and stops with an error if the metadata disagrees with the folder
holding the run.

## Layout

```
README.md          this file
tables/
  all_runs.csv            both sweeps, one row per run, drowned runs included
                          with the reason in `status`; scores in the columns
                          described under "Alongshore variation explained"
  observed_target.csv     the target's mean, spread (sd) and flat-line RMSE
  wave_height_sweep.csv   the wave_height rows of all_runs.csv
  wave_angle_sweep.csv    the wave_angle rows
figures/
  wave_height/   rmse_bias_vs_wave_height_by_offset_scale_1996_2010.png
  wave_angle/    rmse_bias_vs_high_angle_fraction_by_asymmetry_1996_2010.png
  combined/      best_rmse_by_offset_scale_1996_2010.png
                 bias_vs_rmse_all_runs_1996_2010.png
                 alongshore_rate_best_runs_vs_coastsat_1996_2010.png
                 variance_explained_by_offset_scale_1996_2010.png
                 spread_vs_correlation_all_runs_1996_2010.png
  alongshore_sensitivity/<scale>/
                 rate_and_position_change_by_<parameter>_<scale>_1996_2010.png
                 parameter = wave_height, high_angle_fraction, asymmetry;
                 9 figures, one per scale x parameter
  <each>/supporting/   CAPTIONS.md, PDFs, the CSV behind each summary figure
logs/
  wave_height/<scale>_<source>/Hs<h>_asymmetry<a>_highangle<f>.log
  wave_angle/<scale>_<source>/Hs<h>_asymmetry<a>_highangle<f>.log
  drivers/        the console output of each sweep launch
runs_wave_height/<scale>_<source>/1996_2010/zeroBE/<run_name>/
runs_wave_angle/<scale>_<source>/1996_2010/zeroBE/<run_name>/
```

Each run folder holds its parameters yaml, `_run_metadata.json` / `.txt`, the
shoreline matrix (`.npy`), `tables/shoreline_change_rate.csv` (the modelled
rate for each domain) and two per-run figures. The model state (`.npz`) was
not saved.

**What is in git and what is not.** The run folders (`runs_wave_height/`,
`runs_wave_angle/`, `checks/`) are on disk only: their full paths reach 339
characters, past Windows' 260, and git will not index them (`.gitignore` here
says so). Committed: this README, `tables/`, and under each
`figures/<folder>/supporting/` the PDF, caption and data CSV of every figure
(PNGs are ignored under `output/raw_runs/` repo-wide). `tables/all_runs.csv`
has every setting and score of every run; the per-domain rates the figures are
drawn from are in the run folders only, so redrawing the figures needs the runs
on disk (or a re-run: every setting is in the table and the log name).

## Results

The score is the interior RMSE and mean bias (GIS 2–89) of the modelled LRR
rate against the CoastSat LRR target (LOESS, 10 domains). **Bias is the score
to read, not r.** Every tuned value here is fitted on the window it is scored
on, so read the best settings as a band, not a calibrated value.

**Wave height** (dune-line source; the shoreline source is within 0.3 m/yr
everywhere), interior RMSE, m/yr:

| Hs (m) | ÷10 | metres | metres, trend removed |
|---|---|---|---|
| 0.5 | drowned yr 4 | drowned yr 4 | drowned yr 4 |
| 0.6 | drowned yr 8 | drowned yr 6 | drowned yr 6 |
| 0.65 | 2.80 | 3.98 | 3.42 |
| 0.75 | 1.07 | 1.63 | 1.85 |
| 1.0 | **1.05** | 2.08 | 2.38 |
| 1.5 | 1.05 | 3.69 | 4.14 |
| 2.0 | 1.05 | 6.46 | 6.54 |
| 2.5 | 1.07 | 9.92 | 8.91 |
| 3.0 | 1.14 | 13.73 | 11.02 |

**Wave angles** (Hs 1.0 m, dune-line source), interior RMSE / mean bias, m/yr:

| high-angle fraction | ÷10 (any asymmetry) | metres, asym 0.8 | metres, trend removed, asym 0.5 |
|---|---|---|---|
| 0.1 | 1.05 / +0.01 | 2.02 / −0.75 | 2.31 / −0.66 |
| 0.2 | 1.07 / +0.02 | 1.67 / −0.58 | 1.93 / −0.49 |
| 0.3 | 1.09 / +0.04 | 1.32 / −0.39 | 1.54 / −0.31 |
| 0.4 | 1.11 / +0.05 | 1.11 / −0.18 | **1.24 / −0.13** |
| 0.45 | — | **1.07 / −0.04** | — |
| 0.5 | — | 1.11 / +0.02 | — |

**Lower Hs, and Hs combined with the high-angle fraction** (added the same
day). The barrier drowns at Hs 0.6 at every scale (years 6–8) and survives at
0.65, but at 0.65 every scale is far worse (bias about −2 m/yr, RMSE 2.5–4):
close to drowning, the whole island erodes. So in Hs alone the full-scale
minimum is at 0.75, now bracketed. Metres at Hs 0.75 with asymmetry 0.8:

| high-angle fraction | RMSE / bias at Hs 0.75 | RMSE / bias at Hs 1.0 | variation explained, 0.75 / 1.0 |
|---|---|---|---|
| 0.3 | 1.19 / −0.21 | 1.32 / −0.39 | −4% / −27% |
| 0.4 | 1.11 / −0.07 | 1.11 / −0.18 | 10% / 11% |
| 0.45 | 1.10 / +0.02 | **1.07 / −0.04** | 12% / **17%** |
| 0.5 | 1.12 / +0.06 | 1.11 / +0.02 | 9% / 10% |

The two levers do not add up: with the high-angle fraction at its best,
lowering Hs to 0.75 only damps the pattern further (spread 0.30× observed
against 0.37× at Hs 1.0) and scores no better.

What this shows:

1. **The ÷10 runs are insensitive to the wave climate.** RMSE stays within
   1.05–1.14 across Hs 1–3 and across every wave-angle setting. With a nearly
   straight shoreline, the wave climate only nudges the mean.
2. **At full scale, the high-angle fraction is the lever, not Hs.** Tuning Hs
   gets metres no lower than 1.63 before the barrier drowns (Hs sets the
   shoreface depth, 8.9·Hs). Raising the high-angle fraction from 0.1 to 0.4
   takes metres from 2.02 to **1.11**, within 0.06 of the ÷10 runs, and
   brings the bias from −0.75 to −0.18. High-angle waves are anti-diffusive,
   so a larger share slows BRIE's smoothing of the island's real curvature.
3. **For metres the best high-angle fraction is 0.45.** At asymmetry 0.7
   and 0.8, RMSE bottoms out there (1.07 at 0.8) and rises again at 0.5;
   bias crosses zero between 0.45 and 0.5. At asymmetry 0.5 and 0.6 it is
   still falling at 0.5 (1.12, 1.09). Above about 0.5 most of the coast is
   anti-diffusive and BRIE clamps its diffusivity at zero. Trend-removed was
   not extended past 0.4.
4. **Asymmetry matters less, and in opposite directions.** Metres improves
   with more asymmetry (0.8 best); trend-removed improves with less (0.5
   best). That fits asymmetry acting on the net drift, which is what the lean
   drives in metres and what trend removal takes away.
5. **Metres beats trend-removed once the angles are tuned** (1.07 vs 1.24).
6. **The two scales fail in different ways along the island** (alongshore
   figure, best run of each). The ÷10 run is smooth and flattens the
   observed pattern. The full-scale runs recover some of it: the
   trend-removed run has accretion peaks near GIS 18–20 and 30–32 where
   CoastSat has them, and both reach the erosion at Rodanthe (GIS 76–78),
   which the ÷10 run underestimates. But they add sharp local dips
   CoastSat does not have (GIS 6, 35, 81), and metres drops at GIS 90. None
   of them has the accretion near GIS 40–44 or the erosion north of GIS 80.

### Against a flat line: read this before the RMSE numbers

A model that predicted the observed mean rate at every domain would score
RMSE **1.174** (the spread of the CoastSat target along the island, sd 1.17 m/yr).

| run | interior RMSE | spread along the island, sd (m/yr) | r with CoastSat |
|---|---|---|---|
| flat line at the observed mean | 1.174 | 0 | — |
| best ÷10 (Hs 2.0) | 1.046 | 0.58 | 0.47 |
| best metres (asymmetry 0.8, high-angle 0.45) | 1.070 | 0.43 | 0.42 |
| best trend-removed (asymmetry 0.5, high-angle 0.4) | 1.240 | 1.11 | 0.42 |

- **Every run beats the flat line by 11% at most.** The ÷10 runs and the best
  metres run score well mainly because their rates vary little along the island
  (sd 0.4–0.6 against the observed 1.17), not because they reproduce the pattern.
- **Raising the high-angle fraction improves metres largely by flattening it.**
  It slows BRIE's smoothing, so the full-scale curvature drives less change;
  the best metres run is nearly as flat as the ÷10 run (alongshore figure).
- **The trend-removed run is the only one with the observed amount of
  alongshore variation** (sd 1.11), but it puts some of it in the wrong places,
  so its RMSE is worse than the flat line.
- So an RMSE close to ÷10 does not mean the full-scale offset is as good.
  Both are close to the flat-line baseline, and r is about 0.4–0.5 for all
  three. Earlier smoothing-scale tests found r in this setup does not clear
  its null, so r does not separate them either.

### Alongshore variation explained

`score` adds, for every run, how much of the observed alongshore variation
(GIS 2–89) it explains, computed from the run's own rate table against the
same target the runner uses (each run's RMSE is recomputed from it and must
match the runner's):

| column | meaning |
|---|---|
| `variance_explained` | 1 − Σ(model − obs)² / Σ(obs − obs mean)². 100% perfect, 0% = the flat line, negative = worse than it. Bias counts against it. Equals 1 − (RMSE / 1.174)². |
| `pattern_variance_explained` | the same with each series' own mean removed first: the pattern alone |
| `r_alongshore` | correlation of the modelled and observed alongshore series |
| `model_sd_m_yr`, `sd_ratio` | the modelled spread along the island, and that ÷ the observed 1.17 m/yr |

Best run of each scale (dune-line source):

| scale | best run | explained | pattern only | r | spread ÷ observed |
|---|---|---|---|---|---|
| ÷10 | Hs 2.0 | **21%** | 22% | 0.47 | 0.50 |
| metres | asymmetry 0.8, high-angle 0.45 | **17%** | 17% | 0.42 | 0.37 |
| metres, trend removed | asymmetry 0.5, high-angle 0.4 | **−12%** | −10% | 0.42 | 0.95 |

- **No run explains more than about a fifth of the alongshore variation.**
  Every ÷10 run beats the flat line; 28% of metres runs do; no trend-removed
  run does. Bias hardly matters here: the pattern-only score is within 2
  points of the full one at every best run.
- **The two ways to miss are visible in the spread-vs-correlation figure.**
  The ÷10 runs are too flat (spread 0.18–0.87× observed) with r 0.38–0.51. The
  full-scale runs trade one against the other: in metres, the more the
  shoreline moves, the better placed the pattern (r rises from 0.31 at
  spread 0.2× to 0.63 at 3–11×), but no run has both good placement and the
  observed amount of spread. The wave-angle trend-removed runs have about
  the observed spread (1–2×) but r only 0.42–0.48; at higher Hs their
  spread grows to 9× and r falls to 0.31.
- **The best-placed pattern of any run is metres at high Hs** (r 0.61–0.63,
  Hs 1.5–3.0), but it is 3–11× too strong and biased by −1.5 to −5 m/yr, so
  it explains far less than nothing. Treat r cautiously: earlier tests found
  r in this setup does not clear its null.

### Alongshore sensitivity: rate and position change

`figures/alongshore_sensitivity/<scale>/` has one figure per offset scale and
parameter. Each moves one parameter and holds the other two at their
defaults: Hs with asymmetry 0.7 and high-angle 0.1; high-angle fraction and
asymmetry at Hs 1.0 with the other angle parameter at its default. So the
best tuned combinations (e.g. metres at asymmetry 0.8, high-angle 0.45) are
not drawn there; they are in `combined/alongshore_rate_best_runs_vs_coastsat`.
Each figure has two panels:

- **(a) Rate**: the modelled LRR against the CoastSat LRR target.
- **(b) Position change**: the modelled shoreline position at the end of 2010
  minus the start of 1996 (`change_rate_m_yr × 14`), against the observed
  CoastSat change, which is the mean position over calendar 2010 minus the
  mean over calendar 1996 (`5-scr/3-rates/coastsat/total_change/1996_2010/`,
  smoothed at 10 domains). No rate in the observed side. Seaward positive.

Each run's line is labelled with its share of the rate variation explained and
its bias; a drowned run is listed in the legend with no line.

What they show:
- **÷10 (any parameter):** the same broad shape each time, with erosion at
  GIS 35–38 and 70–85; Hs mostly scales how deep it goes. It misses the
  observed accretion near GIS 17, 29, 42 and 70 and the erosion at 7–12
  and 22–24.
- **metres, high-angle fraction:** the erosion hotspots stay at the same
  places (GIS 6, 35–38, 76–78, 90) and shrink as the fraction rises; the
  pattern does not move toward CoastSat, it only gets weaker. GIS 35–38 and
  76–78 coincide with observed erosion, but the model's is several times too
  deep at low fractions, and GIS 6 and 90 are not in the observations.

Runs with no score: every Hs 0.5 run (barrier drowned by height and width in
year 4, all scales) and every Hs 0.6 run (years 6–8), so a floor on Hs
between 0.6 and 0.65, not an offset effect; metres-detrended
shoreline at Hs 3.0 (drowned by height, year 1). NC-12 drowned in up to 4
domains in the metres runs at Hs 2.5 and 3.0; in no other run.

Regression check: the ÷10 dune-line run at Hs 1.0, asymmetry 0.7,
high-angle 0.1 exists in both sweeps and gives identical scores (RMSE
1.052757, bias +0.005503).

## Metres or trend removed: what trend removal does, and why it was rejected

**Leading choice (Hannah, 2026-09-24): metres, no trend removal.** Metres is
correct for the model, and it scores better than trend removal on every
measure here. Trend removal is recorded as an alternative that was tested and
rejected, not a candidate.

**It is not a units choice.** `metres` and `metres-detrended` are both in
metres, the unit BRIE needs. The ÷10 is the units error; trend removal is a
separate question about the island's overall tilt.

**What it does.** The offset gives each domain's cross-shore position against
the model's straight alongshore axis. Along GIS 1–90 that has two parts: a
straight tilt of about 6 km end to end (about 7°), and about 2 km of bends on
top of it (the capes and bays). Trend removal fits a straight line through the
90 values and subtracts it, which keeps the bends and drops the tilt. In
effect it turns the model's axis to run along the island's overall direction
instead of along the baseline the offsets were measured from. Figure
`wave_height/…` panel (a) draws the three versions.

**Why it was tested.**

1. *The tilt depends on the reference frame.* It exists because the offsets
   were measured from a baseline that is not parallel to the island. A
   different baseline changes the tilt; the coast is the same. So it was a
   fair question whether the tilt is physics or bookkeeping.
2. *BRIE treats the tilt as real.* It compares each domain's shoreline angle
   with the wave directions, both against its own axis (`brie.py` ~820 and
   ~1297, `coast_diff[90 − θ]`). A 7° tilt makes every wave arrive 7° off
   and drives net alongshore transport in one direction along the whole
   island. In the 2026-08-21 test (1984–2004, end rates imposed) the tilt
   alone gave 19× the observed net drift, against 2.6× with the trend
   removed.
3. *BRIE's domain wraps around.* The last domain joins back to the first, so
   the 6 km tilt has to be closed through the invented buffer domains, which
   forces a steep artificial shoreline there (up to 29° in `metres`). Without
   the tilt the closure is easy.

**Why it was rejected.**

| | metres | trend removed |
|---|---|---|
| best interior RMSE (m/yr) | **1.07** | 1.24 |
| alongshore variation explained | **17%** | −12% (worse than a flat line) |
| runs that beat the flat line | 28% | none |

The wave-angle sweep shows why the tilt is not a problem in practice.
Asymmetry, the share of waves from each side, is the parameter that sets
net drift, and it is tuned anyway: the best metres run uses asymmetry 0.8 and
the best trend-removed run 0.5. Tuning asymmetry absorbs the lean, so the
problem trend removal was meant to solve is already handled, without changing
the measured geometry. That also fits the project's preference for changing
input data as little as possible.

**Still to watch.** The steep wrap-around buffer in `metres` (29°) stays below
the ~42° where BRIE's shoreline becomes unstable (the sign change in
1.2 sin²θ − cos²θ), but it may be behind the metres dip at GIS 90 in the
alongshore figures. Check it when the end domains are solved.

## Not yet done

- The end domains solved (edgeBE) at the best metres setting (Hs 1.0,
  asymmetry 0.8, high-angle fraction 0.45), for a comparison with the
  calibrated ÷10 edgeBE runs; watch GIS 90 there (see above).

## Adopted: metres is the default (2026-09-24)

After this study the runner's default became `offset_mode: metres`
(`HAT_hindcast_config.py`, `hat_run.yaml`), and every offset build was
re-padded as **v2** with the model's own smooth wrap-around
(`cascade_pipeline.hindcast.pad_offset_ring`), so the padded file is exactly
what metres mode hands Cascade. v2's real domains are identical to v1's; only
the buffer domains changed. The groin sweep now takes an offset mode too
(`HAT_SWEEP_OFFSET_MODE`, default the runner's) and files a non-asrun sweep
under an `offset<mode>` suffix.

`checks/defaults-after-metres-switch/` is one run of the runner with nothing
set: all defaults, v2 offset. It reproduces this study's
`runs_wave_height/metres_duneline/…_offsetmetres_road_bdm_nogroin` (Hs 2.5)
exactly (interior RMSE 9.920179382773128, bias −3.895997977190597), on
`duneline/v2` where the study ran `duneline/v1`. So no run in this study is
changed by the switch. The `asrun` runs here were made on v1 and reproduce
only with `HAT_OFFSET_VERSION_<year>=v1`.

## History

- The wave-height sweep ran first as `experiments/2026-09-24-offset-units-hs/`,
  the wave-angle sweep as `experiments/2026-09-24-offset-units-waveangle/`.
  Both were moved here the same day. `div10` was `asrun`, and
  `metres-detrended` was `detrended`. Each run's metadata tag was rewritten to
  its new path and the run index rebuilt. The run count was checked before
  and after the move (33 + 48 scored runs with metadata).
- The two drivers and the plot script of those folders were replaced by the
  two scripts below and deleted.

## Reproduce

```
python scripts/hatteras_ms/experiments/HAT_offset_scale_wave_tuning.py run wave_height
python scripts/hatteras_ms/experiments/HAT_offset_scale_wave_tuning.py run wave_height --scales metres metres-detrended --hs 0.5 0.75
python scripts/hatteras_ms/experiments/HAT_offset_scale_wave_tuning.py run wave_height --scales div10 --hs 0.5
python scripts/hatteras_ms/experiments/HAT_offset_scale_wave_tuning.py run wave_height --scales div10 --hs 0.75
python scripts/hatteras_ms/experiments/HAT_offset_scale_wave_tuning.py run wave_height --hs 0.6 0.65
python scripts/hatteras_ms/experiments/HAT_offset_scale_wave_tuning.py run wave_angle --scales metres --hs 0.75 --asymmetry 0.8 --high-fraction 0.3 0.4 0.45 0.5
python scripts/hatteras_ms/experiments/HAT_offset_scale_wave_tuning.py run wave_angle
python scripts/hatteras_ms/experiments/HAT_offset_scale_wave_tuning.py run wave_angle --scales metres --high-fraction 0.45 0.5
python scripts/hatteras_ms/experiments/HAT_offset_scale_wave_tuning.py score
python scripts/hatteras_ms/experiments/HAT_plot_offset_scale_wave_tuning.py
```

`run` skips any setting whose log shows a clean finish; `--overwrite` redoes it.
