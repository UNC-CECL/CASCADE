# Storm length selection: which duration rule matches the observed overwash? (2026-09-28)

**Question (Hannah).** Which storm series should the hindcast use? The first priority is that the model's overwash matches the observed record.

**Candidates.** Every candidate keeps every event, including Isabel, March 2018 and Florence. `trimL` cuts any event longer than L hours to the L hours around its peak total water level. The control is `drop72`, the committed series, which drops events over 72 h; its runs are the matrix runs. `full` trims nothing: the longest event in 1996–2024 is 193 h, so it equals a 240 h limit.

- **Built.** The builder's own functions produce the series, written to `storms/`. `drop72` rebuilds the committed series exactly in both windows.
- **Runs.** The unchanged runner, with the storm file swapped in its own process. Barrier3D is the current code (49fd069).
- **Settings.** Edge rates are the matrix values, solved on `drop72`. Stage 2 re-solves them.

Driver: `scripts/hatteras_ms/experiments/HAT_storm_length_selection.py`.

## Stage 1: overwash against the imagery (32 runs, all ran to the end)

**How it is scored.**

- **Matching storms to images.** Each image is compared with the model's overwash since the previous image. A model year's overwash is shared among that year's storms that reached a dune gap, and each storm is dated by its end less the 7-day grace.
- **The gaps.** They are recomputed exactly: the pre-storm crest comes from `DuneGrowth` and the gaps from the model's own `DuneGaps`. The crest matches the replay's exactly.
- **Checking the sharing.** Replaying sampled domain-years (`validate`): the storm that carried the most overwash always received a share, and the median volume credited to the wrong storm was 16%.
- **Scores.** POD (hit rate), POFD (false-alarm rate), PSS = POD − POFD, plus two correlations: per-image counts (timing) and per-domain frequency (space).

Managed runs (full_management; the imagery shows the managed island), overwash threshold 0:

| series | PSS 1996–2010 | PSS 2010–2024 | timing r 96 / 10 | interior RMSE 96 / 10 (m/yr, matrix ends) |
|---|---|---|---|---|
| drop72 (current) | **−0.07** | **0.37** | −0.18 / 0.66 | 1.05 / 2.25 |
| trim24 | 0.36 | 0.05 | 0.42 / 0.34 | 1.11 / 2.26 |
| trim36 | 0.35 | 0.14 | 0.43 / 0.35 | 1.27 / 2.70 |
| trim48 | 0.35 | 0.08 | 0.43 / 0.34 | 1.38 / 3.02 |
| trim72 | 0.34 | 0.08 | 0.42 / 0.35 | 1.58 / 3.41 |
| trim96 to full | 0.34–0.37 | 0.03–0.08 | 0.43–0.47 / 0.32–0.34 | 1.74–1.96 / 3.65–4.00 |

At thresholds of 1 and 5 m³/m the ranking is the same: 1996 PSS 0.37–0.45 for the trims against −0.07 to −0.10 for `drop72`; 2010 PSS 0.26 and 0.17 for `drop72`, −0.09 to 0.06 for the trims.

**What it shows.**

1. **Keeping the dropped storms fixes 1996–2010.** PSS goes from −0.07 to about 0.35, and the timing correlation from −0.18 to about 0.43. Isabel's 2004 image is the difference: 48 domains observed, 12 modelled on `drop72`, 90 on every trim.
2. **In 2010–2024 the same storms make overwash agreement worse.** PSS falls from 0.37 to 0.03–0.14 because the false-alarm rate rises from 0.44 to about 0.76. The added storms overwash 40–90 domains in the 2014, 2017, 2018 and 2019 image windows, where the imagery shows 0–10.
3. **Trim length hardly changes where and when overwash happens.** Whether a domain overwashes depends on whether the storm is in the series and how high Rhigh gets, not on how long it lasts. What length changes is how much sand moves, and so shoreline retreat: interior RMSE rises steadily with trim length in every run.
4. **The biggest mismatch is not storm length.** On every series, `drop72` included, the model overwashes far more domains than the imagery shows in most windows. For example, the 2008-06 image shows 0 domains while the model has 85–90, and 2005-09 shows 0 while the model has 10–70 (`figures/stage1_per_image.png`). That excess is what limits every score.

Possible causes, not tested here:

- the 2004-start v1 dune windows, which clip the crest in 20 domains and feed the 2010–2024 runs;
- run-up computed with a beach slope of 0.06;
- washover that faded, or was cleared off NC-12, before the image.

**Stage 1 reading.** Keep every storm and trim short. `trim24` gives the best or near-best PSS in 1996–2010 at every threshold, has the smallest shoreline cost (RMSE 1.11 against 1.05; 2.26 against 2.25, before the ends are re-solved), and has fewer total storm-hours than `drop72` in both windows. `trim36` is the runner-up and does best among the trims in 2010–2024 at threshold 0.

**Stage 2 is on hold** until Hannah decides whether to run it, and on which candidates.

Files: `storms/`, `runs/<variant>_<scenario>/`, `logs/`, `tables/stage1_scores.csv`, `tables/stage1_cells.csv`, `tables/sharing_validation.csv`, `figures/stage1_scores.png`, `figures/stage1_per_image.png`.
