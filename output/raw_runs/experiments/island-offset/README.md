# island-offset

How the island's planform (the BRIE island offset) is set: its units, and whether it comes from the dune line or the CoastSat shoreline.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-09-22-div10-offset-shoreline-trial-original`](2026-09-22-div10-offset-shoreline-trial-original/NOTE.md) | Does a shoreline-derived offset change 1996–2010? | Could not say: the runs used the ÷10 offset, so the change was a tenth of the real difference. | superseded by `2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8` |
| [`2026-09-24-div10-vs-metres-wave-sweep`](2026-09-24-div10-vs-metres-wave-sweep/README.md) | Offset in metres (as BRIE reads it) or ÷10? Can waves re-tune the metres runs? | Metres is correct; the ÷10 was a units bug. Metres became the default on 09-24. | **current** (the units decision) |
| [`2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8`](2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8/README.md) | Dune line or shoreline as the offset, in metres (Hs 1 / Tp 8 / asym 0.8, high-angle 0.3–0.55)? | Managed: no difference (+23–24% either way). Natural: the shoreline offset scores higher (+22% vs +14% at 0.45). | superseded by `2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a` (same question at option A) |
| [`2026-09-28-div10-offset-duneline-vs-shoreline-rebuild`](2026-09-28-div10-offset-duneline-vs-shoreline-rebuild/README.md) | A REBUILD, not run at the time: the 09-25 design (natural + full management, zeroBE, both offsets) on the ÷10 offset and ÷10-era waves (Hs 2.5 / Tp 8 / asym 0.7 / high-angle 0.1), current Barrier3D. | Kept to set ÷10 beside metres in the same design; the model's profile is smooth and the two offsets differ by at most ~4 m. The record of what ran on ÷10 is `2026-09-22`. | record |
| [`2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a`](2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/README.md) | Dune line or shoreline as the offset at the option A waves, zeroBE; then each graded on its own feature, both periods. | Shoreline offset slightly ahead at the domain scale (+3–4 points raw, level when smoothed); the offset source is a small lever. Relation to `comparisons/target_comparison/` at the end of its README. | superseded by `2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup` (same design on the adopted setup) |
| [`2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup`](2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup/README.md) | Does shoreline offset v2 (±1 yr of the DEM's lidar) change anything against v1 (calendar mean)? And does dune line vs shoreline hold on the adopted setup (dune ceilings, split12 storms, cap fix), both periods, natural and managed? | v1 → v2 is not a lever: scores within 0.015 raw / 1.5 points smoothed, model change within 5.7 m (mean 0.1 m). The shoreline leads the dune line in 1996–2010 (+2–3 raw; +7–8 smoothed managed); 2010–2024 fails in every arm. | **current** |

**Figures.** Every study that asks the dune-line-or-shoreline question draws the
same four figures through one function (`house_figures` in
`scripts/hatteras_ms/experiments/HAT_offset_source_comparison.py`), so they compare
directly: `duneline_offset_vs_duneline_change_*`, `shoreline_offset_vs_coastsat_total_change_*`,
`shoreline_offset_vs_coastsat_projected_change_*`, `total_change_difference_shoreline_minus_duneline_*`;
net change in metres, observations LOWESS 7 domains, model unsmoothed, fills marked, no
scores on the figures. `2026-09-24-div10-vs-metres-wave-sweep` asks a different question
(offset scale against wave tuning) and keeps its sweep figures.

**Status** — **current**: its answer is in use now. **superseded**: a later study
re-asked it; follow the pointer. **record**: a finished check or a result from an
earlier set-up (÷10 offset, Hs 2.5 calibration), kept so the number can be traced.

Back to [the map](../README.md).
