# smoothed_lowess7_with_cascade — the model against both candidate targets

> **Redrawn 2026-09-29 after the dune-cap fix** (the zeroBE matrix run, rerun with the fix; see `../README.md`).
>
> **Redrawn 2026-09-28 on the adopted model** (the rebuilt zeroBE matrix run: Barrier3D `hatteras/adopted` (the three overwash fixes and per-cell dune ceilings), storms `v3_trim24` (every event kept, trimmed to 24 h around its peak); metres offset, Hs 2.0 / Tp 7.5 / asym 0.6 / high-angle 0.5). The numbers in the text below are from earlier versions; the current ones are in the tables and captions.

> **Redrawn 2026-09-27 on option A** (the option A zeroBE matrix run; metres offset, Hs 2.0 / Tp 7.5 / asym 0.6 / high-angle 0.5). The numbers in the text below are the /10-era ones; the current interior means are in `PROVENANCE.md`: model minus dune line +10.7 m (1996-2010) and -10.7 m (2010-2024); model minus shoreline -6.7 / -11.0 m (projected) and +0.7 / -24.7 m (total). The /10 version is in `output/archive/2026-09-27_target-comparison-div10/smoothed_loess7_with_cascade/`.

Built 2026-09-22 (Hannah, by interview). The observations-only pair in
`data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/smoothed_lowess7/`
with the CASCADE hindcast drawn over it in **dark green**.

Two sheets, 1996–2010 above 2010–2024, one fixed ±100 m axis. They differ from
their observations-only twins in exactly one thing — the green line — because
everything about the two observed curves is imported from that script rather
than re-implemented.

| sheet | shoreline target |
|---|---|
| `lowess7_projected_vs_duneline_with_cascade_...` | the **same** 1996–2024 LRR × 14 yr in both panels |
| `lowess7_total_change_vs_duneline_with_cascade_...` | **each panel's own** LRR × its own 14 yr |

The dune-line target is the same either way and always follows the sub-period.

## Why it is here and not beside its twin

`data/hatteras_init/5-scr/4-comparisons/` is observations only. A figure
carrying model output is a product, and `ORGANIZATION.md` rule 1 puts products
in `output/`. `target_comparison/` already holds this exact kind of figure —
both candidate targets with the hindcast over them — so this is its smoothed,
two-panel sibling.

## The run

**zeroBE, full management, groin off**, one per period:
`HAT_1996_2010_zeroBE_road_bdm_nogroin` and
`HAT_2010_2024_zeroBE_road_bdm_nourish_nogroin` (`runs_used.csv`).

It carries **no source/sink term in any domain**, the two ends included, so
all 90 domains are the model's own response and **neither candidate target was
fitted anywhere in it**. That is what makes it readable against both.

Chosen over the headline edgeBE matrix run for that reason: edgeBE solves GIS
1 and 90 against the CoastSat target, which would leave the model partly
fitted to one of the two things it is being compared against. The script
checks the run name for `zeroBE` before drawing — a caption claiming nothing
was fitted, over a solved run, would be false.

## Interior means, GIS 2–89

| shoreline target | window | model − shoreline | model − dune line |
|---|---|---|---|
| projected | 1996–2010 | −10.8 m | +6.6 m |
| projected | 2010–2024 | −14.8 m | −14.5 m |
| total | 1996–2010 | −3.4 m | +6.6 m |
| total | 2010–2024 | −28.6 m | −14.5 m |

The model-vs-dune-line column is the same in both sheets, as it must be — only
the shoreline target changes between them. What changes is how far the model
sits from the *shoreline* target: in 2010–2024 it is −14.8 m from the
long-term projection but −28.6 m from that window's own much more accretional
rate.

## Reading it

All three curves are smoothed at 7 domains (3.5 km) so no pair is a smoothed
quantity against an unsmoothed one. Two things are stated rather than hidden:

- **GIS 1–10 are unsmoothed on every curve** (the Oregon Inlet boundary
  treatment the scoring target uses), so nothing there is the smoother's.
- The observed targets are smoothed at **transect** resolution (~10 CoastSat
  and 5 dune transects per domain) while the model exists only **per domain**,
  so its LOWESS runs over 90 points rather than ~900. At this width that barely
  changes the fitted curve, but the two are not literally the same operation.

**Correlations are deliberately not quoted.** A symmetric smoother inflates r
on both sides regardless of whether the curves genuinely agree better at this
width — `../smoothing_scale/PROVENANCE.md` established that over 96 rows.
Read the bias.

Rebuild:
`python scripts/analyze_output/compare_runs/smoothed_lowess7_with_cascade.py`
(`--window N` for a different width).
