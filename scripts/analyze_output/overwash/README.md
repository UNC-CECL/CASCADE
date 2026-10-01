# overwash - a run's overwash against the imagery record

One script: `compare_overwash_figures.py` draws a CASCADE run's overwash flux
(Qow) against the overwash mapped from aerial imagery. Nothing here runs the
model.

| figure | file in `output/comparisons/overwash/` |
|---|---|
| 1, stacked: model above, observed below | `compare_stacked.png` |
| 2, contingency: hit / miss / false alarm per year and domain | `compare_contingency_T<threshold>.png` |
| 3, spatial profile, one panel | `spatial_profile_normalised.png` |
| 4, spatial profile, two panels | `spatial_profile_dual_panel.png` |

**Is this the live comparison?** Not for current runs. It reads the archived
1984-2004 calibBE run. The run-vs-imagery comparison for the current managed
runs (1996-2010 and 2010-2024, each image against the storms since the
previous one) is
`scripts/input_prep/8-overwash-analysis/4-vs-model/overwash_vs_model.py`.
The observed record itself, and its own figures, live in
`scripts/input_prep/8-overwash-analysis/`.

## compare_overwash_figures.py


It reads one run's `.npz` and the observation workbook (`Overwash_Matrix`
sheet, resolved by `site_layer/hat_overwash.py`) and writes to
`output/comparisons/overwash/`.

| figure | switch | shows |
|---|---|---|
| 1, stacked | `PLOT_STACKED` | model Qow heatmap above (every model year, imagery years highlighted), observed overwash below (binary, imagery years; hatched where there is no image); shared domain axis and island-section bar |
| 2, contingency | `PLOT_CONTINGENCY` | every model year x domain classified hit / miss / false alarm / correct rejection, years without imagery greyed; summary statistics (POD, CSI) printed |
| 3 and 4, spatial | `PLOT_SPATIAL` | overwash frequency alongshore, model (mean with IQR shading) against observed (bars): normalised on one panel, and as two panels |

**The run:** `HAT_1984_2004_calibBE_road_bdm_groin`, 1984-2004, calibBE,
resolved through `cascade_pipeline.run_registry`. *Fixed 2026-09-30:* the run
moved to `raw_runs/archive/2026-09-24-pre-metres/` when the matrix was
archived on 2026-09-24, and the scripts crashed on import until they were
pointed at it there (`RUN_KIND`, `RUN_TAG`). Whether to compare a current
matrix run instead is an open choice.

**Threshold:** a model cell counts as overwash if Qow > `QOW_THRESHOLD`
(dam³/yr). Start at 0 (any non-zero flux) and raise it if there are too many
false alarms.

**Imagery-year remap:** some imagery post-dates the storm it
captures, so `YEAR_REMAP` maps an image year to the model year that best
represents the event visible in it. The May 2004 imagery shows Hurricane
Isabel's (September 2003) deposits, so 2004 is remapped to 2003. 1996 imagery
is flagged as poor quality (`POOR_QUALITY_YEARS`).

**Palette:** muted coastal earth and ocean tones, with
consistent warmth across all elements; the section bar uses desaturated
versions of the warm and cool families in the contingency legend. Section
bar: warm linen #CABB9E for villages, dusty maritime blue #8DAFC2 between
them. Contingency: hit deep maritime green #2E7D5A, miss deep maritime blue
#2A5F8F, false alarm warm brick #A85840, correct rejection warm cream #F2EFE9,
no imagery warm greige #BEB9B4. Observed overwash uses the same brick as a
false alarm. The look was modelled on AGU / Nature Geoscience figures.

**To use on another run:** set `RUN_NAME`, `RUN_PERIOD`, `RUN_PRESET` (and
`RUN_KIND`/`RUN_TAG` if it is not in the matrix), then adjust `START_YEAR`,
`END_YEAR`, the domain constants and `SECTIONS` to match it, and choose
`QOW_THRESHOLD`. If the run is not where you say, `find_run_dir` raises naming
where it *is* on disk.

**History:** the paths used to be absolute literals under a folder spelling
that no longer exists (`input_preperation`), so the scripts could not read
their observations and crashed on the first save; they are anchored on the
repo root now. The workbook left `scripts/input_prep/8-overwash-analysis/`
on 2026-09-10 and is resolved by `site_layer/hat_overwash.py` since
2026-09-18. The `.npz` loader substitutes stand-in classes for anything the
pickled model references but cannot import, and is duplicated from
`plot_overwash.py` (deleted 2026-10-01, see below) so the script stands alone.

**Contingency classes** (threshold T = `QOW_THRESHOLD`):

| class | rule | colour |
|---|---|---|
| hit | Qow > T and observed | deep maritime green |
| miss | Qow <= T and observed | deep maritime blue |
| false alarm | Qow > T and not observed | warm brick |
| correct rejection | Qow <= T and not observed | warm cream |
| not assessed | no image, or the cell is blank in the workbook | warm greige |

<details><summary>Function notes (the original docstrings)</summary>

**`compute_spatial_data()`**

```
Returns per-domain summary arrays used by both spatial profile figures.

qow_mean  : mean annual Qow per domain across all model years (dam³/yr)
qow_p25   : 25th percentile of annual Qow per domain
qow_p75   : 75th percentile of annual Qow per domain
obs_freq  : fraction of assessed imagery years with overwash per domain
            (NaN where no assessed imagery exists for that domain)
bar_colors: colour per domain based on section type (village vs inter-village)
```

</details>

## Deleted 2026-10-01

Hannah chose to keep one script and delete the rest outright, since git keeps
them. This departs from ORGANIZATION.md rule 4 (retired code is parked, not
deleted), as the 5-scr reorganization did on 2026-09-22.

| deleted | why |
|---|---|
| `compare_overwash_observed.py` | made the same stacked and contingency figures, writing the SAME file names into the same folder, so whichever ran last won. It had no imagery-year remap: with model years 1984-2003 it silently dropped the May 2004 image, which is where Hurricane Isabel (Sep 2003) shows. Its contingency figure kept imagery years only. Nothing it drew is missing from `compare_overwash_figures.py`. |
| `superseded_20260918/plot_overwash.py` | the first overwash plot of a hindcast run, read from an absolute path into a run (`HAT_1984_2004_basestorms_Hs2p0`) that no longer exists. Its `.npz` loader lives on, copied, in `compare_overwash_figures.py`. |
| `superseded_20260918/overwash_pea_early.py` | an early copy of Roya's Pea Island overwash script (was `file_from_roya/overwash_pea.py`). The maintained copy is `scripts/other_ms/pea_island_ms/overwash_pea.py`. |
| `superseded_20260918/WHY.md` | the note for the two above; its content is this table. |

To recover any of them:

```
git log --diff-filter=D --oneline -- scripts/analyze_output/overwash/<path>
git show <commit>^:scripts/analyze_output/overwash/<path>
```

The `compare_stacked.png` already in `output/comparisons/overwash/` may have
been drawn by either script; re-run `compare_overwash_figures.py` to be sure.
