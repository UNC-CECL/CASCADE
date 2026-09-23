# Captions — v1 (built as v2; renumbered 2026-09-19)

`superseded_v1` in the comparison caption is `superseded_20260919_pre-redigitized/v1`; `v1` is this build.

> **Files renamed 2026-09-23 (Hannah: the files in a `v1` folder should say v1).** `offset_2010_v1_vs_v2.*` is now `offset_2010_superseded_v1_vs_v1.*` (CSV columns `model_v1_m/model_v2_m` etc. are now `model_superseded_v1_m/model_v1_m` etc.), and every `_buffer_diagnostic_v2` is now `_buffer_diagnostic_1to1`: that `_v2` was the second drawing of the figure, not a build. The step-output logs below keep the names the scripts printed. The comparison figure was redrawn the same day under its new name, so its legend and caption say `superseded_v1` and `v1`: `HAT_compare_offset_versions.py --year 2010 --a superseded_20260919_pre-redigitized/v1 --label-a superseded_v1 --b v1` with the same two raw files; the CSV came out identical.


Written by the figure scripts through `hat_figure_style.caption()`; the images carry no titles or footnotes, this file does.

**`offset_2010_superseded_v1_vs_v1.png`.** Two builds of the 2010 island offset. (a) The unpadded offset each build hands the model, superseded_v1 in red and v1 in blue, each zeroed on its own most seaward domain. (b) The change in the dune line itself, v1 minus superseded_v1, measured from the shared offshore datum along the 100 m transects and averaged per 500 m domain; positive is landward. 50 of 90 domains differ by 0.5 m or more (mean over all domains -6.7 m; the largest, -65.6 m, at GIS 2). Both raw files were produced by the same shapely intersection, so the metre-scale station convention of the earlier ArcGIS export is not part of the difference.

**`Island_Dune_Offsets_2010_buffer_diagnostic_1to1.png`.** The padded 2010 dune-line offset the model reads for geometry base, build 2010/duneline/v1 (GIS 1 to 90 plus 15 buffer domains each side), at 1:1 scale: one metre alongshore is one metre cross-shore, as a map draws it. The offset is the distance from a north-south datum line east of the island to the dune line, so the slope is the coast's bearing relative to north, not its curvature; the mean bearing is 8 degrees over the reach and the steepest domain-to-domain angle is 26 degrees. The compressed planform every calibrated run uses divides these offsets by ten.
