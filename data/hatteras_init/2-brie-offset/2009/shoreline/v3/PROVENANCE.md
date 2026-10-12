# 2009 island offsets (shoreline), v3

> **Test build, not CURRENT (2026-10-09).** Built for `output/raw_runs/experiments/target-end-window/2026-10-09-6mo-starts/`, which runs the +/-6-month start against the +/-1-yr one (v2). `SHORELINE_WINDOW_FOR_YEAR` was switched to the +/-6-month window for the build only and then restored, so the table still pairs this start with the v2 window. The mean shoreline was built with **min 5 positions per transect, not 10**, because the 6-month window is thin (at 10, 570 transects and 69 domains; GIS 1-18 and 20-22 missing). Runs read this build only with `HAT_OFFSET_VERSION_2009_SHORELINE=v3`.

Built 2026-10-09 by `scripts/input_prep/2-brie-offset/1-produce/build_island_offset.py --year 2009 --source shoreline --no-current` from `shoreline_mean_2009-02-17_2010-02-17.geojson` (the CoastSat mean over 2009-02-17 to 2010-02-17, the 2009 USACE midpoint 2009-08-17 +/-6 months; 905 of 906 transects, median 11 positions, transect SE median 4.5 m).

**The driver stopped at its comparison step**, so this file was written by hand and `../../PROVENANCE.md` got no Builds row from it. `HAT_compare_offset_versions.py` reads the column named for the year, and v2's unpadded file still calls it `2010` (built before the 10-05 key rename): `KeyError: '2009'`. Steps 1 and 2 had finished; the files below are complete.

## The intersection (step 1)

`step1_line_to_raw_offsets.py` against `transects/transects_100m.geojson`: 450 transects in GIS 1-90, 0 with no crossing, 0 crossed more than once. Raw: `raw_offsets/2009-02-17_2010-02-17_shoreline_offset_raw.csv`, copied here.

## The model input (step 2)

`step2_raw_offsets_to_island_offset.py --year 2009 --version v3 --source shoreline`: zeroed on GIS 76 (the same zero domain as v2), padded to 120 (real span 0-6145.1 m; largest angle 26.3 deg real, 27.5 deg buffer). Files: `Island_Shoreline_Offsets_2009_PADDED_120.csv` (read by the model), `_CASCADE_Input.csv`, `_CASCADE_Input_unpadded.csv`, `_buffer_diagnostic.png`.

## Against v2 (by hand)

Model frame (each zeroed on its own minimum), v3 minus v2: mean +2.96 m, range -3.31 to +9.87 m, 19 domains beyond 5 m. In this frame a positive difference is the line moved seaward relative to GIS 76.

## Rebuild

```
python scripts/input_prep/2-brie-offset/1-produce/step2_raw_offsets_to_island_offset.py --year 2009 --version v3 --source shoreline --raw-file data/hatteras_init/2-brie-offset/2009/shoreline/v3/2009-02-17_2010-02-17_shoreline_offset_raw.csv
```
