# hurdat2 — NHC Atlantic best tracks

`hurdat2-1851-2025-091226.txt`, downloaded 2026-09-29 from
https://www.nhc.noaa.gov/data/hurdat/ (the newest file listed that day;
Landsea & Franklin 2013 format: a header line per storm, then 6-hourly fixes
with status, position and maximum wind).

**Used by** `scripts/input_prep/3-env-forcings/3-storms/storm_figures.py` to
type each event of the model's storm series: tropical when a TD/TS/HU/SD/SS
fix lies within 500 km of Cape Hatteras (35.25 N, 75.53 W) inside 24 h of the
event's peak, otherwise "other high-water event" (mostly nor'easters, but not
identified as such). A distant-cyclone / swell class was tried and removed on
2026-09-29. The per-event result is
`3-storms/figures/1996_2024/supporting/storm_types_1996_2024.csv`.

Figures only: nothing the model reads depends on this file.
