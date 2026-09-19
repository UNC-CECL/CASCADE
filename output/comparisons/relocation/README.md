# relocation - does CASCADE relocate NC-12 where and when history did?

Arm A (the roadway module deciding on its own) against arm B (the recorded
relocations prescribed), one set per pair of runs, and the readers that put
those sets side by side. One tree since 2026-09-17; before that the same
material was three top-level folders (`relocation_1984_2004/`,
`relocation_1996_2010/`, `relocation_periods/`) plus
`relocation_standard_setback/`. The numbers and the reading are in
`scripts/hatteras_ms/experiments/RELOCATION_COMPARISON_RESULTS.md`.

```
<start>_<end>/<version>/<preset>[_groin]/   one comparison set: a window, the
                                            dune-TOPOGRAPHY version BOTH arms ran
                                            on (not the island-offset version)
                                            (read from their metadata, never
                                            typed), a source/sink preset,
                                            groin on or off
<start>_<end>/<version>/dune_position_check/
                                            HAT_relocation_dune_position_check.py,
                                            one PNG per preset (1984-2004 only)
events/<year>/<version>/                    one event, two windows, read side
                                            by side
versions/v2_vs_v3/                          one window, two topography versions
standard_setback/                           the 20 m rebuild clearance at GIS 11
```

Every set has the same shape, so a reader finds the same thing in the same
place whichever window, version or preset they open:

```
<set>/
    report.txt              the console output, headed by the identity of the
                            two runs (name, date, topo version, commit)
    tables/                 confusion.csv (recall and false positives per
                            tolerance), first_relocation_year.csv,
                            near_miss_margin.csv, setback_by_year.csv,
                            setback_summary.csv, indexing_check.csv,
                            road_outcomes.csv
    1-island/               one folder per place, two animations in each:
    2-event-1999_GIS9-14/       topography.gif    Barrier3D's grids year by
    3-event-1989_GIS84-87/                        year, NC-12 on them
                                dune-and-road.gif the dune line and the road
                                                  as lines
```

## The windows

| window | version | sets | runs | events in the window |
|---|---|---|---|---|
| `1984_2004/` | `v2/` (CURRENT, the re-pick base) | calibBE_groin (09-09); zeroBE, edgeBE (09-15, for the cross-window report) | `raw_runs/versions/version-pair/v2/` | 1989 Pea Island (GIS 84-87), 1999 (GIS 9-14) |
| `1984_2004/` | `v3/` (the 1984 reconstruction built from v2) | calibBE_groin (09-09) | `raw_runs/versions/version-pair/v3/` | the same two |
| `1996_2010/` | `v2/` (dune topography v2; island offsets v3, the 09-18 re-digitized 1997 line) | zeroBE, edgeBE (rebuilt 09-18 on the re-run arms; the 09-15 sets were on offsets v2) | `raw_runs/matrix/1996_2010/` | 1999 only |

`v1/` (six sets on the original 2026-08-27 pick set, superseded 09-04) was
deleted 2026-09-17; its numbers are in the results document.

**1996-2010 differs from 1984-2004 in four ways.** One event, not two: the
1989 relocation precedes the start and is already applied in the derived
1996 setback file, so there is no `3-event-1989` place and the historical set
is six domains, not ten. Three model years before the event, not fifteen.
The position cross-check is mid-window: the only surveyed road position is
the 2004 one, year 8 of 14, so `setback_summary.csv` carries `check_year`,
`free_at_check_m` and `prescribed_at_check_m` beside `measured_2004_m`.
calibBE is absent: no source/sink field has been solved for 1996-2010, and at
GIS 9-14 zeroBE and edgeBE are the same forcing.

## The readers

**`events/1999/<version>/`** reads `1984_2004/<version>/<preset>/` and
`1996_2010/<version>/<preset>/` for the six domains the 1999 relocation
moved and re-scores nothing: retreat accumulated before 1999 in each window,
whether and when the free arm fired, the 2004 position check in both, and one
figure per preset with both windows' setback trajectories on one axis per
domain. `report.txt`, `tables/event_{domains,recall,outcomes}_<preset>.csv`,
`setback_trajectories_1999_<preset>.png/.pdf`, `CAPTIONS.md`. Written by
`scripts/hatteras_ms/experiments/HAT_relocation_period_compare.py`.

**`versions/v2_vs_v3/`** reads `1984_2004/v2/calibBE_groin/` and
`1984_2004/v3/calibBE_groin/` side by side, the like-for-like pair (same
scenario, code and day; nothing but the topography and its setback CSV
different). `report.txt` and `tables/` are every section of a per-version
report with a v2 column, a v3 column and their difference, plus what v3
changed at the road before the model ran, from the v3 footprint audit
(`HAT_version_pair_report.py`). `<scenario>/<place>/` holds the animations
with v2 in the left panel and v3 in the right, for the island, the two event
blocks, Pea Island (rows added) and Avon to Tri-Village (rows removed)
(`HAT_version_pair_gif.py`). Both scripts live with the reconstruction,
`scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/6-result/`,
and the island-wide v2-against-v3 figure is in that step's data folder.

**`standard_setback/`** is the GIS 11 drowning figure behind the 20 m rebuild
clearance in `HAT_hindcast_config.py`: at 10 m a relocation lands at 97 m of
setback and `bulldoze` drowns the road; at 20 m it lands at 87 m and clears
the wet-row criterion by one cell. `GIS11_profiles.npz` is the per-domain
extract taken from the reloc-arm runs before the 2026-08-28 archive was
deleted; it is the ONE input here, untracked (`.npz`), and regenerates only by
re-running the arms at `relocation_setback_m: measured`. Drawn by
`scripts/figure_making/model_output/HAT_gis11_relocation_drown_figure.py`.

## Regenerate

```
# a set: the output folder follows from the arms, <start>_<end>/<version>/<preset>[_groin]/
python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --preset <preset>                  # 1984-2004 nogroin pair
python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --period 1996 --preset <preset>    # 1996-2010
python scripts/hatteras_ms/experiments/HAT_relocation_comparison.py --preset calibBE \
    --arm-a output/raw_runs/versions/version-pair/v3/1984_2004/calibBE/HAT_1984_2004_calibBE_road_bdm_groin \
    --arm-b output/raw_runs/versions/version-pair/v3/1984_2004/calibBE/HAT_1984_2004_calibBE_road_reloc_bdm_groin
# the readers, after the sets they read exist
python scripts/hatteras_ms/experiments/HAT_relocation_period_compare.py --presets zeroBE edgeBE --version v2
python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/6-result/HAT_version_pair_report.py
python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/6-result/HAT_version_pair_gif.py
```

A pair of arms on different versions is refused. Pass `--out` only to put a
set somewhere else on purpose. The comparison needs each run's `.npz` model
state (the roadway managers' time series live only there), so the runs must
be made with `save_model_state: true`.
