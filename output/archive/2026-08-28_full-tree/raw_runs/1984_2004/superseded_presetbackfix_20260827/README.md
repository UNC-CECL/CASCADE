# SUPERSEDED — 1984-2004 zeroBE runs made BEFORE the road setbacks were re-measured

Archived 2026-08-27. **Do not use these in analysis.** Kept because they are
the only directly comparable "before" for the setback change, and unlike the
sibling archive they still have their model states.

## These ARE on the 1996-mosaic topography

Unlike `superseded_pre1996mosaic_20260827/`, these six runs started from
`1-barrier3d-domains/1984-start/dune-topo/v1` — the correct product for this
period. The topography is not what makes them stale.

## Why they are superseded

They were run 2026-08-26 14:30-14:45. `RoadSetback_1984_dunestart.csv` was
regenerated at **16:54 that same day** — the first version measured against
`1984-start` interior row 0 rather than `2004-start`'s. These six consumed the
older file.

The change lands squarely on the domains the relocation analysis scores
(1989 event: GIS 84-87; 1999 event: GIS 9-14):

```
GIS      9    10    11    12    13    14  ...   84    85    86    87
old     40     0     0     0     0    40         0     0     0    30
new     40    30    30    40    40    60        30     0     0    60
```

This is not a cosmetic shift. `HAT_relocation_comparison.py`'s headline caveat
was that six of the ten historical domains start at a setback of 0.0 m, so a
relocation puts the road back where it already was and the trigger re-fires the
next year the dune line moves landward. Under the new file **only two of ten**
start at zero. Relocation timing, the ratcheting behaviour, and the hit/miss
recall all move.

`RoadElevation.csv` was also regenerated (17:00), but only GIS 78 and 79
changed, by ~0.01 m.

## v1 vs v2 topography is NOT why these are archived

`1984-start` has since advanced to v2 (CURRENT and the extractor agree). v2 is
v1 plus 114 unsurveyed cells bridged by linear interpolation in 22 holes across
GIS 4-7. Per `dune-topo/v2/BRIDGE_MANIFEST.txt`, interior shapes and interior
row 0 are **identical to v1**, and GIS 4-7 carry no road — so setbacks measured
on v1 remain valid on v2, and this archive differs from the replacement runs
essentially only in the setback file.

## What is in here

Six zeroBE, groin-off runs — the four scenario cells plus the two relocation
arms:

    HAT_1984_2004_zeroBE_noroad_nobdm_nogroin      (natural)
    HAT_1984_2004_zeroBE_noroad_bdm_nogroin        (beachdune_only)
    HAT_1984_2004_zeroBE_road_nobdm_nogroin        (roadway_only)
    HAT_1984_2004_zeroBE_road_bdm_nogroin          (full_management)       <- relocation arm A
    HAT_1984_2004_zeroBE_road_reloc_nobdm_nogroin  (roadway_only + reloc)
    HAT_1984_2004_zeroBE_road_reloc_bdm_nogroin    (full_management + reloc) <- relocation arm B

`run_index_snapshot.csv` holds their six index rows as they stood before
archiving; they have been removed from the live `run_index.csv`.

## The model states are INTACT

Each directory still holds its ~300 MB `.npz`, so this archive can be fed
straight back into the comparison to reproduce the pre-setback-fix answer:

    python scripts/hatteras_ms/HAT_relocation_comparison.py \
      --arm-a output/raw_runs/1984_2004/superseded_presetbackfix_20260827/zeroBE/HAT_1984_2004_zeroBE_road_bdm_nogroin \
      --arm-b output/raw_runs/1984_2004/superseded_presetbackfix_20260827/zeroBE/HAT_1984_2004_zeroBE_road_reloc_bdm_nogroin \
      --out output/comparisons/relocation_1984_2004/zeroBE_presetbackfix

Both arms must come from the SAME archive — the script derives every number as
a difference between the arms, so pairing an archived arm with a current one
would compare two forcings and report it as a relocation result.

## Note

`git_dirty = True` on all six: the working tree carried uncommitted road-asset
changes when they ran, so the recorded commit does not reproduce them.
