# NC-12 road setback audit

Generated 2026-08-28T11:19:29 by `HAT_road_setback_audit.py`.

**Diagnostic record.** No setback was modified, no setback file was written, no topography was touched.

| | |
|---|---|
| Topography (1984-start) | `C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\1-barrier3d-domains\1984-start\dune-topo\v1\topography` (version v1) |
| Topography (2004-start) | `C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\1-barrier3d-domains\2004-start\dune-topo\v1\topography` (version v1) |
| Road width | 20 m |
| Width rule | `barrier3d.FindWidths` — first cell `<= SL`, verbatim |
| Drown rule | `roadway_manager.bulldoze` — >20% of profiles bordering water |
| Off-island rule | ≥20% of the road's **own** cells in water |

The runner uses **one topography for both hindcast periods**, so the 1984 setbacks and the 2004 setbacks are both spent against a 2009 grid. For 1984 that is a 25-year anachronism, and it is the reason this audit exists.

## What `DROWNS` actually costs

`roadway_drown` is not a warning. Traced through the code, one drown does this:

1. **`bulldoze()` mutates `barrier3d.InteriorDomain` in place before it checks anything.** `xyz_interior_grid[road_start:road_end, :] = new_road_domain` — the grid is passed by reference, so the write lands even though the function returns early afterwards. On GIS 52 at the 1984 setback, 53 of 100 cells in the road rows were water (`-0.3 dam`) and all 100 came out at `0.145 dam`: the model gains a 20 m ribbon of 1.45 m land across open water.
2. `RoadwayManager` sets `_drown_break = 1` and returns.
3. `cascade.py` (~line 625) sees `drown_break` on **every later year** and never calls `update()` again. `_road_break[iB3D] = 1`, dune growth rates reset to natural.

So from that year on there is **no overwash removal, no dune rebuilding, and no road relocation**, and `_road_overwash_volume` / `_dunes_rebuilt_TS` / `_rebuild_dune_volume_TS` stay at zero for the rest of the run.

> A domain flagged `DROWNS` is an **unmanaged barrier wearing a road label** for the whole hindcast. If any result contrasts managed against unmanaged shoreline, those domains are silently in the unmanaged group from year zero.

`group_roadway_abandonment` is `None` in the runner, so this stays per-domain. Set it, and one drowned domain abandons its entire group.

## How the relocation events were corrected

The events used to carry `post_relocation_setback_m`, an **absolute** setback built as `(1984 setback) + (1978→1997 road displacement)`. That double-counted the shoreline retreat: the 1984 setback is measured from the **1984 dune line**, but it is spent against a grid whose row 0 is the **2009** dune crest, which already sits landward by the 1984→2009 retreat *R*. Adding the physical displacement on top counted *R* twice, and the road was placed *R* metres behind where NCDOT actually put it.

Checked against `RoadSetback_2004.csv` — a same-year (2004 road vs 2004 dune) measurement taken **after both events**, so it already is the post-relocation position:

| GIS | Old absolute | 2004 measured | Overshoot |
|---:|---:|---:|---:|
| 11 | 129 | 81 | +48 |
| 12 | 126 | 89 | +37 |
| 13 | 125 | 87 | +38 |
| 14 | 126 | 93 | +33 |
| 15 | 106 | 71 | +35 |
| 84 | 163 | 50 | +113 |
| 85 | 165 | 85 | +80 |
| 86 | 205 | 88 | +117 |
| 87 | 113 | 40 | +73 |

Median overshoot **+35 m** in the 1999 block, **+96 m** at Pea Island — worst exactly where retreat was worst. At GIS 84 the measured 2004 setback (50 m) is *smaller* than the 1984 setback (93 m) even though the road was moved 70 m landward in 1989: the shoreline overtook the road faster than NCDOT could move it. That is the real Pea Island story, and the absolute value erased it.

**The fix** (`HAT_hindcast_1984_2024.py`, `HISTORICAL_ROAD_EVENTS`): the events now carry `relocation_displacement_m` — how far the road physically moved, and nothing else — applied to the setback the model is *currently* carrying:

```python
new_sb = rm._road_setback + relocation_displacement_m[gis_id]
```

CASCADE already decrements `road_setback` by `dune_migrated` every year, so by the event year the model's setback has absorbed the modelled retreat on its own. The retreat is now counted exactly once, from the model's own dune migration.

### What the correction does and does not fix

The corrected setback depends on how much the model has retreated by the event year. Evaluated across plausible retreat *R*:

- **1999 block (GIS 9–15, event at year 15)** converges toward the 2004 measurement at realistic rates. GIS 11 goes 129 → 89 m at *R* = 40 m, against a measured 81 m. The correction largely works here.
- **1989 block (GIS 84–87, event at year 5)** does **not** converge. Only a third of the retreat has accumulated by year 5, so GIS 86 is still ~178 m against a measured 88 m.

The reason is structural, and the correction cannot reach it: the road is placed **correctly relative to the dune**, but the island *behind* that dune is the 2009 island — 25 years too narrow. At Pea Island, where retreat is fastest on the whole island, that width deficit is larger than the relocation itself. Fixing the double-count removes the avoidable error; the remaining error is the width anachronism, and it needs a 1984 DEM that does not exist.

## The road's own cells are never checked

`bulldoze()` tests the rows **flanking** the road — `road_end + 1` on the bay side, `road_start - 1` on the sea side — and never looks at the road footprint itself. It also skips `road_end`: the road occupies `[road_start:road_end]` and the bayside test is at `road_end + 1`, so the cell immediately behind the asphalt is only used for `np.size()`.

That is why `pct_road_cells_water` is reported separately below. A road can be substantially in the bay while the flanking test is still under threshold.

## 1984 initial

82 domains — `FITS` 82. **0 with the road off the island.**

No domain flagged.

## 2004 initial

82 domains — `FITS` 82. **0 with the road off the island.**

No domain flagged.

## 1989 relocation (GIS 84-87) - zero-retreat bound

4 domains — `FITS` 4. **0 with the road off the island.**

No domain flagged.

### Guards the prescribed-relocation path skips

CASCADE has two guards on moving a road landward, and the runner's `HISTORICAL_ROAD_EVENTS` path takes **neither** — it assigns `rm._road_setback = new_sb` directly, leaves `road_ele` at its SLR-decremented value rather than rebuilding at grade, and lets the next `bulldoze()` flatten whatever is there.

| GIS | Setback (m) | Rebuilt road elevation (m) | `get_road_relocation_elevation` | Average barrier width (m) | `road_relocation_checks` |
|---:|---:|---:|---|---:|---|
| 84 | 65 | 1.38 | allows | 263 | allows |
| 85 | 109 | 1.06 | allows | 361 | allows |
| 86 | 94 | 1.08 | allows | 286 | allows |
| 87 | 78 | 1.21 | allows | 205 | allows |

**0 of 4 prescribed relocations would be refused by CASCADE's own logic.** Evaluated on the t=0 interior; the real grid at the event year will have evolved, so treat these as the initialisation-time verdict.

## 1999 relocation (GIS 9-14) - zero-retreat bound

6 domains — `FITS` 6. **0 with the road off the island.**

No domain flagged.

### Guards the prescribed-relocation path skips

CASCADE has two guards on moving a road landward, and the runner's `HISTORICAL_ROAD_EVENTS` path takes **neither** — it assigns `rm._road_setback = new_sb` directly, leaves `road_ele` at its SLR-decremented value rather than rebuilding at grade, and lets the next `bulldoze()` flatten whatever is there.

| GIS | Setback (m) | Rebuilt road elevation (m) | `get_road_relocation_elevation` | Average barrier width (m) | `road_relocation_checks` |
|---:|---:|---:|---|---:|---|
| 9 | 58 | 1.50 | allows | 205 | allows |
| 10 | 67 | 1.44 | allows | 172 | allows |
| 11 | 87 | 1.23 | allows | 161 | allows |
| 12 | 88 | 1.39 | allows | 165 | allows |
| 13 | 81 | 1.46 | allows | 179 | allows |
| 14 | 85 | 1.52 | allows | 166 | allows |

**0 of 6 prescribed relocations would be refused by CASCADE's own logic.** Evaluated on the t=0 interior; the real grid at the event year will have evolved, so treat these as the initialisation-time verdict.

## Where Barrier3D's two guards disagree

The model does not use one definition of "island":

| Routine | What it reads | Used for |
|---|---|---|
| `FindWidths` | Stops at the **first** water cell — land behind a gap is not part of the island | `InteriorWidth_Avg`, i.e. the `average_barrier_width` that `road_relocation_checks()` uses |
| `bulldoze` | The **literal cell** at `road_end + 1` — land behind a gap counts normally | The `roadway_drown` test |

So the road can be founded on ground the relocation logic believes is off the back of the barrier: the drown test passes while the relocation test sees a much narrower island. Counted, not corrected — resolving it would mean changing a Barrier3D definition.

| Scenario | GIS | % of profiles where the road's border cell is land but `FindWidths` has already ended the island | Narrowness |
|---|---:|---:|---|
| 1984 initial | 73 | 22.0 | `HOLED` |
| 2004 initial | 73 | 22.0 | `HOLED` |
| 1984 initial | 74 | 20.0 | `HOLED` |
| 2004 initial | 74 | 20.0 | `HOLED` |
| 1984 initial | 77 | 18.0 | `MIXED` |
| 2004 initial | 77 | 18.0 | `MIXED` |
| 1984 initial | 78 | 12.0 | `MIXED` |
| 1984 initial | 71 | 10.0 | `HOLED` |
| 2004 initial | 78 | 10.0 | `MIXED` |
| 1984 initial | 68 | 8.0 | `MIXED` |
| 1984 initial | 79 | 8.0 | `TRUE_END` |
| 2004 initial | 26 | 8.0 | `HOLED` |
| 2004 initial | 68 | 8.0 | `MIXED` |
| 2004 initial | 71 | 8.0 | `HOLED` |
| 2004 initial | 79 | 8.0 | `TRUE_END` |
| 1984 initial | 19 | 4.0 | `HOLED` |
| 1984 initial | 26 | 4.0 | `HOLED` |
| 1984 initial | 35 | 4.0 | `MIXED` |
| 2004 initial | 19 | 4.0 | `HOLED` |
| 1984 initial | 67 | 2.0 | `MIXED` |
| 1984 initial | 81 | 2.0 | `HOLED` |
| 2004 initial | 35 | 2.0 | `MIXED` |

## Narrowness labels, and what they do not tell you

The source DEM was pre-masked to the island: across all 131 domain arrays there are **822,315 cells of exactly `-10.0` (NoData) and zero real cells below the water clamp**. There is no bathymetry, so a water cell and a LiDAR dropout are the same number and cannot be distinguished by value.

| Label | Meaning |
|---|---|
| `TRUE_END` | Nothing behind the first water cell. The interior really ends where `FindWidths` says. |
| `HOLED` | Land resumes in >50% of profiles, usually after a 1–3 cell gap. The reported width is a truncation. Whether that gap is a NoData hole or a real marsh channel **is not resolved here** — Hatteras has both at 10–30 m scale. |
| `MIXED` | Between the two. |

Barrier3D stops at the first water cell either way, so the reported `InteriorWidth` **is** the width the model uses. The label marks which domains carry the doubt; no data was changed to resolve it.

## Scope of any road-management claim

Some domains put the road off the island even at a correctly measured same-year setback. The 2009 barrier genuinely cannot hold NC-12 where it was, and no change to the setback fixes that — only a contemporaneous DEM would.

These are left **unmodified and allowed to drown**. The consequence is the one described at the top: road management stops permanently for them. So any managed-versus-unmanaged conclusion must **exclude the domains listed in the per-scenario tables above**, by name, rather than averaging them in — they are unmanaged barriers from the year they drown, and including them biases every management statistic toward the unmanaged case.

### Period weighting

The runner uses one topography (2009) for both hindcast periods, so 1984–2024 runs its first 25 years on end-of-period geometry while 2004–2024 brackets its topography. **Both periods are reported as primary results**, with this limitation stated once in the methods.

What that means concretely, and what belongs in the caveat: in period 1 the road is placed correctly *relative to the dune* but the island behind that dune is 25 years too narrow, so the room available for management is systematically understated — most severely at Pea Island, where retreat is fastest. Period 2 does not carry that deficit to the same degree, and its setbacks are same-year measurements.

## Unguarded index in `bulldoze()`

`roadway_manager.bulldoze()` checks the two rows flanking the road. The **seaside** check is guarded:

```python
if road_start > 0:
    seaside_water_cells = ... xyz_interior_grid[road_start - 1, :] ...
else:
    seaside_water_cells = 0
```

The **bayside** check is not:

```python
bayside_water_cells = np.count_nonzero(
    (xyz_interior_grid[road_end + 1, :] * dz) <= drown_threshold
) / number_border_cells          # IndexError if road_end + 1 >= DomainWidth
```

The asymmetry reads as an oversight rather than a design choice: one edge of the road is bounds-checked and the other is not. When the barrier has narrowed past the road, the bayside index runs off the array and the run dies with `IndexError` instead of reporting a drowned road.

**Status: not patched.** CASCADE is run as published. This project previously carried a monkey-patch in the run script that supplied the missing guard (`bayside_water_cells = 1.0` when off-grid, which trips the drown test at any threshold); it has been removed so the model matches the package. The exposure is managed upstream of the model instead — this audit refuses any setback whose border row would exceed `DomainWidth` (`BLOCK_OVERRUN`), which is currently zero across all scenarios.

**Residual risk.** This audit tests **t = 0 only**. The interior is not fixed-size: `barrier3d.py` grows it (`np.vstack`, :648) and shrinks it (`np.delete(..., 0, axis=0)`, :724) as the barrier migrates. Erosion alone should be safe, because `road_setback` is decremented by `dune_migrated` in the same step, so the road and the grid shift together. The exposed path is a **relocation**, which resets `road_setback` to a value that does not account for accumulated narrowing — which is exactly the case that originally motivated the patch (GIS 11, 1999 event). That path is now defended by correcting the relocation to a displacement and by warning when CASCADE's own relocation guards object, but it is not proven impossible.

## Known assumptions

1. Road and dune lines are the same vintage in each setback file, so the shoreline retreat cancels in the subtraction. Verified for 1984 and 2004.
2. A setback is transferable between dune definitions: measured to the digitised dune line of its own year, spent against interior row 0 = the 2009 DEM `argmax` dune crest. **Unverified**, and the largest single assumption here.
3. Interior row 0 of profile *i* is `dune_loc[i] + 1` (`USE_CONST_INTERIOR = False`), so the dune search window sets both the setback origin and how much land remains behind it.
4. Cross-shore distance carries a 1/cos θ inflation (1% at 8°, 11% at 26°); interior width carries the same factor, so the ratio is consistent even though absolute distances are long.
5. `bulldoze` compares against `drown_threshold = 0` described as *m MSL*, but the arrays are **MHW-relative** (`MHW_M = 0.36` subtracted first). The effective test is "at or below MHW". Pre-existing in the model; reported, not changed.
6. `ROAD_ELEVATION = 1.45` in the runner is ambiguous between NAVD88 and MHW-relative — flagged in `write_elevation_csv()`, still unresolved. If NAVD88 was meant, it should be 1.09 MHW-relative.
7. The raw `1984_road_offset_raw.csv` behind `RoadSetback_1984.csv` is **not present in the repo**, so those values cannot be re-derived or verified at source. Audited as given.
8. Relocation scenarios are evaluated on the **t=0** interior. The real grid in 1989/1999 will have evolved under storms and SLR, so those verdicts are initialisation-time, not event-time.

## Full results

| Scenario | GIS | Setback (m) | Verdict | Off-island | % road cells water | % bayside water | DomainWidth | IW min/p20/med/max | Largest fitting (m) | Narrowness |
|---|---:|---:|---|---|---:|---:|---:|---|---:|---|
| 1984 initial | 9 | 40 | `FITS` |  | 0.0 | 0.0 | 163 | 11/14/18/73 | 1540 | `HOLED` |
| 1984 initial | 10 | 20 | `FITS` |  | 0.0 | 0.0 | 154 | 14/15/17/22 | 130 | `HOLED` |
| 1984 initial | 11 | 10 | `FITS` |  | 0.0 | 0.0 | 163 | 14/14/16/22 | 110 | `HOLED` |
| 1984 initial | 12 | 20 | `FITS` |  | 0.0 | 0.0 | 170 | 13/14/17/21 | 110 | `TRUE_END` |
| 1984 initial | 13 | 30 | `FITS` |  | 0.0 | 0.0 | 170 | 15/17/18/20 | 140 | `TRUE_END` |
| 1984 initial | 14 | 60 | `FITS` |  | 0.0 | 0.0 | 35 | 10/13/16/25 | 100 | `TRUE_END` |
| 1984 initial | 15 | 60 | `FITS` |  | 0.0 | 0.0 | 32 | 12/13/14/19 | 100 | `TRUE_END` |
| 1984 initial | 16 | 40 | `FITS` |  | 0.0 | 0.0 | 63 | 13/17/20/23 | 140 | `MIXED` |
| 1984 initial | 17 | 90 | `FITS` |  | 0.0 | 0.0 | 71 | 18/22/28/47 | 210 | `HOLED` |
| 1984 initial | 18 | 120 | `FITS` |  | 0.0 | 0.0 | 87 | 16/29/35/54 | 390 | `HOLED` |
| 1984 initial | 19 | 130 | `FITS` |  | 0.0 | 2.0 | 97 | 8/18/26/47 | 270 | `HOLED` |
| 1984 initial | 20 | 60 | `FITS` |  | 0.0 | 0.0 | 103 | 28/32/43/51 | 350 | `MIXED` |
| 1984 initial | 21 | 130 | `FITS` |  | 0.0 | 0.0 | 114 | 26/42/45/71 | 600 | `HOLED` |
| 1984 initial | 22 | 190 | `FITS` |  | 0.0 | 0.0 | 128 | 34/38/43/100 | 880 | `HOLED` |
| 1984 initial | 23 | 230 | `FITS` |  | 0.0 | 0.0 | 134 | 43/51/68/94 | 900 | `HOLED` |
| 1984 initial | 24 | 210 | `FITS` |  | 0.0 | 0.0 | 141 | 51/59/90/111 | 930 | `HOLED` |
| 1984 initial | 25 | 200 | `FITS` |  | 0.0 | 0.0 | 138 | 42/52/80/116 | 1050 | `HOLED` |
| 1984 initial | 26 | 230 | `FITS` |  | 0.0 | 0.0 | 148 | 18/47/71/117 | 1090 | `HOLED` |
| 1984 initial | 27 | 260 | `FITS` |  | 0.0 | 0.0 | 161 | 40/71/81/121 | 1080 | `HOLED` |
| 1984 initial | 28 | 270 | `FITS` |  | 0.0 | 0.0 | 168 | 46/48/79/114 | 1150 | `HOLED` |
| 1984 initial | 29 | 260 | `FITS` |  | 0.0 | 0.0 | 154 | 37/47/87/121 | 1010 | `HOLED` |
| 1984 initial | 30 | 250 | `FITS` |  | 0.0 | 0.0 | 140 | 29/50/59/81 | 720 | `HOLED` |
| 1984 initial | 31 | 230 | `FITS` |  | 0.0 | 0.0 | 139 | 30/32/41/65 | 340 | `HOLED` |
| 1984 initial | 32 | 210 | `FITS` |  | 0.0 | 0.0 | 141 | 32/44/46/57 | 420 | `TRUE_END` |
| 1984 initial | 33 | 250 | `FITS` |  | 0.0 | 0.0 | 148 | 30/43/46/52 | 430 | `MIXED` |
| 1984 initial | 34 | 310 | `FITS` |  | 0.0 | 0.0 | 149 | 41/43/46/50 | 400 | `TRUE_END` |
| 1984 initial | 35 | 400 | `FITS` |  | 2.0 | 0.0 | 137 | 39/49/51/58 | 470 | `MIXED` |
| 1984 initial | 36 | 370 | `FITS` |  | 0.0 | 0.0 | 136 | 56/57/58/63 | 540 | `TRUE_END` |
| 1984 initial | 37 | 335 | `FITS` |  | 0.0 | 0.0 | 159 | 55/58/80/89 | 740 | `HOLED` |
| 1984 initial | 38 | 300 | `FITS` |  | 0.0 | 0.0 | 166 | 70/80/83/89 | 790 | `TRUE_END` |
| 1984 initial | 39 | 230 | `FITS` |  | 0.0 | 0.0 | 170 | 68/69/71/78 | 660 | `TRUE_END` |
| 1984 initial | 40 | 210 | `FITS` |  | 0.0 | 0.0 | 172 | 54/63/71/77 | 640 | `MIXED` |
| 1984 initial | 41 | 200 | `FITS` |  | 0.0 | 0.0 | 174 | 43/52/55/59 | 500 | `HOLED` |
| 1984 initial | 42 | 180 | `FITS` |  | 0.0 | 0.0 | 180 | 42/43/60/67 | 560 | `HOLED` |
| 1984 initial | 43 | 170 | `FITS` |  | 0.0 | 0.0 | 179 | 54/62/64/66 | 600 | `TRUE_END` |
| 1984 initial | 44 | 160 | `FITS` |  | 0.0 | 0.0 | 169 | 46/50/53/66 | 530 | `HOLED` |
| 1984 initial | 45 | 155 | `FITS` |  | 0.0 | 0.0 | 155 | 32/39/42/57 | 380 | `HOLED` |
| 1984 initial | 46 | 150 | `FITS` |  | 0.0 | 0.0 | 142 | 35/39/40/49 | 360 | `HOLED` |
| 1984 initial | 47 | 150 | `FITS` |  | 0.0 | 0.0 | 127 | 30/32/35/42 | 300 | `HOLED` |
| 1984 initial | 48 | 150 | `FITS` |  | 0.0 | 0.0 | 118 | 25/27/29/31 | 240 | `TRUE_END` |
| 1984 initial | 49 | 130 | `FITS` |  | 0.0 | 0.0 | 98 | 20/23/25/28 | 200 | `TRUE_END` |
| 1984 initial | 50 | 150 | `FITS` |  | 0.0 | 0.0 | 88 | 18/23/24/29 | 200 | `TRUE_END` |
| 1984 initial | 51 | 160 | `FITS` |  | 0.0 | 0.0 | 58 | 19/25/27/34 | 230 | `TRUE_END` |
| 1984 initial | 52 | 170 | `FITS` |  | 0.0 | 0.0 | 66 | 22/26/28/32 | 240 | `TRUE_END` |
| 1984 initial | 53 | 170 | `FITS` |  | 0.0 | 0.0 | 56 | 24/27/29/38 | 240 | `MIXED` |
| 1984 initial | 54 | 160 | `FITS` |  | 0.0 | 0.0 | 114 | 25/27/31/46 | 250 | `HOLED` |
| 1984 initial | 55 | 170 | `FITS` |  | 0.0 | 0.0 | 107 | 24/34/37/49 | 340 | `HOLED` |
| 1984 initial | 56 | 180 | `FITS` |  | 0.0 | 0.0 | 77 | 24/28/30/34 | 250 | `MIXED` |
| 1984 initial | 57 | 210 | `FITS` |  | 0.0 | 0.0 | 68 | 25/30/37/45 | 270 | `TRUE_END` |
| 1984 initial | 58 | 240 | `FITS` |  | 0.0 | 0.0 | 69 | 32/41/45/53 | 400 | `MIXED` |
| 1984 initial | 59 | 270 | `FITS` |  | 0.0 | 0.0 | 78 | 42/44/48/61 | 440 | `HOLED` |
| 1984 initial | 60 | 310 | `FITS` |  | 0.0 | 0.0 | 102 | 47/56/62/79 | 590 | `HOLED` |
| 1984 initial | 61 | 330 | `FITS` |  | 0.0 | 0.0 | 116 | 47/56/65/91 | 900 | `HOLED` |
| 1984 initial | 62 | 320 | `FITS` |  | 0.0 | 0.0 | 135 | 55/61/69/103 | 650 | `HOLED` |
| 1984 initial | 63 | 330 | `FITS` |  | 0.0 | 0.0 | 139 | 59/62/65/73 | 600 | `HOLED` |
| 1984 initial | 64 | 320 | `FITS` |  | 0.0 | 0.0 | 133 | 57/60/64/92 | 880 | `HOLED` |
| 1984 initial | 65 | 320 | `FITS` |  | 0.0 | 0.0 | 138 | 57/65/73/109 | 780 | `HOLED` |
| 1984 initial | 66 | 340 | `FITS` |  | 0.0 | 0.0 | 131 | 64/66/69/86 | 630 | `HOLED` |
| 1984 initial | 67 | 430 | `FITS` |  | 0.0 | 0.0 | 148 | 44/66/69/76 | 640 | `MIXED` |
| 1984 initial | 68 | 565 | `FITS` |  | 1.0 | 2.0 | 103 | 45/62/73/90 | 620 | `MIXED` |
| 1984 initial | 69 | 570 | `FITS` |  | 0.0 | 0.0 | 106 | 66/70/74/90 | 700 | `HOLED` |
| 1984 initial | 70 | 480 | `FITS` |  | 0.0 | 0.0 | 141 | 54/70/73/78 | 670 | `HOLED` |
| 1984 initial | 71 | 420 | `FITS` |  | 5.0 | 0.0 | 137 | 39/55/73/83 | 680 | `HOLED` |
| 1984 initial | 72 | 485 | `FITS` |  | 0.0 | 0.0 | 123 | 68/77/81/93 | 740 | `MIXED` |
| 1984 initial | 73 | 550 | `FITS` |  | 0.0 | 0.0 | 132 | 21/24/82/99 | 830 | `HOLED` |
| 1984 initial | 74 | 590 | `FITS` |  | 0.0 | 0.0 | 131 | 22/61/74/87 | 700 | `HOLED` |
| 1984 initial | 75 | 530 | `FITS` |  | 0.0 | 0.0 | 103 | 68/72/75/82 | 690 | `TRUE_END` |
| 1984 initial | 76 | 500 | `FITS` |  | 0.0 | 2.0 | 140 | 52/62/65/74 | 600 | `TRUE_END` |
| 1984 initial | 77 | 530 | `FITS` |  | 3.0 | 6.0 | 146 | 39/55/66/81 | 610 | `MIXED` |
| 1984 initial | 78 | 550 | `FITS` |  | 0.0 | 2.0 | 137 | 0/59/68/73 | 650 | `MIXED` |
| 1984 initial | 79 | 560 | `FITS` |  | 0.0 | 0.0 | 108 | 25/64/68/71 | 640 | `TRUE_END` |
| 1984 initial | 80 | 540 | `FITS` |  | 0.0 | 0.0 | 99 | 66/67/70/81 | 640 | `TRUE_END` |
| 1984 initial | 81 | 410 | `FITS` |  | 0.0 | 0.0 | 105 | 1/66/68/92 | 630 | `HOLED` |
| 1984 initial | 82 | 155 | `FITS` |  | 0.0 | 0.0 | 106 | 41/46/50/62 | 440 | `TRUE_END` |
| 1984 initial | 83 | 90 | `FITS` |  | 0.0 | 0.0 | 108 | 42/48/50/54 | 450 | `TRUE_END` |
| 1984 initial | 84 | 25 | `FITS` |  | 0.0 | 0.0 | 149 | 18/21/25/51 | 370 | `HOLED` |
| 1984 initial | 85 | 0 | `FITS` |  | 0.0 | 0.0 | 178 | 23/27/37/53 | 360 | `HOLED` |
| 1984 initial | 86 | 0 | `FITS` |  | 0.0 | 0.0 | 189 | 16/17/28/44 | 140 | `HOLED` |
| 1984 initial | 87 | 60 | `FITS` |  | 0.0 | 0.0 | 185 | 15/17/22/26 | 180 | `HOLED` |
| 1984 initial | 88 | 160 | `FITS` |  | 0.0 | 0.0 | 177 | 23/31/36/51 | 280 | `HOLED` |
| 1984 initial | 89 | 160 | `FITS` |  | 0.0 | 0.0 | 179 | 49/56/59/70 | 550 | `MIXED` |
| 1984 initial | 90 | 160 | `FITS` |  | 0.0 | 0.0 | 182 | 27/37/49/81 | 340 | `HOLED` |
| 2004 initial | 9 | 50 | `FITS` |  | 0.0 | 0.0 | 162 | 11/14/17/73 | 1530 | `HOLED` |
| 2004 initial | 10 | 40 | `FITS` |  | 0.0 | 0.0 | 151 | 10/11/14/18 | 90 | `HOLED` |
| 2004 initial | 11 | 20 | `FITS` |  | 0.0 | 0.0 | 157 | 6/8/9/16 | 50 | `HOLED` |
| 2004 initial | 12 | 30 | `FITS` |  | 0.0 | 0.0 | 165 | 7/8/11/15 | 50 | `TRUE_END` |
| 2004 initial | 13 | 50 | `FITS` |  | 0.0 | 0.0 | 166 | 8/10/12/15 | 70 | `TRUE_END` |
| 2004 initial | 14 | 65 | `FITS` |  | 0.0 | 16.0 | 34 | 8/9/12/20 | 60 | `TRUE_END` |
| 2004 initial | 15 | 50 | `FITS` |  | 0.0 | 0.0 | 29 | 8/9/10/18 | 60 | `TRUE_END` |
| 2004 initial | 16 | 60 | `FITS` |  | 0.0 | 0.0 | 64 | 15/18/20/22 | 150 | `TRUE_END` |
| 2004 initial | 17 | 90 | `FITS` |  | 0.0 | 0.0 | 71 | 17/21/26/47 | 200 | `HOLED` |
| 2004 initial | 18 | 120 | `FITS` |  | 0.0 | 0.0 | 87 | 16/29/35/54 | 390 | `HOLED` |
| 2004 initial | 19 | 130 | `FITS` |  | 0.0 | 2.0 | 97 | 8/18/25/48 | 270 | `HOLED` |
| 2004 initial | 20 | 80 | `FITS` |  | 0.0 | 0.0 | 106 | 30/33/45/55 | 350 | `MIXED` |
| 2004 initial | 21 | 125 | `FITS` |  | 0.0 | 0.0 | 113 | 25/42/44/71 | 610 | `HOLED` |
| 2004 initial | 22 | 190 | `FITS` |  | 0.0 | 0.0 | 128 | 34/38/44/100 | 890 | `HOLED` |
| 2004 initial | 23 | 190 | `FITS` |  | 0.0 | 0.0 | 133 | 41/48/64/93 | 890 | `HOLED` |
| 2004 initial | 24 | 170 | `FITS` |  | 0.0 | 0.0 | 137 | 46/55/86/106 | 920 | `HOLED` |
| 2004 initial | 25 | 160 | `FITS` |  | 0.0 | 0.0 | 134 | 38/48/76/112 | 1010 | `HOLED` |
| 2004 initial | 26 | 205 | `FITS` |  | 2.0 | 0.0 | 147 | 17/45/67/116 | 1050 | `HOLED` |
| 2004 initial | 27 | 260 | `FITS` |  | 0.0 | 0.0 | 161 | 40/72/83/121 | 1080 | `HOLED` |
| 2004 initial | 28 | 270 | `FITS` |  | 0.0 | 0.0 | 166 | 46/48/79/114 | 1150 | `HOLED` |
| 2004 initial | 29 | 260 | `FITS` |  | 0.0 | 0.0 | 155 | 38/47/87/121 | 990 | `HOLED` |
| 2004 initial | 30 | 270 | `FITS` |  | 0.0 | 0.0 | 142 | 33/52/62/83 | 740 | `HOLED` |
| 2004 initial | 31 | 230 | `FITS` |  | 0.0 | 0.0 | 140 | 30/32/41/67 | 340 | `HOLED` |
| 2004 initial | 32 | 210 | `FITS` |  | 0.0 | 0.0 | 144 | 32/44/46/50 | 420 | `TRUE_END` |
| 2004 initial | 33 | 250 | `FITS` |  | 0.0 | 0.0 | 149 | 30/44/46/53 | 430 | `MIXED` |
| 2004 initial | 34 | 320 | `FITS` |  | 0.0 | 0.0 | 150 | 41/43/46/52 | 400 | `TRUE_END` |
| 2004 initial | 35 | 320 | `FITS` |  | 0.0 | 2.0 | 129 | 33/42/44/52 | 400 | `MIXED` |
| 2004 initial | 36 | 330 | `FITS` |  | 0.0 | 0.0 | 133 | 51/52/54/60 | 490 | `TRUE_END` |
| 2004 initial | 37 | 320 | `FITS` |  | 0.0 | 0.0 | 157 | 53/57/78/86 | 740 | `HOLED` |
| 2004 initial | 38 | 280 | `FITS` |  | 0.0 | 0.0 | 165 | 69/78/81/88 | 770 | `TRUE_END` |
| 2004 initial | 39 | 230 | `FITS` |  | 0.0 | 0.0 | 170 | 68/69/71/77 | 660 | `TRUE_END` |
| 2004 initial | 40 | 210 | `FITS` |  | 0.0 | 0.0 | 172 | 53/63/70/75 | 640 | `MIXED` |
| 2004 initial | 41 | 190 | `FITS` |  | 0.0 | 0.0 | 174 | 43/52/54/58 | 490 | `HOLED` |
| 2004 initial | 42 | 180 | `FITS` |  | 0.0 | 0.0 | 176 | 41/42/60/65 | 560 | `HOLED` |
| 2004 initial | 43 | 170 | `FITS` |  | 0.0 | 0.0 | 179 | 53/61/63/66 | 590 | `TRUE_END` |
| 2004 initial | 44 | 160 | `FITS` |  | 0.0 | 0.0 | 169 | 46/50/53/66 | 530 | `HOLED` |
| 2004 initial | 45 | 150 | `FITS` |  | 0.0 | 0.0 | 155 | 32/38/41/56 | 360 | `HOLED` |
| 2004 initial | 46 | 140 | `FITS` |  | 0.0 | 0.0 | 140 | 33/38/40/46 | 360 | `HOLED` |
| 2004 initial | 47 | 140 | `FITS` |  | 0.0 | 0.0 | 126 | 29/32/34/41 | 300 | `HOLED` |
| 2004 initial | 48 | 120 | `FITS` |  | 0.0 | 0.0 | 116 | 21/23/25/30 | 200 | `TRUE_END` |
| 2004 initial | 49 | 130 | `FITS` |  | 0.0 | 0.0 | 95 | 20/22/23/27 | 190 | `TRUE_END` |
| 2004 initial | 50 | 110 | `FITS` |  | 0.0 | 0.0 | 86 | 15/19/21/24 | 160 | `TRUE_END` |
| 2004 initial | 51 | 120 | `FITS` |  | 0.0 | 0.0 | 54 | 18/22/24/28 | 190 | `TRUE_END` |
| 2004 initial | 52 | 95 | `FITS` |  | 0.0 | 0.0 | 59 | 15/18/21/27 | 160 | `TRUE_END` |
| 2004 initial | 53 | 140 | `FITS` |  | 0.0 | 0.0 | 53 | 18/23/26/33 | 210 | `MIXED` |
| 2004 initial | 54 | 150 | `FITS` |  | 0.0 | 0.0 | 111 | 22/25/28/45 | 230 | `HOLED` |
| 2004 initial | 55 | 160 | `FITS` |  | 0.0 | 0.0 | 106 | 23/33/36/49 | 310 | `HOLED` |
| 2004 initial | 56 | 180 | `FITS` |  | 0.0 | 0.0 | 77 | 24/26/29/34 | 240 | `MIXED` |
| 2004 initial | 57 | 200 | `FITS` |  | 0.0 | 0.0 | 64 | 24/29/36/45 | 270 | `TRUE_END` |
| 2004 initial | 58 | 240 | `FITS` |  | 0.0 | 0.0 | 68 | 32/40/45/52 | 390 | `MIXED` |
| 2004 initial | 59 | 270 | `FITS` |  | 0.0 | 0.0 | 78 | 42/44/48/61 | 430 | `HOLED` |
| 2004 initial | 60 | 310 | `FITS` |  | 0.0 | 0.0 | 102 | 48/56/62/73 | 590 | `HOLED` |
| 2004 initial | 61 | 330 | `FITS` |  | 0.0 | 0.0 | 117 | 48/53/62/90 | 900 | `HOLED` |
| 2004 initial | 62 | 330 | `FITS` |  | 0.0 | 0.0 | 136 | 56/62/70/103 | 650 | `HOLED` |
| 2004 initial | 63 | 330 | `FITS` |  | 0.0 | 0.0 | 139 | 58/62/65/75 | 610 | `HOLED` |
| 2004 initial | 64 | 370 | `FITS` |  | 0.0 | 0.0 | 138 | 62/65/70/94 | 910 | `HOLED` |
| 2004 initial | 65 | 370 | `FITS` |  | 0.0 | 0.0 | 139 | 60/69/77/114 | 810 | `HOLED` |
| 2004 initial | 66 | 380 | `FITS` |  | 0.0 | 0.0 | 131 | 67/70/72/85 | 670 | `HOLED` |
| 2004 initial | 67 | 420 | `FITS` |  | 0.0 | 0.0 | 144 | 46/67/68/71 | 640 | `MIXED` |
| 2004 initial | 68 | 555 | `FITS` |  | 2.0 | 2.0 | 99 | 46/61/72/90 | 610 | `MIXED` |
| 2004 initial | 69 | 560 | `FITS` |  | 0.0 | 0.0 | 104 | 68/70/75/90 | 720 | `HOLED` |
| 2004 initial | 70 | 480 | `FITS` |  | 0.0 | 0.0 | 140 | 54/69/72/78 | 670 | `HOLED` |
| 2004 initial | 71 | 410 | `FITS` |  | 2.0 | 2.0 | 137 | 39/55/73/82 | 670 | `HOLED` |
| 2004 initial | 72 | 490 | `FITS` |  | 0.0 | 0.0 | 123 | 70/79/82/93 | 760 | `MIXED` |
| 2004 initial | 73 | 550 | `FITS` |  | 0.0 | 0.0 | 132 | 22/25/84/99 | 850 | `HOLED` |
| 2004 initial | 74 | 600 | `FITS` |  | 0.0 | 0.0 | 132 | 23/62/74/88 | 710 | `HOLED` |
| 2004 initial | 75 | 540 | `FITS` |  | 0.0 | 0.0 | 103 | 68/73/76/82 | 700 | `TRUE_END` |
| 2004 initial | 76 | 505 | `FITS` |  | 0.0 | 2.0 | 140 | 52/62/65/74 | 600 | `TRUE_END` |
| 2004 initial | 77 | 500 | `FITS` |  | 3.0 | 6.0 | 143 | 35/52/61/79 | 560 | `MIXED` |
| 2004 initial | 78 | 490 | `FITS` |  | 1.0 | 2.0 | 130 | 35/52/61/67 | 570 | `MIXED` |
| 2004 initial | 79 | 500 | `FITS` |  | 0.0 | 0.0 | 101 | 18/61/63/66 | 590 | `TRUE_END` |
| 2004 initial | 80 | 490 | `FITS` |  | 0.0 | 0.0 | 90 | 62/63/64/72 | 600 | `TRUE_END` |
| 2004 initial | 81 | 365 | `FITS` |  | 0.0 | 0.0 | 100 | 61/63/64/84 | 600 | `HOLED` |
| 2004 initial | 82 | 120 | `FITS` |  | 0.0 | 0.0 | 104 | 39/44/48/61 | 420 | `TRUE_END` |
| 2004 initial | 83 | 70 | `FITS` |  | 0.0 | 0.0 | 106 | 41/47/48/51 | 440 | `TRUE_END` |
| 2004 initial | 84 | 30 | `FITS` |  | 0.0 | 2.0 | 147 | 5/18/21/48 | 340 | `HOLED` |
| 2004 initial | 85 | 70 | `FITS` |  | 0.0 | 0.0 | 176 | 20/23/33/50 | 300 | `HOLED` |
| 2004 initial | 86 | 50 | `FITS` |  | 0.0 | 0.0 | 185 | 12/13/22/40 | 190 | `HOLED` |
| 2004 initial | 87 | 40 | `FITS` |  | 0.0 | 0.0 | 181 | 12/13/19/24 | 140 | `HOLED` |
| 2004 initial | 88 | 120 | `FITS` |  | 0.0 | 0.0 | 174 | 20/28/32/49 | 250 | `HOLED` |
| 2004 initial | 89 | 160 | `FITS` |  | 0.0 | 0.0 | 177 | 47/55/58/69 | 540 | `MIXED` |
| 2004 initial | 90 | 150 | `FITS` |  | 0.0 | 0.0 | 181 | 26/37/49/81 | 350 | `HOLED` |
| 1989 relocation (GIS 84-87) - zero-retreat bound | 84 | 65 | `FITS` |  | 0.0 | 0.0 | 149 | 18/21/25/51 | 370 | `HOLED` |
| 1989 relocation (GIS 84-87) - zero-retreat bound | 85 | 109 | `FITS` |  | 0.0 | 0.0 | 178 | 23/27/37/53 | 360 | `HOLED` |
| 1989 relocation (GIS 84-87) - zero-retreat bound | 86 | 94 | `FITS` |  | 0.0 | 0.0 | 189 | 16/17/28/44 | 140 | `HOLED` |
| 1989 relocation (GIS 84-87) - zero-retreat bound | 87 | 78 | `FITS` |  | 0.0 | 0.0 | 185 | 15/17/22/26 | 180 | `HOLED` |
| 1999 relocation (GIS 9-14) - zero-retreat bound | 9 | 58 | `FITS` |  | 0.0 | 0.0 | 163 | 11/14/18/73 | 1540 | `HOLED` |
| 1999 relocation (GIS 9-14) - zero-retreat bound | 10 | 67 | `FITS` |  | 0.0 | 0.0 | 154 | 14/15/17/22 | 130 | `HOLED` |
| 1999 relocation (GIS 9-14) - zero-retreat bound | 11 | 87 | `FITS` |  | 0.0 | 0.0 | 163 | 14/14/16/22 | 110 | `HOLED` |
| 1999 relocation (GIS 9-14) - zero-retreat bound | 12 | 88 | `FITS` |  | 0.0 | 0.0 | 170 | 13/14/17/21 | 110 | `TRUE_END` |
| 1999 relocation (GIS 9-14) - zero-retreat bound | 13 | 81 | `FITS` |  | 0.0 | 0.0 | 170 | 15/17/18/20 | 140 | `TRUE_END` |
| 1999 relocation (GIS 9-14) - zero-retreat bound | 14 | 85 | `FITS` |  | 0.0 | 6.0 | 35 | 10/13/16/25 | 100 | `TRUE_END` |
