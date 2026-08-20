# NC-12 road setback audit

Generated 2026-08-19T13:35:33 by `HAT_road_setback_audit.py`.

**Diagnostic record.** No setback was modified, no setback file was written, no topography was touched.

| | |
|---|---|
| Topography | `C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\1-barrier3d-domains\2009-dune-topo\2009_v4\topography` |
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

## 1999 relocation (inter-village south, GIS 9-15) - zero-retreat bound

7 domains — `DROWNS` 5, `FITS` 2. **5 with the road off the island.**

| GIS | Setback (m) | % of road's own cells in water | % of profiles fully submerged | % water at bayside border | Largest that fits (m) | Verdict | Narrowness |
|---:|---:|---:|---:|---:|---:|---|---|
| 11 | 129 | 89.0 | 88.0 | 92.0 | 50 | `DROWNS` **OFF-ISLAND** | `HOLED` |
| 12 | 126 | 80.0 | 72.0 | 100.0 | 40 | `DROWNS` **OFF-ISLAND** | `TRUE_END` |
| 13 | 125 | 71.0 | 62.0 | 98.0 | 70 | `DROWNS` **OFF-ISLAND** | `TRUE_END` |
| 14 | 126 | 52.0 | 52.0 | 66.0 | 50 | `DROWNS` **OFF-ISLAND** | `TRUE_END` |
| 15 | 106 | 38.0 | 22.0 | 64.0 | 60 | `DROWNS` **OFF-ISLAND** | `TRUE_END` |

### Guards the prescribed-relocation path skips

CASCADE has two guards on moving a road landward, and the runner's `HISTORICAL_ROAD_EVENTS` path takes **neither** — it assigns `rm._road_setback = new_sb` directly, leaves `road_ele` at its SLR-decremented value rather than rebuilding at grade, and lets the next `bulldoze()` flatten whatever is there.

| GIS | Setback (m) | Rebuilt road elevation (m) | `get_road_relocation_elevation` | Average barrier width (m) | `road_relocation_checks` |
|---:|---:|---:|---|---:|---|
| 9 | 73 | 1.27 | allows | 199 | allows |
| 10 | 97 | 0.49 | allows | 150 | allows |
| 11 | 129 | -1.94 | **refuses** — road ≤ 0 m MSL | 95 | **refuses** — island too narrow |
| 12 | 126 | -1.59 | **refuses** — road ≤ 0 m MSL | 97 | **refuses** — island too narrow |
| 13 | 125 | -0.78 | **refuses** — road ≤ 0 m MSL | 112 | **refuses** — island too narrow |
| 14 | 126 | -0.74 | **refuses** — road ≤ 0 m MSL | 123 | **refuses** — island too narrow |
| 15 | 106 | 0.50 | allows | 115 | **refuses** — island too narrow |

**5 of 7 prescribed relocations would be refused by CASCADE's own logic.** Evaluated on the t=0 interior; the real grid at the event year will have evolved, so treat these as the initialisation-time verdict.

## 1989 relocation (Pea Island, GIS 84-87) - zero-retreat bound

4 domains — `DROWNS` 3, `FITS` 1. **2 with the road off the island.**

| GIS | Setback (m) | % of road's own cells in water | % of profiles fully submerged | % water at bayside border | Largest that fits (m) | Verdict | Narrowness |
|---:|---:|---:|---:|---:|---:|---|---|
| 84 | 163 | 31.0 | 30.0 | 34.0 | 320 | `DROWNS` **OFF-ISLAND** | `HOLED` |
| 86 | 205 | 22.0 | 18.0 | 38.0 | 170 | `DROWNS` **OFF-ISLAND** | `HOLED` |
| 87 | 113 | 2.0 | 0.0 | 38.0 | 130 | `DROWNS` | `HOLED` |

### Guards the prescribed-relocation path skips

CASCADE has two guards on moving a road landward, and the runner's `HISTORICAL_ROAD_EVENTS` path takes **neither** — it assigns `rm._road_setback = new_sb` directly, leaves `road_ele` at its SLR-decremented value rather than rebuilding at grade, and lets the next `bulldoze()` flatten whatever is there.

| GIS | Setback (m) | Rebuilt road elevation (m) | `get_road_relocation_elevation` | Average barrier width (m) | `road_relocation_checks` |
|---:|---:|---:|---|---:|---|
| 84 | 163 | 0.21 | allows | 202 | **refuses** — island too narrow |
| 85 | 165 | 0.27 | allows | 300 | allows |
| 86 | 205 | 0.17 | allows | 223 | **refuses** — island too narrow |
| 87 | 113 | 0.38 | allows | 166 | allows |

**2 of 4 prescribed relocations would be refused by CASCADE's own logic.** Evaluated on the t=0 interior; the real grid at the event year will have evolved, so treat these as the initialisation-time verdict.

## Where Barrier3D's two guards disagree

The model does not use one definition of "island":

| Routine | What it reads | Used for |
|---|---|---|
| `FindWidths` | Stops at the **first** water cell — land behind a gap is not part of the island | `InteriorWidth_Avg`, i.e. the `average_barrier_width` that `road_relocation_checks()` uses |
| `bulldoze` | The **literal cell** at `road_end + 1` — land behind a gap counts normally | The `roadway_drown` test |

So the road can be founded on ground the relocation logic believes is off the back of the barrier: the drown test passes while the relocation test sees a much narrower island. Counted, not corrected — resolving it would mean changing a Barrier3D definition.

| Scenario | GIS | % of profiles where the road's border cell is land but `FindWidths` has already ended the island | Narrowness |
|---|---:|---:|---|
| 1984 initial | 73 | 24.0 | `MIXED` |
| 2004 initial | 73 | 24.0 | `MIXED` |
| 1984 initial | 74 | 20.0 | `MIXED` |
| 2004 initial | 74 | 20.0 | `MIXED` |
| 1984 initial | 77 | 18.0 | `MIXED` |
| 2004 initial | 77 | 18.0 | `MIXED` |
| 1989 relocation (Pea Island, GIS 84-87) - zero-retreat bound | 86 | 14.0 | `HOLED` |
| 1989 relocation (Pea Island, GIS 84-87) - zero-retreat bound | 84 | 12.0 | `HOLED` |
| 1984 initial | 78 | 10.0 | `MIXED` |
| 2004 initial | 78 | 10.0 | `MIXED` |
| 1984 initial | 26 | 8.0 | `HOLED` |
| 1984 initial | 68 | 8.0 | `MIXED` |
| 1984 initial | 71 | 8.0 | `HOLED` |
| 1984 initial | 79 | 8.0 | `TRUE_END` |
| 2004 initial | 26 | 8.0 | `HOLED` |
| 2004 initial | 68 | 8.0 | `MIXED` |
| 2004 initial | 71 | 8.0 | `HOLED` |
| 2004 initial | 79 | 8.0 | `TRUE_END` |
| 1984 initial | 19 | 4.0 | `HOLED` |
| 2004 initial | 19 | 4.0 | `HOLED` |
| 1984 initial | 35 | 2.0 | `MIXED` |
| 1984 initial | 69 | 2.0 | `HOLED` |
| 2004 initial | 35 | 2.0 | `MIXED` |
| 2004 initial | 69 | 2.0 | `HOLED` |
| 1989 relocation (Pea Island, GIS 84-87) - zero-retreat bound | 87 | 2.0 | `HOLED` |

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
| 1984 initial | 9 | 40 | `FITS` |  | 0.0 | 0.0 | 100 | 11/14/17/73 | 120 | `HOLED` |
| 1984 initial | 10 | 20 | `FITS` |  | 10.0 | 10.0 | 64 | 0/14/16/21 | 110 | `HOLED` |
| 1984 initial | 11 | 0 | `FITS` |  | 0.0 | 0.0 | 62 | 6/8/9/16 | 50 | `HOLED` |
| 1984 initial | 12 | 0 | `FITS` |  | 0.0 | 0.0 | 58 | 6/7/10/14 | 40 | `TRUE_END` |
| 1984 initial | 13 | 0 | `FITS` |  | 0.0 | 0.0 | 16 | 8/9/11/15 | 70 | `TRUE_END` |
| 1984 initial | 14 | 30 | `FITS` |  | 0.0 | 0.0 | 22 | 7/8/11/19 | 50 | `TRUE_END` |
| 1984 initial | 15 | 50 | `FITS` |  | 0.0 | 0.0 | 21 | 9/9/10/19 | 60 | `TRUE_END` |
| 1984 initial | 16 | 75 | `FITS` |  | 0.0 | 0.0 | 25 | 17/19/21/24 | 170 | `TRUE_END` |
| 1984 initial | 17 | 90 | `FITS` |  | 0.0 | 0.0 | 47 | 17/21/26/46 | 200 | `HOLED` |
| 1984 initial | 18 | 120 | `FITS` |  | 0.0 | 0.0 | 62 | 16/29/35/54 | 390 | `HOLED` |
| 1984 initial | 19 | 130 | `FITS` |  | 0.0 | 2.0 | 55 | 8/18/25/48 | 270 | `HOLED` |
| 1984 initial | 20 | 80 | `FITS` |  | 0.0 | 0.0 | 58 | 29/33/45/55 | 390 | `MIXED` |
| 1984 initial | 21 | 120 | `FITS` |  | 0.0 | 0.0 | 92 | 25/42/44/71 | 610 | `HOLED` |
| 1984 initial | 22 | 185 | `FITS` |  | 0.0 | 0.0 | 100 | 33/38/43/97 | 870 | `HOLED` |
| 1984 initial | 23 | 220 | `FITS` |  | 0.0 | 0.0 | 101 | 42/49/66/93 | 880 | `HOLED` |
| 1984 initial | 24 | 170 | `FITS` |  | 0.0 | 0.0 | 99 | 46/56/86/97 | 920 | `HOLED` |
| 1984 initial | 25 | 160 | `FITS` |  | 0.0 | 0.0 | 102 | 38/48/76/98 | 950 | `HOLED` |
| 1984 initial | 26 | 205 | `FITS` |  | 2.0 | 0.0 | 102 | 17/45/67/100 | 950 | `HOLED` |
| 1984 initial | 27 | 260 | `FITS` |  | 0.0 | 0.0 | 104 | 44/72/85/103 | 980 | `HOLED` |
| 1984 initial | 28 | 270 | `FITS` |  | 0.0 | 0.0 | 107 | 46/48/79/105 | 1010 | `HOLED` |
| 1984 initial | 29 | 260 | `FITS` |  | 0.0 | 0.0 | 110 | 38/47/87/108 | 990 | `HOLED` |
| 1984 initial | 30 | 250 | `FITS` |  | 0.0 | 0.0 | 112 | 32/51/60/83 | 720 | `HOLED` |
| 1984 initial | 31 | 250 | `FITS` |  | 0.0 | 0.0 | 69 | 32/34/43/67 | 320 | `HOLED` |
| 1984 initial | 32 | 230 | `FITS` |  | 0.0 | 0.0 | 52 | 33/46/48/51 | 430 | `TRUE_END` |
| 1984 initial | 33 | 260 | `FITS` |  | 0.0 | 0.0 | 57 | 31/45/48/53 | 450 | `MIXED` |
| 1984 initial | 34 | 320 | `FITS` |  | 0.0 | 0.0 | 54 | 41/43/46/52 | 400 | `TRUE_END` |
| 1984 initial | 35 | 320 | `FITS` |  | 0.0 | 2.0 | 54 | 33/43/44/52 | 400 | `MIXED` |
| 1984 initial | 36 | 330 | `FITS` |  | 0.0 | 0.0 | 69 | 51/52/54/60 | 490 | `TRUE_END` |
| 1984 initial | 37 | 320 | `FITS` |  | 0.0 | 0.0 | 88 | 53/57/78/86 | 740 | `HOLED` |
| 1984 initial | 38 | 280 | `FITS` |  | 0.0 | 0.0 | 90 | 69/78/81/88 | 770 | `TRUE_END` |
| 1984 initial | 39 | 230 | `FITS` |  | 0.0 | 0.0 | 78 | 68/69/71/77 | 660 | `TRUE_END` |
| 1984 initial | 40 | 210 | `FITS` |  | 0.0 | 0.0 | 77 | 53/63/70/75 | 640 | `MIXED` |
| 1984 initial | 41 | 190 | `FITS` |  | 0.0 | 0.0 | 65 | 43/52/54/58 | 490 | `HOLED` |
| 1984 initial | 42 | 190 | `FITS` |  | 0.0 | 0.0 | 66 | 42/43/61/65 | 570 | `HOLED` |
| 1984 initial | 43 | 170 | `FITS` |  | 0.0 | 0.0 | 68 | 53/61/63/66 | 590 | `TRUE_END` |
| 1984 initial | 44 | 160 | `FITS` |  | 0.0 | 0.0 | 66 | 46/50/53/65 | 530 | `HOLED` |
| 1984 initial | 45 | 150 | `FITS` |  | 0.0 | 0.0 | 64 | 31/38/41/56 | 370 | `HOLED` |
| 1984 initial | 46 | 140 | `FITS` |  | 0.0 | 0.0 | 57 | 30/36/39/46 | 340 | `HOLED` |
| 1984 initial | 47 | 140 | `FITS` |  | 0.0 | 0.0 | 50 | 29/32/34/41 | 300 | `HOLED` |
| 1984 initial | 48 | 120 | `FITS` |  | 0.0 | 0.0 | 36 | 21/23/25/30 | 200 | `TRUE_END` |
| 1984 initial | 49 | 130 | `FITS` |  | 0.0 | 0.0 | 29 | 20/22/23/27 | 190 | `TRUE_END` |
| 1984 initial | 50 | 110 | `FITS` |  | 0.0 | 0.0 | 26 | 15/19/21/24 | 160 | `TRUE_END` |
| 1984 initial | 51 | 120 | `FITS` |  | 0.0 | 0.0 | 30 | 18/22/24/28 | 190 | `TRUE_END` |
| 1984 initial | 52 | 95 | `FITS` |  | 0.0 | 0.0 | 32 | 15/18/21/27 | 160 | `TRUE_END` |
| 1984 initial | 53 | 140 | `FITS` |  | 0.0 | 0.0 | 37 | 21/24/26/33 | 210 | `MIXED` |
| 1984 initial | 54 | 150 | `FITS` |  | 0.0 | 0.0 | 45 | 22/25/28/44 | 230 | `HOLED` |
| 1984 initial | 55 | 160 | `FITS` |  | 0.0 | 0.0 | 52 | 23/33/36/49 | 310 | `HOLED` |
| 1984 initial | 56 | 180 | `FITS` |  | 0.0 | 0.0 | 52 | 24/27/29/34 | 250 | `MIXED` |
| 1984 initial | 57 | 200 | `FITS` |  | 0.0 | 0.0 | 47 | 25/29/36/45 | 270 | `TRUE_END` |
| 1984 initial | 58 | 240 | `FITS` |  | 0.0 | 0.0 | 55 | 32/41/45/52 | 400 | `MIXED` |
| 1984 initial | 59 | 270 | `FITS` |  | 0.0 | 0.0 | 58 | 42/44/48/57 | 430 | `HOLED` |
| 1984 initial | 60 | 310 | `FITS` |  | 0.0 | 0.0 | 60 | 48/56/58/59 | 540 | `TRUE_END` |
| 1984 initial | 61 | 330 | `FITS` |  | 0.0 | 0.0 | 63 | 48/53/59/62 | 560 | `MIXED` |
| 1984 initial | 62 | 350 | `FITS` |  | 0.0 | 0.0 | 68 | 56/62/64/67 | 590 | `MIXED` |
| 1984 initial | 63 | 350 | `FITS` |  | 0.0 | 0.0 | 70 | 60/63/65/69 | 620 | `MIXED` |
| 1984 initial | 64 | 370 | `FITS` |  | 0.0 | 0.0 | 72 | 62/65/69/70 | 660 | `MIXED` |
| 1984 initial | 65 | 370 | `FITS` |  | 0.0 | 0.0 | 71 | 60/67/68/70 | 640 | `TRUE_END` |
| 1984 initial | 66 | 370 | `FITS` |  | 0.0 | 0.0 | 72 | 65/67/68/71 | 640 | `TRUE_END` |
| 1984 initial | 67 | 420 | `FITS` |  | 0.0 | 0.0 | 74 | 46/66/68/71 | 640 | `TRUE_END` |
| 1984 initial | 68 | 555 | `FITS` |  | 2.0 | 2.0 | 78 | 46/61/72/76 | 610 | `MIXED` |
| 1984 initial | 69 | 560 | `FITS` |  | 0.0 | 0.0 | 82 | 37/70/75/79 | 720 | `HOLED` |
| 1984 initial | 70 | 480 | `FITS` |  | 0.0 | 0.0 | 80 | 54/69/72/78 | 670 | `HOLED` |
| 1984 initial | 71 | 410 | `FITS` |  | 2.0 | 2.0 | 83 | 39/55/73/82 | 670 | `HOLED` |
| 1984 initial | 72 | 490 | `FITS` |  | 0.0 | 0.0 | 82 | 70/79/80/81 | 760 | `TRUE_END` |
| 1984 initial | 73 | 550 | `FITS` |  | 0.0 | 0.0 | 80 | 22/25/77/79 | 730 | `MIXED` |
| 1984 initial | 74 | 600 | `FITS` |  | 0.0 | 0.0 | 76 | 22/62/74/75 | 710 | `MIXED` |
| 1984 initial | 75 | 540 | `FITS` |  | 0.0 | 0.0 | 74 | 68/69/71/73 | 670 | `TRUE_END` |
| 1984 initial | 76 | 505 | `FITS` |  | 0.0 | 2.0 | 68 | 52/62/64/67 | 600 | `TRUE_END` |
| 1984 initial | 77 | 500 | `FITS` |  | 3.0 | 6.0 | 63 | 35/52/55/62 | 520 | `MIXED` |
| 1984 initial | 78 | 470 | `FITS` |  | 0.0 | 16.0 | 68 | 35/49/51/67 | 470 | `MIXED` |
| 1984 initial | 79 | 400 | `FITS` |  | 0.0 | 20.0 | 67 | 18/42/62/66 | 400 | `TRUE_END` |
| 1984 initial | 80 | 440 | `FITS` |  | 0.0 | 0.0 | 68 | 47/47/63/67 | 440 | `TRUE_END` |
| 1984 initial | 81 | 350 | `FITS` |  | 0.0 | 0.0 | 58 | 50/50/51/57 | 470 | `TRUE_END` |
| 1984 initial | 82 | 115 | `FITS` |  | 0.0 | 0.0 | 57 | 39/44/48/56 | 420 | `TRUE_END` |
| 1984 initial | 83 | 60 | `FITS` |  | 0.0 | 0.0 | 52 | 39/45/46/51 | 420 | `TRUE_END` |
| 1984 initial | 84 | 0 | `FITS` |  | 0.0 | 0.0 | 51 | 3/15/19/44 | 320 | `HOLED` |
| 1984 initial | 85 | 0 | `FITS` |  | 0.0 | 0.0 | 50 | 18/20/30/47 | 290 | `HOLED` |
| 1984 initial | 86 | 0 | `FITS` |  | 0.0 | 0.0 | 40 | 10/11/20/39 | 170 | `HOLED` |
| 1984 initial | 87 | 20 | `FITS` |  | 0.0 | 0.0 | 40 | 11/12/18/23 | 130 | `HOLED` |
| 1984 initial | 88 | 110 | `FITS` |  | 0.0 | 0.0 | 56 | 20/28/32/49 | 250 | `HOLED` |
| 1984 initial | 89 | 160 | `FITS` |  | 0.0 | 0.0 | 73 | 47/55/58/69 | 540 | `MIXED` |
| 1984 initial | 90 | 150 | `FITS` |  | 0.0 | 0.0 | 76 | 27/37/49/74 | 350 | `HOLED` |
| 2004 initial | 9 | 50 | `FITS` |  | 0.0 | 0.0 | 100 | 11/14/17/73 | 120 | `HOLED` |
| 2004 initial | 10 | 70 | `FITS` |  | 10.0 | 10.0 | 64 | 0/14/16/21 | 110 | `HOLED` |
| 2004 initial | 11 | 20 | `FITS` |  | 0.0 | 0.0 | 62 | 6/8/9/16 | 50 | `HOLED` |
| 2004 initial | 12 | 20 | `FITS` |  | 0.0 | 0.0 | 58 | 6/7/10/14 | 40 | `TRUE_END` |
| 2004 initial | 13 | 40 | `FITS` |  | 0.0 | 0.0 | 16 | 8/9/11/15 | 70 | `TRUE_END` |
| 2004 initial | 14 | 55 | `FITS` |  | 0.0 | 4.0 | 22 | 7/8/11/19 | 50 | `TRUE_END` |
| 2004 initial | 15 | 50 | `FITS` |  | 0.0 | 0.0 | 21 | 9/9/10/19 | 60 | `TRUE_END` |
| 2004 initial | 16 | 70 | `FITS` |  | 0.0 | 0.0 | 25 | 17/19/21/24 | 170 | `TRUE_END` |
| 2004 initial | 17 | 90 | `FITS` |  | 0.0 | 0.0 | 47 | 17/21/26/46 | 200 | `HOLED` |
| 2004 initial | 18 | 120 | `FITS` |  | 0.0 | 0.0 | 62 | 16/29/35/54 | 390 | `HOLED` |
| 2004 initial | 19 | 130 | `FITS` |  | 0.0 | 2.0 | 55 | 8/18/25/48 | 270 | `HOLED` |
| 2004 initial | 20 | 80 | `FITS` |  | 0.0 | 0.0 | 58 | 29/33/45/55 | 390 | `MIXED` |
| 2004 initial | 21 | 125 | `FITS` |  | 0.0 | 0.0 | 92 | 25/42/44/71 | 610 | `HOLED` |
| 2004 initial | 22 | 185 | `FITS` |  | 0.0 | 0.0 | 100 | 33/38/43/97 | 870 | `HOLED` |
| 2004 initial | 23 | 220 | `FITS` |  | 0.0 | 0.0 | 101 | 42/49/66/93 | 880 | `HOLED` |
| 2004 initial | 24 | 170 | `FITS` |  | 0.0 | 0.0 | 99 | 46/56/86/97 | 920 | `HOLED` |
| 2004 initial | 25 | 160 | `FITS` |  | 0.0 | 0.0 | 102 | 38/48/76/98 | 950 | `HOLED` |
| 2004 initial | 26 | 205 | `FITS` |  | 2.0 | 0.0 | 102 | 17/45/67/100 | 950 | `HOLED` |
| 2004 initial | 27 | 260 | `FITS` |  | 0.0 | 0.0 | 104 | 44/72/85/103 | 980 | `HOLED` |
| 2004 initial | 28 | 270 | `FITS` |  | 0.0 | 0.0 | 107 | 46/48/79/105 | 1010 | `HOLED` |
| 2004 initial | 29 | 260 | `FITS` |  | 0.0 | 0.0 | 110 | 38/47/87/108 | 990 | `HOLED` |
| 2004 initial | 30 | 250 | `FITS` |  | 0.0 | 0.0 | 112 | 32/51/60/83 | 720 | `HOLED` |
| 2004 initial | 31 | 250 | `FITS` |  | 0.0 | 0.0 | 69 | 32/34/43/67 | 320 | `HOLED` |
| 2004 initial | 32 | 230 | `FITS` |  | 0.0 | 0.0 | 52 | 33/46/48/51 | 430 | `TRUE_END` |
| 2004 initial | 33 | 260 | `FITS` |  | 0.0 | 0.0 | 57 | 31/45/48/53 | 450 | `MIXED` |
| 2004 initial | 34 | 320 | `FITS` |  | 0.0 | 0.0 | 54 | 41/43/46/52 | 400 | `TRUE_END` |
| 2004 initial | 35 | 320 | `FITS` |  | 0.0 | 2.0 | 54 | 33/43/44/52 | 400 | `MIXED` |
| 2004 initial | 36 | 330 | `FITS` |  | 0.0 | 0.0 | 69 | 51/52/54/60 | 490 | `TRUE_END` |
| 2004 initial | 37 | 320 | `FITS` |  | 0.0 | 0.0 | 88 | 53/57/78/86 | 740 | `HOLED` |
| 2004 initial | 38 | 280 | `FITS` |  | 0.0 | 0.0 | 90 | 69/78/81/88 | 770 | `TRUE_END` |
| 2004 initial | 39 | 230 | `FITS` |  | 0.0 | 0.0 | 78 | 68/69/71/77 | 660 | `TRUE_END` |
| 2004 initial | 40 | 210 | `FITS` |  | 0.0 | 0.0 | 77 | 53/63/70/75 | 640 | `MIXED` |
| 2004 initial | 41 | 190 | `FITS` |  | 0.0 | 0.0 | 65 | 43/52/54/58 | 490 | `HOLED` |
| 2004 initial | 42 | 190 | `FITS` |  | 0.0 | 0.0 | 66 | 42/43/61/65 | 570 | `HOLED` |
| 2004 initial | 43 | 170 | `FITS` |  | 0.0 | 0.0 | 68 | 53/61/63/66 | 590 | `TRUE_END` |
| 2004 initial | 44 | 160 | `FITS` |  | 0.0 | 0.0 | 66 | 46/50/53/65 | 530 | `HOLED` |
| 2004 initial | 45 | 150 | `FITS` |  | 0.0 | 0.0 | 64 | 31/38/41/56 | 370 | `HOLED` |
| 2004 initial | 46 | 140 | `FITS` |  | 0.0 | 0.0 | 57 | 30/36/39/46 | 340 | `HOLED` |
| 2004 initial | 47 | 140 | `FITS` |  | 0.0 | 0.0 | 50 | 29/32/34/41 | 300 | `HOLED` |
| 2004 initial | 48 | 120 | `FITS` |  | 0.0 | 0.0 | 36 | 21/23/25/30 | 200 | `TRUE_END` |
| 2004 initial | 49 | 130 | `FITS` |  | 0.0 | 0.0 | 29 | 20/22/23/27 | 190 | `TRUE_END` |
| 2004 initial | 50 | 110 | `FITS` |  | 0.0 | 0.0 | 26 | 15/19/21/24 | 160 | `TRUE_END` |
| 2004 initial | 51 | 120 | `FITS` |  | 0.0 | 0.0 | 30 | 18/22/24/28 | 190 | `TRUE_END` |
| 2004 initial | 52 | 95 | `FITS` |  | 0.0 | 0.0 | 32 | 15/18/21/27 | 160 | `TRUE_END` |
| 2004 initial | 53 | 140 | `FITS` |  | 0.0 | 0.0 | 37 | 21/24/26/33 | 210 | `MIXED` |
| 2004 initial | 54 | 150 | `FITS` |  | 0.0 | 0.0 | 45 | 22/25/28/44 | 230 | `HOLED` |
| 2004 initial | 55 | 160 | `FITS` |  | 0.0 | 0.0 | 52 | 23/33/36/49 | 310 | `HOLED` |
| 2004 initial | 56 | 180 | `FITS` |  | 0.0 | 0.0 | 52 | 24/27/29/34 | 250 | `MIXED` |
| 2004 initial | 57 | 200 | `FITS` |  | 0.0 | 0.0 | 47 | 25/29/36/45 | 270 | `TRUE_END` |
| 2004 initial | 58 | 240 | `FITS` |  | 0.0 | 0.0 | 55 | 32/41/45/52 | 400 | `MIXED` |
| 2004 initial | 59 | 270 | `FITS` |  | 0.0 | 0.0 | 58 | 42/44/48/57 | 430 | `HOLED` |
| 2004 initial | 60 | 310 | `FITS` |  | 0.0 | 0.0 | 60 | 48/56/58/59 | 540 | `TRUE_END` |
| 2004 initial | 61 | 330 | `FITS` |  | 0.0 | 0.0 | 63 | 48/53/59/62 | 560 | `MIXED` |
| 2004 initial | 62 | 350 | `FITS` |  | 0.0 | 0.0 | 68 | 56/62/64/67 | 590 | `MIXED` |
| 2004 initial | 63 | 350 | `FITS` |  | 0.0 | 0.0 | 70 | 60/63/65/69 | 620 | `MIXED` |
| 2004 initial | 64 | 370 | `FITS` |  | 0.0 | 0.0 | 72 | 62/65/69/70 | 660 | `MIXED` |
| 2004 initial | 65 | 370 | `FITS` |  | 0.0 | 0.0 | 71 | 60/67/68/70 | 640 | `TRUE_END` |
| 2004 initial | 66 | 370 | `FITS` |  | 0.0 | 0.0 | 72 | 65/67/68/71 | 640 | `TRUE_END` |
| 2004 initial | 67 | 420 | `FITS` |  | 0.0 | 0.0 | 74 | 46/66/68/71 | 640 | `TRUE_END` |
| 2004 initial | 68 | 555 | `FITS` |  | 2.0 | 2.0 | 78 | 46/61/72/76 | 610 | `MIXED` |
| 2004 initial | 69 | 560 | `FITS` |  | 0.0 | 0.0 | 82 | 37/70/75/79 | 720 | `HOLED` |
| 2004 initial | 70 | 480 | `FITS` |  | 0.0 | 0.0 | 80 | 54/69/72/78 | 670 | `HOLED` |
| 2004 initial | 71 | 410 | `FITS` |  | 2.0 | 2.0 | 83 | 39/55/73/82 | 670 | `HOLED` |
| 2004 initial | 72 | 490 | `FITS` |  | 0.0 | 0.0 | 82 | 70/79/80/81 | 760 | `TRUE_END` |
| 2004 initial | 73 | 550 | `FITS` |  | 0.0 | 0.0 | 80 | 22/25/77/79 | 730 | `MIXED` |
| 2004 initial | 74 | 600 | `FITS` |  | 0.0 | 0.0 | 76 | 22/62/74/75 | 710 | `MIXED` |
| 2004 initial | 75 | 540 | `FITS` |  | 0.0 | 0.0 | 74 | 68/69/71/73 | 670 | `TRUE_END` |
| 2004 initial | 76 | 505 | `FITS` |  | 0.0 | 2.0 | 68 | 52/62/64/67 | 600 | `TRUE_END` |
| 2004 initial | 77 | 500 | `FITS` |  | 3.0 | 6.0 | 63 | 35/52/55/62 | 520 | `MIXED` |
| 2004 initial | 78 | 470 | `FITS` |  | 0.0 | 16.0 | 68 | 35/49/51/67 | 470 | `MIXED` |
| 2004 initial | 79 | 400 | `FITS` |  | 0.0 | 20.0 | 67 | 18/42/62/66 | 400 | `TRUE_END` |
| 2004 initial | 80 | 440 | `FITS` |  | 0.0 | 0.0 | 68 | 47/47/63/67 | 440 | `TRUE_END` |
| 2004 initial | 81 | 350 | `FITS` |  | 0.0 | 0.0 | 58 | 50/50/51/57 | 470 | `TRUE_END` |
| 2004 initial | 82 | 115 | `FITS` |  | 0.0 | 0.0 | 57 | 39/44/48/56 | 420 | `TRUE_END` |
| 2004 initial | 83 | 60 | `FITS` |  | 0.0 | 0.0 | 52 | 39/45/46/51 | 420 | `TRUE_END` |
| 2004 initial | 84 | 0 | `FITS` |  | 0.0 | 0.0 | 51 | 3/15/19/44 | 320 | `HOLED` |
| 2004 initial | 85 | 50 | `FITS` |  | 0.0 | 0.0 | 50 | 18/20/30/47 | 290 | `HOLED` |
| 2004 initial | 86 | 30 | `FITS` |  | 0.0 | 0.0 | 40 | 10/11/20/39 | 170 | `HOLED` |
| 2004 initial | 87 | 30 | `FITS` |  | 0.0 | 0.0 | 40 | 11/12/18/23 | 130 | `HOLED` |
| 2004 initial | 88 | 110 | `FITS` |  | 0.0 | 0.0 | 56 | 20/28/32/49 | 250 | `HOLED` |
| 2004 initial | 89 | 160 | `FITS` |  | 0.0 | 0.0 | 73 | 47/55/58/69 | 540 | `MIXED` |
| 2004 initial | 90 | 150 | `FITS` |  | 0.0 | 0.0 | 76 | 27/37/49/74 | 350 | `HOLED` |
| 1999 relocation (inter-village south, GIS 9-15) - zero-retreat bound | 9 | 73 | `FITS` |  | 0.0 | 0.0 | 100 | 11/14/17/73 | 120 | `HOLED` |
| 1999 relocation (inter-village south, GIS 9-15) - zero-retreat bound | 10 | 97 | `FITS` |  | 10.0 | 10.0 | 64 | 0/14/16/21 | 110 | `HOLED` |
| 1999 relocation (inter-village south, GIS 9-15) - zero-retreat bound | 11 | 129 | `DROWNS` | yes | 89.0 | 92.0 | 62 | 6/8/9/16 | 50 | `HOLED` |
| 1999 relocation (inter-village south, GIS 9-15) - zero-retreat bound | 12 | 126 | `DROWNS` | yes | 80.0 | 100.0 | 58 | 6/7/10/14 | 40 | `TRUE_END` |
| 1999 relocation (inter-village south, GIS 9-15) - zero-retreat bound | 13 | 125 | `DROWNS` | yes | 71.0 | 98.0 | 16 | 8/9/11/15 | 70 | `TRUE_END` |
| 1999 relocation (inter-village south, GIS 9-15) - zero-retreat bound | 14 | 126 | `DROWNS` | yes | 52.0 | 66.0 | 22 | 7/8/11/19 | 50 | `TRUE_END` |
| 1999 relocation (inter-village south, GIS 9-15) - zero-retreat bound | 15 | 106 | `DROWNS` | yes | 38.0 | 64.0 | 21 | 9/9/10/19 | 60 | `TRUE_END` |
| 1989 relocation (Pea Island, GIS 84-87) - zero-retreat bound | 84 | 163 | `DROWNS` | yes | 31.0 | 34.0 | 51 | 3/15/19/44 | 320 | `HOLED` |
| 1989 relocation (Pea Island, GIS 84-87) - zero-retreat bound | 85 | 165 | `FITS` |  | 0.0 | 4.0 | 50 | 18/20/30/47 | 290 | `HOLED` |
| 1989 relocation (Pea Island, GIS 84-87) - zero-retreat bound | 86 | 205 | `DROWNS` | yes | 22.0 | 38.0 | 40 | 10/11/20/39 | 170 | `HOLED` |
| 1989 relocation (Pea Island, GIS 84-87) - zero-retreat bound | 87 | 113 | `DROWNS` |  | 2.0 | 38.0 | 40 | 11/12/18/23 | 130 | `HOLED` |
