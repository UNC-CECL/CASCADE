# 2026-09-22 — the island offset built from the CoastSat shoreline

Does BRIE's 1996 island offset have to come from the digitised dune line?
This runs 1996–2010 on the **shoreline-derived** offset instead
(`2-brie-offset/1996/shoreline/v1/`, the 1995–1997 CoastSat window mean) and
compares against the matrix runs, which are identical in every other respect.

## Arms

| run | offset source | matrix control |
|---|---|---|
| `1996_2010/edgeBE/HAT_1996_2010_edgeBE_road_bdm_nogroin` | `shoreline/v1` | `matrix/1996_2010/edgeBE/` same name |
| `1996_2010/edgeBE/HAT_1996_2010_edgeBE_road_reloc_bdm_nogroin` | `shoreline/v1` | `matrix/1996_2010/edgeBE/` same name |

Both are full management, no groin, edgeBE, base geometry, calibration wave
climate. The only thing changed is the island offset.

`island_offset_version` in each run's metadata reads `shoreline/v1`; the
controls read `v1`, from before the sources were split on 2026-09-22 (which
means `duneline/v1` — see `2-brie-offset/1996/PROVENANCE.md`).

## READ THIS FIRST: both runs used `offset_mode: asrun`, which divides by ten

`build_island_offset(mode="asrun")` returns `offset_m / 10` — its own
docstring calls it "reproducing the historical unit error exactly", and it is
the config default, so the matrix controls use it too.

So this experiment changed the model input by a **tenth** of what the two
offset files actually differ by. Verified against the runs' initial `x_s`,
which match to the second decimal on all four statistics:

| shoreline − duneline, padded 120 domains | mean | sd | min | max |
|---|---|---|---|---|
| the offset FILES, metres | +2.85 | 23.63 | −58.45 | +28.86 |
| what `asrun` hands the model | +0.29 | 2.36 | −5.84 | +2.89 |
| measured initial `x_s` difference | +0.29 | 2.36 | −5.84 | +2.89 |

**The skill numbers below are therefore not a test of the shoreline offset.**
They are a test of one tenth of it, and the near-identical scores are what
that predicts. A real test needs `offset_mode: metres` (or `detrended`) for
BOTH arms, which means a new control as well — the matrix has none at that
mode.

## Result at one-tenth scale: the dune-derived offset still scores better, barely

Skill against the CoastSat LOESS 10-domain target, interior GIS 2–89, LRR
estimator — **bias is the score here, not r²**:

| | bias (m/yr) | RMSE (m/yr) | LRR r² median | r² below floor |
|---|---|---|---|---|
| duneline (control) | **+0.0229** | **1.1284** | 0.799 | 32 of 90 |
| shoreline | +0.0350 | 1.1560 | 0.811 | 29 of 90 |

The shoreline offset moves the bias **away** from zero, +0.023 → +0.035 m/yr,
and RMSE up by 0.028. It improves the LRR fit quality (r² median up, four
fewer domains below the floor), but that is not the score — see
[[cascade-smoothing-scale-null]].

**Keep the size in view.** 0.012 m/yr of bias over a 14-year window is
**0.17 m** — a sixtieth of a Barrier3D cell. Given the input differed by a
mean 0.29 m and sd 2.36 m after `asrun` scaling, an effect this small is
exactly what should be expected. It says nothing about whether the model is
sensitive to the offset at full scale.

Against the 2010 dune line, the shoreline arm gives observed +16.2 m, modeled
+4.8 m, misfit −11.3 m, RMSE 29.0 m. The control's dune-line misfit is not
recorded — the runner prints it but does not store it — so that row is not a
comparison and should not be read as one.

## The two arms are the same run

Both report **0 relocation events across 0 domains**, and their skill differs
in the fourth decimal. No historical relocation falls inside 1996–2010 as the
events are dated, so `road_bdm` and `road_reloc_bdm` simulate the same thing
here. That matches what the 2010 period showed
([[cascade-2010-period-solved]]).

## Reproduce

```
HAT_ISLAND_OFFSET_SOURCE=shoreline \
HAT_RUN_KIND=experiment HAT_RUN_TAG=2026-09-22-shoreline-offset \
HAT_SOURCE_SINK_PRESET=edgeBE HAT_START_YEAR=1996 \
HAT_SCENARIO=full_management [HAT_RELOCATIONS=true] \
python scripts/hatteras_ms/HAT_hindcast_1984_2024.py
```

`HAT_ISLAND_OFFSET_SOURCE` was added on 2026-09-22 and is read by
`hatteras_site_config`. The runner **refuses** it with `HAT_RUN_KIND=matrix`:
the run name carries no offset token, so a matrix run on a non-default source
would overwrite the control it is being compared against.

`HAT_SCENARIO=full_management` alone gives the `road_bdm` arm, not
`road_reloc_bdm` — `relocations: null` in `hat_run.yaml` inherits from the
scenario and every named scenario sets it False. `HAT_RELOCATIONS=true` is
what adds the `reloc` token.

## The beach-width question, settled

The concern was that the 1996 interior topography is extracted in the
dune-line frame, so swapping the offset alone might count the beach width
twice. Traced through the code, it is real, and it is better named a **frame
mismatch** than double counting:

1. `brie_coupler.offset_shoreline` adds the offset to BOTH `x_t` and `x_s`,
   so `LShoreface = x_s − x_t` is unchanged — the offset is a rigid
   cross-shore translation of the shoreface, not a change to its shape.
   Confirmed: shoreface length is identical in both runs, 2639.4–2640.4 m.
2. Barrier3D hangs its interior grid off `x_s`, and the extractor defines
   **interior row 0 as the cell one landward of the dune crest**
   (`HAT_dune_topo_extractor.interior_row0_line`). Its own plan view places
   "each domain's topo row 0 (ocean side) on canvas row = offset_cells".

So the offset positions a DUNE-referenced row. A dune-derived offset is
therefore the consistent input, and a shoreline-derived one places each
domain's dune at its shoreline's alongshore position — displacing it seaward
by that domain's own beach width, relative to the minimum.

The damage is the alongshore VARIATION, not the mean: on the 90 real domains
the two offsets differ by sd 17.2 m, range −51 to +29 m, which at full scale
is up to five Barrier3D cells of spurious dune displacement.

**The tension is structural, not a bug.** One number, `x_s`, does two jobs:
it is BRIE's shoreline, which its alongshore diffusion acts on and which
should carry the shoreline's shape, and it is Barrier3D's interior anchor,
which should carry the dune's. The two differ by the beach width. Whichever
source is chosen, one of the two roles is wrong by the alongshore variation
in beach width. Using the shoreline offset consistently would mean
re-referencing the extraction to the shoreline — not just a `CURRENT` swap.

## What this does NOT license

Nothing here argues for making the shoreline source `CURRENT` for 1996. Before
that could happen, the 1996 interior topography is extracted in the dune-line
frame, so swapping the offset alone risks counting the beach width twice —
see `2-brie-offset/1996/shoreline/PROVENANCE.md`.
