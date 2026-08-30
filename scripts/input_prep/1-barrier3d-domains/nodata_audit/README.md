# nodata_audit - what the unsurveyed cells do to a CASCADE run

The DEM chain leaves cells that **no survey ever saw**. Barrier3D has no
representation for "unknown", so every one of them is written to the water
sentinel and read as an elevation of exactly **-3.0 m MHW**. The
`<stem>_nodata.npy` mask that records which cells those were is a sidecar:
`hat_topo_version.domain_arrays()` hands `Cascade()` the topography and dune
paths only, so nothing in a run ever opens it.

These five scripts find those cells, measure what they cost, and test whether
each one is a real pond or a lidar dropout.

## Why they are here and not in `0-elevation/`

`0-elevation/3-figures/HAT_plot_dem_holes.py` asks the same question of the
**DEM**. These ask it of the **extracted Barrier3D domains** - after the beach
and dune are clipped off, profiles are sheared straight, water rows are
trimmed and everything below -3.0 m is clamped. Different stage, different
answer, so they sit beside the extractor whose output they read.

## Run order

| | script | writes |
|---|---|---|
| 1 | `HAT_plot_island_nodata.py` | `HAT_island_nodata_<ver>_<year>_padded.png`, a `_D1-8.png` zoom, and **`bracketed_hole_cells.csv`** |
| 2 | `HAT_plot_topo_retained.py` | `HAT_topo_retained_<ver>.png` |
| 3 | `HAT_test_hole_pond_or_dropout.py` | `hole_verdicts.csv`, `dropout_mask/domain_<N>.npy` |
| 4 | `HAT_hole_aerial_chips.py` | `aerial_1996_conflicts/sheet_*.png`, `aerial_review.csv` |
| 5 | `HAT_hole_aerial_picker.py` | fills `aerial_review.csv`, interactively |
| 6 | `HAT_test_hole_pond_or_dropout.py` again | re-decides with the review folded in |
| 7 | `HAT_bridge_dropouts.py` | a NEW dune-topo version with the cleared cells filled |

**None of this is currently live.** Steps 3-7 have not been run against the
`v1` any run now loads, and step 7 has not been run at all since 2026-08-27.
The chain is here to be rerun if a bridged surface is ever wanted; the decision
as of 2026-08-28 is that one is not. See "Re-audit of the current v1" below.

**1 must run before 3**, which needs `bracketed_hole_cells.csv`. 2 is
independent. 4 and 5 both need 3's verdicts; 5 supersedes 4 - the contact
sheets are for reading away from a machine, the picker is for deciding.
**3 must run a second time after 5**: it reads `aerial_review.csv` at load, so
a review pass finished afterwards changes nothing until it is re-run. That
caught us once - the test printed byte-identical output after a completed
review and looked like the review had not saved.

## Everything lands in one place

```
data/hatteras_init/1-barrier3d-domains/<product>/dune-topo/<version>/nodata-audit/
```

`audit_dir()` in each script is the only place that name is spelled. Nothing
here writes into the extraction's own `figures/`, `topography/` or `dunes/`
folders - this is an audit *of* that run, not part of it, and mixing the two
is how a reader loses track of which script owns which file.

One consequence worth knowing: `HAT_island_nodata_*_padded.png` is built on
exactly the canvas `HAT_dune_topo_island_planview_*_padded.png` uses - same
offsets, same padding, same dune row - so the two overlay cell for cell. They
now live one directory apart. That is the price of a clean tree, and the
alternative was leaving audit output loose in the run folder.

## The question the chain answers

An unsurveyed cell read as -3.0 m water is not cosmetic. `barrier3d.FindWidths`
measures the island as the run of land from interior row 0 to **the first water
cell**, and land behind that cell is invisible to the model. So one unsurveyed
cell can delete hundreds of metres of measured barrier from the width Barrier3D
uses.

For the **live** `1984-start/v1` (extracted 2026-08-27T21:12, measured
2026-08-28): **365 of 4,500 profiles** are truncated by an unsurveyed cell, and
**117** of those are holes with measured land on *both* sides - the only ones
where anything is actually lost. They carry 818 cells and hide 46,190 m,
concentrated in D6 (33), D22 (15), D7 (13), D24 (12) and D1 (7).

Earlier revisions of this file quoted 362 / 99 / 731 / 42,210 m. Those describe
the `v1` that was **deleted on 2026-08-27**, not the surface any run now loads.
The two are the same DEM and the same unsurveyed cells - see
"Re-audit of the current v1" below for why the counts moved anyway.

## How a hole is judged

Two automatic references, and a third that needs your eyes:

| | reference | independent of |
|---|---|---|
| A | 2014 NCFMP hydro-flattening stamp | the lidar returns |
| C | shape of the nodata blob | time, and A |
| B | 1996 aerial imagery | both, and contemporaneous |

A is disqualified as an *elevation* source - 94.6% of its coverage in the gap
is two stamped constants - which is exactly what makes it authoritative as a
*water mask*. It is read for whether a value repeats, never for what it means,
so its unknown vertical datum does not matter.

The rule is deliberately asymmetric: a cell is cleared as a dropout only on
**affirmative agreement**, and unknown or conflicting evidence leaves it as
water, which is current behaviour. The DEM changes only where there is evidence
it is wrong. The cost is a one-directional bias - the island stays too narrow
wherever the test cannot resolve a hole - and that belongs in any methods
paragraph built on this.

**A and C conflict on 58 of 99 holes.** With two references and no tiebreaker
the default, not the evidence, decides most of the outcome. That is why B
exists and why step 5 is not optional.

## What it concluded on the first pass — HISTORY, against a deleted tree

Run 2026-08-26 against the `v1` that no longer exists. Kept because the
reference-by-reference verdict below is the reason `aerial_review.csv` is
trusted, and because it is the only place the three references have ever been
scored against each other. Of the 99 bracketed holes:

| reference | POND | DROPOUT | unresolved |
|---|---:|---:|---:|
| A, NCFMP stamp | 41 | 57 | 1 |
| C, blob shape | 70 | 29 | - |
| B, 1996 aerial review | 44 | 8 | 6 unclear |
| **agreed** | **77** | **22** | - |

**B decided it.** A and C conflicted on 58 of 99, so on the majority of holes
the automated pair cancelled and the manual review carried the decision:
52 holes `B breaks tie`, 40 `A+C agree`, 6 `B unclear`, 1 no tiebreak.

**A was the reference that failed.** NCFMP called DROPOUT 57 times and the
imagery agreed with 8 of them. Its modal-fraction test reads ordinary terrain
variation as "not a stamped water polygon", which turns out to say nothing
about whether the ground is dry. C tracked the imagery far better. If this is
written up, C earned its place and A did not.

`HAT_bridge_dropouts.py` then wrote a **`v2`** - that `v1` plus 114 cells filled
by linear interpolation across the 22 cleared holes, in domains 4-7. Mean island
width gained: D4 +13 m, D5 +9 m, D6 +62 m, D7 +54 m, recovering 6,750 m of the
42,210 m those 99 holes were hiding from `FindWidths`. **Both that `v1` and that
`v2` were deleted on 2026-08-27** (see `1-barrier3d-domains/LINEAGE.md`). No
bridged surface exists today and none is planned - see below.

The headline conclusion is close to the null that was flagged as possible from
the start: **the sub-MHW cells in this DEM are overwhelmingly real ponds.** A
cell is unsurveyed here only because 1996 ALACE, 2009 USACE *and* 2014
Post-Sandy all failed at it, and three surveys failing in one spot is what
water looks like. The genuine dropouts are a small, localised set in D4-D7,
none of which carries a road - the 1984 setbacks cover domains 9-90 - so no
roadway, drowning or relocation result depends on any of this.

## Re-audit of the current v1 — 2026-08-28

The live `v1` was re-extracted on 2026-08-27 from the **same DEM and the same
`npy-arrays/`** (unchanged since 2026-08-26 10:12) against a **new pick set**,
`picks/HAT_dune_search_windows_v1.json`. Only the cross-shore window origin
moved. The question was whether the 2026-08-26 verdicts still cover it.

**The 99 old holes carry over exactly.** All 99 `(domain, profile)` keys are
still bracketed holes, all 99 with identical cell counts, and the subtotal
reproduces the old pass to the cell: 731 cells, 42,210 m. `aerial_review.csv`
transfers with no re-review, exactly as `aerial-review/README.md` predicted.

**18 holes are new** - not new ground, but profiles where the moved window
changed which cell `FindWidths` stops at, so the old pass never saw them:

| | cells | hidden |
|---|---:|---:|
| the 99 carried-over holes | 731 | 42,210 m |
| the 18 newly exposed holes | 87 | 3,980 m |
| **live v1 total** | **818** | **46,190 m** |

16 of the 18 sit in domains 9-90, so the first pass's "none of which carries a
road" exit line does not transfer on location and had to be re-established on
substance.

### The neighbour-consistency test, and why none of the 18 was bridged

A cheap fourth check settles them without a manual review. For each hole
profile, compare its `FindWidths` width to the **median width of the profiles
within ±4 alongshore that are truncated by *measured* water** - neighbours
whose stopping cell is real, surveyed, sub-MHW ground. A dropout punching a
false hole in a wide barrier reads far narrower than that median. A hole lying
inside a genuine water body reads the same as it.

**15 of 18 came back consistent, ratios 0.91-1.35.** They are unsurveyed cells
sitting *inside* surveyed water, not holes in a barrier.

The clearest case is D73/D74, which on hidden-metres alone looked like the worst
of the set - four holes of 1-2 cells apiece, hiding 560-640 m each, in the road
reach. D73 profiles 39-49 are *all* truncated at 220-260 m and only 44, 45 and
46 stop at an unsurveyed cell; the other nine stop at measured water of -0.16 to
-1.40 m at the same row. D74 profiles 0-9 are all truncated at 230-290 m, only
prof 7 by nodata, flanked by measured water at -0.98 and -1.03. The 220 m width
is the island's real shape there, and the "hidden" land is the far shore of a
bay that is equally invisible behind genuine water on every neighbour. Bridging
those cells would interpolate a land surface across a measured bay.

The 3 that flagged anomalous do not survive inspection either. D3/2 and D22/41
have **no** measured-water neighbours to form a median, so the test returns
nothing rather than a verdict. D8/46 sits in a bimodal domain and matches its
wide mode (1360 m against neighbours of 1370 and 1250), so it truncates nothing
of consequence. That leaves **D22 prof 41, 420 m on one profile**, as the only
genuinely unresolved hole in the 18 - well inside the uncertainty the first
pass already accepted with its 6 UNCLEAR holes.

### Decision: the live v1 is not bridged, and stays that way

No `v2`. Steps 3-7 were not run and `HAT_bridge_dropouts.py` was not run. The
whole 18-hole delta is a **pick artefact** - same DEM, same cells, a window
origin that moved - and the only surface CASCADE loads is the measured one.

This keeps the same one-directional bias the asymmetric rule was chosen for:
**~46,190 m of measured barrier across 117 profiles is hidden from `FindWidths`
by unsurveyed cells in the live `1984-start/v1`, a known and quantified
residual, not an oversight.** Almost all of it is land that genuine measured
water hides on the neighbouring profiles anyway. The bias runs toward a
narrower island, which makes the barrier look more vulnerable rather than less.

The neighbour-consistency rule is stated above in full and has **no script in
this folder** - it was run ad hoc on 2026-08-28. Reproducing it needs only the
topography and nodata arrays.

## The iterative residual — HISTORY, and now moot

This described the deleted `v2` and is kept only because it is the reason
bridging is known to be a multi-pass correction rather than a one-shot fix.
**Nothing in the live `v1` carries this residual, because the live `v1` is not
bridged.** Its residual is the 46,190 m figure above.

Bridging exposed a second layer, because the analysis only ever examined the
FIRST water cell on each profile. Fill that cell and `FindWidths` advances to
the next one, which on most of these profiles is also unsurveyed and was
therefore invisible to the original pass:

| of the 22 bridged profiles | |
|---|---:|
| now stop at genuine water - fully fixed | 6 |
| now stop at another unsurveyed cell | 16 |
| of those, bracketed by measured land, so bridgeable | 8 |
| further land still hidden behind them | **950 m** |

That is why truncated profiles fell only 362 to 356 rather than by 22. The
width gain is real regardless - the new stopping cell is further landward -
but the correction is ITERATIVE and one pass does not finish it.

**Left as is, on purpose.** A second round would need another full
three-reference pass, including a manual review, over 8 holes, to recover
950 m against the 6,750 m already recovered - about 14% more. The residual is
smaller than the uncertainty already carried by the 6 UNCLEAR holes and the
conservative default, and it runs in the same direction as every other
unresolved cell: the island stays slightly too narrow, which makes the barrier
look more vulnerable rather than less.

That 950 m applied to the deleted `v2` only. **Do not quote it for the live
`v1`** - the figure to quote is 46,190 m across 117 profiles, from the re-audit
above. The machinery to close a residual iteratively is still in this folder if
a bridged surface is ever wanted: rerun the chain from step 1 against whatever
version step 7 writes.

## Nothing here is committed

`.gitignore` covers `*.png` (line 158) and `*.npz`. The CSVs and this README are
tracked; the figures and the chip cache are not. Everything is regenerable from
the arrays and the D: drive sources.

## Do not hardcode the product or version

Both resolve through `scripts/hat_topo_version.py`. Set `TOPO_PRODUCT` at the
top of a script and leave `VERSION_OVERRIDE = None`. The header of
`hat_topo_version.py` records the incident that rule comes from: four road
scripts pinned a version string and kept reading a stale interior for 18
domains without erroring.
