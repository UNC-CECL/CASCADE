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
| 1 | `HAT_bracketed_hole_cells.py` | **`bracketed_hole_cells.csv`** |
| 1 | `HAT_plot_island_nodata.py` | `HAT_island_nodata_<ver>_<year>_padded.png`, a `_D1-8.png` zoom |
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

**`HAT_bracketed_hole_cells.py` must run before 3**, which needs
`bracketed_hole_cells.csv`; the two step-1 scripts and 2 are independent.
Until 2026-10-01 this table credited the CSV to `HAT_plot_island_nodata.py`,
which never wrote it: the original writer was never committed (see
"HAT_bracketed_hole_cells.py" below). 4 and 5 both need 3's verdicts; 5 supersedes 4 - the contact
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

**Correction, 2026-10-01.** The 18 are not an effect of the moved window.
`HAT_bracketed_hole_cells.py`, run on this same live `v1`, finds the strict
bracketed set (the cell past the hole is measured land) to be exactly the
original 99, at the same npy cells, and the 18 are profiles where the
unsurveyed run ends in *measured water* with land somewhere further back. They
are a looser count of "land hidden behind", not new holes. That agrees with the
neighbour test below, which found them to be unsurveyed cells inside surveyed
water. The dune-topo `v2` gives the same 99 and the same 18.

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

Both resolve through `scripts/site_layer/hat_topo_version.py`. Set `TOPO_PRODUCT` at the
top of a script and leave `VERSION_OVERRIDE = None`. The header of
`hat_topo_version.py` records the incident that rule comes from: four road
scripts pinned a version string and kept reading a stale interior for 18
domains without erroring.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### HAT_bracketed_hole_cells.py

The unsurveyed holes that stop FindWidths with measured land on both sides, cell by cell.

Written 2026-10-01 to replace a writer that was never committed. Four scripts
read its output (`HAT_test_hole_pond_or_dropout.py`, `HAT_hole_aerial_chips.py`,
`HAT_hole_aerial_picker.py`, `HAT_bridge_dropouts.py`). Until then the only copy
was `aerial-review/bracketed_hole_cells_v1.csv`, kept from the deleted tree.

**What counts as a bracketed hole.** On each profile, find the first cell at or
below 0 m MHW, the one `FindWidths` stops at. If that cell is unsurveyed, the
hole is the unbroken run of unsurveyed cells starting there. It is bracketed
when the cell before the run is land (interior row ≥ 1) and the cell after it is
*measured land* (above 0 m MHW, surveyed). Those are the holes
`HAT_bridge_dropouts.py` can fill, since it interpolates between exactly those
two cells. One hole per profile, the first, as before.

Profiles whose run ends in measured water, with land further back, are counted
and printed ("hide land, not bracketed") but not written. They are the "18 new
holes" of the 08-28 re-audit (see the correction there).

**Where each cell is on the raw grid.** `npy_row`/`npy_col` index the product's
`npy-arrays/domain_<N>.npy`. The script imports the extractor and rebuilds each
domain from those arrays and the version's own picks file. The steps are
`load_profiles`, `find_dunes`, `build_interior`, `remove_water_rows` and
`interior_row0_line`. It then **refuses** unless the rebuilt interior equals the
saved topography and nodata mask. A version whose arrays were edited after
extraction (a row insert, a bridge) therefore stops the script rather than
getting wrong coordinates. Saved row `r` on profile `p` is processed column
`row0[p] + r`, raw ocean-first column `+ c0 + shear[p]`, and so
`npy_col = W - 1 - that` and `npy_row = 49 - p` (ocean on the right, alongshore
flip). `utm_x`/`utm_y` are the cell centre, `origin + 10 (col + ½)`, with the
origin taken from the elevation product's `2-resampled-10m/resample_audit.csv`.

**Checked against the surviving file, 2026-10-01.** On the live `1984-start/v1`
the output matches `aerial-review/bracketed_hole_cells_v1.csv` row for row:
99 holes, 731 cells, 42,210 m hidden, with every `domain`, `profile`, `npy_row`,
`npy_col`, `utm_x` and `utm_y` identical. Only `interior_row` differs, by a
constant 1-3 rows per profile in D1, D3, D4, D6, D7 and D24, which is the
moved window origin. Every interior cell of v1 (586,840) and v2 (581,185)
maps to the raw cell it was extracted from. The only exceptions are cells the
shear shifted in and cells past the profile end, neither of which can be in a
bracketed hole.

**The review diff is automatic.** `aerial-review/README.md` asks for the
new cells to be diffed against the reviewed ones before old verdicts are
trusted. The script prints that diff: holes unchanged, moved, new and gone.
For `v2` it is 99 unchanged and none moved, new or gone, so `aerial_review.csv`
applies to `v2` as is.

### HAT_bridge_dropouts.py

Fill the unsurveyed cells that three references agreed are survey dropouts, as a new dune-topo version.

From the script's original header:

```text
Fill the unsurveyed cells that three references agreed are survey dropouts, and
write the result as a NEW dune-topo version.

WHAT IT CHANGES, AND WHAT IT REFUSES TO CHANGE
Only the cells that HAT_test_hole_pond_or_dropout.py cleared as DROPOUT. Each
such hole is bracketed by MEASURED land on both sides along its own profile, so
the fill is a straight linear interpolation between two measurements. No value
is invented beyond what those two measurements already imply, and nothing that
was judged POND or left UNCLEAR is touched - those keep the -3.0 m sentinel.

    interior row:   r-1      r    r+1   r+2    r+3
    v1 (dam):      0.111   -0.3  -0.3  -0.3   0.116     <- three unsurveyed
    v2 (dam):      0.111   0.112 0.113 0.114  0.116     <- interpolated

THE NODATA MASK IS DELIBERATELY LEFT ALONE
A bridged cell is still a cell no survey ever saw. Its value is now an
interpolation, not a measurement, and any consumer asking "was this measured?"
must still get NO. So `<stem>_nodata.npy` is copied through unchanged and a
SECOND mask, `<stem>_bridged.npy`, records which cells were filled. Two masks,
two different questions:

    nodata   True = no survey saw this cell            (unchanged from v1)
    bridged  True = and its value is now interpolated  (new in v2)

Collapsing them into one would destroy the only record that these values are
inferred, which is the same conflation that put three roadways underwater at
t = 0 in an earlier product.

WHY A NEW VERSION RATHER THAN AN EDIT
v1 is what the extractor produced from the DEM. v2 is v1 with a documented,
evidence-backed repair applied on top. Keeping both means the repair can be
audited, reverted, or re-derived, and any figure or run can say which it used.

    v1  extraction, untouched
    v2  v1 + bridged dropouts

Interior SHAPES and interior ROW 0 are identical between them - this only
rewrites values inside existing arrays. That matters because every road setback
is measured from interior row 0, so the 1984 setback CSVs stay valid across the
version bump and do not need re-measuring.

TWO THINGS HAVE TO MOVE FOR v2 TO TAKE EFFECT
hat_topo_version.resolve_version puts the EXTRACTOR'S `VERSION` literal ABOVE
the CURRENT file. Writing CURRENT alone is silently ineffective. So this script
writes both, and says so:

    dune-topo/CURRENT                        -> v2
    HAT_dune_topo_extractor.py  VERSION      -> "v2"

THE CLOBBER RISK, STATED PLAINLY
Because the extractor now says VERSION = "v2", re-running it would overwrite
these arrays with a fresh, UNBRIDGED extraction. That is recoverable in one
command - this script reads v1 and rebuilds v2 - but it is silent, so
BRIDGE_MANIFEST.txt in the run folder says it too.

INPUT   dune-topo/<src>/topography, dunes
        dune-topo/<src>/nodata-audit/hole_verdicts.csv
        dune-topo/<src>/nodata-audit/bracketed_hole_cells.csv

OUTPUT  dune-topo/<dst>/topography/domain_<N>_{topography,nodata,bridged}.npy
        dune-topo/<dst>/dunes/domain_<N>_dune.npy      (copied unchanged)
        dune-topo/<dst>/BRIDGE_MANIFEST.txt
```

Notes that were in the code:

```text
--- activate ------------------------------------------------------------
Picks FIRST. The moment the extractor's VERSION literal moves, every
caller that resolves a window file through it looks for this version's
pick set. Writing it after the bump would leave a window -- however
short -- in which the tree resolves to a path that does not exist.
```

<details><summary>Function notes (the original docstrings)</summary>

**`carry_picks_forward()`**

```text
Give the new version its own pick set, seeded from the source version.

WHY THIS IS PART OF THE VERSION BUMP AND NOT AN AFTERTHOUGHT

The extractor derives its window file from the VERSION literal:

    PICK_SET    = RUN_NAME = VERSION
    WINDOW_JSON = PICKS_DIR / f"HAT_dune_search_windows_{PICK_SET}.json"

So bumping VERSION without writing that file leaves every consumer that
resolves picks through the extractor pointing at a path that does not
exist. That is not hypothetical: v1 -> v2 did exactly this on 2026-08-26,
and `HAT_road_offset_from_dune_start.py` -- which re-derives interior row 0
from these windows in order to measure the road setbacks -- could not be
re-run at all afterwards. The setbacks already on disk stayed correct; they
simply could not be regenerated.

COPYING IS THE RIGHT ANSWER HERE, NOT SHARING

The extractor's own guard prescribes this procedure verbatim: "set
PICK_SET = RUN_NAME and copy the old file to
HAT_dune_search_windows_<RUN_NAME>.json first". The reason is that
`save_windows()` writes back to WINDOW_JSON after EVERY domain during a
pick pass, so pointing a re-pick at a shared file destroys the earlier
version's picks one domain at a time. A per-version copy makes that
impossible.

And the copy is semantically honest for a BRIDGED version: bridging fills
unsurveyed cells inside existing arrays. It does not move a dune, so it
cannot move a dune window. v2's windows ARE v1's windows.

PROVENANCE IS RECORDED RATHER THAN IMPLIED

A bare copy would leave a file that looks like it was picked for this
version. `_meta.inherited_from` says otherwise, so a later reader can tell
an inherited pick set from a re-picked one. It goes INSIDE `_meta` on
purpose: consumers count domains as "every key that is not `_meta`", so a
second top-level underscore key would be counted as a domain.

Returns:
    (destination path, "written" | "exists") or None if the source has no
    pick file to carry.
```

**`bridge_profile()`**

```text
Linear fill of `rows` in one profile. Returns False if not bracketed.

Refuses rather than guesses. A hole touching either end of the array has no
measured value on one side, so there is nothing to interpolate between -
that is an extrapolation, which is the thing this whole exercise exists to
avoid.
```

</details>

### HAT_hole_aerial_chips.py

Reference B: 1996 aerial chips for the holes where references A and C disagree.

From the script's original header:

```text
Reference B: 1996 aerial chips for the holes where references A and C disagree.

WHY ONLY THE CONFLICTS
HAT_test_hole_pond_or_dropout.py runs two references over the 99 unsurveyed
holes that truncate the model island:

    A  the 2014 NCFMP hydro-flattening stamp
    C  the shape of the nodata blob

They agree on 40 holes and CONFLICT on 58. With only two references and no
tiebreaker, the conservative default - not the evidence - decides those 58, so
the test asserts "pond" on cells where its own NCFMP vote says otherwise. That
is what this fixes, and it is why the aerial pass is worth its cost after all:
58 chips, not the 99 the full pass would have needed.

WHY THE AERIAL IS THE RIGHT TIEBREAKER
It is the only reference contemporaneous with the survey. A is 2014 standing in
for 1996 across 18 years of marsh change; C has no date at all. The 1996 frames
were flown the same year as the ALACE survey that this DEM's beach comes from,
so a pond visible in them is a pond the lidar would have been looking at.

WHAT YOU ARE JUDGING
Each chip is the 1996 imagery around one hole, with the unsurveyed cells drawn
as an outline. The question is only:

    is there standing water inside the outline?

    yes            -> POND      the -3.0 m sentinel is right, leave it
    no, it is land -> DROPOUT   the lidar failed over ground, bridge it
    cannot tell    -> UNCLEAR   no vote; the conservative default keeps it water

Do not judge the ring, and do not try to reconcile it with the other two votes -
the whole point is an independent third opinion.

COORDINATES
The frames are NAD83 / North Carolina State Plane in US SURVEY FEET, 3-band
RGB, 1 ft pixels. The domain rasters are EPSG 3725, UTM 18N, metres. Every
transform goes through pyproj from the frame's own CRS, never a hardcoded
factor - a foot is not 0.3048 m in this projection, it is 0.304800609601219 m,
and the difference over 3 million feet of easting is metres.

Frames overlap. The one chosen for a hole is the one giving the cleanest chip -
these are scanned frames with a black surround baked into the raster, so bounds
margin is not a guide. See pick_frame.

INPUT   dune-topo/<version>/hole_verdicts.csv          which holes conflict
        dune-topo/<version>/bracketed_hole_cells.csv   the cells to outline
        D:\Hatteras_GIS\Aerial\1996_henderson\1996_georef_TIF\*.tif

OUTPUT  dune-topo/<version>/figures/aerial_1996_conflicts/
            sheet_D<a>-<b>.png        contact sheets, 12 chips each
            aerial_review.csv         one row per hole, blank verdict column
```

Notes that were in the code:

```text
Every output of this folder lands under one directory beside the extraction it
describes, rather than being scattered through the run folder it did not
produce. audit_dir() is the only place that name is spelled.
```

```text
Headless only when this file is the program. HAT_hole_aerial_picker.py
imports the helpers below and needs a real interactive backend, so the
module must not claim Agg at import time.
```

```text
--- review sheet --------------------------------------------------------
The PICKER owns this file. Writing a fresh template over a reviewed one
would silently destroy the pass it took to fill in, so an existing file
with any verdict in it is left alone. Delete it to start over.
```

<details><summary>Function notes (the original docstrings)</summary>

**`pick_frame()`**

```text
The covering frame that gives the CLEANEST chip, not the widest margin.

Choosing by distance from the frame's bounds was the first version and it
put a black wedge across a third of the D6 chips. These are scanned aerial
frames: each carries a black surround baked INTO the raster, so a point can
sit far inside the bounds and still land on unexposed film. Frames overlap,
so the fix is to read the candidate window from each and keep the one with
the least black - which costs a few extra reads for 58 chips and nothing
else.
```

</details>

### HAT_hole_aerial_picker.py

Review the aerial chips one hole at a time: a keystroke per hole records pond, dropout or unclear.

From the script's original header:

```text
Interactive review of the 1996 aerial chips: reference B, one keystroke per hole.

WHAT YOU ARE DECIDING
For each unsurveyed hole where the NCFMP stamp (A) and the blob shape (C)
disagree, one question:

    is there standing water INSIDE the yellow outline, in 1996?

    w   POND      water. The -3.0 m sentinel is right; leave the DEM alone.
    g   DROPOUT   ground. The lidar failed over land; the cell is bridgeable.
    u   UNCLEAR   no vote. The conservative default keeps it water.

UNCLEAR is a real answer and costs nothing. It is there so you are never forced
to guess, which is the failure mode that would quietly turn a review into a
coin-flip dressed as evidence.

Judge only what is inside the outline. Do not try to reconcile it with the two
votes shown in the title - B is worth having precisely because it is an
independent opinion, and a tiebreaker that has already read the other votes is
not one.

THE TRAP IN THE 1996 IMAGERY, AND WHAT TO DO ABOUT IT
The 1996 frames have a narrow tonal range. Wet sand, damp marsh and shallow
standing water all land on the same mid-grey, and a dark patch is as often
shadow or dense vegetation as it is water. Texture and a closed edge separate
water better than darkness does.

When 1996 will not resolve, press z to cycle the other sources on the drive:

    2019 NGS      0.30 m colour, covers every hole - the most legible
    2008 IOCM     0.50 m colour
    2014          0.35 m colour
    2018 NAIP     0.60 m colour AND near-infrared, plus an NDWI view
    2004          colour, closest in time to 1996 after 1996 itself

Only 2018 carries a verified NIR band; 2014 and 2016 ship four bands but the
fourth is not infrared, so they get no NDWI. See verify_nir().

NDWI is (green - NIR) / (green + NIR): BLUE is water, RED is land. Water
absorbs near-infrared almost completely, so it separates open water from wet
sand in a way no visible-band image can. It is the view to reach for on exactly
the chips that are hard - but it exists only for 2018.

The catch, and it is a real one: only 1996 is contemporaneous with the survey
whose dropouts are in question. Everything else answers "is there a pond here
NOW", which is strong evidence for a pond that has sat in one place for
decades and weak evidence for a marsh pool that migrates. Let the later
imagery break a tie; do not let it overrule a clear 1996 view. The 2018 NIR
also covers only domains 1-8, which happens to be where 48 of the 57 conflicts
are.

KEYS
----
    w / g / u     verdict, then auto-advance
    left / right  move to another hole without deciding
    n             jump to the next hole with no verdict
    z / x         cycle imagery forward / back for this hole
    r             clear this hole's verdict
    q             quit

Every verdict is written to aerial_review.csv IMMEDIATELY, the same way
HAT_dune_topo_extractor.save_windows writes after every domain. Quit whenever
you like and re-run to resume; 58 holes is more than one sitting.

Any write that would REDUCE the number of verdicts on disk copies the old file
to aerial_review.<timestamp>.bak.csv first. That guard exists because the file
was once blanked by a helper that had checked it was empty earlier in the
session and did not re-check before overwriting. A hand-entered review cannot
be regenerated from anything, so it does not get overwritten silently.

CACHE
Rendering a chip means reading a window from up to 33 scanned frames to find
the one with the least black surround, which is slow enough to feel in an
interactive loop. So chips are rendered once into figures/aerial_1996_conflicts
/chip_cache/ as PNGs plus an index of the outline geometry, and the picker
reads those. Delete the folder to force a rebuild. The PNGs are covered by
.gitignore's *.png rule, like every other figure here.

INPUT   dune-topo/<version>/hole_verdicts.csv
        dune-topo/<version>/bracketed_hole_cells.csv
        the 1996 frames, via HAT_hole_aerial_chips

OUTPUT  dune-topo/<version>/figures/aerial_1996_conflicts/aerial_review.csv
            the aerial_verdict column, filled in
```

Notes that were in the code:

```text
Chip half-widths are PER HOLE, not fixed. The truncating holes run from 1 to
56 cells, so a fixed 90 m view fits a 3-cell hole in D6 and cuts a 560 m hole
in D1 clean off the edge - which is what the first version did. Each hole gets
a close view sized to its own footprint and a wide view for context.
```

```text
--- IMAGERY SOURCES ----------------------------------------------------
1996 is the only CONTEMPORANEOUS evidence and stays the primary view: it was
flown the same year as the ALACE survey whose dropouts are in question. The
rest are CORROBORATION, and they answer a slightly different question - "is
there a pond here in 2019" rather than "in 1996". For a pond that has sat in
the same place for decades that is strong support; for a marsh pool that
migrates it is weaker. Cycle to them when 1996 is ambiguous, which is often,
because 1996 is a scanned panchromatic-ish frame where wet sand, damp marsh
and shallow water all land on the same mid-grey.

Near-infrared is the band that actually settles a hard chip. Water absorbs it
almost completely, so a pond is near-black in NIR while wet sand stays bright
- the discrimination the visible bands cannot make. A source with real NIR
gets an extra NDWI view, (green - NIR) / (green + NIR), positive over water.

WHICH BAND IS NIR IS VERIFIED, NOT CONFIGURED. Three of these datasets ship
four bands, and only one of them is genuinely RGB+NIR:

2018 NAIP   band 4 veg/water 3.33 against 2.06 for the best visible  ->  NIR
2014        band 4 veg/water 0.86 - BRIGHTER over water              ->  not NIR
2016        band 4 veg/water 1.16                                    ->  not NIR
2019 NGS    band 4 has 12 distinct values                            ->  alpha mask

Assuming band 4 was NIR produced an NDWI panel for 2014 that called open
water "land" while its own photo showed a pond. So `nir` below is a CANDIDATE
index, and verify_nir() has to agree before any NDWI view is written.

NOT Google Earth: its imagery is licensed, bulk tile extraction breaches its
terms, and it is lower resolution here than the 2019 NGS tiles and has no NIR
band at all. Nothing it offers is missing from this list.
```

```text
Every output of this folder lands under one directory beside the extraction it
describes, rather than being scattered through the run folder it did not
produce. audit_dir() is the only place that name is spelled.
```

```text
Backends that deliver key_press_event to a live window. Anything else - Agg,
or PyCharm's SciView (module://backend_interagg), which renders a figure as a
static image in a tool pane - makes plt.show() return immediately and the
script exit with nothing reviewed and no error. That is not a failure mode a
user should have to diagnose, so it is checked and fixed here.
```

```text
'q' is left alone: matplotlib's quit closes the window, which is what the
picker wants it to do anyway.
```

```text
Clip to the image. Without this a polygon reaching past the chip
autoscales the axes and the imagery shrinks into a white field.
```

```text
A scale bar of a round length near a quarter of the view, so it stays
useful whether the chip is 180 m or 1.5 km across.
```

<details><summary>Function notes (the original docstrings)</summary>

**`free_the_keys()`**

```text
Drop matplotlib's own bindings for the keys this picker uses.

Found by collision, not by reading the docs: 'g' toggles a grid over the
imagery, 'r' resets the view, and left/right walk matplotlib's own view
history. All three fire alongside the picker's handler, so a verdict could
also silently rescale the chip you were judging.
```

**`verify_nir()`**

```text
Is band `cand` really near-infrared? (ok, ratio, best_visible_ratio).

NIR is bright over vegetation and near-black over water. So classify from
the VISIBLE bands only - vegetation as green-dominant and above-median
brightness, water as the darkest decile - then ask which band separates
them best. A real NIR band wins by a clear margin; a mislabelled fourth
band does not, and an alpha mask is not even monotonic.
```

**`ImagerySource()`**

```text
One aerial dataset, sampled in the domain CRS regardless of its own.

Each dataset carries its own projection and linear unit - the 1996 frames
are State Plane in US survey feet, the rest are UTM in metres - so the
transform and the metres-to-unit factor are read from each file rather
than assumed anywhere.
```

**`build_cache()`**

```text
Render every conflict hole in every source that covers it.

Was 1996 only. Extended because the 1996 frames cannot separate wet sand
from shallow water, which is most of what makes a hole hard to call, and
the drive already holds colour at 0.30 m and near-infrared at 0.35 m.
```

**`save_review()`**

```text
Write the review file, backing it up first if this would LOSE verdicts.

A review pass is hand-entered judgement that cannot be regenerated from
anything. This file has already been destroyed once - blanked by a helper
that had checked it was empty earlier in the same session and did not
re-check before overwriting - so any write that reduces the verdict count
now leaves a timestamped copy behind and says so.

The check is on the count rather than on content because that is the only
thing that matters here: a write that keeps or adds verdicts is the normal
path, and a write that drops them is either a deliberate reset or a bug,
and both deserve a copy on disk.
```

**`chip()`**

```text
(rgb, ndwi|None, fx, fy, half, tile) or None if nothing covers it.

Picks the covering tile with the least black, for the reason in
HAT_hole_aerial_chips.pick_frame: scanned frames carry an unexposed
surround baked into the raster, so bounds margin is not a guide.
```

</details>

### HAT_plot_island_nodata.py

The island plan view reduced to one question: where is the unsurveyed ground in what CASCADE runs on?

From the script's original header:

```text
The island plan view, stripped to one question: where is the unsurveyed ground
in what CASCADE actually runs on?

WHY A SECOND PLAN VIEW
HAT_dune_topo_island_planview_<run>_<year>_padded.png shows elevation, and at
that colour scale an unsurveyed cell is indistinguishable from water: both sit
at the -3.0 m sentinel and both render as the same dark blue. That is not a
flaw in the figure - it is the honest consequence of Barrier3D having no
representation for "unknown" - but it means the elevation view cannot answer
"is any of this no-data affecting my model".

So this draws the same canvas, at the same offsets, with the same padding, and
throws the elevation away. Land is one flat grey, water another, and the only
thing with a colour is the no-data.

THE CANVAS IS THE SAME ONE, DELIBERATELY
Every geometric rule here is copied from _build_island_canvas() in
HAT_dune_topo_extractor.py so the two figures overlay cell for cell:

    offsets     2-brie-offset/<year>/Island_Dune_Offsets_*.csv,
                metres, seaward positive, row 0 = domain 1 (Cape Point).
                A 120-row file is stripped of its 15 buffer domains per end.
    origin      round(offset_m / 10) - the canvas row interior row 0 lands on
    padding     every domain padded landward to ISLAND_PAD_ROWS = 200 cells,
                or cropped to it
    dune        written into canvas row origin - 1, one row, matching
                ISLAND_INCLUDE_DUNE
    columns     domains concatenated in ascending order, 50 profiles each,
                no per-domain flip - the arrays already run south to north

If those constants move in the extractor they must move here. The alternative -
importing the extractor - drags in its interactive picker and its own
TOPO_PRODUCT literal, which is how the figure scripts came to disagree with the
road scripts before.

TWO SHADES OF RED, AND THE DIFFERENCE MATTERS
    unsurveyed              a cell CASCADE reads as -3.0 m water that was
                            never measured
    unsurveyed, truncating   the same, AND it is the first water cell on its
                            profile, so barrier3d.FindWidths stops there

The second is the one with a demonstrable effect. FindWidths measures the
island as the run of land from interior row 0 to the first water cell, and land
behind that cell is invisible to the model. A truncating unsurveyed cell
therefore deletes every real, measured cell behind it from the island width
Barrier3D uses. The rest of the red is inside the barrier and may or may not
matter, depending on what the run does with it.

Panel (b) counts both per domain, so nothing is missed at 45 km: a single
unsurveyed cell is a third of a pixel wide in panel (a) and can be invisible
there while still being a real bar below.

INPUT   <product>/dune-topo/<version>/topography/domain_<N>_topography.npy  dam
        <product>/dune-topo/<version>/topography/domain_<N>_nodata.npy     bool
        <product>/dune-topo/<version>/dunes/domain_<N>_dune.npy            dam
        2-brie-offset/<year>/Island_Dune_Offsets_*.csv            m

        Product and version resolve through scripts/site_layer/hat_topo_version.py.

OUTPUT  <product>/dune-topo/<version>/HAT_dune_topo_island_nodata_<version>_<year>_padded.png
        Written beside the elevation plan view it is meant to be compared with.
```

Notes that were in the code:

```text
Every output of this folder lands under one directory beside the extraction it
describes, rather than being scattered through the run folder it did not
produce. audit_dir() is the only place that name is spelled.
```

```text
The CURRENT build (2026-09-18). This took the first sorted match under
2-brie-offset/, which for 1984 and 2004 is superseded_20260915_flat/ --
a build that differs from the current one.
```

```text
Measured land behind the truncating cell. Flagged ONLY on
profiles an unsurveyed cell truncated: land behind a genuine
water cell is also invisible to FindWidths, but that is a real
bay, not a data artefact, and colouring it here would blame the
survey for the island's actual shape.
```

```text
Truncated profiles are a count of profiles, not of cells, so they get
their own axis rather than being stacked onto a bar they do not belong on.
```

<details><summary>Function notes (the original docstrings)</summary>

**`pad_or_crop()`**

```text
Pad landward to ISLAND_PAD_ROWS with the sentinel, or crop to it.

Returns the padded arrays and how many real land cells the crop discarded,
which the extractor also reports - a domain wider than 2000 m loses its bay
margin to this figure's frame, not to the model.
```

**`classify()`**

```text
Category codes for one domain, plus per-profile truncation flags.

Two things are separated here, and the distinction is the whole point of
the figure:

INSIDE the island envelope - row 0 up to that profile's last cell above
MHW - an unsurveyed cell is a hole in the barrier. Beyond it, the same cell
is open sound the survey never flew over, which is the expected state of a
lidar return over water and changes nothing about the barrier.

first_water is barrier3d.FindWidths' stopping point: the first cell at or
below sea level, scanning landward from interior row 0. Sea level is 0 in
the Lagrangian frame, and these arrays are MHW-relative. When that cell is
unsurveyed, the profile's island is truncated there and every measured cell
behind it is invisible to the model.
```

**`draw_zoom()`**

```text
The same canvas, cropped to a few domains, at one pixel per cell.

The FindWidths boundary is drawn on top as a step line. Above it, on a
truncated profile, is measured land the model cannot see - which is the
whole reason this zoom exists. The step is drawn per profile rather than
smoothed: it moves by whole cells, and interpolating it would suggest a
precision the 10 m grid does not have.
```

</details>

### HAT_plot_topo_retained.py

What the DEM's unsurveyed ground becomes in a CASCADE input domain, and what Barrier3D does with it at t = 0.

From the script's original header:

```text
What the DEM's unsurveyed ground becomes once it is a CASCADE input domain, and
what Barrier3D does with it at t = 0.

THE ANSWER TO "ARE THERE CELLS THAT STAY NO-DATA AT MODEL START"
No. There is no such state to stay in. Barrier3D has one float per cell and no
representation for "unknown", so the extractor writes every unsurveyed cell to
SENTINEL_WATER_M and CASCADE reads it as an elevation of exactly -3.0 m MHW.
The `<stem>_nodata.npy` mask that records which cells those were is a sidecar:
hat_topo_version.domain_arrays() hands Cascade() the topography and dune paths
only, so nothing in the run ever opens it.

So the question is not whether no-data survives. It is what the model believes
instead, and the answer is: open water, at the bottom of the clamp, in the
middle of the barrier.

WHY THAT IS NOT COSMETIC
barrier3d.FindWidths - transcribed in roadway.interior_widths - measures the
island as the run of land from interior row 0 to THE FIRST WATER CELL. Land
behind a water cell is invisible to it. One unsurveyed cell at row k therefore
truncates that profile's island at row k-1 and discards every real, measured
cell behind it.

Panel (c) is that cost. It is an UPPER BOUND on the damage, not an estimate:
it compares the width Barrier3D actually sees against the width it would see if
every unsurveyed cell turned out to be land. Nobody knows that they are - that
is what unsurveyed means - so the true loss is somewhere between zero and the
bar drawn. The bound is still worth having, because it is the number that would
have to be small for the truncation not to matter.

The same conflation is what drowned three roadways at t = 0 in an earlier
product: roadway_manager.bulldoze drowns a road when more than 20% of the cells
flanking it sit at or below 0 m MHW, and an unsurveyed cell passes that test.
predict_drowning() is run here against this product's own setbacks and the
verdict is printed.

THE PANELS
a  The CASCADE input domain for all 90 domains: the 2 dune rows Barrier3D
   builds from dunes/domain_<N>_dune.npy, then the interior rows from
   topography/. Ocean at the bottom. This is the whole stack the model starts
   from, in the order it starts from it.

b  The seaward 300 m of the same stack, so the dune rows and the first interior
   rows are actually resolvable. At island scale two 10 m rows are one pixel.

c  Island width Barrier3D sees, and the upper bound on what unsurveyed cells
   cost it.

d  Unsurveyed cells per domain: in the DEM, and still there in the input domain.

DUNE ROWS
dunes/domain_<N>_dune.npy is one height above berm per profile, (50,). Barrier3D
runs DuneWidth = 2, and row 1 is a copy of row 0 - see the dune-rows note in
the extractor. Both rows are drawn. Their elevation is BERM_ELEV + height; the
berm is 1.7 m NAVD88, so 1.34 m MHW.

No dune cell in this product is unsurveyed: all 4500 carry a measured height.
That is checked at run time, not assumed, and the count is printed.

INPUT   <product>/npy-arrays/domain_<N>.npy                    m NAVD88
        <product>/dune-topo/<version>/topography/domain_<N>_topography.npy   dam
        <product>/dune-topo/<version>/topography/domain_<N>_nodata.npy      bool
        <product>/dune-topo/<version>/dunes/domain_<N>_dune.npy             dam

        Product and version resolve through scripts/site_layer/hat_topo_version.py.
        Do not hardcode either.

OUTPUT  <product>/dune-topo/<version>/figures/HAT_topo_retained_<version>.png
```

Notes that were in the code:

```text
Every output of this folder lands under one directory beside the extraction it
describes, rather than being scattered through the run folder it did not
produce. audit_dir() is the only place that name is spelled.
```

<details><summary>Function notes (the original docstrings)</summary>

**`classify_input()`**

```text
Raw ocean-first m NAVD88 -> (unsurveyed, sub-MHW) counts in the envelope.

The envelope rule is the one HAT_plot_dem_holes.py uses, so panel (d)'s
'in the DEM' bars mean the same thing that figure's red does.
```

**`stack_domain()`**

```text
The CASCADE input domain as one categorical array, ocean-first.

Rows 0..DUNE_ROWS-1 are the dune, then the interior. Returns the codes and
the per-profile last-land row of the INTERIOR part, in interior numbering.
```

**`widths_and_bound()`**

```text
(width Barrier3D sees, width if every unsurveyed cell were land), cells.

The first is roadway.interior_widths verbatim - barrier3d.FindWidths, which
stops at the first cell at or below sea level. The second re-runs it on a
copy with the unsurveyed cells lifted above the threshold, which is the
most land those cells could possibly be hiding.
```

**`read_setbacks()`**

```text
1984 road setback in metres per GIS domain, or None if unavailable.

Only used for the t = 0 drowning check, which is a printout. A missing file
downgrades that check rather than failing the figure.
```

</details>

### HAT_test_hole_pond_or_dropout.py

Are the unsurveyed holes that truncate the model island real ponds, or lidar dropouts?

From the script's original header:

```text
Are the unsurveyed holes that truncate the model island real ponds, or lidar
dropouts?

THE QUESTION AND WHY IT IS NOT OBVIOUS
99 profiles in 1984-start have their island truncated by an unsurveyed cell:
barrier3d.FindWidths stops at the first cell at or below sea level, and an
unsurveyed cell is written to the -3.0 m water sentinel, so the scan stops there
and every measured cell behind it is discarded. 731 cells are involved.

Whether that is wrong depends on what those cells are:

    a POND      the survey was right to return nothing, water is water, and the
                truncation is the island's real shape. Change nothing.
    a DROPOUT   the survey failed over dry ground, and the model is running on
                an island several hundred metres narrower than the data.

An earlier diagnostic claimed the holes are dropouts because they are bracketed
by measured land at ~0.9 m MHW on both sides. That argument is wrong and is
recorded here so it is not made again: a pond is BY DEFINITION surrounded by dry
ground, so "dry on both sides" is equally consistent with either. The elevation
test discriminates a survey stopping at a false shoreline from a real bay
margin, which is a different question about a different set of cells.

Evidence pointing the other way, which any result has to beat: a cell is
unsurveyed in this product only if 1996 ALACE, 2009 USACE AND 2014 NOAA
Post-Sandy all failed at it. Three independent surveys over 18 years failing at
one spot is what persistent water looks like. "Mostly ponds, change nothing" is
a live outcome, not a failed test.

TWO INDEPENDENT REFERENCES, BOTH MUST AGREE TO OVERTURN
A  2014 NCFMP hydro-flattening stamp.
   NCFMP is DISQUALIFIED as an elevation source - 94.6% of its coverage in the
   gap is two stamped constants, -0.762 and -0.914 m, i.e. -2.5 and -3.0 ft.
   That is exactly what makes it authoritative here. A hydro-flattening
   compiler delineates water-body polygons and stamps a flat surface inside
   them, so the constant IS a water classification, made independently of
   whether the lidar returned anything. This reads its water mask, never its
   elevation, so its "unknown" vertical datum and foot-derived values do not
   matter: the test is whether the value is CONSTANT, not what it means.

C  Blob shape.
   Marsh ponds are compact and roughly convex. Lidar dropouts follow flight
   lines and scan geometry, so they come out elongated and sparse in their
   bounding box. Computed on the raw nodata mask, which no other reference
   touches, and independent of time - which matters because A is 2014 standing
   in for 1996.

A third reference, the 1996 aerial frames, was considered and dropped: they are
scanned historical imagery with per-frame exposure variation, so classifying
them needs either a manual pass over 99 chips or a threshold that would not
survive review. Dropping it leaves the temporal assumption resting on A alone.
That is a real weakness of this design and belongs in the methods paragraph.

THE DECISION RULE IS DELIBERATELY ASYMMETRIC
    A and C agree             ->  their verdict
    A and C conflict, B votes ->  B decides
    otherwise                 ->  WATER, and the DEM is left alone

B is the 1996 aerial review from HAT_hole_aerial_picker.py. Chips were rendered
only where A and C conflict, so B exists exactly where the automated pair
cancels out - which makes "2 of 3 agree" and "B breaks the tie" the same rule.
B is also the only CONTEMPORANEOUS reference, so where it has an opinion it
deserves to carry the decision rather than be outvoted by two proxies.

Unknown, absent or conflicting evidence defaults to water, which is the
current behaviour. The DEM changes only where there is affirmative evidence it
is wrong. The cost of that choice is that the island stays too narrow wherever
the test cannot resolve a hole - a known, stated, one-directional bias.

INPUT   dune-topo/<version>/bracketed_hole_cells.csv   from HAT_plot_island_nodata
        <product>/npy-arrays/domain_<N>.npy            raw nodata, for shape
        D:\Hatteras_GIS\Elevation\Polygons\2014\2014_NCFMP_*\*.tif

OUTPUT  dune-topo/<version>/hole_verdicts.csv          per hole, both votes
        dune-topo/<version>/dropout_mask/domain_<N>.npy  bool, cleared cells only
```

Notes that were in the code:

```text
Every output of this folder lands under one directory beside the extraction it
describes, rather than being scattered through the run folder it did not
produce. audit_dir() is the only place that name is spelled.
```

```text
THE RULE, with reference B folded in.

Chips were rendered only where A and C conflict, so B exists exactly
where the automated pair cancels out. That makes "2 of 3 agree" and
"B breaks the tie" the same rule, not two:

A and C agree            -> their verdict
A and C conflict, B votes -> B decides (B + one of A/C = 2 of 3)
otherwise                 -> POND, the conservative default

UNCLEAR is not a vote for water. It is the absence of a vote, and it
falls through to the default for the same reason anything unresolved
does: the DEM changes only on affirmative evidence that it is wrong.
```

<details><summary>Function notes (the original docstrings)</summary>

**`ring_points()`**

```text
Measured-ground cell centres in a ring around the hole footprint.

Excludes the hole itself and every other unsurveyed cell, so the ring is
ground some survey actually saw. Without that the ring could be more
no-data, and the depression test would compare two stamps.
```

**`verdict_a()`**

```text
NCFMP stamp test -> (verdict, modal fraction, depression).

The statistic is the MODAL FRACTION - what share of the NCFMP pixels under
the hole carry a single repeated value. That is deliberately the same
metric HAT_survey_dem_coverage.py uses to disqualify a DEM as void-filled
("genuine surveys show their most common value in ~0.2% of cells"), applied
here at one hole instead of island-wide.

It replaced a max-minus-min range test, which failed in both directions and
is worth recording:
  * one NCFMP pixel per 10 m cell made a single-cell hole trivially flat,
    range exactly 0, and returned UNKNOWN for 45 of 99 holes;
  * the whole footprint fixed that but picked up the EDGE of a stamped
    water polygon, where values vary by construction, so a pond-edge cell
    read as varying and therefore as a dropout.
A modal fraction survives both: a pond cell is mostly one stamped value
even when its footprint clips the polygon edge, and dry ground is not.
```

**`read_aerial()`**

```text
Reference B: {(domain, profile): POND|DROPOUT|UNCLEAR} from the review.

Absent file, or an unreviewed row, simply yields no vote - which the
decision rule below treats as no evidence, not as evidence of water.
```

**`blob_shape()`**

```text
Shape of the connected nodata blob this hole belongs to.

8-connectivity, per domain. A blob straddling a domain seam is measured
only within its own 500 m tile; that under-measures a few blobs and is
noted rather than corrected, because the domains are the unit everything
else here is expressed in.
```

**`sample()`**

```text
Every NCFMP pixel inside the 10 m cell footprint around each point.

Point-sampling one NCFMP pixel per 10 m cell was the first version and
it broke the flatness test: NCFMP is 1.52 m, so a 10 m cell contains
~43 of its pixels, and taking one made a single-cell hole trivially
"flat" with a range of exactly zero. 45 of 99 holes returned UNKNOWN
for that reason alone. Reading the footprint measures flatness WITHIN
one 10 m cell, which is what the stamp test actually needs.
```

</details>
