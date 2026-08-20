# Old road-setback method — how it works and where it strains

Generated 2026-08-20T11:53:20 by `HAT_old_method_figures.py`. Read-only: no forcing file is written here.

## The calculation

Per domain, from two ArcGIS point exports that share one transect system:

```
setback = min(road ORIG_LEN)  −  min(dune ORIG_LEN)
```

`ORIG_LEN` is the **index of a point along its transect**, not a measured length — always an integer, always with `ORIG_SEQ == ORIG_LEN + 1`. A difference of indices is a count of 1 m points, so the result is metres and the output label is correct. Two independent confirmations:

- every crossing is **3–4 consecutive** indices in both layers;
- both layers were buffered at a 1.5 m half-width — road `BUFF_DIST` 4.92125 (US survey feet = 1.5 m), dune `BUFF_DIST` 1.5 (m). A 3 m buffer sampled at 1 m spacing yields exactly 3–4 points.

- **1984**: recomputing from the raw exports reproduces exactly (`RoadSetback_1984.csv`).
- **2004**: recomputing from the raw exports reproduces exactly (`RoadSetback_2004.csv`).

> This matters for the record: the legacy setbacks **are reproducible** from the CSVs in the repo. What is missing is the upstream `duneline_<year>.geojson`, not the offsets derived from it.

## The shortcut — two minima, taken independently

Both minima are taken over the whole domain, each within its own layer, so they need not fall on the same transect. The method can therefore measure a road at one alongshore position against a dune at another.

This is avoidable. `LineID` is a **shared transect key** — 410 of 410 road transects are also dune transects across GIS 9–90 — so a per-transect difference, then a median over transects, is available and is the like-for-like comparison.

| Year | Median legacy − paired | p10 | p90 | Max abs | Domains > 25 m |
|---:|---:|---:|---:|---:|---:|
| 1984 | +1.0 m | -17.7 | +24.9 | 55 | 15 |
| 2004 | +1.0 m | -15.8 | +23.0 | 54 | 14 |

Small in the median, not small in the tail. The median is ~1 m — which is why this was recorded as a minor shortcut — but the disagreement exceeds 25 m in roughly a sixth of the road and reaches ~55 m. Those domains are flagged `DISAGREE` in the per-domain CSV.

Sample size compounds it. The transect layer is regular — **every** domain in GIS 9–90 carries exactly 5 shared transects, one per 100 m of a 500 m domain — so both the legacy minimum and the paired median rest on 5 observations. The spread those 5 span is a median of 52 m for the road and 66 m for the dune, reaching 172 m and 131 m. Figure 2(a).

## Road source vintages — checked, and deliberate

The road inputs do not carry the vintages their filenames suggest:

| Year label | Road source layer | Dune source layer |
|---|---|---|
| 1984 | `NC12_1978_buffer` | `dune_1984_buffer` |
| 2004 | `NC12_2008_buffer` | `dune_2004_buffer` |

The dune layers match their labels; the road layers do not. Git records `nc12_1984.csv` as a rename of `1978_road_offset_raw.csv`, and the 2004 export carries `NC12_2008_buffer`.

**This is deliberate and it is sound.** NC-12 did not move between 1978 and 1984, or between 2004 and 2008, so the 1978 line stands in for 1984 and the 2008 line for 2004. The digitised lines that exist are the ones that were flown; substituting an unchanged neighbour beats not measuring.

That is checkable rather than assumed, and the script re-checks it every run against `HATTERAS_ROAD_EVENTS` — the same list the runner uses:

| Substitution | Road events in between | Verdict |
|---|---|---|
| 1978 → 1984 | none | **valid** |
| 2008 → 2004 | none | **valid** |

The three events on record — 1989 Pea Island relocation, 1999 inter-village relocation, 2022 Jug Handle Bridge — all fall outside both intervals.

Consequence: because the road is unchanged across both substitutions, **panel (c) of figure 1 is a genuine 20-year change** in the forcing, not a vintage artifact. A new method fed the same road CSVs inherits the same, valid, substitution.

## A note on precision

Figure 1(a) is the point: the setback is a **difference of two large numbers**. Both distances run to several thousand metres from the transect origin, and the answer is a few hundred. Nothing here is wrong — the indices are exact integers, so there is no floating-point erosion — but it does mean the result inherits any error in either baseline in full, and a transect origin that shifted between the two layers would propagate straight through.

## Per-domain results

Full table in `RoadSetback_oldmethod_domains.csv`. Flagged domains only:

| Year | GIS | Legacy | Paired | Difference | Shared transects | Flags |
|---:|---:|---:|---:|---:|---:|---|
| 1984 | 17 | 141 | 113 | +28 | 5 | `DISAGREE(+28)` |
| 1984 | 21 | 163 | 136 | +27 | 5 | `DISAGREE(+27)` |
| 1984 | 33 | 307 | 273 | +34 | 5 | `DISAGREE(+34)` |
| 1984 | 67 | 454 | 402 | +52 | 5 | `DISAGREE(+52)` |
| 1984 | 70 | 448 | 495 | -47 | 5 | `DISAGREE(-47)` |
| 1984 | 71 | 483 | 454 | +29 | 5 | `DISAGREE(+29)` |
| 1984 | 72 | 582 | 548 | +34 | 5 | `DISAGREE(+34)` |
| 1984 | 74 | 616 | 648 | -32 | 5 | `DISAGREE(-32)` |
| 1984 | 77 | 575 | 547 | +28 | 5 | `DISAGREE(+28)` |
| 1984 | 79 | 580 | 617 | -37 | 5 | `DISAGREE(-37)` |
| 1984 | 81 | 409 | 461 | -52 | 5 | `DISAGREE(-52)` |
| 1984 | 83 | 176 | 121 | +55 | 5 | `DISAGREE(+55)` |
| 1984 | 85 | 40 | 71 | -31 | 5 | `DISAGREE(-31)` |
| 1984 | 87 | 55 | 99 | -44 | 5 | `DISAGREE(-44)` |
| 1984 | 88 | 158 | 194 | -36 | 5 | `DISAGREE(-36)` |
| 2004 | 21 | 168 | 140 | +28 | 5 | `DISAGREE(+28)` |
| 2004 | 33 | 367 | 339 | +28 | 5 | `DISAGREE(+28)` |
| 2004 | 37 | 346 | 373 | -27 | 5 | `DISAGREE(-27)` |
| 2004 | 38 | 269 | 296 | -27 | 5 | `DISAGREE(-27)` |
| 2004 | 55 | 177 | 123 | +54 | 5 | `DISAGREE(+54)` |
| 2004 | 67 | 453 | 410 | +43 | 5 | `DISAGREE(+43)` |
| 2004 | 69 | 555 | 582 | -27 | 5 | `DISAGREE(-27)` |
| 2004 | 70 | 437 | 486 | -49 | 5 | `DISAGREE(-49)` |
| 2004 | 71 | 473 | 438 | +35 | 5 | `DISAGREE(+35)` |
| 2004 | 79 | 526 | 552 | -26 | 5 | `DISAGREE(-26)` |
| 2004 | 81 | 354 | 382 | -28 | 5 | `DISAGREE(-28)` |
| 2004 | 82 | 196 | 165 | +31 | 5 | `DISAGREE(+31)` |
| 2004 | 83 | 115 | 89 | +26 | 5 | `DISAGREE(+26)` |
| 2004 | 88 | 121 | 159 | -38 | 5 | `DISAGREE(-38)` |
