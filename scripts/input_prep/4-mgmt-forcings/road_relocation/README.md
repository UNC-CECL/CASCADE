# road_relocation

How far NC-12 moved between two digitised vintages, per Barrier3D domain.

Port of `roya_files/road_relocation_dis.py` (Pea Island, 1992→1996) onto
Hatteras. Same measurement — sample the old road inside each domain, measure
each sample to the whole new road — with the Hatteras road files, the
90-polygon domain file, and a signed direction added.

```
HAT_road_relocation_distance.py
```

```
python HAT_road_relocation_distance.py
HAT_RELOC_FROM=1984 HAT_RELOC_TO=2004 python HAT_road_relocation_distance.py
```

Writes to `data/hatteras_init/4-mgmt-forcing/road_relocation/<from>_<to>/`:
a per-domain CSV, the sample points as GeoJSON, and three figures (below).

**This is not a forcing.** `road_setback` comes from `../road_offset/`,
measured against interior row 0 of the model grid. This is a GIS-frame
observation of two digitised lines, and nothing reads it. Its use is to say
whether an observed relocation happened and where, so a *modelled* relocation
can be checked against one.

## The zeros are not measurements

**70.5% of the 1984 line's vertices are identical to a 2004 vertex to the
millimetre.** The 2004 road was digitised by editing a copy of the 1984 one,
and only the stretches that visibly moved were re-drawn. Every shared vertex
measures exactly 0.000 m by construction.

So a 0 m domain here means *nobody edited that stretch* — not that the road
held still. The script measures the shared fraction itself and prints it
before any result.

A second artefact sits behind the first: some stretches *were* re-drawn, but
only **re-traced** — the new centreline wanders a metre or two and never
leaves the old road's own footprint. That is one road digitised twice, not a
road that moved.

Every domain therefore carries a `classification`:

| | |
|---|---|
| `no_edit` | `coincident_fraction ≥ 0.90` — the line was copied through unedited |
| `redigitized` | largest displacement anywhere `< REDIGITIZE_MAX_M` (5 m) |
| `relocated` | the road actually moved — **the only measured subset** |

Nothing is dropped from the CSV; only `relocated` domains are summarised and
drawn. **Never average across the whole table** — it averages real movement
against copy and re-tracing alike, which is why no island-wide mean is printed.

### Why 5 m

About half a road width: NC-12 is ~8 m wide, and two digitisings of one road
off different photos disagree by that much from georeferencing alone. The
Hatteras data splits cleanly there and the exact value does not matter —
re-traced domains top out at **3.25 m**, real relocations start at **12.2 m**,
so anything from ~3.5 to ~12 m gives the same answer.

Direction corroborates it independently. The re-traced domains are internally
coherent (`sign_agreement = 1.00`) but flip sign arbitrarily between
neighbours — 16–18 seaward, 19–22 landward, 24–27 seaward. That is what a
smoothed line does. A road that moved does not change its mind every few
hundred metres.

## What the columns mean

| Column | |
|---|---|
| `mean_relocation_m` | unsigned nearest-distance, old road → new road. Magnitude only, and a *lower bound* on cross-shore movement: where the lines cross obliquely the nearest point is diagonally alongshore |
| `mean_signed_landward_m` | same displacements, signed **+ landward / − seaward**, from which side of the northward-oriented old road the new road falls on (`OCEAN_ON_RIGHT`) |
| `coincident_fraction` | share of samples on unedited shared geometry |
| `sign_agreement` | share of *moved* samples agreeing on direction; `ROADS_CROSS` below 0.90 means the domain mean is averaging both ways |
| `oblique_fraction` | share of samples where the road runs >60° off north; over half and the domain is flagged `OBLIQUE_SIGN` |
| `classification` | `no_edit` / `redigitized` / `relocated` — see above |

The sign has one blind spot. Around the Cape Point bend the island turns
east–west, the ocean stops being on the right of a northward road, and the
landward/seaward convention breaks. Those domains are flagged `OBLIQUE_SIGN` —
on the 1984→2004 pair that is **domain 8 and nothing else**. Magnitude is
unaffected; it never used the tangent.

## The three figures

They answer three questions at three incompatible scales, so they are three
files rather than one page.

| File | Question |
|---|---|
| `..._alongshore.png` | **Where** along the island did the road move? |
| `..._sites.png` | **What** did the move look like? |
| `..._domain_map.png` | **Which** domains carry a relocation? |

The island is 8 km wide and 45 km long — aspect 0.18 — so any true-scale map of
the whole thing is a hairline in a column of white space. `sites` sidesteps
that by only ever drawing ~2 km; `domain_map` rotates the island 90° clockwise
so it runs south (left) → north (right), matching the domain axis of
`alongshore`. Rotation preserves distance, so its scale bar is still honest.

**alongshore** — signed bars per domain, with a pale bar behind for the largest
displacement anywhere in that domain (the gap between them is how much of the
domain moved). Artefact domains shaded, town spans along the foot so a domain
number means a place, each site bracketed and named.

**sites** — one true-scale zoom per relocation site. Sites are found from the
data (contiguous runs of `relocated`), not named in the script, so this still
works on a different pair of vintages. **All panels share one scale**: Buxton
spans 8 domains and Rodanthe 4, and framing each on its own extent would render
them at different metres-per-inch, making the shorter site look like the bigger
relocation. The old road is a dashed spine drawn **on top of** its own sample
points — underneath, the points hide it and the two alignments no longer
visibly pull apart — plus an OCEAN label, because the signed metric depends
entirely on which side that is.

**domain_map** — every domain in its real place, filled by what it carries:
white where the road never reaches, light grey never-edited, darker grey
re-traced, and the relocated ones outlined in black and filled by distance on
the same colour scale as `sites`.

## Result, 1984 → 2004

Of 83 road-carrying domains: **56 never edited, 15 re-traced, 12 relocated.**

The 12 are two coherent stretches, **every one of them landward**, with
nothing at all in between:

| Domains | Where | Largest domain mean |
|---|---|---|
| 8–15 | Buxton | **77 m landward** (domain 11) |
| 84–87 | Rodanthe | **109 m landward** (domain 85) |

Across those 12: mean 45.1 m, median 40.1 m, 0 seaward. Domains 8 and 15 are
flagged `ROADS_CROSS` — they are the tapered ends of the Buxton relocation,
where the new alignment rejoins the old and the domain mean mixes both
directions. Domain 8 also carries `OBLIQUE_SIGN` (Cape Point bend), so read its
magnitude and ignore its direction.

Both stretches sit where the hindcast already expects trouble; Rodanthe is the
S-curves. Treat the magnitudes as the movement of a *digitised centreline*
between a 1978 and a 2008 photo, which is the interval those two files really
span — the 1984/2004 labels are the hindcast periods they stand in for.
