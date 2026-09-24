# 1996 island offsets from the SHORELINE

The 1996 start's island offset, measured from the **CoastSat satellite
shoreline** instead of the digitised dune line. Added 2026-09-22 (Hannah, by
interview).

## Why this is a source and not a version

`../v1/` is built from `duneline_1997.geojson` — a line traced on aerial
imagery flown on one day. This folder is built from the mean satellite
shoreline over calendar 1995–1997. Those are **two different features on the
island**, not two readings of the same one, so this is a separate source with
its own `v1`, never a `v2` of the dune build
([[feedback-version-numbering-restarts]]).

The dune build keeps the flat `1996/v1/` layout it has always had, so nothing
the runner resolves moved. `hat_topo_version.offset_version(1996)` still
returns the dune build; the shoreline arm is reached explicitly with
`source="shoreline"`.

## What the model reads

**Nothing, yet.** This is an alternative arm, built to be compared against the
dune profile. `hatteras_site_config._island_offset_file(1996)` still resolves
`1996/v1/Island_Dune_Offsets_1996_PADDED_120.csv`.

**Before this could ever become the 1996 offset**, one thing has to be checked
that has not been: the Barrier3D interior topography for 1996 is extracted in
the *dune-line* frame. Swapping the offset alone, without re-examining that,
risks counting the beach width twice.

## Layout

```
CURRENT           the build every reader of this source takes -> v1
v1/               the build from the 1995-1997 mean shoreline, padded with the model's wrap-around
superseded_20260924_pre-metres/v1/  the same build with the old slope-and-bridge buffers
```

Resolved by `hat_topo_version.offset_version(1996, "shoreline")` and
`offset_file(1996, source="shoreline")`. Overridden for one run by
`HAT_OFFSET_VERSION_1996_SHORELINE` in the environment — a separate key from
the dune build's `HAT_OFFSET_VERSION_1996`, so overriding this arm cannot
silently move what the runner reads.

## Renumbered (2026-09-24)

`v1/` is the build now in `superseded_20260924_pre-metres/v1/` re-padded with the model's smooth wrap-around (`cascade_pipeline.hindcast.pad_offset_ring`); the real domains are identical and the padded file is exactly what offset_mode `metres` hands Cascade. `CURRENT` = `v1`. Built as `v2` and renumbered the same day (numbering restarted). See `v1/PROVENANCE.md`.
