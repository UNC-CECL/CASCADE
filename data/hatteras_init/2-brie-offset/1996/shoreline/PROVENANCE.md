# 1996 island offsets from the SHORELINE

The 1996 start's island offset, measured from the **CoastSat satellite
shoreline** instead of the digitised dune line. Added 2026-09-22 (Hannah, by
interview).

## Why this is a source and not a version

`../v1/` is built from `duneline_1997.geojson` — a line traced on aerial
imagery flown on one day. This folder is built from the mean satellite
shoreline over a window (`v2`, CURRENT: 1995-10-12 to 1997-10-12, ±1 yr of
the 1996 ALACE lidar; `v1`: calendar 1995–1997). Those are **two different features on the
island**, not two readings of the same one, so this is a separate source with
its own `v1`, never a `v2` of the dune build
([[feedback-version-numbering-restarts]]).

The dune build sits beside it at `../duneline/` (flat at `1996/v1/` until
2026-09-22). `hat_topo_version.offset_version(1996)` with no `source` returns
the dune build; this one is `source="shoreline"`.

## What the model reads

**This build, since 2026-10-05.** Runs read the shoreline source by default
(`hat_topo_version.RUN_OFFSET_SOURCE`), so
`hatteras_site_config._island_offset_file(1996)` resolves
`1996/shoreline/v2/Island_Shoreline_Offsets_1996_PADDED_120.csv`. Run metadata
records `island_offset_version: "shoreline/v2"`. Matrix runs before 2026-10-05
read the dune build; `HAT_ISLAND_OFFSET_SOURCE=duneline` reproduces them.

**Open, and not recorded as checked:** the Barrier3D interior topography is extracted in the *dune-line* frame, so swapping the offset alone, without re-examining that, risks counting the beach width twice. This was the stated condition for the shoreline becoming a default; the default moved on 2026-10-05 and nothing in the repository records the check.

## Layout

```
CURRENT           the build every reader of this source takes -> v2 (since 2026-09-29; was v1)
v1/               the build from the 1995-1997 mean shoreline, padded with the model's wrap-around
v2/               the build from the 1995-10-12 to 1997-10-12 mean shoreline, ±1 yr of the 1996 ALACE flights (2026-09-29; CURRENT)
superseded_20260924_pre-metres/v1/  the same build with the old slope-and-bridge buffers
v<n>/comparisons/duneline_vs_shoreline/  the dune line against that build (compare_offset_sources.py --shoreline-version v<n>)
```

Resolved by `hat_topo_version.offset_version(1996, "shoreline")` and
`offset_file(1996, source="shoreline")`. Overridden for one run by
`HAT_OFFSET_VERSION_1996_SHORELINE` in the environment — a separate key from
the dune build's `HAT_OFFSET_VERSION_1996`, so overriding this arm cannot
silently move what the runner reads.

## Renumbered (2026-09-24)

`v1/` is the build now in `superseded_20260924_pre-metres/v1/` re-padded with the model's smooth wrap-around (`cascade_pipeline.hindcast.pad_offset_ring`); the real domains are identical and the padded file is exactly what offset_mode `metres` hands Cascade. `CURRENT` = `v1`. Built as `v2` and renumbered the same day (numbering restarted). See `v1/PROVENANCE.md`.
