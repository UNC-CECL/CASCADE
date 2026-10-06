# 2010 island offsets — which SOURCE, then which version

> **Renamed 2026-10-05: this folder was `2010/`.** The period now starts in the year of the DEM it
> sits on (USACE 2009, flown 2009-08-10 to 08-24). The contents are unchanged. In the current builds
> only the file names went from `2010` to `2009`; superseded builds keep their old names.
> Text below this note still says 2010.

This folder holds **no build of its own**. Every build sits under the feature
it was measured from, so a folder listing says what it is:

| source | what it was measured from | CURRENT | also here |
|---|---|---|---|
| `duneline/` | a dune line digitised from aerial imagery | `v1` | `superseded_20260919_pre-redigitized/` |
| `shoreline/` | the mean CoastSat shoreline (v2: 2008-08-17 to 2010-08-17, ±1 yr of the 2009 USACE lidar; v1: calendar 2009-2011) | `v2` (since 2026-09-29) | |



## Resolving a build

```python
from site_layer.hat_topo_version import offset_file
offset_file(2009, source="shoreline")    # the build runs read since 2026-10-05 (shoreline/v2)
offset_file(2009)                        # the dune build; HAT_ISLAND_OFFSET_SOURCE=duneline
```

`CURRENT` inside each source folder names the build that source's readers
take; `HAT_OFFSET_VERSION_2010` in the environment outranks it for one run
(a non-default source uses `HAT_OFFSET_VERSION_2010_<SOURCE>`, so overriding
an arm cannot silently move what the runner reads). **Never join these paths
by hand** — `hat_topo_version` owns the layout, and it changed on 2026-09-22.

## What changed on 2026-09-22

Dune builds used to sit flat at `2010/v<n>/`, from when they were the only
kind. A listing could not then say what `v1` was measured from, and the
shoreline source added in the same week sat one level deeper than the dune
one. Everything now nests under its source. Run metadata records
`island_offset_version` as `"duneline/v1"` rather than `"v1"`; a token with no
`/` is from before the split, and every build then was dune-derived.
