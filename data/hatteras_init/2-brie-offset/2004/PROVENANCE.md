# 2004 island offsets — which SOURCE, then which version

This folder holds **no build of its own**. Every build sits under the feature
it was measured from, so a folder listing says what it is:

| source | what it was measured from | CURRENT | also here |
|---|---|---|---|
| `duneline/` | a dune line digitised from aerial imagery | `v1` | `superseded_20260915_flat/` |

`shoreline/` — not built for 2004.



## Resolving a build

```python
from site_layer.hat_topo_version import offset_file
offset_file(2004)                        # the dune build the model reads
offset_file(2004, source="shoreline")    # the shoreline arm
```

`CURRENT` inside each source folder names the build that source's readers
take; `HAT_OFFSET_VERSION_2004` in the environment outranks it for one run
(a non-default source uses `HAT_OFFSET_VERSION_2004_<SOURCE>`, so overriding
an arm cannot silently move what the runner reads). **Never join these paths
by hand** — `hat_topo_version` owns the layout, and it changed on 2026-09-22.

## What changed on 2026-09-22

Dune builds used to sit flat at `2004/v<n>/`, from when they were the only
kind. A listing could not then say what `v1` was measured from, and the
shoreline source added in the same week sat one level deeper than the dune
one. Everything now nests under its source. Run metadata records
`island_offset_version` as `"duneline/v1"` rather than `"v1"`; a token with no
`/` is from before the split, and every build then was dune-derived.
