# 1996-start

Barrier3D domains for the **1996 run** (calibration, 1996-2009). Named
`1984-start/` until 2026-10-19, when the 1984 start was retired; nothing in
the arrays changed.

| | |
|---|---|
| DEM | `0-elevation/2009-2014-1996/` - the 2009-2014 baseline plus the 1996 ALACE survey grafted wherever ALACE has data, so the beach and dunes are 1996 |
| npy arrays | 90 domains, exported 2026-08-26 from the no-boundary DEM |
| **what loads** | **`v2`** - `CURRENT` says so |
| picks | `1-extraction/picks/HAT_dune_search_windows_v2.json` (not regenerable; `_v1.json` kept for the offloaded v1) |

## Layout

```
1-extraction/
    npy-arrays/          domain_<N>.npy   m NAVD88, -10 nodata   <- extractor INPUT
    npy-arrays_survey/   domain_<N>.npy   provenance codes (0/1996/2009/2014)
    picks/               the dune search windows that define each version
    aerial-review/       58 holes of 1996 imagery review, keyed on (domain, profile) - LIVE
dune-topo/
    CURRENT              v2
    v2/                  extraction on the 2026-09-02 re-pick, NC-12 visible;
                         what the road tree measures against
```

`v2` is an unmodified extraction: no rows added, no cell edited.

## Offloaded 2026-10-19

To `D:\archive\CASCADE_offload\data\hatteras_init\1-barrier3d-domains\1984-start\`,
SHA-256 verified, under the old folder name:

- `dune-topo/v1/` - extraction on the original pick set (2026-08-27); the
  2026-09-01 calibration tree and the `pea1989base` arm ran on it
- `dune-topo/v3/` - `v2` + the 1984 footprint behind the road
- `2-domain-reconstruction-1984/` - the 1984 footprint work in six steps,
  including `duneline_retreat_1984_1997.csv`, the measurement any insert layer
  is built from

`dune-topo/README.md` describes all three versions as they stood before the
offload. Earlier clears (2026-08-27, 09-03, 09-07) are in `../LINEAGE.md`.
