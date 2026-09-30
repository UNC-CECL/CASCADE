# 2026-09-28 (late): 2010 ends re-solved after the dune-cap fix

Hannah, 2026-09-28: "go with option 1, clip only the bulldozed sand". CASCADE's beach/dune manager used to clip whole dune cells to 4 m above the berm every year; now its 4 m cap limits only the overwash sand it adds (`cascade/beach_dune_manager.py`, `DUNE_CAP_APPLIES_TO`). The evidence is in `output/comparisons/adoption_2026-09-28/README.md`. The change moves the managed runs, so the ends were checked with the edgeBE full-management runs at the adopted ends:

| window | ends | residuals GIS 1 / GIS 90 |
|---|---|---|
| 1996-2010 | +4.3509 / +19.0935 | −0.006 / +0.007: holds, not re-solved |
| 2010-2024 | +8.0 / +22.4937 | +0.002 / **+0.185**: re-solved |

GIS 89 was one of the clipped domains in 2010, which is why GIS 90 moved.

The protocol is as `2026-09-28-ends-resolved-adopted/`: option A waves, full management, edgeBE, no groin, CoastSat LRR target (GIS 1 raw, GIS 90 LOWESS-7). The driver was:

    HAT_resolve_ends_metres.py --periods 2010 --hs 2.0 --tp 7.5 --asym 0.6 --ahf 0.5
        --accept 0.03 --tag end-domain-boundaries/2026-09-28-ends-resolved-dunecap
        --seed "2010=8.0,21.57"

| step | ends | residuals |
|---|---|---|
| 1 | +8.000 / +21.570 | +0.002 / +0.088 |
| 2 | +8.000 / +19.413 | +0.002 / −0.316 |
| 3 | +8.000 / +21.099 | +0.002 / −0.045 |
| 4 | **+8.000 / +21.258** | **+0.002 / −0.008** |

**Adopted in `HATTERAS_BE_EDGE_ONLY`: 2010 (+8.0, +21.2582).** GIS 1 did not need to move. `tables/ends.json` is the record, and `runs/` holds the probes.
