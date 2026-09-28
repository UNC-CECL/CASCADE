# island-offset

How the island's planform (the BRIE island offset) is set: its units, and whether it comes from the dune line or the CoastSat shoreline.

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| [`2026-09-22-shoreline-offset-at-div10-scale`](2026-09-22-shoreline-offset-at-div10-scale/NOTE.md) | Does a shoreline-derived offset change 1996–2010? | Could not say: the runs used the ÷10 offset, so the change was a tenth of the real difference. | superseded by `2026-09-25-offset-source-duneline-vs-shoreline` |
| [`2026-09-24-metres-1-offset-units`](2026-09-24-metres-1-offset-units/README.md) | Offset in metres (as BRIE reads it) or ÷10? Can waves re-tune the metres runs? | Metres is correct; the ÷10 was a units bug. Metres became the default on 09-24. | **current** (the units decision) |
| [`2026-09-25-offset-source-duneline-vs-shoreline`](2026-09-25-offset-source-duneline-vs-shoreline/README.md) | Dune line or shoreline as the offset, in metres (Hs 1 / Tp 8 / asym 0.8, high-angle 0.3–0.55)? | Managed: no difference (+23–24% either way). Natural: the shoreline offset scores higher (+22% vs +14% at 0.45). | **current**; not re-run at the recommended waves |

**Status** — **current**: its answer is in use now. **superseded**: a later study
re-asked it; follow the pointer. **record**: a finished check or a result from an
earlier set-up (÷10 offset, Hs 2.5 calibration), kept so the number can be traced.

Back to [the map](../README.md).
