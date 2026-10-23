# Is the model's fillet the right shape? (2026-10-14)

**Question.** In the 1996-2009 matrix runs the groin helps only at GIS 5-6, and GIS 7-14 erode in the record but not in the model. Does the groin's effect not reach far enough north? Checked on the one period where only the groin acted: its build-up, 1967-1978.

**Setup.** The wet/dry photo lines of 1974 and 1978 against the 1967-start "full life" runs of `../2026-10-13-pregroin-planform/` (no groin, 0.16, and the whole-life 0.25). Change since 1967 at GIS 2-12, seaward positive, with each line's straight line between its end domains taken off so only the shape the groins made is compared.

    python fillet_shape.py   # -> fillet_shape.csv, figures/fillet_shape_1967_1978.png

## Findings

Change since 1967 (m, seaward +, end-to-end line removed):

| | GIS 3 | 4 | 5 | groin | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| 1974 observed | +6 | -62 | -60 | | +7 | +3 | -4 |
| 1974 model 0.25 | -4 | -44 | -68 | | +38 | -2 | -5 |
| 1978 observed | -43 | -129 | -114 | | +7 | +28 | +5 |
| 1978 model 0.25 | -12 | -59 | -78 | | +32 | 0 | -5 |
| 1978 model, no groin | -5 | -38 | -30 | | -26 | -13 | -9 |

1. **By 1974 the model has the shape about right downdrift** (GIS 4-5, -44 to -68 against -60 to -62 m). It over-builds the updrift side: +38 m at GIS 6 against +7.
2. **By 1978 the observed fillet is almost all downdrift, and wider and deeper than the model's.** GIS 3-5 have retreated 43, 129 and 114 m; the model's 0.25 run 12, 59 and 78 m. The model's downdrift trough is half as deep and misses GIS 3.
3. **Updrift, the observed accretion is small and sits one domain further north**: +7 m at GIS 6 and +28 m at GIS 7, against the model's +32 m at GIS 6 and nothing at GIS 7. So the model's updrift build-up is in the right amount but one domain too far south, packed against the face. In that narrow sense the groin's updrift effect does not reach far enough north.
4. **The model splits the groin's effect about evenly between the two sides** (1978, groin minus no groin: GIS 4-5 -21 and -48 m, GIS 6 +58 m); the record puts most of it downdrift. The corrected groin conserves sand across the face (what it takes from GIS 5 it puts at GIS 6), so it cannot make a downdrift deficit larger than the updrift gain; the record's asymmetry needs sand to leave the downdrift side by some other route, such as offshore or round the cape.
5. **This does not explain the 1996-2009 erosion at GIS 7-14.** In the build-up the coast north of the groin advanced; a groin reaching further north would add accretion there. That erosion is a separate signal (BE set 1's domain sources and sinks are the place it would show).

**Caveats.** One photo line per date; the 1978 line has no GIS 12, so its end-to-end line runs to GIS 11. 1968-1983 storms replay 1984-1999. The runs use the 1996 topography outside GIS 2-12.
