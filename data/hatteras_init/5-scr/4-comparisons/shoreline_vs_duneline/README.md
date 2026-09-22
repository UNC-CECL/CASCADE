# shoreline_vs_duneline — does the dune line move with the shoreline?

> **Lost?** [`FIGURES.md`](../../../../../FIGURES.md) is the one-page index of which figure answers which question. [`WINDOWS.md`](../../WINDOWS.md) says which window is which — two chains, and one context window nothing is graded against.

Every folder here compares a **CoastSat shoreline** reading against the
**dune line**. The dune side is the same measurement in all of them:

> the end digitized dune line minus the start line, along the 100 m transects,
> averaged over the ~5 in each 500 m domain. **Measured, never fitted** — read
> from `3-rates/duneline/endpoint/<window>/`, never recomputed here.

**So the only thing that differs between these folders is how the SHORELINE is
read.** That is why each folder is named for both sides (Hannah, 2026-09-22):
the name alone says what was measured and how, with nothing to look up. Before
that they were `endpoint_net_change/` and `total_change/`, which named neither
side — and `total_change/` collided with a folder of the same name in
`3-rates/` meaning a different comparison.

```
coastsat_endpoint_vs_duneline_endpoint/     BOTH sides are two snapshots
    <window>/          _alongshore.png  the two changes and the gap
                       _scatter.png     dune rate against shoreline rate
    all_windows_stacked/                 was chains/
    ..._four_windows_alongshore.png      the four model windows on one axis

coastsat_total_change_vs_duneline_endpoint/ shoreline = the window's OWN trend
    <window>/          _two_panel.png, _shaded_gap.png, _overlay.png,
                       domain_comparison.csv, PROVENANCE.md
    all_windows_stacked/  ..._stacked.png         1996–2024 above its halves
                          ..._halves_overlay.png  the two halves only, in the
                                                  same form as the projected
                                                  sheet, for side-by-side
    change_between_periods/              was difference/

coastsat_projected_vs_duneline_endpoint/    shoreline = the LONG-TERM trend,
    1996_2010/, 2010_2024/               carried onto a window it was NOT
                                         fitted on (built 2026-09-22)
    all_windows_stacked/  ..._halves_overlay.png  the two halves on one sheet:
                                                  ONE prediction, two dune-line
                                                  outcomes

smoothing_test/                          BOTH sheets above with both
                                         curves LOESS-smoothed at 7
                                         domains; a test, not a product

superseded_20260921/                     the pre-rename build; see its README
```

## How each shoreline side is read

| folder | shoreline side | rate fitted on | × |
|---|---|---|---|
| `coastsat_endpoint_vs_duneline_endpoint` | mean CoastSat position ±6 months about each dune-line date, differenced | nothing — no rate | — |
| `coastsat_total_change_vs_duneline_endpoint` | `3-rates/coastsat/lrr/<window>` | **that same window** | its calendar span (14, 14, 28 yr) |
| `coastsat_projected_vs_duneline_endpoint` | `3-rates/coastsat/lrr/1996_2024` | **the full period** | 14 yr |

The middle one is **total** shoreline change — the rate is evaluated over the
window it was fitted on, so nothing is extrapolated. The last is **projected**
— the rate is carried onto a window it was *not* fitted on. See
[`3-rates/README.md`](../../3-rates/README.md) for the vocabulary.

`coastsat_projected_vs_duneline_endpoint` has **no `1996_2024`**: over the full
period the rate window *is* the change window, so that case is
`coastsat_total_change_vs_duneline_endpoint/1996_2024`. It also has no
`change_between_periods/`: both halves carry the same 1996–2024 rate over the
same 14 yr, so the shoreline side is identical in them and the difference
between halves is zero by construction on that side.

**The two `_halves_overlay.png` sheets are the pair to read together.** They
are the same figure in the same form, differing only in where the shoreline
rate came from: `coastsat_projected_.../all_windows_stacked/` holds the
shoreline side FIXED at the 1996–2024 rate, while
`coastsat_total_change_.../all_windows_stacked/` lets each half use its own.
Between them they separate what the long-term trend *predicts* from what it
was *fitted on*.

It **does** have `all_windows_stacked/`, and that identical shoreline side is
exactly why it is worth drawing: `..._1996_2010_2024_stacked.png` puts the two
halves one above the other, so **one prediction** is read against **two
different dune-line outcomes**. The method is stated once in the header rather
than twice, and each panel title carries only what differs — the window and
that half's dune-line dates.

## Which interval each side spans

In `coastsat_endpoint_vs_duneline_endpoint`, both sides are measured between
the same two dates (the dune-line survey dates; the 2023 one is assumed to be
1 July), so the gap is beach-width change and nothing else.

The two rate-based folders are the exception, by choice (Hannah, 2026-09-21:
the year is what a period means, so the fit and the interval are both
calendar). The shoreline is carried over the full calendar span — 14, 14 and
28 yr — while the dune line spans 11.63, 14.09 and 25.72 yr, so its gap also
holds that much shoreline drift. It moves the island-mean beach width by 0.8 m
in the first half, 0.1 m in the second and 0.4 m over the whole. Every table
carries the interval-matched value beside the headline one
(`*_dune_interval_m`) rather than correcting it, **and every figure states both
spans in a header line above the panels**, so the mismatch cannot be missed.

## The house form

Every alongshore figure is the same: **net change in metres, seaward positive,
±100 m fixed axis**, the shoreline as the blue/red sign fill with the dune line
in black (or both as lines, shoreline blue and dune red), the gap between them
shaded as beach-width change — solid grey where the beach widened, hatched
where it narrowed — with the village, groin, pier, shoal and fill marks.
Anything beyond ±100 m is marked with a triangle at the edge and named in the
caption. The scatter figures stay in m/yr; everything alongshore has been in
metres since 2026-09-19.

Beach-width change is shoreline change minus dune-line change. The two sides
use different transects (CoastSat ~10 per domain, dune line ~5), so they meet
as domain means.

Paths resolve through `hat_observed_rates.SHORELINE_VS_DUNELINE`,
`COASTSAT_ENDPOINT_VS_DUNELINE`, `TOTAL_CHANGE_VS_DUNELINE`,
`PROJECTED_VS_DUNELINE_ENDPOINT` and `NET_CHANGE_1996_2024`.

| folder | script (`scripts/input_prep/5-scr/`) |
|---|---|
| `coastsat_endpoint_vs_duneline_endpoint/<window>/` | `4-comparisons/shoreline_vs_duneline/coastsat_vs_duneline.py --start-year S --end-year E` (then `--grid`) |
| `coastsat_endpoint_vs_duneline_endpoint/all_windows_stacked/` | `4-comparisons/shoreline_vs_duneline/net_change_vs_duneline.py` |
| `coastsat_total_change_vs_duneline_endpoint/` | `4-comparisons/shoreline_vs_duneline/total_change_vs_duneline.py` |
| `coastsat_projected_vs_duneline_endpoint/` | `4-comparisons/shoreline_vs_duneline/total_change_vs_duneline.py --product projected` |
