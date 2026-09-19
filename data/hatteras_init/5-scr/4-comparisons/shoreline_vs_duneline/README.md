# shoreline_vs_duneline — does the dune line move with the shoreline?

One question folder since 2026-09-19 (Hannah). It replaced three sibling
folders that asked the same question: `duneline_vs_coastsat/` (now
`net_change/<window>/`), `net_change_1996_2024/` (now `net_change/chains/`)
and `projected_vs_duneline/` (now `projected/`). Every alongshore figure is in
the same form: **net change in metres, seaward positive, shoreline blue and
dune line red, the gap between them shaded as beach-width change** (solid grey
where the beach widened, hatched where it narrowed), with the village, groin,
pier, shoal and fill marks.

```
net_change/                      the shoreline as OBSERVED net change
    <window>/                    1984_2004, 1996_2010, 2004_2024, 2010_2024, 1996_2024
        alongshore_dune_vs_coastsat.png   the two net changes and the gap
        scatter_dune_vs_coastsat.png      dune rate against shoreline rate, per domain
        supporting/                       domain_comparison.csv, PROVENANCE.md, CAPTIONS.md
    alongshore_four_windows.png  the four model windows stacked on one y axis
    chains/
        net_change_chain_1996_2010_2024.png   1997–2023 above its two halves
                                              (1997–2009, 2009–2023); also published
                                              to output/figures/shoreline/
        supporting/  domain_comparison.csv, island_summary.csv
projected/                       the shoreline as the LRR PROJECTED over the dune interval
    1996_2024/   two_panel, shaded_gap, overlay; domain_comparison.csv, PROVENANCE.md
```

| folder | shoreline side | dune side | script (`scripts/input_prep/5-scr/`) |
|---|---|---|---|
| `net_change/<window>/` | `3-rates/coastsat/endpoint`: mean position ±6 months about each dune-line date, end minus start | `3-rates/duneline/endpoint` | `duneline_vs_coastsat/duneline_vs_coastsat.py --start-year S --end-year E` (then `--grid` for the stack) |
| `net_change/chains/` | as above, 1996_2024 and its halves | as above | `net_change/net_change_1996_2024.py` |
| `projected/1996_2024/` | `3-rates/coastsat/lrr/1996_2024` LRR x the dune interval (25.7 yr) | as above | `projected_vs_duneline/projected_vs_duneline.py` |

Both sides are always measured between the same two dates (the dune-line
survey dates; the 2023 one is assumed to be 1 July). Beach-width change is
shoreline change minus dune-line change. The two use different transects
(CoastSat ~10 per domain, dune line ~5), so they meet as domain means.

The scatter figures stay in m/yr; everything alongshore is in metres since
2026-09-19 (before that the per-window alongshore figures were m/yr with no
gap). Paths resolve through `hat_observed_rates.SHORELINE_VS_DUNELINE`,
`DUNELINE_VS_COASTSAT`, `NET_CHANGE_1996_2024` and `PROJECTED_VS_DUNELINE`.
