# 8-overwash-analysis

The observed overwash record on Hatteras Island, 1984–2024, and the figures
of it. Scripts: `scripts/input_prep/8-overwash-analysis/`.

```
1-observations/                     THE RECORD, and the record reshaped
    Hatteras_Overwash_Data.xlsx     One row per image assessed (Overwash_Matrix),
                                    1 / 0 / blank per domain; the storm reference
                                    sheet; edit here.
    overwash_observations.csv       long form: Obs_ID, date, year, season, domain, overwash
    storms_by_image.csv             per storm: last day, first image after it, days between
2-record/                           what the record shows
    heatmaps/
        overwash_heatmap_period1.png            1984–2004
        overwash_heatmap_period2.png            2004–2024
        overwash_heatmap_combined.png           1984–2024, period bars on the left
        overwash_heatmap_combined_uniform.png   same, one row per year (no gap collapse)
    map/
        overwash_map_period1.png                the island 1984–2004, domains shaded by
        overwash_map_period2.png                images with overwash; same shade scale in both
                                                (overwash_map_periods.png, side by side, with --both)
3-vs-footprint/                     the record against the 1984 footprint
    overwash_vs_footprint_alongshore.png    images, footprint bars and flags, by domain
    overwash_vs_footprint_summary.png       shift by overwash status; share overwashed per action
    overwash_vs_footprint_map.png           three alongshore sections: overwash, footprint, reading
    overwash_vs_footprint_map_island.png    the whole island, the same three layers
    tables/
        overwash_vs_footprint_by_domain.csv     the joined table: action, rows, shift, overwash, flags
        overwash_vs_footprint_contingency.csv   overwashed x add/none/remove, all and unflagged
        overwash_vs_footprint_summary.txt       the two readings in words, with the domain lists
CAPTIONS.md                         one entry per figure, headed by its folder; the
                                    figures carry no in-image text
archive/superseded_20260910/
    figures-2026-05-18/             the three May 2026 heatmaps the current ones replace
    Figures_1984_2004-2026-05-12/   earlier Period 1 drafts and model-comparison panels,
                                    with the workbook as it was on 2026-05-11
```

Grouped by job since 2026-09-18. Before that the folders were split by file
type (`observations/`, `figures/`, `tables/`), which put the footprint
comparison's figures in `figures/vs-footprint/` and its tables in a separate
top-level `vs-footprint/`. Resolve every path through
`scripts/site_layer/hat_overwash.py`; do not type them. The two
model-comparison scripts in `scripts/analyze_output/overwash/` read the
workbook through it too.

Period 1 rows are the Hapke and Henderson (2007) delineations; Period 2
rows are read from Google Earth imagery. 2004 belongs to both periods.
The domain axis is the CASCADE domain, 1 at Cape Point to 90 at Pea
Island; the map draws the boxes from `D:/Hatteras_GIS/domains.geojson`.

The footprint comparison (2026-09-10) takes the nine images strictly between
the frames the two dune lines were digitised on (19 Sep 1984 out, 12 Oct 1997
in); overwash is expected with rows added; a rows-added domain without overwash
is listed as unexplained, not counted against the footprint, because the images
are years apart and washover fades.

The PNGs are not tracked by git (`*.png` is ignored repository-wide); the
PDFs beside them are. Rerun the scripts to rebuild the PNGs.
