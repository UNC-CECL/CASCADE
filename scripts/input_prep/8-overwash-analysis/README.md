# 8-overwash-analysis — scripts

Figures of the observed overwash record on Hatteras Island, per image and
per CASCADE domain, for the two model periods. The record itself is the
workbook in `data/hatteras_init/8-overwash-analysis/observations/`; every
output lands under `data/hatteras_init/8-overwash-analysis/` (see the README
there).

| file | what it does |
|---|---|
| `overwash_data.py` | Reads the workbook (observation matrix, storm reference sheet) and decides, from dates, which image first shows each storm. No plotting. Imported by both scripts. |
| `overwash_heatmap_multiperiod.py` | The heatmap figure, one per period (`period1`, `period2`, `combined`), plus the two tables and `CAPTIONS.md`. |
| `overwash_map_periods.py` | The island map, one figure per period, with the domains shaded by images with overwash and the date strip beside it; `--both` adds the side-by-side figure. Needs `D:/Hatteras_GIS` (domain boxes, coastline). |
| `HAT_overwash_vs_footprint.py` | Sets the overwash seen between the two dune-line frames (Aug 1985 to Oct 1997, nine images) against the rows the 1984 reconstruction adds and removes, per domain; contingency, lists, flags and three figures (alongshore, summary, three-panel island map). Reads the footprint, the road relocation table and the shoreline rates as inputs. |
| `superseded_20260910/` | The May 2026 single-period script. Kept for the record; its paths are dead. |

## Run

```
python overwash_heatmap_multiperiod.py                     # all three periods
python overwash_heatmap_multiperiod.py period2             # one period
python overwash_heatmap_multiperiod.py combined --uniform-years
python overwash_map_periods.py                             # one figure per period
python overwash_map_periods.py --both                      # plus side by side
python HAT_overwash_vs_footprint.py                        # overwash vs the 1984 footprint
```

Each script replaces only its own entries in `CAPTIONS.md` (`upsert_caption`
in `overwash_data.py`), so the run order does not matter.

## The storm-to-image rule

An image shows a storm if it is the first image taken on or after the
storm's last listed day (`Storm_Reference`, column `Search_GE_After`), with
a 7-day grace because the sheet's date ranges run to dissipation. The one
hand override is the May 2022 nor'easter, routed to the October 2023 image
on Hannah's own note. Storms not in the sheet (Ida 2009, Debby 2024) are
listed in `EXTRA_STORMS`. The decisions are written out to
`tables/storms_by_image.csv`, so a change in the sheet is visible there
rather than only in the figure.

Style: `scripts/hat_figure_style.py`. No in-image titles or footnotes; the
words are in `CAPTIONS.md`.
