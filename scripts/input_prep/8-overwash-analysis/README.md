# 8-overwash-analysis — scripts

Figures of the observed overwash record on Hatteras Island, per image and
per CASCADE domain, for the two model periods. The record itself is the
workbook in `data/hatteras_init/8-overwash-analysis/1-observations/`; every
output lands under `data/hatteras_init/8-overwash-analysis/` (see the README
there).

Split into steps 2026-09-22 to match the data tree, which was already
`1-observations/ 2-record/ 3-vs-footprint/` while the code sat flat.

| file | what it does |
|---|---|
| `1-observations/overwash_data.py` | Reads the workbook (observation matrix, storm reference sheet) and decides, from dates, which image first shows each storm. No plotting. Imported by all three of the others, which is why it sits in the step it builds rather than in a folder of its own. |
| `2-record/overwash_heatmap_multiperiod.py` | The heatmap figure, one per period (`period1`, `period2`, `combined`), plus the two tables and `CAPTIONS.md`. |
| `2-record/overwash_map_periods.py` | The island map, one figure per period, with the domains shaded by images with overwash and the date strip beside it; `--both` adds the side-by-side figure. Needs `D:/Hatteras_GIS` (domain boxes, coastline). |
| `3-vs-footprint/overwash_vs_footprint.py` | Sets the overwash seen between the two dune-line frames (Aug 1985 to Oct 1997, nine images) against the rows the 1984 reconstruction adds and removes, per domain; contingency, lists, flags and three figures (alongshore, summary, three-panel island map). Reads the footprint, the road relocation table and the shoreline rates as inputs. |
| `superseded_20260910/` | The May 2026 single-period script. Kept for the record; its paths are dead. |

## Run

```
python 2-record/overwash_heatmap_multiperiod.py                 # all three periods
python 2-record/overwash_heatmap_multiperiod.py period2         # one period
python 2-record/overwash_heatmap_multiperiod.py combined --uniform-years
python 2-record/overwash_map_periods.py                         # one per period
python 2-record/overwash_map_periods.py --both                  # plus side by side
python 3-vs-footprint/overwash_vs_footprint.py                  # vs the 1984 footprint
```

Each script replaces only its own entries in `CAPTIONS.md` (`upsert_caption`
in `1-observations/overwash_data.py`), so the run order does not matter.

## The storm-to-image rule

An image shows a storm if it is the first image taken on or after the
storm's last listed day (`Storm_Reference`, column `Search_GE_After`), with
a 7-day grace because the sheet's date ranges run to dissipation. The one
hand override is the May 2022 nor'easter, routed to the October 2023 image
on Hannah's own note. Storms not in the sheet (Ida 2009, Debby 2024) are
listed in `EXTRA_STORMS`. The decisions are written out to
`1-observations/storms_by_image.csv`, so a change in the sheet is visible there
rather than only in the figure.

Style: `scripts/site_layer/hat_figure_style.py`. No in-image titles or footnotes; the
words are in `CAPTIONS.md`.

## Naming and the sibling imports

No `HAT_` prefix; `overwash_vs_footprint.py` lost its on 2026-09-22, when it
was the only prefixed file of the four. See `../README.md` for which stages
are bare and which are not.

Three scripts import `overwash_data` by name, and `overwash_vs_footprint.py`
also imports the two in `2-record/` for their shared colours and geometry.
Each puts the folders it needs on `sys.path`, anchored on the `pyproject.toml`
found by searching upward — never counted from the file's own depth, so a
script can change step without breaking (rule 5).
