# 8-overwash-analysis

The observed overwash record on Hatteras Island, 1984–2024, and the figures
of it (the record figures cover the model windows, 1996–2025). Scripts: `scripts/input_prep/8-overwash-analysis/`.

```
1-observations/                     THE RECORD, and the record reshaped
    Hatteras_Overwash_Data.xlsx     One row per image assessed (Overwash_Matrix),
                                    1 / 0 / blank per domain; the storm reference
                                    sheet; edit here.
    overwash_observations.csv       long form: Obs_ID, date, year, season, domain, overwash
    storms_by_image.csv             per storm: last day, first image after it, days between
2-record/                           what the record shows (reorganised 2026-10-19)
    1-where/                        which parts of the island overwash most, all images together
        overwash_where_record_1996_2025.png       domains shaded by images with overwash,
                                                  19 images taken 1997-2025, three hotspot zooms
    2-when/                         which images showed it, domain by domain
        overwash_by_image_1997_2009.png           image x domain grid + storms, calibration window
        overwash_by_image_2010_2025.png           the same, test window
        overwash_by_image_1996_2025.png           the same, full model window, window bars on the left
        overwash_map_by_image_1997_2009.png       the island beside one column per image; same
        overwash_map_by_image_2010_2025.png       shade scale in both
        overwash_map_by_image_1997_2009_2010_2025.png   the two side by side (--both)
3-vs-footprint/                     the record against the 1984 footprint
    overwash_vs_footprint_alongshore.png    images, footprint bars and flags, by domain
    overwash_vs_footprint_summary.png       shift by overwash status; share overwashed per action
    overwash_vs_footprint_map.png           three alongshore sections: overwash, footprint, reading
    overwash_vs_footprint_map_island.png    the whole island, the same three layers
    tables/
        overwash_vs_footprint_by_domain.csv     the joined table: action, rows, shift, overwash, flags
        overwash_vs_footprint_contingency.csv   overwashed x add/none/remove, all and unflagged
        overwash_vs_footprint_summary.txt       the two readings in words, with the domain lists
4-vs-model/                         the record against modelled overwash (calendar-year runs since 2026-10-16)
    overwash_vs_model_1997_2009.png         each image vs the model since the previous image
    overwash_vs_model_2010_2025.png         (managed runs, groin on; 7-day grace as the storm table)
    storms_vs_overwash_1997_2025.png        skill, storms vs model dunes, per image and per domain
    overwash_timeline_1997_2025.png         plain view: each image paired with the model since the previous image
    tables/
        overwash_vs_model_cells.csv         every image x domain x run x threshold
        overwash_vs_model_summary.csv       both/observed-only/model-only counts, hit rates
CAPTIONS.md                         one entry per figure, headed by its folder; the
                                    figures carry no in-image text
```

`archive/superseded_20260910/` (the May 2026 heatmaps and Period 1 drafts, with
the workbook as it was on 2026-05-11) was deleted on 2026-10-01. The workbook
copy and the WHY.md come back from git
(`git log --diff-filter=D --oneline -- data/hatteras_init/8-overwash-analysis/archive/superseded_20260910`,
then `git show <commit>^:<path>`); the PNGs were git-ignored and are gone.

Grouped by job since 2026-09-18. Before that the folders were split by file
type (`observations/`, `figures/`, `tables/`), which put the footprint
comparison's figures in `figures/vs-footprint/` and its tables in a separate
top-level `vs-footprint/`. Resolve every path through
`scripts/site_layer/hat_overwash.py`; do not type them. The
model-comparison script in `scripts/analyze_output/overwash/` reads the
workbook through it too.

The windows are the model's: calibration (storms 1997-2009), test (storms
2010-2025) and the full record, labelled 1996-2025 like its run (storms
1997-2025). An image belongs to a window when it was taken in one of the
window's storm years, the same images 4-vs-model compares, so the 1984-1996
images appear only in the vs-footprint figures. Before 2026-10-19 the figures
used 1984-2004 / 2004-2024 / 1984-2024 and sat in heatmaps/ and map/; their
PDFs are in git history. A same-layout map of the overwash_events polygons in the
ArcGIS project (1984-2004, year field only) was made and removed the same day:
its years sit outside the model windows. Its script was never committed. Images to 1998 are the Hapke and Henderson (2007)
delineations; from 2004 they are read from Google Earth imagery.
The domain axis is the CASCADE domain, 1 at Cape Point to 90 at Pea
Island; the map draws the boxes from `D:/Hatteras_GIS/domains.geojson`.

The footprint comparison (2026-09-10) takes the nine images strictly between
the frames the two dune lines were digitised on (19 Sep 1984 out, 12 Oct 1997
in); overwash is expected with rows added; a rows-added domain without overwash
is listed as unexplained, not counted against the footprint, because the images
are years apart and washover fades.

The PNGs are not tracked by git (`*.png` is ignored repository-wide); the
PDFs beside them are. Rerun the scripts to rebuild the PNGs.
