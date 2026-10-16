# 1-observations — what the groin did to the shoreline

Measured or digitized. Nothing here is fitted to the model, and nothing here changes when the model does.
Split into `data/` (what was recorded) and `analysis/` (what was measured from it) on 2026-10-12. On
2026-10-14 the gap analyses were split by source (photos, CoastSat, the two compared) and the wide and
map views renumbered 4 and 5.

```
data/                          what was recorded or digitized; no scripts
    structure_history.md           dates, repairs, decay, the fills, press sources; start here
    gis/                           GIS layers: groins_hatteras, wet_dry_groin, transects_100m, the
                                   offshore datum line, the CASCADE domains and study area; its README
                                   says which are copies and who reads each; ArcGIS originals in shapefiles/
analysis/                      one folder per study, in reading order: scripts at the top,
                               tables/, figures/ (PDFs and CAPTIONS.md in figures/supporting/)
    1-aerial_photo_gap/            PHOTOS ALONE: the 24 dated wet/dry surveys, 1967-2023, as change since
                                   1967 per domain (GIS 2-12); the gap across the groins from them. Its
                                   Change_from_wetdry table is the observed gap every fit reads,
                                   and stays at the folder's top level because code reads it by path
    2-coastsat_gap/                COASTSAT ALONE, narrow view: the gap across the groins each year,
                                   1984-2025, and when it stopped widening (1995)
    3-photos_vs_coastsat/          THE COMPARISON: what the gap is in each record, where they disagree
                                   (the 2004 survey), and the groin's state in 1996 from both
    4-shoreline_rates_by_era/      WIDE VIEW, ±60 km: change rates on every transect before, during and
                                   after the groin's working life, by decade, how far the effect reaches
    5-shoreline_maps/              MAP VIEW: the CoastSat shoreline each year, 1984-2025, as true-scale
                                   maps of the groin window and the whole coast (GIFs and one static
                                   figure); split out of the rates folder on 2026-10-13
```

## What the observations show

- **The gap, not new beach.** Both sides eroded; the north (updrift) side eroded less. By 2004 the
  south side had retreated 175 m since 1967 and the north side 25 m (`data/structure_history.md`).
- **Built fast, then slowly.** The gap grew about +11 m/yr in 1967-78, then about +1.3 m/yr to its
  1995 peak of 155 m (photos). CoastSat gives +4.5 m/yr over 1984-95.
- **Stopped widening in 1995**, the last repair, not at Isabel (2003): the break is 1995 (90% range
  1992-98) and holds under every check; a 2002-05 break gets 1-2% of the bootstraps.
- **Then the south side recovered, and the gap went with the fill.** After 1995 the gap narrowed about
  1.4-2.4 m/yr, but that is GIS 5 regaining what Gordon took (1996-2003) while GIS 6 held; from 2004
  both sides retreat with the coast north of the field (`analysis/2-coastsat_gap/` figure 6). There is
  no separate drop at 2004. It fell about 40 m around the 2017 fill, as the reported structural
  losses (2015-2021) began.
- **One date disagrees.** The photos and CoastSat agree within 13 m except the 2004 survey, 42 m
  above CoastSat. It is the only support for "held until 2004".
- **Reach.** In the 1970s the effect reached about 3 km north and 2 km south; the wide analysis
  cannot separate the groin from the 1971/1973 and 2017/2022 fills.

**The groin in 1996, the start of the hindcast** (`analysis/3-photos_vs_coastsat/`, figure 5): a
full-sized gap (135 m since 1967 in the photos, 87% of the peak), a year after the last repair, no
longer widening (−2.2 ± 1.8 m/yr in 1996-2003, not distinguishable from zero), and nothing marking a
collapse before the 2015-2021 losses. The 1995 peak is partly Hurricane Gordon's (1994) loss on the
south side.

## Who reads these

- `analysis/1-aerial_photo_gap/Change_from_wetdry_1967_D2_D12.csv`: every hindcast fit, the
  figures, `scripts/hatteras_ms/groin-sweep/HAT_groin_sweep_config.py`.
- `analysis/2-coastsat_gap/`: the groin fits from the 2026-10-08 schedule refit on
  (`../4-hindcast/2-blocking-1996-2025/` and `../4-hindcast/3-wave-driven-1967-2025/`), which import its CoastSat gap.
- `analysis/4-shoreline_rates_by_era/tables/groin_analysis_chainage_all.csv`: the dipole sweep's
  target (`scripts/hatteras_ms/groin-sweep/`), and the input to `analysis/5-shoreline_maps/`.
- `analysis/5-shoreline_maps/shoreline_plan_view_gif.py`: its map layers feed the corrected groin's
  plan-view figure (`../4-hindcast/3-wave-driven-1967-2025/2026-10-09-corrected-groin/plan_view_figure.py`).
- `data/gis/groins_hatteras.geojson`: `scripts/input_prep/4-mgmt-forcings/beach_nourishment.py`.
- `data/gis/transects_100m.geojson`: `analysis/4-shoreline_rates_by_era/`; the 2-brie-offset stage
  reads its byte-identical twin (`data/hatteras_init/2-brie-offset/transects/README.md`).
