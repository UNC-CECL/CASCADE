# Hatteras figure style

Written 2026-10-02 18:28 by `scripts/site_layer/hat_figure_style.py` (`write_style_sheet()`); the
module is the source, this page is its rendering. `HAT_figure_style_sheet.png`
beside it shows every colour, the elevation classes, a chart and a map drawn
under the rules.

## How a script uses it

```python
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import apply_style, C, C_1984, C_1997, INK, INK_MUTED, _title, _scalebar, _north_arrow, caption
apply_style()                 # before any figure is made
```

A map also takes the map elements: `scale_bar_km`, `north_dart`, `water_label`,
`place_label`, `map_label`, `groin_marker`, `groin_handle`, `top_row`, `MAP_TEXT`.
`figure_making/island/study_area_figures.py` (`fig_study_area`) is the worked example.

Every figure script under `scripts/` that draws for this project calls
`apply_style()` first. `0-elevation/3-figures/HAT_plot_duneline_offset.py` also
re-exports these names, so `import HAT_plot_duneline_offset as off` still gives
`off.INK`, `off._title()` and so on to the scripts that take their map loaders
from it.

## The rules

| | |
|---|---|
| size | drawn at the printed width: `figsize("single")` = 3.54 in (90 mm), `figsize("double")` = 7.48 in (190 mm), height from `aspect`; never wider, so the type below is the type on the page |
| typeface | Arial, Helvetica, Liberation Sans, DejaVu Sans: the first one installed; 9 pt body, 10 pt panel titles, 8 pt ticks and legends |
| alongshore axis | one label, `DOMAIN_AXIS_LABEL` = "GIS domain (south → north)"; villages as light bands named once by `town_bands(ax)`; the endpoints (1 at Cape Point, 90 at Pea Island) go in the caption |
| ink | text and axes `0.15`, secondary text and rulers `0.42`, grid `0.88`; axes 0.6 pt |
| panels | a bold letter at the left of the title, the title centred (`_title(ax, i, text)`); inside the corner when the title is wide (`_letter_inside`) |
| maps | closed frame (`spines_for_image`), no coordinate ticks. Since 2026-10-02 the house map is the study-area map: a segmented scale bar, `scale_bar_km(ax)` (black and white segments with a white keyline, 0 / 5 / 10 km below it with a clear gap; `unit="m"` under 1 km), and the split-dart north arrow `north_dart(ax, c, north=...)`, half black and half white, pointing to true north in the map's own frame (a rotated map passes its rotation's north). `_scalebar` (says the cell count under 1 km) and `_north_arrow` remain for the model-grid panels. A labelled UTM frame needs neither |
| map labels | one rule, `MAP_TEXT`: white 8 pt type with a thin dark halo (1.2 pt). Water bodies upright in letter-spaced capitals (`water_label`, `spaced_caps`), villages and places italic (`place_label`), domain numbers and scale figures plain (`map_label`). On a pale map (an outline on a light canvas) the same rule in ink with a thin white halo, `MAP_TEXT_DARK`. A name never sits on land: a short leader takes it out over water |
| map corners | with an inset in the upper right, `top_row(ax, ax_in, handles, north=...)` puts the legend (`MAP_LEGEND`: translucent white, no frame, 8 pt, tight) top-aligned immediately left of the inset, the panel letter at the left edge on the same top line, and the north arrow just left of the legend; the scale bar sits alone in the lower left |
| map symbols | the Buxton groin field is `groin_marker` (red `#d7191c` bar with a white edge; `groin_handle()` for the legend, labelled "Buxton groins"); the study area on a locator inset is teal `#1b7f8c`, never red, so it cannot be read as the groins; NC-12 on imagery is road yellow `#ffd23f`; a locator inset carries a light graticule without labels, its spacing in the caption |
| charts | top and right spines off (`open_frame`), hairline grid on the value axis only when it helps |
| legends | frameless (a faint white backing when inside), outside the axes where the layout allows: `fig.legend(handles, loc="outside lower center", ncol=n, frameon=False)` under `constrained_layout` |
| vintages | the earlier line or surface is red `#b2182b`, the later blue `#2166ac`, everywhere the two are drawn together; the light fills `#f4a582` / `#92c5de` are the band between them |
| semantic colours | `C["BASE"]` #7f7f7f unmodified input · `C["ACCENT"]` #7b3294 the modification under test · `C["ROAD"]` #1a1a1a NC-12 · `C["ADDED"]` #c8880f fabricated ground · `C["WATER"]` #a8c8e0 · `C["REF"]` #2c6e49 a reference value |
| elevation | terrain, everywhere domain topography is drawn (since 2026-10-01; classes before): `elevation_cmap()` gives water below 0 m MHW one colour, terrain's navy `#333399`, and land the terrain ramp from 0 to 4 m, higher ground its top colour; colorbar ticks at 0, 0.5, 1, 1.5, 2, 3, 4 m |
| error surfaces | greyscale, no hue: `error_cmap()` (dark is worse; `reverse=True` where high is better). A scalar error or cost over a parameter grid is BACKGROUND, and all colour is reserved for what is marked on top of it -- the best cell, the chosen pair, a constraint, an iso-product curve |
| the canvas | no title sentences, statistics lines or footnote paragraphs on the image. That text goes in `supporting/CAPTIONS.md` beside the figure. ONE exception since 2026-09-22: a figure comparing two different measurements may carry a `compare_header()` line saying WHAT each side is and over what dates - a naming line, never a result. `caption(fig, text)` writes it there on the figure's next `savefig`; scripts with their own captions file (dune-line offset, footprint, road relocation) write it themselves |
| the folder | a figure folder shows figures: PNGs at the top level and nothing else. The PDFs, `CAPTIONS.md`, any table or `PROVENANCE.md` a figure script writes go under `supporting/` (`save()` and `record_caption()` do this; a script's own files use `support_dir(folder)`). Since 2026-09-15 |
| legend wording | no working vocabulary: not "today's setback", "v2"/"v3", "as placed", "blank". Say what the thing is: "setback measured on the 1996 surface", "1984 setback (model input)", "rows inserted landward of NC-12", "centreline unchanged between surveys" |
| output | `save(fig, path)`: a 300 dpi PNG and, under `supporting/`, a PDF with the same stem for anything drawn with lines and bars (`vector=False` for image-only panels); white background; `bbox_inches="tight"` only when nothing is positioned absolutely |
| semantic accent | `C["ACCENT"]` is purple since 2026-09-10; it was a red indistinguishable from the 1984 vintage red, so "the change under test" and "1984" read as one colour |

## Where it came from

The 2026-09-04 restyle of the dune-line figures (Hannah: "more academic /
professionally styled so they are informative and look good to present") set
these rules; the older `hat_figure_style.py` of the row-insert work carried the
colour semantics and the elevation classes. The two were merged here on
2026-09-10, with the 09-04 rules winning wherever they disagreed (typeface
order, legend frames, captions on the canvas), so that one style can be applied
across every figure.
