# 9-figures - the house style, and the shared map layers

Not a stage. What every figure in the project shares.

```
STYLE.md                  the written style: the palette pair, letters, what
                          never goes on the canvas, captions in CAPTIONS.md
HAT_figure_style_sheet.png  the style rendered, so it can be looked at
supporting/               its PDF, as every figure folder keeps its PDFs
map_elements/             the domain polygons and the island outline, with
                          every shapefile sidecar
```

A figure folder shows figures: PNGs at the top level, and the PDFs,
`CAPTIONS.md`, tables and `PROVENANCE.md` under `supporting/` (Hannah,
2026-09-15). `save()` and `caption()` place theirs; a script's own files
go through `support_dir(folder)`. Folders written before that date still
hold PDFs and captions at the top level until their script is re-run;
move the old copies out when it is, or the captions file is duplicated.

The style is applied through `scripts/hat_figure_style.py`. Import it rather
than copying values out of the sheet: the sheet is a picture of the style, not
the style.
