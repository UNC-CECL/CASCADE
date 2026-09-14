# 9-figures - the house style, and the shared map layers

Not a stage. What every figure in the project shares.

```
STYLE.md                  the written style: the palette pair, letters, what
                          never goes on the canvas, captions in CAPTIONS.md
HAT_figure_style_sheet.*  the style rendered, so it can be looked at
map_elements/             the domain polygons and the island outline, with
                          every shapefile sidecar
```

The style is applied through `scripts/hat_figure_style.py`. Import it rather
than copying values out of the sheet: the sheet is a picture of the style, not
the style.
