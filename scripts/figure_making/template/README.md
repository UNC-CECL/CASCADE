# template — a figure in the house style, to share

One standalone script that draws a two-panel figure the way this project does:
at print size, in the house colours, with the caption kept off the image.
Meant to be **copied out of this repository** and handed to a colleague.

```
figure_template.py   draw() -> figures/<name>.png, with the PDF and
                     CAPTIONS.md under figures/supporting/
```

It imports nothing from this repository; it needs `numpy` and `matplotlib`.
Change CONFIG and replace `draw()` with your own panels. Everything else is
the style.

```
python figure_template.py
```

## The rules it follows

| | rule | why |
|---|---|---|
| **size** | drawn at the printed width: 90 mm (`"single"`) or 190 mm (`"double"`) | the 9 pt type stays 9 pt on the page instead of shrinking to 4 pt |
| **type** | Arial (or the first of Helvetica, Liberation Sans, DejaVu Sans), 9 pt body, 10 pt titles, 8 pt ticks and legends | readable in print and consistent across figures |
| **panels** | bold `(a)` at the left, title centred | one way to refer to a panel, in the caption and in the text |
| **charts** | top and right spines off; a hairline grid on the value axis only when it helps read values | less ink that isn't data |
| **two years** | the earlier is red `#b2182b`, the later blue `#2166ac`; light fills `#f4a582` / `#92c5de` for the band between | the same pair everywhere, so a reader learns it once; it survives greyscale and colour-blindness |
| **one meaning per colour** | grey `#7f7f7f` is the unmodified baseline, purple `#7b3294` the change under test, green `#2c6e49` a reference value | once red means "1996" it cannot also mean "erosion" in the same figure, which is why the example's change bars are grey |
| **legend** | no frame, outside the axes at the bottom when there is room | keeps the data area clear |
| **the canvas** | no title sentence, no statistics, no footnotes on the image | that text is the caption, and a caption can be edited without redrawing |
| **the caption** | written to `supporting/CAPTIONS.md`, one entry per figure, replaced on re-run | the figure still explains itself months later, out of its folder |
| **the folder** | PNGs at the top; the PDF and captions under `supporting/` | a figure folder shows figures |
| **output** | PNG at 300 dpi; a PDF for anything drawn with lines or bars | the PDF stays sharp in a manuscript |
| **wording** | labels say what the thing is ("shoreline in 1996"), never working names ("v2", "run B") | the reader never saw your folders |

## Before you trust a figure

Open the PNG at 100% and check it at the size it will be printed. If the type
looks small, the figure is too wide: make it `"single"` and adjust `ASPECT`,
rather than raising the font size.

---

Hannah A. Henry, Coastal Environmental Change Lab, University of North Carolina at Chapel Hill  
hahenry@unc.edu · version 2026-09-30
