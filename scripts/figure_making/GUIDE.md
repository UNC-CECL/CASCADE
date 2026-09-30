# How a figure script is written

How a figure *looks* is [`STYLE.md`](STYLE.md) (generated from
`site_layer/hat_figure_style.py`; don't edit it by hand). How any script is
laid out, commented and signed is [`../STYLE.md`](../STYLE.md). This page is
the part in between: how a figure script in this repo is put together and
where its output goes.

A standalone version for colleagues outside the repo is in
[`template/`](template/).

---

## 1. Where the script lives

Under `figure_making/<subject>/`, by what the figure is **of**: `island/`,
`management/`, `shoreline/`, `model_output/`, `model/`, `pipeline/<step>/`.
A figure made by the script that makes the data (a 5-scr comparison, say)
stays with that script, and publishes the same way.

## 2. The skeleton

```python
"""
Shoreline position in 1996 and 2010, and the change between them.

    python shoreline_change_figure.py

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1997, DOMAIN_AXIS_LABEL, apply_style, caption, figure_dir,
    figsize, open_frame, save, title, town_bands)

# --- CONFIG ------------------------------------------------------------------
OUT = figure_dir("observations", "shoreline")
PUBLISH = True           # False holds the figure back while its layout is broken
# -----------------------------------------------------------------------------


# Run: style, draw, caption, save
def main() -> None:
    apply_style()
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", 0.42), layout="constrained")
    ...
    title(axes[0], 0, "Shoreline position")
    open_frame(axes[0])
    caption(fig, "(a) ... (b) ...")
    if PUBLISH:
        save(fig, OUT / "shoreline_change_1996_2010", close=True)


if __name__ == "__main__":
    main()
```

The calls that matter:

| call | what it does |
|---|---|
| `apply_style()` | the house type, ink and sizes; first thing, every time |
| `figsize("single" \| "double", aspect)` | the only way a figure gets its size: 90 or 190 mm |
| `title(ax, i, text)` | bold `(a)` at the left, title centred |
| `open_frame(ax)` / `spines_for_image(ax)` | a chart drops two spines; a map keeps all four |
| `town_bands(ax)`, `DOMAIN_AXIS_LABEL` | the village bands and the one alongshore axis label |
| `C_1984` / `C_1997` (`C["EARLY"]` / `C["LATE"]`) | the earlier/later pair; `C["BASE"]`, `C["ACCENT"]`, `C["REF"]` for the other meanings |
| `caption(fig, text)` | the caption, written to `supporting/CAPTIONS.md` on save; nothing on the canvas |
| `figure_dir(subject, ...)` | the output folder; raises on an unknown subject |
| `save(fig, path)` | PNG at 300 dpi, PDF under `supporting/` |

## 3. Where the figure goes

Never beside the script. `figure_dir()` resolves `output/figures/`:

| subject | folder | for |
|---|---|---|
| `site` | `1-site/` | the reach, the domains, one domain |
| `observations` | `2-observations/<...>/` | CoastSat, dune lines, mean shoreline |
| `inputs` | `3-model-inputs/<step>/` | how each model input is built (`INPUT_STEPS`) |
| `mechanics` | `4-model-mechanics/<...>/` | how Barrier3D, BRIE and CASCADE work |
| `results` | `5-results/` | hindcasts and scenarios |
| `talk` | `talk/<same path>` | projector versions |

A cross-run comparison goes to `COMPARISONS_ROOT`; when it is also a finished
manuscript figure, the same run of the same script writes it to `5-results/`
too, so the two copies cannot drift.

The folder shows PNGs only. The PDF, `CAPTIONS.md`, and any table or
`PROVENANCE.md` the script writes go under `supporting/`
(`support_dir(folder)`).

## 4. Once it draws

1. Open the PNG at 100%: type should read at print size, and nothing on the
   canvas should be a sentence.
2. Add the script to `STEPS` in `tools/regenerate_all_figures.py`, so it is
   redrawn with everything else.
3. Run `python scripts/figure_making/tools/figure_index.py` to put it in
   `output/figures/README.md`. A dash under "shows" means no caption was
   recorded.

## 5. Don't

- Put a title sentence, a statistic or a footnote on the image; that is the
  caption.
- Use red or blue for anything but the earlier/later pair in a figure that
  shows two years.
- Write working names in a legend ("v2", "as placed", "run B"): say what the
  thing is.
- Set a figure width by hand, or raise font sizes to rescue a figure that is
  too wide.
