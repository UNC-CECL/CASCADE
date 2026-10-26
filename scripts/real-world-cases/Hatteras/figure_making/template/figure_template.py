"""
A figure in the house style: panels at print size, PNG + PDF, caption off the canvas.

    python figure_template.py

Writes figures/example_figure.png, with the PDF and CAPTIONS.md under
figures/supporting/. Copy it and replace draw() with your own panels; the
rules behind each choice are in README.md. Needs numpy, matplotlib.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- CONFIG ------------------------------------------------------------------
OUT_DIR = Path("figures")
NAME    = "example_figure"
WIDTH   = "double"          # "single" = 90 mm column, "double" = 190 mm page
ASPECT  = 0.42              # height / width
CAPTION = ("(a) Shoreline position along the reach in 1996 (red) and 2010 "
           "(blue). (b) Net change 1996-2010; positive is seaward.")
# -----------------------------------------------------------------------------

# House colours: the earlier/later pair, and one meaning per other colour
INK, INK_MUTED, GRID = "0.15", "0.42", "0.88"
EARLY, LATE = "#b2182b", "#2166ac"
EARLY_FILL, LATE_FILL = "#f4a582", "#92c5de"
BASE, ACCENT, REF = "#7f7f7f", "#7b3294", "#2c6e49"

STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "grid.color": GRID, "grid.linewidth": 0.5,
    "legend.frameon": False, "lines.solid_capstyle": "butt",
    "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.dpi": 300,
}


# Size at the printed width, so 9 pt type is 9 pt on the page
def figsize(width: str | float, aspect: float) -> tuple[float, float]:
    w = {"single": 3.54, "double": 7.48}.get(width, width)
    return float(w), min(float(w) * aspect, 9.4)


# Bold panel letter at the left, title centred
def title(ax, i: int, text: str) -> None:
    ax.set_title(f"({chr(97 + i)})", loc="left", fontweight="bold")
    ax.set_title(text, loc="center")


# Charts drop the top and right spines
def open_frame(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)


# PNG in the folder; PDF and the caption under supporting/
def save(fig, out_dir: Path, name: str, caption: str) -> Path:
    support = out_dir / "supporting"
    support.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{name}.png"
    fig.savefig(png)
    fig.savefig(support / f"{name}.pdf")
    captions = support / "CAPTIONS.md"
    entries = captions.read_text(encoding="utf-8").split("\n\n") if captions.exists() else []
    entries = [e for e in entries if e.strip() and not e.startswith(f"**`{png.name}`")]
    entries.append(f"**`{png.name}`.** {caption}")
    captions.write_text("\n\n".join(entries) + "\n", encoding="utf-8")
    return png


# Example panels -- replace with your own
def draw(fig, axes) -> None:
    rng = np.random.default_rng(0)
    x = np.arange(1, 91)
    early = 60 + 15 * np.sin(x / 14) + rng.normal(0, 2, x.size)
    late = early - 8 + 10 * np.cos(x / 20) + rng.normal(0, 2, x.size)

    ax = axes[0]
    ax.plot(x, early, color=EARLY, lw=1.2, label="1996")
    ax.plot(x, late, color=LATE, lw=1.2, label="2010")
    ax.fill_between(x, early, late, where=late < early, color=EARLY_FILL, lw=0)
    ax.fill_between(x, early, late, where=late >= early, color=LATE_FILL, lw=0)
    ax.set(xlabel="Alongshore position (south → north)", ylabel="Shoreline position (m)")
    title(ax, 0, "Shoreline position")

    ax = axes[1]
    change = late - early
    ax.bar(x, change, width=1.0, lw=0, color=BASE)   # red/blue are taken by the years
    ax.axhline(0, color=INK, lw=0.6)
    ax.set(xlabel="Alongshore position (south → north)", ylabel="Net change (m)")
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    title(ax, 1, "Change, 1996–2010")

    for ax in axes:
        open_frame(ax)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(labels))


# Run: apply the style, draw, save
def main() -> None:
    plt.rcParams.update(STYLE)
    fig, axes = plt.subplots(1, 2, figsize=figsize(WIDTH, ASPECT), layout="constrained")
    draw(fig, axes)
    png = save(fig, OUT_DIR, NAME, CAPTION)
    plt.close(fig)
    print(f"wrote {png.resolve()} (+ supporting/{NAME}.pdf, supporting/CAPTIONS.md)")


if __name__ == "__main__":
    main()
