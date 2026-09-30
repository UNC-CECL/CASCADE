#!/usr/bin/env python3
"""
Where a higher Hs changes the source/sink correction, zone by zone (Hs 2.5 vs 3.0).

    python scripts/sensitivity_analysis/plot_hs_experiment.py

Reads the two pass-0 calibrations under output/calibration/hs/ and writes a
diverging bar figure and its CSV to output/calibration/hs/comparison/. Details: scripts/sensitivity_analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# House style (site_layer/hat_figure_style.py), applied at import
import sys as _sys
from pathlib import Path as _HP
_sys.path.insert(0, str(next(_q for _q in _HP(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import apply_style, figsize  # noqa: E402
apply_style()
import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if str(_HERE.parent) not in sys.path:
    sys.path.insert(0, str(_HERE.parent))

from plot_sensitivity import HOUSE_STYLE, panel_label, tidy  # noqa: E402

plt.rcParams.update(HOUSE_STYLE)

# --- CONFIG ------------------------------------------------------------------
EXPERIMENT = PROJECT_BASE_DIR / "output" / "calibration" / "hs"
CONTROL = EXPERIMENT / "02_zones_Hs2p5" / "be_zone_metrics.csv"
TEST = EXPERIMENT / "03_zones_Hs3" / "be_zone_metrics.csv"
OUT_DIR = EXPERIMENT / "comparison"
# -----------------------------------------------------------------------------

# Encoding detected, not assumed: the analysis writes cp1252 or UTF-8
METRICS_ENCODINGS = ("utf-8", "cp1252")

# Diverging pair: less correction is the cool hue, more the warm one
BETTER, WORSE, NEUTRAL = "#2166AC", "#B2182B", "#8A8F94"
INK, INK_MUTED = "#222222", "#666666"


# One pass-0 calibration, indexed by domain, in whichever encoding decodes cleanly
def load(path):
    frame = None
    for encoding in METRICS_ENCODINGS:
        try:
            frame = pd.read_csv(path, encoding=encoding).set_index("domain")
            break
        except UnicodeDecodeError:
            continue
    if frame is None:
        raise ValueError(f"{path} decoded as none of {METRICS_ENCODINGS}")

    # Repair the en dash the analysis already lost (a literal U+FFFD in the file)
    frame["physical_zone"] = frame["physical_zone"].str.replace(
        "�", "–", regex=False)
    return frame


# Per-zone RMS residual for both arms, ordered south to north
def zone_table(control, test, period):
    col = f"smooth_residual_{period}"
    rows = []
    for zone, group in control.groupby("physical_zone"):
        idx = group.index
        a = float(np.sqrt((control.loc[idx, col] ** 2).mean()))
        b = float(np.sqrt((test.loc[idx, col] ** 2).mean()))
        # Columns read with [], not attributes: first/last are DataFrame methods
        rows.append(dict(zone=zone, domain_from=int(idx.min()),
                         domain_to=int(idx.max()),
                         control=a, test=b, change_pct=100 * (b / a - 1)))
    # Geographic order, not rank order -- see the module docstring.
    return pd.DataFrame(rows).sort_values("domain_from").reset_index(drop=True)


# Run: both arms' zone residuals, the diverging bar figure and its CSV
def main():
    control, test = load(CONTROL), load(TEST)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    periods = (("p1", "1984–2004"), ("p2", "2004–2024"))
    tables = {key: zone_table(control, test, key) for key, _ in periods}
    order = tables["p1"]

    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.31), sharey=True,
                             constrained_layout=True)
    y = np.arange(len(order))[::-1]        # D1 at the top, so south is up

    for col, ((key, label), ax) in enumerate(zip(periods, axes)):
        table = tables[key]
        colors = [BETTER if v < 0 else WORSE for v in table.change_pct]
        ax.barh(y, table.change_pct, color=colors, height=0.62, zorder=3)
        ax.axvline(0, color=NEUTRAL, lw=1.0, zorder=4)

        for yi, (_, row) in zip(y, table.iterrows()):
            offset = 1.6 if row.change_pct >= 0 else -1.6
            ax.annotate(f"{row.change_pct:+.0f}%", (row.change_pct, yi),
                        xytext=(offset, 0), textcoords="offset points",
                        va="center", fontsize=8, color=INK,
                        ha="left" if row.change_pct >= 0 else "right")
        ax.set_title(label, pad=8)
        ax.set_xlim(-45, 40)
        ax.grid(axis="y", visible=False)
        tidy(ax, minor=False)
        panel_label(ax, "ab"[col], dx=-0.02 if col else -0.30)

    axes[0].set_yticks(y)
    axes[0].set_yticklabels(
        [f'{r["zone"]}\nD{r["domain_from"]}–{r["domain_to"]}'
         for _, r in order.iterrows()], fontsize=8)
    fig.supxlabel("Change in RMS source/sink correction required at Hs 3.0 (%)",
                  fontsize=9.5)

    fig.suptitle("Raising Hs to 3.0 m helps the north end and hurts the "
                 "mid-island\nBlue = less source/sink correction needed; "
                 "red = more. Groin held at M = 60, f = 0.6.",
                 fontsize=10)

    path = OUT_DIR / "zone_change_Hs2p5_vs_Hs3.png"
    fig.savefig(path)
    plt.close(fig)

    for key, label in periods:
        tables[key].insert(0, "period", label)
    pd.concat([tables[k] for k, _ in periods]).to_csv(
        OUT_DIR / "zone_change.csv", index=False)
    print(f"wrote {path.name} and zone_change.csv")
    print(pd.concat([tables[k] for k, _ in periods]).round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
