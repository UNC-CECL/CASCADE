#!/usr/bin/env python3
"""Where a higher Hs changes the source/sink correction, zone by zone.

WHAT THIS ANSWERS
    The totals say the required correction field shrinks ~6% at Hs 3.0. They
    also hide the finding: the reaches move in OPPOSITE directions, and the
    ones that improve are not the ones carrying most of the correction. A
    single number for the island would report a modest win and conceal that
    the mid-island got worse.

WHY DIVERGING, AND WHY ORDERED SOUTH TO NORTH
    The quantity is a signed change either side of "no difference", which is a
    polarity encoding: two hues with a neutral midpoint, never a sequential
    ramp. And the zones are laid out by their domain range rather than sorted
    by value, so the panel reads as an alongshore profile -- which is what
    exposes that the improvement is concentrated at one end of the island.

INPUT
    The two pass-0 calibrations under output/hs_experiment/, produced by
    HAT_be_zone_residual_fit.py with HAT_BE_OUTPUT_DIR redirected. Nothing
    here reads or writes the production calibration.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# HOUSE STYLE: one typeface and one palette across every figure in this
# project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
# found by searching upward (ORGANIZATION.md rule 5). This file drew in
# matplotlib's defaults until 2026-09-17 -- it never called apply_style().
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

from HAT_plot_sensitivity import HOUSE_STYLE, panel_label, tidy  # noqa: E402

plt.rcParams.update(HOUSE_STYLE)

EXPERIMENT = PROJECT_BASE_DIR / "output" / "hs_experiment"
CONTROL = EXPERIMENT / "02_zones_Hs2p5" / "be_zone_metrics.csv"
TEST = EXPERIMENT / "03_zones_Hs3" / "be_zone_metrics.csv"
OUT_DIR = EXPERIMENT / "comparison"

# THE ENCODING IS NOT FIXED, so it is detected rather than assumed. The
# analysis writes through whatever encoding stdout has: run at a Windows
# console it emits cp1252 (the case HAT_be_apply_fit_to_config.py's RATES_ENCODING documents),
# run with stdout redirected to a file it emits UTF-8. Assuming either one
# turns the en-dash in "Buxton-Avon Transition" into mojibake in the zone
# labels -- a replacement character one way, "a-EUR-quote" the other.
METRICS_ENCODINGS = ("utf-8", "cp1252")

# A diverging pair with a neutral midpoint: less correction needed is the good
# direction and gets the cool hue, more correction the warm one. Deliberately
# NOT the sensitivity figures' sequential ramp -- that encodes magnitude along
# one hue, and this quantity has a sign.
BETTER, WORSE, NEUTRAL = "#2166AC", "#B2182B", "#8A8F94"
INK, INK_MUTED = "#222222", "#666666"


def load(path):
    """One pass-0 calibration, indexed by domain.

    Tries each encoding in turn and takes the first that decodes cleanly. A
    wrong guess does not raise -- cp1252 decodes any byte -- so UTF-8 is tried
    first and only a genuine decode failure falls through to it.
    """
    frame = None
    for encoding in METRICS_ENCODINGS:
        try:
            frame = pd.read_csv(path, encoding=encoding).set_index("domain")
            break
        except UnicodeDecodeError:
            continue
    if frame is None:
        raise ValueError(f"{path} decoded as none of {METRICS_ENCODINGS}")

    # The en-dash in "Buxton-Avon Transition" is ALREADY LOST in the file: the
    # analysis wrote a literal U+FFFD because its own output encoding could not
    # represent the character. No read encoding recovers it, so it is repaired
    # here for display. Fixing it at the source would mean the analysis writing
    # UTF-8 explicitly, which is a change to a script this experiment is
    # deliberately not modifying.
    frame["physical_zone"] = frame["physical_zone"].str.replace(
        "�", "–", regex=False)
    return frame


def zone_table(control, test, period):
    """Per-zone RMS residual for both arms, ordered south to north."""
    col = f"smooth_residual_{period}"
    rows = []
    for zone, group in control.groupby("physical_zone"):
        idx = group.index
        a = float(np.sqrt((control.loc[idx, col] ** 2).mean()))
        b = float(np.sqrt((test.loc[idx, col] ** 2).mean()))
        # NOT `first`/`last`: those are DataFrame methods, so a column of either
        # name is reachable by [] but NOT by attribute -- row.first silently
        # returns a bound method, which is how the zone labels came out as
        # "<bound method NDFrame.first of zone ...>" on the first draft.
        rows.append(dict(zone=zone, domain_from=int(idx.min()),
                         domain_to=int(idx.max()),
                         control=a, test=b, change_pct=100 * (b / a - 1)))
    # Geographic order, not rank order -- see the module docstring.
    return pd.DataFrame(rows).sort_values("domain_from").reset_index(drop=True)


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
