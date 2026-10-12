#!/usr/bin/env python3
"""
The two start lines side by side: island offset from the ±1-yr mean (v2) and the ±6-month mean (v3).

    python start_offsets_1yr_vs_6mo.py  ->  figures/start_offsets_1yr_vs_6mo.png

One column per start year (1996, 2009). Top: the unpadded offset each build hands the model,
each zeroed on its own most seaward domain. Bottom: ±6 months minus ±1 yr in that frame,
on one shared axis so the two years compare directly.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-10
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C_1984, C_1997, DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style, figsize,
    record_caption, save, town_bands)

# --- CONFIG ------------------------------------------------------------------
# year -> {version: (label, first, last of the CoastSat window)}
BUILDS = {
    1996: {"v2": ("±1 yr", "12 Oct 1995", "12 Oct 1997"),
           "v3": ("±6 months", "12 Apr 1996", "12 Apr 1997")},
    2009: {"v2": ("±1 yr", "17 Aug 2008", "17 Aug 2010"),
           "v3": ("±6 months", "17 Feb 2009", "17 Feb 2010")},
}
COLOURS = {"v2": C_1984, "v3": C_1997}           # red ±1 yr, blue ±6 months, in both columns
# The lines differ by metres on a 6 km profile: wide red under dashed blue so both show
STYLE = {"v2": dict(lw=3.0, alpha=0.6), "v3": dict(lw=1.0, ls="--")}
PNG = HERE / "figures" / "start_offsets_1yr_vs_6mo.png"
# -----------------------------------------------------------------------------


def unpadded(year, version):
    p = tv.offset_file(year, "unpadded", version=version, source="shoreline")
    return pd.read_csv(p).set_index("Domain_ID").iloc[:, 0]   # 2009 v2 heads its column "2010"


def main():
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=figsize("double", aspect=0.62), sharex=True,
                             sharey="row", constrained_layout=True,
                             gridspec_kw=dict(height_ratios=[1, 1]))
    stats = []
    for j, (year, builds) in enumerate(BUILDS.items()):
        off = {v: unpadded(year, v) for v in builds}
        x = off["v2"].index.to_numpy()
        top, bot = axes[0, j], axes[1, j]
        for v, (lab, a, b) in builds.items():
            top.plot(x, off[v], color=COLOURS[v], **STYLE[v], label=f"{lab} ({a} to {b})")
        top.legend(loc="lower left", fontsize=7, frameon=False)
        _title(top, j, f"{year} start")
        d = off["v3"] - off["v2"]
        bot.axhline(0, color=INK_MUTED, lw=0.6)
        bot.bar(x, d, color=C_1997, width=0.8)
        _title(bot, 2 + j, f"{year}: ±6 months − ±1 yr (mean {d.mean():+.1f} m)")
        bot.set_xlabel(DOMAIN_AXIS_LABEL)
        bot.set_xlim(x.min() - 0.5, x.max() + 0.5)
        for a_ in (top, bot):
            town_bands(a_, where="top", strip=0.08, label=(a_ is top))
        stats.append(f"{year}: mean {d.mean():+.1f} m, range {d.min():+.1f} to {d.max():+.1f} m")
        print(f"{year} v3 - v2: mean {d.mean():+.2f}, range {d.min():+.2f} .. {d.max():+.2f}")
    axes[0, 0].set_ylabel("Island offset (m)")
    axes[1, 0].set_ylabel("Difference (m)")
    save(fig, PNG, dpi=300, close=True)
    record_caption(PNG, (
        "The two candidate start lines for each period, as island offsets the model reads. "
        "Top: the unpadded offset built from the CoastSat mean shoreline over ±1 yr of the DEM "
        "date (v2, red; the current input) and over ±6 months (v3, blue), each zeroed on its "
        "own most seaward domain; dates in the legends. At this scale the two coincide. Bottom: ±6 months minus ±1 yr in that "
        "frame, shared axis; positive is a larger offset in the ±6-month build. "
        + "; ".join(stats) + ". The ±6-month lines were built with a minimum of 5 positions "
        "per transect, against 10 for ±1 yr."))
    print(f"wrote {PNG}")


if __name__ == "__main__":
    main()
