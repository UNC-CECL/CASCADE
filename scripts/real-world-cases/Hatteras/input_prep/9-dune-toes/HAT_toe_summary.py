"""
Summary figure of the 1996 dune toes (2026-11-05): the new domain toe along the island against the model's current
one, with the villages marked.

    python scripts/real-world-cases/Hatteras/input_prep/9-dune-toes/HAT_toe_summary.py

Reads (D:, 1-transects-1996/): 1-toes/domain_toes_1996.csv (HAT_toe_method_line.py, on your crest picks).
Writes 8-summary/toe_alongshore_1996.png (caption in supporting/CAPTIONS.md).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from HAT_toe_common import BERM_M, CREST, MIN_KEPT, OUT, TOE, TOES, YEAR  # noqa: E402

SUMMARY = OUT / "8-summary"
VILLAGES = {"Buxton": (7, 8), "Avon": (21, 31), "Salvo": (68, 71), "Waves": (73, 75), "Rodanthe": (78, 83)}


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    fs.apply_style()
    SUMMARY.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(TOES / f"domain_toes_{YEAR}.csv").sort_values("gis")
    g, new, cur = d.gis.to_numpy(), d.toe_m_navd88.to_numpy(), d.current_toe_m_navd88.to_numpy()
    filled = (d.source != "median").to_numpy()

    fig, ax = plt.subplots(figsize=(fs.figsize("double")[0], 3.4), layout="constrained")
    for name, (a, b) in VILLAGES.items():
        ax.axvspan(a - 0.5, b + 0.5, color="#EFE7DA", lw=0, zorder=0)
        ax.text((a + b) / 2, 1.01, name, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7.5,
                color="0.3")
    ax.axhline(BERM_M, color=CREST, lw=0.6, ls=(0, (2, 2)), zorder=1)
    ax.text(90.4, BERM_M + 0.04, "Berm", fontsize=7, color=CREST, ha="right", va="bottom")
    ax.plot(g, cur, color="0.55", lw=1.0, ls=(0, (4, 2)), zorder=2, label="Current toe")
    ax.plot(g, new, color=TOE, lw=1.8, zorder=3, label="New toe (your crests)")
    ax.plot(g[~filled], new[~filled], "o", ms=3, color=TOE, zorder=4)
    if filled.any():
        ax.plot(g[filled], new[filled], "o", ms=4.5, mfc="white", mec=TOE, mew=1.2, zorder=5, lw=0,
                label="Filled from neighbours (" + ", ".join(f"GIS {x}" for x in g[filled]) + ")")
    ax.set_xlim(0.5, 90.5)
    ax.set_ylim(BERM_M - 0.15, np.nanmax([new, cur]) + 0.4)
    ax.set_xlabel("GIS domain (south → north)")
    ax.set_ylabel("Dune toe (m NAVD88)")
    ax.legend(fontsize=7.5, frameon=False, loc="upper center", ncol=3)
    fs.open_frame(ax)

    png = SUMMARY / f"toe_alongshore_{YEAR}.png"
    cap = SUMMARY / "supporting" / "CAPTIONS.md"
    if cap.exists():                                             # one figure here: its caption alone
        cap.unlink()
    fs.save(fig, png, close=True)
    dd = new - cur
    fs.record_caption(png, (
        f"Dune toe along Hatteras Island from the {YEAR} lidar transects. Solid, the new toe: per domain the median "
        f"of the transects' toes (the point furthest below the line from the shoreline to the front-dune crest, no "
        f"lower than the berm), on crests checked by hand; dashed, the model's current toe; shaded, the villages; "
        f"open, a domain with fewer than {MIN_KEPT} transects kept, filled from its neighbours; dotted, the berm ({BERM_M:g} m NAVD88). New minus current: mean {np.nanmean(dd):+.2f} m, "
        f"{int((np.abs(dd) > 0.5).sum())} of {np.isfinite(dd).sum()} domains differ by more than 0.5 m."))
    print(f"wrote {png}")


if __name__ == "__main__":
    main()
