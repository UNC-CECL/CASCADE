"""
How the dune toe is picked, on four real profiles: the cape, Buxton-Avon, the Avon-Salvo road and Rodanthe.

    python scripts/real-world-cases/Hatteras/input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toe_examples.py

For one cell of each site (the cell nearest its domain's median beach width), on the start's own arrays:
the profile, the beach start (first cell above the beach threshold), the straight line from it to the
crest, the toe HAT_dune_toes.py picks (the point furthest below that line) and the fixed toe the model
uses for that cell (the median over the surveys, dune_toe_cells.csv). Shows the clean case (Rodanthe)
beside the two ways one profile's pick goes wrong: a berm and runnel puts the toe in the runnel (cape,
Avon-Salvo), and a near-straight face puts it high on the face (Buxton-Avon); the survey median and the
screens correct both. Reads <start>/dune-toes/dune_toe_cells.csv (run HAT_dune_toes.py --start first).
Draws output/figures/3-model-inputs/1-domains/dune_toes/dune_toes_four_sites.png.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-29
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
for _path in (_HAT_DIR, _HERE.parent, _HERE.parents[1] / "2-survey-crests"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
import HAT_dune_toes as T  # noqa: E402
import site_layer.hat_figure_style as fs  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
START = "1996-start"
# (title, GIS): one site per setting
SITES = (("Cape", 3), ("Buxton-Avon road", 12), ("Avon-Salvo road", 45), ("Rodanthe north", 86))
# -----------------------------------------------------------------------------


# The start's own profiles (m MHW) and dune windows for the sites' domains
def start_profiles():
    out = {}
    for survey, g, z, i0, i1 in T.domain_profiles(START, [g for _, g in SITES]):
        if survey != "start":
            break
        out[g] = (z, i0, i1)
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fs.apply_style()

    cells = pd.read_csv(T.S.tv.product_dir(START) / "dune-toes" / "dune_toe_cells.csv").set_index(["gis", "cell"])
    profiles = start_profiles()
    mhw = T.MHW_NAVD88_M
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.4))
    for ax, (name, g) in zip(axs.flat, SITES):
        z, i0, i1 = profiles[g]
        r = T.crest_and_toe(z, i0, i1, T.S.X)
        # The cell nearest the domain's median beach width, among cells with a pick
        ok = np.flatnonzero(~np.isnan(r["toe"]))
        w = (r["crest_idx"][ok] - r["beach_start_idx"][ok]) * 10
        cell = int(ok[np.argmin(np.abs(w - np.median(w)))])
        s, c, k = int(r["beach_start_idx"][cell]), int(r["crest_idx"][cell]), int(r["toe_idx"][cell])
        x = np.arange(z.shape[1]) * 10.0
        prof = z[cell] + mhw
        hi = min(z.shape[1], c + 6)
        ax.fill_between(x[:hi], prof[:hi], -1, color="#F2E6C9", lw=0)
        ax.plot(x[:hi], prof[:hi], "-", color="#5A4632", lw=1.8, label="start profile")
        ax.plot([x[s], x[c]], [prof[s], prof[c]], "--", color="0.45", lw=1.2, label="beach start to crest")
        ax.plot(x[s], prof[s], "o", color="#1B9E77", ms=8, label="beach start")
        ax.plot(x[c], prof[c], "^", color="#7570B3", ms=10, label="crest")
        line_k = prof[s] + (prof[c] - prof[s]) * (k - s) / (c - s)
        ax.vlines(x[k], prof[k], line_k, color="#D95F02", lw=1.2, ls=":")
        ax.plot(x[k], prof[k], "o", color="#D95F02", ms=10, mec="white", mew=1.5,
                label="toe on this profile: furthest below the line")
        used = cells.toe_m_navd88.get((g, cell), np.nan)
        if np.isfinite(used):
            ax.axhline(used, color="#D95F02", lw=1.0, ls="-.", alpha=0.8,
                       label=f"toe the model uses: survey median ({used:.2f} m)")
        ax.axhline(mhw, color="#2171B5", lw=0.8)
        ax.text(x[0] + 2, mhw + 0.08, "MHW", color="#2171B5", fontsize=7)
        width = cells.beach_width_m.get((g, cell), np.nan)
        ax.set_title(f"{name}: GIS {g}, cell {cell}   (beach {(c - s) * 10:.0f} m here, median {width:.0f} m)",
                     fontsize=9, loc="left")
        ax.set_xlim(x[0], x[hi - 1])
        ax.set_ylim(-0.5, max(prof[:hi].max(), 3) + 0.8)
        ax.set_xlabel("Distance from the ocean edge of the profile (m)")
        ax.set_ylabel("Elevation (m NAVD88)")
    axs[0, 0].legend(fontsize=7, loc="upper left")
    fig.tight_layout()

    # Write the figure and its caption
    out = fs.figure_dir("inputs", "1-domains", "dune_toes") / "dune_toes_four_sites.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=200, facecolor="white")
    fs.record_caption(out, (
        f"How the dune toe is picked, on four {START} profiles (the cell nearest each domain's median beach width). "
        "Green, the beach start (first cell above the beach threshold); purple, the crest; dashed, the line "
        "between them; orange, the toe on this profile, the point furthest below the line; dash-dot, the toe the "
        "model uses for the cell, the median over the surveys that pass the screens (dune_toe_cells.csv). Rodanthe "
        "is the clean case. At the cape and on the Avon-Salvo road a berm and runnel put this profile's toe in the "
        "runnel; at Buxton-Avon a near-straight face puts it high on the face. The screens and the survey median "
        "correct both. The beach width is the beach start to the crest, so it includes the dune face."))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
