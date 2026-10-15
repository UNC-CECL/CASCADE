#!/usr/bin/env python3
"""
Start, observed end and modelled end shoreline positions around the groin, for each groin version.

    python groin_positions.py   ->  figures/groin_positions_detrended_1996_2009_and_2009_2025.png

Model frame, landward up (as the runner's position GIFs and full_window's positions
figure). Start = the run's t=0 shoreline; CoastSat end = start moved by the
observed net change (coastsat/net_change); modelled ends for no groin, the two
dipoles and the pinned blocking groin, all on the solved edgeBE ends.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-05
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1997, DOMAIN_AXIS_LABEL, INK, _title, apply_style, figsize, open_frame,
    record_caption, save, town_bands)
from site_layer.hat_observed_rates import net_change_domain_csv  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW = REPO / "output" / "raw_runs"
STUDY = RAW / "experiments" / "groin" / "2026-10-05-blocking-fit-dem-to-dem"
PERIODS = ((1996, 2009), (2009, 2025))
PAD = 15                                   # padded index of GIS 1
ZOOM = (1, 20)
C_START = "0.78"
# label, member folder (None = matrix no-groin), colour, line style
SERIES = (
    ("no groin", None, C["BASE"], "-"),
    ("dipole M 60, f 0.6 (old pin)", "dipole_M60_f0.6", C["ADDED"], "--"),
    ("dipole M 12, f 0.3 (best dipole)", "dipole_M12_f0.3", C["ADDED"], ":"),
    ("blocking b 0.6, f 0.6 (new pin)", "b0.60_f0.6", C["ACCENT"], "-"),
)
# -----------------------------------------------------------------------------


def window(p):
    return f"{p[0]}_{p[1]}"


def run_dir(period, member):
    if member is None:
        fill = "_nourish" if period[0] == 2009 else ""
        return RAW / "matrix" / window(period) / "edgeBE" / (
            f"HAT_{window(period)}_edgeBE_offsetmetres_road_bdm{fill}_nogroin")
    return sorted((STUDY / member / window(period) / "edgeBE").glob("*/"))[-1]


def matrix(rd):
    return np.load(next(rd.glob("*_shoreline_matrix.npy")))


def main():
    apply_style()
    gis = np.arange(ZOOM[0], ZOOM[1] + 1)
    idx = PAD + gis - 1
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=7.6), constrained_layout=True)
    rows = []
    for i, (ax, period) in enumerate(zip(axes, PERIODS)):
        start = matrix(run_dir(period, None))[0][idx]
        obs_net = pd.read_csv(net_change_domain_csv(*period), index_col=0)["net_change_m"].loc[gis].values
        obs_end = start - obs_net                   # landward positive frame
        # One straight line through the start, removed from every curve: the cape's tilt goes, metres stay
        trend = np.polyval(np.polyfit(gis, start, 1), gis)
        start, obs_end = start - trend, obs_end - trend
        ax.plot(gis, start, color=C_START, lw=3.2, zorder=2, solid_capstyle="round")
        for label, member, col, ls in SERIES:
            end = matrix(run_dir(period, member))[-1][idx] - trend
            ax.plot(gis, end, color=col, ls=ls, lw=1.5, zorder=4)
            rows.append(dict(period=window(period), run=label,
                             **{f"gis{g}_m": float(v) for g, v in zip(gis, end)}))
        ax.plot(gis, obs_end, color=C_1997, lw=1.3, marker="o", ms=4, mfc=C_1997, mec="white",
                mew=0.6, zorder=6)
        rows.append(dict(period=window(period), run="CoastSat end",
                         **{f"gis{g}_m": float(v) for g, v in zip(gis, obs_end)}))
        rows.append(dict(period=window(period), run="start",
                         **{f"gis{g}_m": float(v) for g, v in zip(gis, start)}))
        ax.set_xlim(*ZOOM)
        ax.set_xticks(gis)
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.12 * (hi - lo))
        ax.axvline(5.5, color=C["GROIN"], lw=0.8, zorder=1)
        ax.text(5.6, lo + 0.03 * (hi - lo), "groin", fontsize=6.5, color=C["GROIN"])
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=(i == 0))
        ax.set_ylabel("Position about the start's\nalongshore trend (m), landward ▲")
        role = "Calibration" if i == 0 else "Test"
        end_label = "1 Jan 2009" if i == 0 else "1 Jan 2025"
        _title(ax, i, f"{role} {period[0]}-{period[1]}: start and end shorelines (model end {end_label})")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    handles = [Line2D([], [], color=C_START, lw=3.2, label="start (DEM-centred CoastSat mean)"),
               Line2D([], [], color=C_1997, lw=1.3, marker="o", ms=4, mec="white",
                      label="CoastSat end (start + observed net change)")]
    handles += [Line2D([], [], color=c, ls=ls, lw=1.5, label=f"model end, {lab}")
                for lab, _, c, ls in SERIES]
    fig.legend(handles=handles, loc="outside upper center", ncol=2, fontsize=6.8, frameon=False)

    (HERE / "tables").mkdir(exist_ok=True)
    pd.DataFrame(rows).round(2).to_csv(HERE / "tables" / "groin_positions.csv", index=False)
    png = HERE / "figures" / "groin_positions_detrended_1996_2009_and_2009_2025.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Shoreline positions around the Buxton groin (red line, between GIS 5 and 6) in the model's "
        "cross-shore frame, landward up and the ocean at the bottom, with one straight line (least "
        "squares through the start shoreline over GIS 1-20) removed from every curve, so the cape's "
        "~1.7 km alongshore tilt does not hide differences of tens of metres; distances between "
        "curves are unchanged metres. Thick grey: the start shoreline "
        "(each run's t = 0, built on the CoastSat mean over +/-1 yr of the start DEM). Blue dots: the "
        "observed end, the start moved by the CoastSat net change (end-window mean minus start-window "
        "mean, domain means; ends centred on 2009-08-17 and 2025-08-17). Lines: the modelled end "
        "shoreline (1 Jan 2009 and 1 Jan 2025) with no groin, the old dipole (M 60, f 0.6), the best "
        "option-A dipole (M 12, f 0.3) and the pinned blocking groin (b 0.6, f 0.6). All runs: full "
        "management, the solved edgeBE ends, relocations off, failure instant from the 2004 step. "
        "(a) Calibration 1996-2009. (b) Test 2009-2025, after the Rodanthe 2014, Buxton 2017 and "
        "2022 fills. Domain values, unsmoothed."))
    print(f"wrote {png.relative_to(REPO)}")


if __name__ == "__main__":
    main()
