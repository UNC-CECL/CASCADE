#!/usr/bin/env python3
"""Animated D2-D12 shoreline over the groin's whole life, 1967-2017.

WHAT THIS ADDS OVER `HAT_groin_zoom_gifs.py`
    That script animates the two hindcast windows and draws the observed change
    as a FIXED endpoint target, because inside 1984-2004 the observation is a
    single endpoint. Over the full life it is not: the wet/dry record carries
    25 dated surveys between 1967 and 2023 across D2-D12, so THE OBSERVATIONS
    CAN ANIMATE TOO. Each frame shows the most recent survey at or before that
    model year, with the surveys already passed left behind as fading ghosts.

    That is the comparison the hindcast windows cannot show -- the module
    building a fillet from a standing start at installation, holding it, and
    then losing it through the deterioration ramp, against a survey record that
    is doing the same thing on its own clock.

ORIENTATION: EROSION IS UP, SO THE PANEL READS AS A PLAN VIEW
    The y axis is LANDWARD-POSITIVE. A retreating shoreline moves UP the panel
    and an accreting one moves down, so the reader is looking down on the
    island with the ocean below the axis and the island above it. This is the
    OPPOSITE of `HAT_groin_zoom_gifs.py`, which is seaward-positive -- the two
    must not be read side by side without noticing.

    Both source arrays are already landward-positive and are therefore NOT
    negated here:
      `Change_from_wetdry_1967_*.csv`  landward-positive (a rising value is
          retreat) -- which is why the fillet is built as downdrift minus
          updrift throughout this directory.
      `shoreline_matrix.npy`           Barrier3D's x_s, landward-positive and
          already in METRES -- no dam conversion. (A dam->m rescale added to
          the zoom gifs on 2026-08-30 made every curve ten times too large;
          anything that rescales this must be re-checked against the cell's own
          shoreline_change_rate.csv.)

    The fillet is therefore reported as D5 - D6, matching GROIN_PLAN.md, rather
    than the D6 - D5 the earlier seaward-positive revision of this script used.

WHAT IS PLOTTED
      model, M = 60 f = 0.6   the rig's groin run
      model, groin OFF        the paired baseline at the same edge calibration
      observed                the most recent survey <= this year, with its own
                              year and coverage; the previous five as ghosts

    NOT DEMEANED, unlike the zoom gifs. The question here is how well the
    module reproduces the MEASURED POSITION, so the level is left in. A uniform
    alongshore offset is owned by the source/sink calibration, not the groin --
    so read a whole-curve offset as that term's business, and read the SHAPE
    around D5/D6 as the groin's.

Usage:
    python HAT_groin_full_life_gif.py

Writes output/groin_sweep/figures/full_life_gif/shoreline_D2-D12_1967_2017.gif

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from HAT_groin_sweep_config import WETDRY_CHANGE_TABLE  # noqa: E402

# The rig pads 11 real domains (D2-D12) with 15 buffer either side: D2 -> 15,
# D5 -> 18, D6 -> 19, D12 -> 25. The RIG's convention, which differs from
# production's -- see HAT_groin_hindcast_1967_2017.py:76.
RIG_BUFFER, RIG_FIRST_GIS, RIG_START_YEAR = 15, 2, 1967
DOMAINS = list(range(2, 13))
PAD = [RIG_BUFFER + (n - RIG_FIRST_GIS) for n in DOMAINS]

RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
GROIN_RUN = "HAT_1967_2018_edge_calibrated_groin"
NO_GROIN_RUN = "HAT_1967_2018_edge_calibrated_no_groin"
OUT_DIR = (PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"
           / "full_life_gif")

EXPECT_M, EXPECT_F = 60.0, 0.6
UPDRIFT_GIS, DOWNDRIFT_GIS = 6, 5
GHOSTS = 5                       # how many past surveys stay on screen
HOLD_FRAMES = 6                  # frames held on the last year, so it can be read

# Okabe-Ito, colour-vision-safe and muted enough to print.
MODEL_C = "#D55E00"              # vermillion -- the groin run
BASE_C = "#8C8C8C"               # grey       -- the groin-OFF baseline
OBS_C = "#111111"                # near-black -- the surveys
MARK_C = "#0072B2"               # blue       -- the structure
GHOST_C = "#CFD3D8"
INK, SUBTLE, GRID = "#1A1A1A", "#6A6F76", "#E3E6EA"

# Documented structure history -- GROIN_PLAN.md section 1. The rig's own
# install year is 1970 (it keeps 1967-69 as a free control window); the
# documented installation is 1969. Both are shown rather than reconciled.
PHASES = [
    (1967, 1969, "before the groin", "#E9ECEF"),
    (1970, 1995, "structure sound", "#CDE3F0"),
    (1996, 2002, "deteriorating", "#F6DFC8"),
    (2003, 2017, "failed — held at M·f", "#EFCFCF"),
]
EVENTS = [(1970, "installed"), (1996, "last repair"), (2003, "storm damage")]


def _matrix(run_name: str) -> Path:
    return RAW_RUNS / run_name / f"{run_name}_shoreline_matrix.npy"


def load_change(run_name: str):
    """Change since 1967 per domain, LANDWARD-positive, metres."""
    path = _matrix(run_name)
    if not path.is_file():
        return None
    matrix = np.load(path)[:, PAD]        # metres already, landward-+ already
    return matrix - matrix[0]


def observed_by_year() -> dict:
    """{survey year: change since 1967 over DOMAINS, LANDWARD-positive}.

    NaNs are preserved: several surveys cover only 8-10 of the 11 domains, and
    interpolating across a gap would invent a shoreline.
    """
    frame = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    columns: dict = {}
    for column in frame.columns:
        # The SECOND year is the survey year; the first is the 1967 datum.
        match = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", column)
        if match:
            columns.setdefault(int(match.group(1)), []).append(column)
    out = {}
    for year, cols in sorted(columns.items()):
        # 1972 was surveyed twice (July and August); average the pair.
        stacked = np.array([[frame.loc[d, c] for d in DOMAINS] for c in cols],
                           dtype=float)
        # An all-NaN column for a domain is a real gap in the survey, not an
        # error -- keep the NaN and let the line break there.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            out[year] = np.nanmean(stacked, axis=0)
    return out


def main() -> None:
    groin = load_change(GROIN_RUN)
    if groin is None:
        raise SystemExit(
            f"missing rig output: {_matrix(GROIN_RUN)}\n"
            "Run HAT_groin_hindcast_1967_2017.py first (run key 'groin').")
    no_groin = load_change(NO_GROIN_RUN)

    diagnostics = pd.read_csv(
        RAW_RUNS / GROIN_RUN / f"{GROIN_RUN}_groin_diagnostics.csv")
    active = diagnostics[diagnostics["groin_active"]]
    applied_M = float(active["trapping_rate_applied_m_yr"].max())
    applied_f = float(active["trapping_rate_applied_m_yr"].min()) / applied_M
    if not (np.isclose(applied_M, EXPECT_M)
            and np.isclose(applied_f, EXPECT_F, atol=0.02)):
        raise SystemExit(
            f"the run directory reports M = {applied_M:g}, f = {applied_f:.2f}, "
            f"not {EXPECT_M:g} / {EXPECT_F:g}. A rig run directory does not name "
            "its own parameters -- re-run HAT_groin_hindcast_1967_2017.py.")
    rate_by_year = dict(zip(diagnostics["model_year"],
                            diagnostics["trapping_rate_applied_m_yr"]))
    # The diagnostics end one year before the shoreline matrix does.
    last_rate = float(diagnostics["trapping_rate_applied_m_yr"].iloc[-1])

    observed = observed_by_year()
    years = RIG_START_YEAR + np.arange(groin.shape[0])
    survey_years = [y for y in observed if years[0] <= y <= years[-1]]

    stack = [groin] + ([no_groin] if no_groin is not None else [])
    finite = np.concatenate([np.asarray(a).ravel() for a in stack]
                            + [v[np.isfinite(v)] for v in observed.values()])
    low, high = np.nanmin(finite), np.nanmax(finite)
    pad = 0.12 * (high - low)
    ylim = (low - pad, high + pad)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10.5,
        "axes.linewidth": 0.8,
        "axes.edgecolor": "#4A4F55",
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 4.0, "ytick.major.size": 4.0,
        "xtick.color": "#4A4F55", "ytick.color": "#4A4F55",
        "text.color": INK, "axes.labelcolor": INK,
        "legend.frameon": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })

    figure = plt.figure(figsize=(10.8, 7.4))
    grid = figure.add_gridspec(
        2, 1, height_ratios=[7.0, 0.75], hspace=0.42,
        left=0.115, right=0.975, top=0.805, bottom=0.135)
    axis = figure.add_subplot(grid[0])
    strip = figure.add_subplot(grid[1])

    # ---- static furniture, drawn once ------------------------------------
    figure.text(0.115, 0.947, "The Buxton groin over its whole life",
                fontsize=16.5, fontweight="bold", color=INK,
                ha="left", va="center")
    figure.text(0.115, 0.905,
                f"Modelled shoreline change against {len(survey_years)} dated "
                f"wet/dry surveys   ·   M = {applied_M:g} m/yr, "
                f"f = {applied_f:.1f}   ·   41-domain 1967 rig",
                fontsize=10.5, color=SUBTLE, ha="left", va="center")

    handles = [
        Line2D([], [], color=OBS_C, marker="s", markersize=6, linestyle="--",
               linewidth=1.6, label="Observed (wet/dry survey)"),
        Line2D([], [], color=MODEL_C, marker="o", markersize=5.5,
               linewidth=2.4, label="Model, groin on"),
        Line2D([], [], color=BASE_C, linestyle=":", linewidth=2.0,
               label="Model, groin off"),
    ]
    figure.legend(handles=handles, loc="upper right",
                  bbox_to_anchor=(0.975, 0.975), fontsize=9.5,
                  handlelength=2.4, labelspacing=0.55)

    figure.text(
        0.115, 0.022,
        "Landward-positive: erosion moves UP, so the panel reads as a plan view "
        "with the ocean below and the island above. Not demeaned — a uniform "
        "alongshore offset belongs to the source/sink calibration, not the "
        "groin, so read the shape around D5/D6 as the groin's. Surveys cover "
        "8–11 of the 11 domains; lines break at gaps rather than interpolating.",
        fontsize=7.9, color=SUBTLE, ha="left", va="bottom", wrap=True)

    def frame(index):
        t = min(index, len(years) - 1)
        year = int(years[t])

        # ---- main panel ---------------------------------------------------
        axis.clear()
        axis.axhline(0.0, color="#B9BEC4", linewidth=0.9, zorder=1)
        axis.axvline(5.5, color=MARK_C, linewidth=1.3, linestyle=(0, (5, 3)),
                     zorder=2, alpha=0.85)
        axis.annotate("Buxton groin field\nD6 updrift  |  D5 downdrift",
                      xy=(5.40, 0.025), xycoords=("data", "axes fraction"),
                      rotation=90, ha="right", va="bottom", fontsize=8.2,
                      color=MARK_C)

        passed = [y for y in survey_years if y < year][-GHOSTS:]
        for rank, past in enumerate(passed):
            alpha = 0.18 + 0.34 * (rank + 1) / max(len(passed), 1)
            axis.plot(DOMAINS, observed[past], "-", color=GHOST_C,
                      linewidth=1.1, alpha=alpha, zorder=2)

        if no_groin is not None:
            axis.plot(DOMAINS, no_groin[t], ":", color=BASE_C, linewidth=2.0,
                      zorder=4)
        axis.plot(DOMAINS, groin[t], "-o", color=MODEL_C, linewidth=2.4,
                  markersize=5.5, zorder=5)

        current = [y for y in survey_years if y <= year]
        if current:
            latest = max(current)
            fresh = latest == year
            axis.plot(DOMAINS, observed[latest], "s--", color=OBS_C,
                      markersize=6.5 if fresh else 5.2,
                      linewidth=2.0 if fresh else 1.4,
                      alpha=1.0 if fresh else 0.78, zorder=6)
            coverage = int(np.isfinite(observed[latest]).sum())
            axis.annotate(
                f"survey {latest}   ({coverage}/{len(DOMAINS)} domains)"
                + ("   ● new this year" if fresh else ""),
                xy=(0.015, 0.965), xycoords="axes fraction", ha="left",
                va="top", fontsize=9.0, color=OBS_C,
                fontweight="bold" if fresh else "normal")

        axis.set_xlim(1.6, 12.4)
        axis.set_ylim(*ylim)
        axis.set_xticks(DOMAINS)
        axis.set_xlabel("Alongshore GIS domain      "
                        "(D2 = toward Cape Point   →   D12 = north)",
                        labelpad=7)
        axis.set_ylabel("Shoreline change since 1967 (m)", labelpad=26)
        axis.grid(axis="y", color=GRID, linewidth=0.7)
        axis.set_axisbelow(True)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)

        # Orientation cues, so "up = erosion" needs no caption to decode.
        axis.annotate("erosion\nlandward", xy=(-0.088, 0.94),
                      xycoords="axes fraction", ha="center", va="top",
                      fontsize=8.4, color=SUBTLE)
        axis.annotate("▲", xy=(-0.088, 0.955), xycoords="axes fraction",
                      ha="center", va="bottom", fontsize=7.5, color=SUBTLE)
        axis.annotate("accretion\nseaward", xy=(-0.088, 0.06),
                      xycoords="axes fraction", ha="center", va="bottom",
                      fontsize=8.4, color=SUBTLE)
        axis.annotate("▼", xy=(-0.088, 0.045), xycoords="axes fraction",
                      ha="center", va="top", fontsize=7.5, color=SUBTLE)

        # Year as a quiet watermark, plus the state of the structure.
        axis.annotate(str(year), xy=(0.988, 0.965), xycoords="axes fraction",
                      ha="right", va="top", fontsize=31, color="#DADEE3",
                      fontweight="bold", zorder=0)
        m_eff = rate_by_year.get(year, last_rate)
        axis.annotate(f"applied M$_{{eff}}$ = {m_eff:.0f} m/yr",
                      xy=(0.988, 0.845), xycoords="axes fraction", ha="right",
                      va="top", fontsize=9.2, color=SUBTLE)

        up_i, down_i = DOMAINS.index(UPDRIFT_GIS), DOMAINS.index(DOWNDRIFT_GIS)
        model_fillet = groin[t][down_i] - groin[t][up_i]
        line = f"fillet  D5 − D6        model  {model_fillet:+.0f} m"
        if current:
            latest = max(current)
            observed_fillet = (observed[latest][down_i]
                               - observed[latest][up_i])
            if np.isfinite(observed_fillet):
                line += f"        observed  {observed_fillet:+.0f} m"
        # Bottom-RIGHT: the structure label occupies the bottom-left, and the
        # accretion half of the panel is empty in every frame.
        axis.annotate(line, xy=(0.988, 0.035), xycoords="axes fraction",
                      ha="right", va="bottom", fontsize=9.4, color=INK)

        # ---- timeline strip -------------------------------------------------
        strip.clear()
        for start, end, label, colour in PHASES:
            strip.axvspan(start, end + 1, color=colour, zorder=1)
            # A three-year segment cannot hold a label; the "installed"
            # event marker already names that stretch.
            if end - start >= 5:
                strip.annotate(label, xy=((start + end + 1) / 2, 0.5),
                               xycoords=("data", "axes fraction"),
                               ha="center", va="center", fontsize=8.0,
                               color="#3E434A")
        for event_year, label in EVENTS:
            strip.axvline(event_year, color="#6A6F76", linewidth=0.9, zorder=3)
            strip.annotate(label, xy=(event_year, -0.30),
                           xycoords=("data", "axes fraction"), ha="center",
                           va="top", fontsize=7.6, color=SUBTLE)
        strip.axvline(year, color=MODEL_C, linewidth=2.4, zorder=5)
        strip.plot([year], [1.0], marker="v", color=MODEL_C, markersize=9,
                   transform=strip.get_xaxis_transform(), clip_on=False,
                   zorder=6)
        strip.set_xlim(years[0], years[-1] + 1)
        strip.set_yticks([])
        strip.set_xticks([1970, 1980, 1990, 2000, 2010])
        strip.tick_params(axis="x", labelsize=9, pad=22)
        for side in ("top", "right", "left"):
            strip.spines[side].set_visible(False)
        return ()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "shoreline_D2-D12_1967_2017.gif"
    animation = FuncAnimation(figure, frame,
                              frames=len(years) + HOLD_FRAMES, blit=False)
    animation.save(out, writer=PillowWriter(fps=3))
    plt.close(figure)
    print(f"wrote {out}")
    print(f"  {len(years)} frames + {HOLD_FRAMES} held, {years[0]}-{years[-1]}, "
          f"{len(survey_years)} surveys inside the window")
    print("  y axis is LANDWARD-positive: erosion up, plan-view orientation")


if __name__ == "__main__":
    main()
