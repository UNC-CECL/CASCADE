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
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline.run_layout import resolve  # noqa: E402
from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              figsize, open_frame, record_caption)
from HAT_groin_sweep_config import WETDRY_CHANGE_TABLE  # noqa: E402

# The rig pads 11 real domains (D2-D12) with 15 buffer either side: D2 -> 15,
# D5 -> 18, D6 -> 19, D12 -> 25. The RIG's convention, which differs from
# production's -- see HAT_groin_hindcast_1967_2017.py:76.
RIG_BUFFER, RIG_FIRST_GIS, RIG_START_YEAR = 15, 2, 1967
DOMAINS = list(range(2, 13))
PAD = [RIG_BUFFER + (n - RIG_FIRST_GIS) for n in DOMAINS]

# The rig lives in output/rig_runs/, not output/raw_runs/ (moved 2026-08-31).
# It is a DIFFERENT GRID -- 41 domains against production's 120 -- and M is
# grid-specific, so mixing the two invited quoting a rig number as a
# production one. raw_runs is production only, and run_index.csv covers it.
RAW_RUNS = PROJECT_BASE_DIR / "output" / "rig_runs"
GROIN_RUN = "HAT_1967_2018_edge_calibrated_groin"
NO_GROIN_RUN = "HAT_1967_2018_edge_calibrated_no_groin"
OUT_DIR = (PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"
           / "full_life_gif")

EXPECT_M, EXPECT_F = 60.0, 0.6
UPDRIFT_GIS, DOWNDRIFT_GIS = 6, 5
GHOSTS = 5                       # how many past surveys stay on screen
HOLD_FRAMES = 6                  # frames held on the last year, so it can be read

# House palette (2026-09-11), replacing an Okabe-Ito set chosen in this file:
# the run under test is the ACCENT, the groin-off baseline BASE grey, the
# surveys INK, and the structure a muted guide line rather than a fourth hue.
MODEL_C = C["ACCENT"]            # the groin run
BASE_C = C["BASE"]               # the groin-OFF baseline
OBS_C = INK                      # the surveys
MARK_C = INK_MUTED               # the structure
GHOST_C = "0.86"                 # surveys already passed

# Documented structure history -- GROIN_PLAN.md section 1. The rig's own
# install year is 1970 (it keeps 1967-69 as a free control window); the
# documented installation is 1969. Both are shown rather than reconciled.
# The strip is a ruler, not data: four greys, light to dark as the structure
# degrades, so the timeline reads in one glance without spending four hues on
# it. It was a grey/blue/orange/red set until 2026-09-11.
PHASES = [
    (1967, 1969, "before the groin", "0.97"),
    (1970, 1995, "structure sound", "0.93"),
    (1996, 2002, "deteriorating", "0.88"),
    (2003, 2017, "failed, held at M·f", "0.82"),
]
EVENTS = [(1970, "installed"), (1996, "last repair"), (2003, "storm damage")]


def _matrix(run_name: str) -> Path:
    return resolve(RAW_RUNS / run_name, "matrix", run_name)


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
        resolve(RAW_RUNS / GROIN_RUN, "groin_csv", GROIN_RUN))
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

    apply_style()

    figure = plt.figure(figsize=figsize("double", aspect=0.66))
    grid = figure.add_gridspec(
        2, 1, height_ratios=[7.0, 0.75], hspace=0.42,
        left=0.115, right=0.975, top=0.805, bottom=0.135)
    axis = figure.add_subplot(grid[0])
    strip = figure.add_subplot(grid[1])

    # ---- static furniture, drawn once ------------------------------------
    # Title and subtitle stay ON the canvas: a GIF is watched standalone, with
    # no caption file beside it in a viewer. House sizes, not a 16.5 pt banner.
    figure.text(0.115, 0.955, "The Buxton groin over its whole life",
                fontsize=11, fontweight="bold", color=INK,
                ha="left", va="center")
    figure.text(0.115, 0.912,
                f"{len(survey_years)} dated wet/dry surveys · "
                f"M {applied_M:g} m/yr, f {applied_f:.1f}",
                fontsize=8, color=INK_MUTED, ha="left", va="center")

    handles = [
        Line2D([], [], color=OBS_C, marker="s", markersize=3.6,
               linestyle="--", linewidth=1.4,
               label="observed, wet/dry survey"),
        Line2D([], [], color=MODEL_C, marker="o", markersize=3.4,
               linewidth=1.8, label="modelled, groin on"),
        Line2D([], [], color=BASE_C, linestyle=":", linewidth=1.4,
               label="modelled, groin off"),
    ]
    figure.legend(handles=handles, loc="upper right",
                  bbox_to_anchor=(0.975, 0.975), fontsize=8,
                  handlelength=2.2, labelspacing=0.45, frameon=False)

    # The footnote paragraph goes to CAPTIONS.md beside the GIF, written after
    # the animation is saved.

    def frame(index):
        t = min(index, len(years) - 1)
        year = int(years[t])

        # ---- main panel ---------------------------------------------------
        axis.clear()
        axis.axhline(0.0, color=INK_MUTED, linewidth=0.8,
                     linestyle=(0, (4, 3)), zorder=1)
        axis.axvline(5.5, color=MARK_C, linewidth=0.8, linestyle=(0, (5, 3)),
                     zorder=2)
        axis.annotate("the groin field: D6 updrift, D5 downdrift",
                      xy=(5.40, 0.025), xycoords=("data", "axes fraction"),
                      rotation=90, ha="right", va="bottom", fontsize=7,
                      color=MARK_C)

        passed = [y for y in survey_years if y < year][-GHOSTS:]
        for rank, past in enumerate(passed):
            alpha = 0.18 + 0.34 * (rank + 1) / max(len(passed), 1)
            axis.plot(DOMAINS, observed[past], "-", color=GHOST_C,
                      linewidth=1.1, alpha=alpha, zorder=2)

        if no_groin is not None:
            axis.plot(DOMAINS, no_groin[t], ":", color=BASE_C, linewidth=1.4,
                      zorder=4)
        axis.plot(DOMAINS, groin[t], "-o", color=MODEL_C, linewidth=1.8,
                  markersize=3.4, zorder=5)

        current = [y for y in survey_years if y <= year]
        if current:
            latest = max(current)
            fresh = latest == year
            axis.plot(DOMAINS, observed[latest], "s--", color=OBS_C,
                      markersize=4.0 if fresh else 3.2,
                      linewidth=1.6 if fresh else 1.1,
                      alpha=1.0 if fresh else 0.78, zorder=6)
            coverage = int(np.isfinite(observed[latest]).sum())
            axis.annotate(
                f"survey {latest}   ({coverage}/{len(DOMAINS)} domains)"
                + ("   ● new this year" if fresh else ""),
                xy=(0.015, 0.965), xycoords="axes fraction", ha="left",
                va="top", fontsize=8, color=OBS_C,
                fontweight="bold" if fresh else "normal")

        axis.set_xlim(1.6, 12.4)
        axis.set_ylim(*ylim)
        axis.set_xticks(DOMAINS)
        axis.set_xlabel("GIS domain (D2 toward Cape Point → D12 north)",
                        labelpad=6)
        # labelpad clears the orientation cues at x = -0.088, which the
        # 24 pt pad of the pre-2026-09-11 layout ran through.
        axis.set_ylabel("shoreline change since 1967 (m)", labelpad=8)
        axis.grid(axis="y")
        axis.set_axisbelow(True)
        open_frame(axis)

        # Orientation cues, so "up = erosion" needs no caption to decode.
        axis.annotate("▲ erosion, landward", xy=(0.012, 0.975),
                      xycoords="axes fraction", ha="left", va="top",
                      fontsize=7.5, color=INK_MUTED,
                      bbox=dict(facecolor="white", alpha=0.85,
                                edgecolor="none",
                                boxstyle="square,pad=0.12"))
        axis.annotate("▼ accretion, seaward", xy=(0.012, 0.025),
                      xycoords="axes fraction", ha="left",
                      va="bottom", fontsize=7.5, color=INK_MUTED,
                      bbox=dict(facecolor="white", alpha=0.85,
                                edgecolor="none",
                                boxstyle="square,pad=0.12"))

        # Year as a quiet watermark, plus the state of the structure.
        axis.annotate(str(year), xy=(0.988, 0.965), xycoords="axes fraction",
                      ha="right", va="top", fontsize=24, color="0.88",
                      fontweight="bold", zorder=0)
        m_eff = rate_by_year.get(year, last_rate)
        axis.annotate(f"applied trapping rate {m_eff:.0f} m/yr",
                      xy=(0.988, 0.845), xycoords="axes fraction", ha="right",
                      va="top", fontsize=8, color=INK_MUTED)

        up_i, down_i = DOMAINS.index(UPDRIFT_GIS), DOMAINS.index(DOWNDRIFT_GIS)
        model_fillet = groin[t][down_i] - groin[t][up_i]
        line = f"fillet, D5 − D6:  modelled {model_fillet:+.0f} m"
        if current:
            latest = max(current)
            observed_fillet = (observed[latest][down_i]
                               - observed[latest][up_i])
            if np.isfinite(observed_fillet):
                line += f",  observed {observed_fillet:+.0f} m"
        # Bottom-RIGHT: the structure label occupies the bottom-left, and the
        # accretion half of the panel is empty in every frame.
        axis.annotate(line, xy=(0.988, 0.035), xycoords="axes fraction",
                      ha="right", va="bottom", fontsize=8, color=INK)

        # ---- timeline strip -------------------------------------------------
        strip.clear()
        for start, end, label, colour in PHASES:
            strip.axvspan(start, end + 1, color=colour, zorder=1)
            # A three-year segment cannot hold a label; the "installed"
            # event marker already names that stretch.
            if end - start >= 5:
                strip.annotate(label, xy=((start + end + 1) / 2, 0.5),
                               xycoords=("data", "axes fraction"),
                               ha="center", va="center", fontsize=7.5,
                               color=INK_MUTED)
        for event_year, label in EVENTS:
            strip.axvline(event_year, color=INK_MUTED, linewidth=0.8, zorder=3)
            strip.annotate(label, xy=(event_year, -0.30),
                           xycoords=("data", "axes fraction"), ha="center",
                           va="top", fontsize=7, color=INK_MUTED)
        strip.axvline(year, color=MODEL_C, linewidth=1.8, zorder=5)
        strip.plot([year], [1.0], marker="v", color=MODEL_C, markersize=6,
                   transform=strip.get_xaxis_transform(), clip_on=False,
                   zorder=6)
        strip.set_xlim(years[0], years[-1] + 1)
        strip.set_yticks([])
        strip.set_xticks([1970, 1980, 1990, 2000, 2010])
        strip.tick_params(axis="x", labelsize=8, pad=20)
        for side in ("top", "right", "left"):
            strip.spines[side].set_visible(False)
        return ()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "shoreline_D2-D12_1967_2017.gif"
    animation = FuncAnimation(figure, frame,
                              frames=len(years) + HOLD_FRAMES, blit=False)
    animation.save(out, writer=PillowWriter(fps=3))
    plt.close(figure)

    # The prose the frames used to carry, in a file beside the GIF.
    record_caption(
        out,
        "The modelled shoreline change of the 41-domain 1967 rig against every "
        "dated wet/dry survey, year by year from {first} to {last}, at "
        "M = {M:g} and f = {f:.1f}. The frames keep their own title, legend and "
        "year clock because a GIF is watched standalone; everything else about "
        "how to read it is here. The y axis is LANDWARD-POSITIVE, so erosion "
        "moves up and the panel reads as a plan view with the ocean below and "
        "the island above. The curves are NOT demeaned: a uniform alongshore "
        "offset belongs to the source/sink calibration rather than to the "
        "groin, so it is the shape around D5 and D6 that is the groin's to get "
        "right. Surveys cover 8 to 11 of the 11 drawn domains and the lines "
        "break at gaps rather than interpolating across them; the {ghosts} "
        "most recent past surveys stay on screen as faint lines so the "
        "trajectory is visible. The strip below the panel is the structure's "
        "own history, light to dark as it degrades, with the rig's 1970 "
        "install marked — the documented installation is 1969, and both are "
        "shown rather than reconciled."
        .format(first=years[0], last=years[-1], M=applied_M, f=applied_f,
                ghosts=GHOSTS))

    print(f"wrote {out}")
    print(f"  {len(years)} frames + {HOLD_FRAMES} held, {years[0]}-{years[-1]}, "
          f"{len(survey_years)} surveys inside the window")
    print("  y axis is LANDWARD-positive: erosion up, plan-view orientation")


if __name__ == "__main__":
    main()
