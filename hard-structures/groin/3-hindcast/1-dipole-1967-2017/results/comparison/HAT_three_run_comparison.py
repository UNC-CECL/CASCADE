"""
Baseline vs nourishment-only vs nourishment + groin: three panels of the 1967-2017 Buxton rig runs against the 2018 survey.

    python HAT_three_run_comparison.py

Make the three runs first (README); set RUN_BASELINE / RUN_NOURISHMENT_ONLY /
RUN_FULL_MODEL to their folder names. Writes two static figures and two GIFs
to three_run_comparison/ beside this file.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- CONFIG ------------------------------------------------------------------
# The three run folder names (README: how each run is made)
RUN_BASELINE          = "HAT_1967_2018_edge_calibrated_no_BN_no_groin"
RUN_NOURISHMENT_ONLY  = "HAT_1967_2018_edge_calibrated_no_groin"
RUN_FULL_MODEL        = "HAT_1967_2018_edge_calibrated_groin"

HERE = Path(__file__).resolve().parent
PROJECT_BASE_DIR = str(next(p for p in HERE.parents if (p / "pyproject.toml").exists()))
OUTPUT_BASE_DIR  = os.path.join(PROJECT_BASE_DIR, "output", "calibration", "groin_rig")

# Geometry, which must match the runs
NUM_REAL_DOMAINS   = 11
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 2
LAST_FILE_NUMBER   = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1     # 12
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS

START_YEAR          = 1967
MODEL_FINAL_YEAR     = 2017   # model's TRUE final simulated year (row 50 of 51)
OBSERVED_FINAL_YEAR  = 2018   # full D2-D12 coverage; 2017 has only 6/11 (README)

WETDRY_CHANGE_TABLE = os.path.join(
    PROJECT_BASE_DIR, "hard-structures", "groin", "1-observations",
    "wetdry_photo_positions", "Change_from_wetdry_1967_D2_D12.csv",
)

# Red / amber / blue: problem -> partial fix -> full fix, on purpose (README)
BASELINE_COLOR = "#D32F2F"   # red -- no groin, no nourishment
NOURISH_COLOR  = "#F9A825"   # amber -- nourishment only
FULL_COLOR     = "#1565C0"   # blue -- nourishment + groin (full model)
GROIN_BOUNDARY_GIS = 5.5
GROIN_BOUNDARY_DISPLAY = GROIN_BOUNDARY_GIS - FIRST_FILE_NUMBER + 1   # 4.5 on the 1-11 scale
GROIN_LABEL_X_OFFSET = 0.15   # keeps the "Buxton groin" text off the dashed line
GROIN_LABEL_Y_FRACTION = 0.75   # 0 bottom, 1 top; 0.75 clears the legend (README)
DOMAIN_TICK_STEP = 2   # every other domain labelled, room for the large fonts

# Font sizes, large on purpose: the figure is shrunk to a third for the page (README)
TICK_LABEL_FONTSIZE  = 18
AXIS_LABEL_FONTSIZE  = 20
TITLE_FONTSIZE       = 20
SUPTITLE_FONTSIZE    = 22
LEGEND_FONTSIZE      = 16
GROIN_LABEL_FONTSIZE = 15
YEAR_TEXT_FONTSIZE   = 20

FIGURE_DPI = 300
GIF_FPS = 4    # frames per second in the animated version
GIF_DPI = 100  # kept lower than FIGURE_DPI to control file size
OUTPUT_DIR = str(HERE / "three_run_comparison")
# -----------------------------------------------------------------------------


# Simple 1-11 index for this figure's x-axis, not the GIS D2-D12 numbering
def _display_axis():
    return np.arange(1, NUM_REAL_DOMAINS + 1)


# GIS domain numbers D2-D12
def _gis_axis():
    return np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)


# Observed 1967->year change per domain, raw sign ('+' = landward), from the wet/dry table
def load_observed_changes(years):
    if not os.path.isfile(WETDRY_CHANGE_TABLE):
        raise FileNotFoundError(f"Wet/dry change table not found:\n  {WETDRY_CHANGE_TABLE}")
    df = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    gis = _gis_axis()
    out = {}
    for year in years:
        col = f"change_from_wetdry_1967_wetdry_{year}_m"
        if col not in df.columns:
            print(f"  WARNING: no observed column for {year} -- skipped.")
            out[year] = None
            continue
        out[year] = np.array([df[col].get(d, np.nan) for d in gis])
    return out


# One run's shoreline matrix
def _load_shoreline(run_name):
    path = os.path.join(OUTPUT_BASE_DIR, run_name, f"{run_name}_shoreline_matrix.npy")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"shoreline matrix not found:\n  {path}")
    m = np.load(path)
    print(f"  Loaded {run_name}: shape {m.shape}")
    return m


# Flip Barrier3D's landward-positive x_s so '+' = seaward/accretion (README)
def _flip(v):
    return -v


# Model row for a calendar year, or None (with a warning) outside the run
def _year_to_row(year, nt, label=""):
    row = year - START_YEAR
    if not (0 <= row < nt):
        print(f"  WARNING [{label}]: year {year} (row {row}) is outside this run's "
              f"{nt} modeled years -- skipped.")
        return None
    return row


# Updrift and downdrift halves shaded either side of the groin
def _updrift_downdrift_shading(ax):
    display = _display_axis()
    ax.axvspan(display[0] - 0.5, GROIN_BOUNDARY_DISPLAY, alpha=0.06, color="firebrick", zorder=0)
    ax.axvspan(GROIN_BOUNDARY_DISPLAY, display[-1] + 0.5, alpha=0.06, color="seagreen", zorder=0)


# The groin line and its label, placed from the final y-limits
def _mark_groin(ax, color):
    ax.axvline(GROIN_BOUNDARY_DISPLAY, color=color, lw=1.5, ls="--", alpha=0.9, zorder=5)
    yl = ax.get_ylim()
    y = yl[0] + GROIN_LABEL_Y_FRACTION * (yl[1] - yl[0])
    ax.text(GROIN_BOUNDARY_DISPLAY + GROIN_LABEL_X_OFFSET, y, "Buxton\ngroin", color=color,
            fontsize=GROIN_LABEL_FONTSIZE, rotation=90, va="top", ha="left", alpha=0.9)


# Static figure: 1967 shoreline, 2017 modelled and 2018 observed, one panel per run
def fig_three_run_comparison(show_title=True):
    display = _display_axis()

    panels = [
        ("No groin, no nourishment", RUN_BASELINE, BASELINE_COLOR),
        ("Nourishment only", RUN_NOURISHMENT_ONLY, NOURISH_COLOR),
        ("Nourishment + groin (best fit)", RUN_FULL_MODEL, FULL_COLOR),
    ]

    observed_changes = load_observed_changes([OBSERVED_FINAL_YEAR])

    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5), sharey=True, constrained_layout=True)

    for ax, (title, run_name, color) in zip(axes, panels):
        m = _load_shoreline(run_name)
        nt = m.shape[0]

        # Flipped, so the observed position is the 1967 planform MINUS the raw change (README)
        raw_1967 = _flip(m[0][START_REAL_INDEX:END_REAL_INDEX])
        ref_mean = np.nanmean(raw_1967)
        planform_1967 = raw_1967 - ref_mean

        ax.plot(display, planform_1967, color="0.5", ls="--", lw=1.8, marker="o", ms=4,
                label=f"{START_YEAR} shoreline (observed)", zorder=3)

        # Modeled final (2017)
        row_final = _year_to_row(MODEL_FINAL_YEAR, nt, label=title)
        if row_final is not None:
            modeled_final = _flip(m[row_final][START_REAL_INDEX:END_REAL_INDEX]) - ref_mean
            ax.plot(display, modeled_final, color=color, ls="-", lw=2.4, marker="D", ms=6,
                    label=f"{MODEL_FINAL_YEAR} modeled", zorder=5)

        # Observed final (2018): the model's 1967 reference minus the observed change
        obs_change_final = observed_changes.get(OBSERVED_FINAL_YEAR)
        if obs_change_final is not None:
            observed_final = planform_1967 - obs_change_final
            ax.plot(display, observed_final, color="black", ls="--", lw=2.2, marker="s", ms=6,
                    label=f"{OBSERVED_FINAL_YEAR} observed (target)", zorder=6)

        _updrift_downdrift_shading(ax)
        ax.set_xticks(np.arange(1, NUM_REAL_DOMAINS + 1, DOMAIN_TICK_STEP))
        ax.set_xlabel("Model Domain ID (1–11)", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
        ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
        if show_title:
            ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=LEGEND_FONTSIZE, loc="best")

    axes[0].set_ylabel("Shoreline position (m)",
                       fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    # Inverted axis: erosion (negative) plots above zero (README)
    for ax in axes:
        ax.set_ylim(ax.get_ylim()[::-1])

    # Groin marked only after the shared y-limits are final (README)
    for ax, (title, run_name, color) in zip(axes, panels):
        _mark_groin(ax, color)

    if show_title:
        fig.suptitle("Effect of Nourishment and Groin on Modeled Shoreline Position",
                     fontsize=SUPTITLE_FONTSIZE, fontweight="bold")

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    tag = "HAT_three_run_comparison" if show_title else "HAT_three_run_comparison_v2"
    fig_out = os.path.join(OUTPUT_DIR, f"{tag}.png")
    fig.savefig(fig_out, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {fig_out}")
    return fig


# Animated version: the three panels evolve year by year, references held fixed
def make_three_run_evolution_gif(show_title=True):
    display = _display_axis()
    panels = [
        ("No groin, no nourishment", RUN_BASELINE, BASELINE_COLOR),
        ("Nourishment only", RUN_NOURISHMENT_ONLY, NOURISH_COLOR),
        ("Nourishment + groin (best fit)", RUN_FULL_MODEL, FULL_COLOR),
    ]
    observed_changes = load_observed_changes([OBSERVED_FINAL_YEAR])
    obs_change_final = observed_changes.get(OBSERVED_FINAL_YEAR)

    loaded = {}
    for title, run_name, color in panels:
        m = _load_shoreline(run_name)
        raw_1967 = _flip(m[0][START_REAL_INDEX:END_REAL_INDEX])
        ref_mean = np.nanmean(raw_1967)
        planform_1967 = raw_1967 - ref_mean
        nt = m.shape[0]
        traj = np.array([_flip(m[t][START_REAL_INDEX:END_REAL_INDEX]) - ref_mean
                          for t in range(nt)])
        observed_final = (planform_1967 - obs_change_final
                           if obs_change_final is not None else None)
        loaded[title] = dict(color=color, planform_1967=planform_1967,
                              traj=traj, observed_final=observed_final, nt=nt)

    nt_common = min(d["nt"] for d in loaded.values())
    if len(set(d["nt"] for d in loaded.values())) > 1:
        print(f"  WARNING: runs have different lengths {[d['nt'] for d in loaded.values()]} "
              f"-- animating only the first {nt_common} shared years.")
    years = START_YEAR + np.arange(nt_common)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6.5), sharey=True, constrained_layout=True)

    lines = {}
    for ax, (title, run_name, color) in zip(axes, panels):
        d = loaded[title]
        ax.plot(display, d["planform_1967"], color="0.5", ls="--", lw=1.8, marker="o", ms=4,
                label=f"{START_YEAR} shoreline (observed)", zorder=3)
        if d["observed_final"] is not None:
            ax.plot(display, d["observed_final"], color="black", ls="--", lw=2.2, marker="s", ms=6,
                    label=f"{OBSERVED_FINAL_YEAR} observed (target)", zorder=6)
        line, = ax.plot(display, d["traj"][0], color=color, ls="-", lw=2.4, marker="D", ms=6,
                         label="Modeled", zorder=5)
        lines[title] = line

        _updrift_downdrift_shading(ax)
        ax.set_xticks(np.arange(1, NUM_REAL_DOMAINS + 1, DOMAIN_TICK_STEP))
        ax.set_xlabel("Model Domain ID (1–11)", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
        ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
        if show_title:
            ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=LEGEND_FONTSIZE, loc="best")

    axes[0].set_ylabel("Shoreline position (m)",
                       fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

    # One y-range for every panel and year, set before inverting (README)
    all_vals = []
    for d in loaded.values():
        all_vals.append(d["planform_1967"])
        if d["observed_final"] is not None:
            all_vals.append(d["observed_final"])
        all_vals.append(d["traj"][:nt_common].ravel())
    all_vals = np.concatenate(all_vals)
    pad = 0.08 * (np.nanmax(all_vals) - np.nanmin(all_vals) + 1e-9)
    for ax in axes:
        ax.set_ylim(np.nanmin(all_vals) - pad, np.nanmax(all_vals) + pad)
        ax.set_ylim(ax.get_ylim()[::-1])   # same flip+invert convention as the static figure

    # Groin marked only after the shared y-range is set
    for ax, (title, run_name, color) in zip(axes, panels):
        _mark_groin(ax, color)

    if show_title:
        fig.suptitle("Effect of Nourishment and Groin on Modeled Shoreline Position Over Time",
                     fontsize=SUPTITLE_FONTSIZE, fontweight="bold")

    year_text = fig.text(0.015, 0.02, "", fontsize=YEAR_TEXT_FONTSIZE, fontweight="bold", zorder=10)

    def update(frame):
        for title, run_name, color in panels:
            lines[title].set_ydata(loaded[title]["traj"][frame])
        year_text.set_text(f"Year: {years[frame]}")
        return list(lines.values()) + [year_text]

    anim = animation.FuncAnimation(fig, update, frames=nt_common, blit=False)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    tag = "HAT_three_run_evolution" if show_title else "HAT_three_run_evolution_v2"
    gif_out = os.path.join(OUTPUT_DIR, f"{tag}.gif")
    anim.save(gif_out, writer=animation.PillowWriter(fps=GIF_FPS), dpi=GIF_DPI)
    plt.close(fig)
    print(f"Saved: {gif_out}")
    return anim


# Run: both static figures, then both GIFs, each titled and untitled
def main():
    fig_three_run_comparison(show_title=True)
    fig_three_run_comparison(show_title=False)
    make_three_run_evolution_gif(show_title=True)
    make_three_run_evolution_gif(show_title=False)


if __name__ == "__main__":
    main()
