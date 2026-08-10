"""
HAT_three_run_comparison.py
=============================
Three-panel comparison figure: baseline (no groin, no nourishment) vs
nourishment-only vs full model (nourishment + best-fit groin from the
sensitivity sweep). Each panel shows the 1967 initial shoreline, the 2017
modeled position, and the observed 2017/2018 target -- so you can see, left
to right, how each piece of added functionality moves the model closer to
(or away from) reality. All three panels are treated identically (no
intermediate-year overlay on any of them, for consistency across panels).

X-AXIS: domain positions are shown as a simple 1-11 index for this figure
(NOT the underlying GIS D2-D12 numbering used everywhere else in the
project) -- easier to read for this specific comparison. GROIN_BOUNDARY_GIS
(5.5, the real D5/D6 interface) is converted to this same 1-11 scale
(4.5) so the groin marker lands in the correct spot either way.

SIGN CONVENTION -- corrected after visual review, verified numerically before
applying (see _flip()'s docstring and fig_three_run_comparison()): raw x_s
increases LANDWARD (Barrier3D's native direction). This script flips it so
'+' = seaward/accretion, '-' = landward/erosion -- matching FLIP_SIGN_MODEL
used in every OTHER script in this project, AND matching the standard
erosion=negative / accretion=positive convention already used throughout
this dissertation's own LRR analysis. With that flip in place, reconstructing
an observed position requires SUBTRACTING the observed wet/dry CHANGE from
the 1967 reference (planform_1967 - observed_change), matching
HAT_groin_effect_comparison.py's formula exactly -- an earlier version of
this script skipped the flip and used addition instead, which made erosion
plot as a positive number, backwards from the rest of the project.

Produces TWO versions each run: HAT_three_run_comparison.png (full title +
subtitle) and HAT_three_run_comparison_v2.png (no figure-level title, for
dropping into a paper with its own caption) -- per-panel titles ("No groin,
no nourishment", etc.) are kept in BOTH versions, since with three visually
distinct scenarios some in-figure labeling is still needed for a reader to
tell panels apart even when a caption describes the figure as a whole.

WORKFLOW -- produce the three runs BEFORE running this script:
  1. Baseline (no groin, no nourishment):
       RUN_MATRIX = ["no_groin"]; ENABLE_HISTORICAL_NOURISHMENT = False
       Give this run a DISTINCT RUN_NAME_SUFFIX (e.g. "..._no_BN") --
       nourishment on/off isn't part of the run_name, so without a distinct
       suffix this run and #2 below would silently overwrite each other.
  2. Nourishment only (no groin):
       RUN_MATRIX = ["no_groin"]; ENABLE_HISTORICAL_NOURISHMENT = True
  3. Full model (nourishment + best-fit groin):
       RUN_MATRIX = ["groin"]; ENABLE_HISTORICAL_NOURISHMENT = True
       GROIN_TRAPPING_RATE_M_YR / GROIN_DETERIORATION_FRACTION = your
       sweep's best combination

Set RUN_BASELINE / RUN_NOURISHMENT_ONLY / RUN_FULL_MODEL below to the three
resulting run folder names, then run this script directly. Produces both the
static comparison figures (titled + v2) AND an animated GIF version --
make_three_run_evolution_gif() -- showing all three panels evolving
year-by-year (1967-2017) simultaneously, same layout/colors/conventions as
the static figure, with each panel's own 1967 reference and 2018 observed
target held fixed while its modeled line moves.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# =============================================================================
# CONFIG -- edit these three to your actual run folder names
# =============================================================================
RUN_BASELINE          = "HAT_1967_2018_edge_calibrated_no_BN_no_groin"
RUN_NOURISHMENT_ONLY  = "HAT_1967_2018_edge_calibrated_no_groin"
RUN_FULL_MODEL        = "HAT_1967_2018_edge_calibrated_groin"

PROJECT_BASE_DIR = r"/"
OUTPUT_BASE_DIR  = os.path.join(PROJECT_BASE_DIR, "shoreline_position_output", "raw_runs")

# --- Geometry (must match the runs) ---
NUM_REAL_DOMAINS   = 11
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 2
LAST_FILE_NUMBER   = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1     # 12
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS

START_YEAR          = 1967
MODEL_FINAL_YEAR     = 2017   # model's TRUE final simulated year (row 50 of 51)
OBSERVED_FINAL_YEAR  = 2018   # full D2-D12 coverage; 2017 itself only has 6/11 -- see
                               # earlier scripts for the same substitution

WETDRY_CHANGE_TABLE = os.path.join(
    PROJECT_BASE_DIR, "scripts", "groin_module", "hindcast_groin_test",
    "input_prep", "shoreline_position", "shoreline_position_output",
    "Change_from_wetdry_1967_D2_D12.csv",
)

# --- Styling ---
BASELINE_COLOR = "#D32F2F"   # red -- no groin, no nourishment
NOURISH_COLOR  = "#F9A825"   # amber -- nourishment only
FULL_COLOR     = "#1565C0"   # blue -- nourishment + groin (full model)
# NOTE: this deliberately differs from the orange/red no-groin/groin
# convention used in every OTHER figure in this project -- three categories
# to distinguish here, not a binary, so a red/amber/blue "problem -> partial
# fix -> full fix" progression reads better for THIS figure specifically.
GROIN_BOUNDARY_GIS = 5.5
GROIN_BOUNDARY_DISPLAY = GROIN_BOUNDARY_GIS - FIRST_FILE_NUMBER + 1   # 4.5 on the 1-11 scale
GROIN_LABEL_X_OFFSET = 0.15   # nudges the "Buxton groin" text right of the dashed
                               # line itself, so the text doesn't sit on top of it
GROIN_LABEL_Y_FRACTION = 0.75   # how far up the y-range the label sits (0=bottom,
                                  # 1=top). Lowered from 0.9 -- at the larger legend
                                  # font size needed for legibility, the legend box
                                  # got physically bigger and the label started
                                  # overlapping it at 0.9. Checked directly via
                                  # bounding-box overlap, not just eyeballed: 0.8 was
                                  # the first value that cleared it, 0.75 adds a
                                  # small safety margin for real data with a
                                  # different y-range than the synthetic test data.
DOMAIN_TICK_STEP = 2   # every OTHER domain labeled (1,3,5,7,9,11) -- deliberately not
                         # every tick, to leave room for the larger font sizes below
                         # without crowding (advisor's suggested fix, not just mine)
# No OCEAN_AT_BOTTOM / axis-inversion flag here, unlike the older flipped-
# convention scripts -- see the y-axis section in fig_three_run_comparison()
# for why this script's raw ('+' = landward/erosion) convention needs no
# inversion at all to get the same visual result.

# --- Font sizes -- deliberately large. This figure gets shrunk to fit a
# paper/dissertation page width (roughly 1/3 of its native 20" rendered
# width), so anything sized to "look right" at native resolution reads as
# illegibly tiny once shrunk. Advisor feedback: must be legible at 100%
# scale viewing -- tick labels, axis labels, AND legend. ---
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
OUTPUT_DIR = os.path.join(PROJECT_BASE_DIR, "scripts", "groin_module",
                           "hindcast_groin_test", "comparison", "three_run_comparison")


# =============================================================================
# HELPERS
# =============================================================================
def _gis_axis():
    return np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)


def _display_axis():
    """Simple 1-11 index for this figure's x-axis, instead of the underlying
    GIS D2-D12 numbering used everywhere else in the project."""
    return np.arange(1, NUM_REAL_DOMAINS + 1)


def _mark_groin(ax, color):
    ax.axvline(GROIN_BOUNDARY_DISPLAY, color=color, lw=1.5, ls="--", alpha=0.9, zorder=5)
    yl = ax.get_ylim()
    y = yl[0] + GROIN_LABEL_Y_FRACTION * (yl[1] - yl[0])
    ax.text(GROIN_BOUNDARY_DISPLAY + GROIN_LABEL_X_OFFSET, y, "Buxton\ngroin", color=color,
            fontsize=GROIN_LABEL_FONTSIZE, rotation=90, va="top", ha="left", alpha=0.9)


def _updrift_downdrift_shading(ax):
    display = _display_axis()
    ax.axvspan(display[0] - 0.5, GROIN_BOUNDARY_DISPLAY, alpha=0.06, color="firebrick", zorder=0)
    ax.axvspan(GROIN_BOUNDARY_DISPLAY, display[-1] + 0.5, alpha=0.06, color="seagreen", zorder=0)


def _flip(v):
    """Barrier3D's raw x_s increases LANDWARD (erosion). Flip so '+' =
    seaward/accretion, '-' = landward/erosion -- matches FLIP_SIGN_MODEL=True
    used in every OTHER script in this project, and matches the standard
    DSAS/shoreline-change convention (erosion = negative, accretion =
    positive) already used throughout this dissertation's own LRR analysis.
    Verified numerically before using this (see module docstring)."""
    return -v


def _load_shoreline(run_name):
    path = os.path.join(OUTPUT_BASE_DIR, run_name, f"{run_name}_shoreline_matrix.npy")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"shoreline matrix not found:\n  {path}")
    m = np.load(path)
    print(f"  Loaded {run_name}: shape {m.shape}")
    return m


def _year_to_row(year, nt, label=""):
    row = year - START_YEAR
    if not (0 <= row < nt):
        print(f"  WARNING [{label}]: year {year} (row {row}) is outside this run's "
              f"{nt} modeled years -- skipped.")
        return None
    return row


def load_observed_changes(years):
    """Observed 1967->year change (RAW sign, '+' = landward) for each year,
    from the validated wet/dry change table."""
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


# =============================================================================
# FIGURE
# =============================================================================
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

        # FLIPPED convention: '+' = seaward/accretion, '-' = landward/erosion
        # (see _flip() docstring). Reconstructing the observed position now
        # means SUBTRACTING the raw observed change (opposite of an earlier,
        # incorrect version of this script that didn't flip and added
        # instead) -- verified numerically before making this change.
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

        # Observed final (2018, full coverage) -- reconstructed from the
        # model's own 1967 reference MINUS the observed CHANGE since 1967
        # (raw sign, '+' = landward -- subtracting it from the flipped
        # planform correctly pushes erosion negative).
        obs_change_final = observed_changes.get(OBSERVED_FINAL_YEAR)
        if obs_change_final is not None:
            observed_final = planform_1967 - obs_change_final
            ax.plot(display, observed_final, color="black", ls="--", lw=2.2, marker="s", ms=6,
                    label=f"{OBSERVED_FINAL_YEAR} observed (target)", zorder=6)

        _updrift_downdrift_shading(ax)
        ax.set_xticks(np.arange(1, NUM_REAL_DOMAINS + 1, DOMAIN_TICK_STEP))
        ax.set_xlabel("Model Domain ID (1\u201311)", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
        ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
        if show_title:
            ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=LEGEND_FONTSIZE, loc="best")

    axes[0].set_ylabel("Shoreline position (m)",
                       fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    # Axis inversion restored: with the FLIPPED convention ('+' = seaward/
    # accretion, '-' = landward/erosion), inverting the axis puts erosion
    # (negative values) above the zero line and accretion (positive) below
    # it -- matching both the standard erosion=negative convention AND the
    # visual "erosion above zero" layout, simultaneously.
    for ax in axes:
        ax.set_ylim(ax.get_ylim()[::-1])

    # _mark_groin() is called HERE, only after the y-limits are fully
    # finalized -- NOT inside the loop above. With sharey=True, each panel's
    # ax.get_ylim() keeps changing as LATER panels get plotted (matplotlib
    # synchronizes the shared range progressively), so calling it mid-loop
    # used a stale, too-narrow range and placed the label far from the top.
    # Verified this directly before fixing it (not just guessed).
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


def make_three_run_evolution_gif(show_title=True):
    """Animated version of fig_three_run_comparison(): all three panels
    evolve simultaneously, year-by-year (1967 through each run's true final
    year), each keeping its OWN static 1967 reference and 2018 observed
    target fixed while its modeled line moves -- same layout, colors, and
    sign convention as the static figure, so the two are a direct visual
    match (just one frozen at 2017, the other animating through it)."""
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
        ax.set_xlabel("Model Domain ID (1\u201311)", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
        ax.tick_params(axis="both", labelsize=TICK_LABEL_FONTSIZE)
        if show_title:
            ax.set_title(title, fontsize=TITLE_FONTSIZE, fontweight="bold")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=LEGEND_FONTSIZE, loc="best")

    axes[0].set_ylabel("Shoreline position (m)",
                       fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

    # Fix y-limits across the WHOLE animation (all three panels, all years,
    # both static reference lines) up front, computed BEFORE inverting --
    # same two-step pattern (set normal range, then reverse) already used
    # in the other GIF-making functions in this project, to avoid any
    # order-of-operations confusion with combining a shared range and an
    # inverted axis.
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

    # _mark_groin() must be called AFTER the shared ylim above is set, not
    # inside the loop -- same ordering bug as the static figure had: calling
    # it earlier uses each axis's original auto-scaled range (from its own
    # plot() calls only), not the explicit all-panel range set here.
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


if __name__ == "__main__":
    fig_three_run_comparison(show_title=True)
    fig_three_run_comparison(show_title=False)
    make_three_run_evolution_gif(show_title=True)
    make_three_run_evolution_gif(show_title=False)
