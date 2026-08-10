"""
CoastSat Reporter-Ready Shoreline Change Plot
=============================================
Reads the transect_lrr_results.csv produced by coastsat_transect_lrr_plot.py
and generates a clean, publication/press-ready figure with geographic labelling
instead of transect IDs.

No domain mapping or lookup tables required — just point it at the CSV,
fill in a place name, and run.

Usage
-----
Edit the CONFIG section below, then run:
    python coastsat_reporter_plot.py
"""

# ============================================================
# CONFIG  –  edit these values before running
# ============================================================

# Path to the transect_lrr_results.csv produced by coastsat_transect_lrr_plot.py
INPUT_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\custom\buxton_2000_2025\transect_lrr_results.csv"

# Place name shown in the title (e.g. "Buxton", "Hatteras Village", "Rodanthe")
PLACE_NAME = "Buxton"

# Date range label shown in the title  (just for display — no filtering done here)
PERIOD_LABEL = "2000–2025"

# Approximate real-world length of the transect stretch in kilometres.
# Used to scale the x-axis into geographic distance.
# Rule of thumb: CoastSat transects are spaced ~50 m apart along the NC coast,
# so for N transects: length_km ≈ N * 0.05
# Leave as None to auto-calculate from transect count * 0.05 km.
STRETCH_LENGTH_KM = None   # e.g. 1.4  or  None for auto

# Which end of your transect range is south?
# "first"  → the first row of the CSV is the southernmost transect (default)
# "last"   → the first row is northernmost (flip the axis)
SOUTH_END = "first"

# ---- Colour thresholds (m/yr) ----
# Bars are coloured into four bins based on these cutoffs:
#   strong erosion  : LRR  <  EROSION_STRONG
#   mild erosion    : EROSION_STRONG  <=  LRR  <  0
#   mild accretion  : 0  <=  LRR  <=  ACCRETION_STRONG
#   strong accretion: LRR  >  ACCRETION_STRONG
EROSION_STRONG   = -0.3   # m/yr  (more negative = stronger erosion threshold)
ACCRETION_STRONG = +0.5   # m/yr  (more positive = stronger accretion threshold)

# ---- Colours for each bin ----
COLOR_STRONG_EROSION   = "#c0392b"   # deep red
COLOR_MILD_EROSION     = "#e8a090"   # light salmon
COLOR_MILD_ACCRETION   = "#a8c8e8"   # light blue
COLOR_STRONG_ACCRETION = "#1a6ea8"   # deep blue

# ---- Annotation callouts ----
# The script auto-finds the most eroding and most accreting transect and
# draws a callout arrow for each. Set to False to suppress them.
SHOW_CALLOUTS = True

# ---- Output ----
OUTPUT_DIR  = r"/scripts/input_prep/5-scr/CoastSat\custom_range&dates\buxton_2000_2025"
OUTPUT_NAME = "shoreline_change_SIMPLE.png"   # filename for the saved figure
DPI         = 180


# ============================================================
# IMPORTS
# ============================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")


# ============================================================
# HELPERS
# ============================================================

def load_results(csv_path: str) -> pd.DataFrame:
    """Load and validate the transect LRR CSV."""
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"Input CSV not found:\n  {csv_path}")
    df = pd.read_csv(csv_path)
    if "lrr_m_yr" not in df.columns:
        raise ValueError("CSV does not contain an 'lrr_m_yr' column. "
                         "Check that INPUT_CSV points to a transect_lrr_results.csv "
                         "from coastsat_transect_lrr_plot.py.")
    df_valid = df.dropna(subset=["lrr_m_yr"]).reset_index(drop=True)
    n_skipped = len(df) - len(df_valid)
    if n_skipped:
        print(f"  ⚠️  {n_skipped} transect(s) had no valid LRR and were skipped.")
    return df_valid


def assign_colors(lrr_series: pd.Series) -> list:
    """Map LRR values to the four-bin colour scheme."""
    colors = []
    for v in lrr_series:
        if v < EROSION_STRONG:
            colors.append(COLOR_STRONG_EROSION)
        elif v < 0:
            colors.append(COLOR_MILD_EROSION)
        elif v <= ACCRETION_STRONG:
            colors.append(COLOR_MILD_ACCRETION)
        else:
            colors.append(COLOR_STRONG_ACCRETION)
    return colors


def build_xticks(dist_km: np.ndarray, n_ticks: int = 5) -> tuple:
    """Generate evenly spaced geographic x-tick positions and labels."""
    positions = np.linspace(dist_km[0], dist_km[-1], n_ticks)
    labels = []
    for i, p in enumerate(positions):
        if i == 0:
            labels.append("South\nend" if SOUTH_END == "first" else "North\nend")
        elif i == n_ticks - 1:
            labels.append("North\nend" if SOUTH_END == "first" else "South\nend")
        else:
            labels.append(f"↑\n~{p:.2f} km")
    return positions, labels


# ============================================================
# MAIN PLOT FUNCTION
# ============================================================

def make_reporter_plot(df: pd.DataFrame) -> None:

    lrr    = df["lrr_m_yr"].values
    n      = len(lrr)

    if SOUTH_END == "last":
        lrr = lrr[::-1]

    # ---- X-axis as geographic distance ----
    stretch = STRETCH_LENGTH_KM if STRETCH_LENGTH_KM else round(n * 0.05, 2)
    dist_km = np.linspace(0, stretch, n)
    bar_width = stretch / n * 0.85

    # ---- Colours ----
    colors = assign_colors(lrr)

    # ---- Figure ----
    fig, ax = plt.subplots(figsize=(11, 6))
    fig.patch.set_facecolor("white")

    ax.bar(dist_km, lrr, width=bar_width,
           color=colors, edgecolor="white", linewidth=0.4)
    ax.axhline(0, color="black", lw=1.2)


    # ---- X-axis: geographic distance ----
    tick_pos, tick_labels = build_xticks(dist_km)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=9)
    ax.tick_params(axis="x", length=0)

    # ---- Legend ----
    legend_elements = [
        mpatches.Patch(facecolor=COLOR_STRONG_EROSION,
                       label=f"Strong erosion  (< {EROSION_STRONG:+.1f} m/yr)"),
        mpatches.Patch(facecolor=COLOR_MILD_EROSION,
                       label=f"Mild erosion  ({EROSION_STRONG:+.1f} to 0 m/yr)"),
        mpatches.Patch(facecolor=COLOR_MILD_ACCRETION,
                       label=f"Mild accretion  (0 to {ACCRETION_STRONG:+.1f} m/yr)"),
        mpatches.Patch(facecolor=COLOR_STRONG_ACCRETION,
                       label=f"Strong accretion  (> {ACCRETION_STRONG:+.1f} m/yr)"),
    ]
    ax.legend(handles=legend_elements, loc="upper left",
              fontsize=8.5, framealpha=0.9, edgecolor="lightgray")

    # ---- Labels and title ----
    ax.set_ylabel("Shoreline change rate (m/yr)", fontsize=11)
    ax.set_xlabel(f"Location along {PLACE_NAME} shoreline  (south → north)",
                  fontsize=11)
    ax.set_title(f"{PLACE_NAME} Shoreline Change Rate  |  {PERIOD_LABEL}",
                 fontsize=14, fontweight="bold", pad=14)

    fig.text(
        0.5, 0.92,
        "Satellite-derived rates from CoastSat  •  "
        "Negative = shoreline moving landward  •  "
        "Positive = shoreline moving seaward",
        ha="center", fontsize=8, color="gray"
    )

    # ---- Summary stats text box ----
    n_eroding   = (lrr < 0).sum()
    n_accreting = (lrr >= 0).sum()
    mean_lrr    = lrr.mean()
    stats_text  = (f"n = {n}  |  mean = {mean_lrr:+.2f} m/yr\n"
                   f"eroding: {n_eroding}  |  accreting: {n_accreting}")
    ax.text(0.99, 0.97, stats_text, transform=ax.transAxes,
            va="top", ha="right", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.75,
                      edgecolor="lightgray"))

    # ---- Axis formatting ----
    pad = stretch * 0.03
    ax.set_xlim(-pad, stretch + pad)
    y_margin = max(abs(lrr.min()), abs(lrr.max())) * 0.25
    ax.set_ylim(lrr.min() - y_margin, lrr.max() + y_margin)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{v:+.1f}")
    )
    ax.grid(True, axis="y", alpha=0.25, lw=0.7)
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.92])

    # ---- Save ----
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    out_path = os.path.join(OUTPUT_DIR, OUTPUT_NAME)
    plt.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.show()
    print(f"\nFigure saved: {out_path}")


# ============================================================
# MAIN
# ============================================================

def main():
    print(f"Loading: {INPUT_CSV}")
    df = load_results(INPUT_CSV)
    print(f"  {len(df)} transects with valid LRR\n")

    print(f"Building reporter plot for {PLACE_NAME}, {PERIOD_LABEL}...")
    make_reporter_plot(df)


if __name__ == "__main__":
    main()
