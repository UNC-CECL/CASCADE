#!/usr/bin/env python3
"""
HATTERAS ISLAND: Management Investigation — Plot from Saved Runs
Natural vs Roadway Management vs Roadway Management + Historical Beach Nourishment

This script loads previously saved CASCADE NPZ files and produces the
management comparison plots without re-running any simulations.

USAGE
-----
1. Run your three scenarios separately using HAT_hindcast_1984_2024_old version.py
   (or any other script) with the appropriate management flags.
2. Fill in the RUN_PATHS dict below with the paths to each saved run folder.
3. Run this script -- it loads the cascade objects and plots.

Each RUN_PATHS entry:
    "Label shown on plot": r"C:/path/to/saved/run/folder"

The folder must contain a cascade.npz file (written by cascade.save()).
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

# =============================================================================
# SECTION 1: CONFIGURE PATHS AND LABELS  ← edit this section
# =============================================================================

# Output folder for plots
OUTPUT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\figure_making\management_investigation\output"

# Path to the CoastSat transect CSV for the active period.
# Columns expected: domain_number (GIS 1–90), lrr_m_yr (m/yr per transect).
COASTSAT_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\transect_lrr_full.csv"
COASTSAT_LABEL = "CoastSat LRR (1984–2004)"   # label shown in legend

# Period covered by these runs — used for axis titles and comparison filenames.
START_YEAR = 1984
END_YEAR   = 2004

# ── RUN PATHS ────────────────────────────────────────────────────────────────
# Keys   = labels shown on the plot (keep them short)
# Values = path to the saved run folder (must contain cascade.npz)
#
# Order determines plotting order (first = bottom of legend).
# Add or remove entries freely — the script handles 2, 3, or more scenarios.

RUN_PATHS = {
    "Natural": (
        r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs"
        r"\HAT_1984_2004_natural_Hs2p5"
    ),
    "Roadway": (
        r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs"
        r"\HAT_1984_2004_roadway_Hs2p5"
    ),
    "Roadway + Historical BN": (
        r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs"
        r"\HAT_1984_2004_roadway_histBN_Hs2p5"
    ),
}

# Historical BN volume schedule — used for the BN bar panel on yearly plots.
# Only needed for Period 2 (2004–2024). For Period 1 leave as empty dict {}.
# Format: {calendar_year: [volume_m3 per padded domain index, length=120]}
# If you leave this empty the bar panel will show zeros (correct for Period 1).
HIST_BN_VOL_BY_YEAR = {}   # populated automatically below if ENABLE_BN = True

# Set True for Period 2 to inject BN volumes into the bar panel.
ENABLE_BN = False   # ← change to True for 2004–2024

# =============================================================================
# SECTION 2: DOMAIN CONSTANTS  (match your Hatteras CASCADE setup exactly)
# =============================================================================

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 1
LAST_FILE_NUMBER   = 90
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120
START_REAL_INDEX   = NUM_BUFFER_DOMAINS        # 15
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS   # 105

DOMAIN_LENGTH_M  = 500.0
DOMAIN_SPACING_M = 500.0
DOMAIN_TICK_STEP = 5

FLIP_SIGN_MODEL  = True   # True → positive = accretion (matches hindcast convention)
TO_METERS        = True

MAKE_YEARLY_GIF      = True
GIF_DURATION_SECONDS = 4

# =============================================================================
# SECTION 3: HISTORICAL BN SCHEDULE (Period 2 only)
# =============================================================================

_CY_TO_M3 = 0.764555

HAT_BN_YEARS = [2014, 2022]

HAT_BN_VOLUME_BY_DOMAIN = {
    # Buxton: GIS 6–15 (10 domains, 1,200,000 cy in 2022)
     6: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
     7: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
     8: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
     9: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    10: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    11: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    12: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    13: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    14: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    15: [0, round(1_200_000 / 10 * _CY_TO_M3, 1)],
    # Avon: GIS 23–26 (4 domains, 2,200,000 cy in 2022)
    23: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    24: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    25: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    26: [0, round(2_200_000 / 4 * _CY_TO_M3, 1)],
    # Rodanthe: GIS 85–88 (4 domains, 1,620,000 cy in 2014)
    85: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
    86: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
    87: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
    88: [round(1_620_000 / 4 * _CY_TO_M3, 1), 0],
}

# =============================================================================
# SECTION 4: STYLING
# =============================================================================

# Add or change colors here — keys must match RUN_PATHS keys exactly.
SCENARIO_COLORS = {
    "Natural":                 "#2ca02c",   # green
    "Roadway":                 "#1f77b4",   # blue
    "Roadway + Historical BN": "#d62728",   # red
}
SCENARIO_LW = 2.2   # line width for all scenario lines

COASTSAT_COLOR = "#FF8C00"   # orange — matches hindcast script
COASTSAT_LW    = 2.0

# Geographic annotation colors
COMM_COLORS = {
    "Cape Point":  "#d6e4f0",
    "Buxton":      "#d5e8d4",
    "Avon":        "#d5e8d4",
    "Tri-Village": "#d5e8d4",
}

# =============================================================================
# STARTUP
# =============================================================================

os.makedirs(OUTPUT_DIR, exist_ok=True)
print("=" * 80)
print("HATTERAS ISLAND — MANAGEMENT INVESTIGATION — PLOT FROM SAVED RUNS")
print("=" * 80)
print(f"Period:     {START_YEAR}–{END_YEAR}")
print(f"Output dir: {OUTPUT_DIR}")
print(f"Runs to load: {list(RUN_PATHS.keys())}")
print("=" * 80 + "\n")


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _gis_to_pad(gis_id):
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


def load_cascade(run_folder):
    """
    Load a saved CASCADE object from a run folder.

    CASCADE's save() method writes a cascade.npz file containing the full
    pickled Cascade object under the key 'cascade.npy'.

    Parameters
    ----------
    run_folder : str
        Path to the saved run folder (the argument passed to cascade.save()).

    Returns
    -------
    cascade object or None if loading fails
    """
    npz_path = os.path.join(run_folder, "cascade.npz")
    if not os.path.exists(npz_path):
        print(f"  ❌ cascade.npz not found in: {run_folder}")
        return None
    try:
        data = np.load(npz_path, allow_pickle=True)
        cascade = data["cascade.npy"].item()
        print(f"  ✓ Loaded: {os.path.basename(run_folder)}")
        return cascade
    except Exception as e:
        print(f"  ❌ Failed to load {npz_path}: {e}")
        return None


def get_x_s_TS(b3d):
    for attr in ("x_s_TS", "_x_s_TS"):
        if hasattr(b3d, attr):
            return np.asarray(getattr(b3d, attr), dtype=float)
    raise AttributeError("No shoreline time series found on Barrier3D object.")


def build_shoreline_matrix(cascade, to_meters=True):
    b3d_list = cascade.barrier3d
    nt = len(get_x_s_TS(b3d_list[0]))
    shoreline = np.zeros((nt, TOTAL_DOMAINS), dtype=float)
    for j, b3d in enumerate(b3d_list):
        shoreline[:, j] = get_x_s_TS(b3d)
    if to_meters:
        shoreline *= 10.0
    return shoreline


def build_relative_shoreline_change_matrix(cascade, to_meters=True, flip_sign=True):
    sm = build_shoreline_matrix(cascade, to_meters=to_meters)
    sc = sm - sm[0, :]
    if flip_sign:
        sc *= -1.0
    return sc


def build_bn_arrays():
    """Build per-year BN volume arrays for the bar panel (Period 2 only)."""
    vol_by_year = {yr: np.zeros(TOTAL_DOMAINS)
                   for yr in range(START_YEAR, END_YEAR + 1)}
    if not ENABLE_BN:
        return vol_by_year
    for gis_id, volumes_m3 in HAT_BN_VOLUME_BY_DOMAIN.items():
        pad = _gis_to_pad(gis_id)
        if not (0 <= pad < TOTAL_DOMAINS):
            continue
        for year, total_m3 in zip(HAT_BN_YEARS, volumes_m3):
            if total_m3 == 0 or year < START_YEAR or year > END_YEAR:
                continue
            vol_by_year[year][pad] = float(total_m3)   # m³ total per domain
    return vol_by_year


def load_coastsat():
    """Load CoastSat transect CSV and aggregate to per-domain means."""
    if not COASTSAT_CSV or not os.path.exists(COASTSAT_CSV):
        print(f"  ⚠️  CoastSat CSV not found: {COASTSAT_CSV}")
        return None, None
    try:
        cs_df = pd.read_csv(COASTSAT_CSV)
        cs_df["domain_number"] = pd.to_numeric(cs_df["domain_number"], errors="coerce")
        cs_df["lrr_m_yr"]      = pd.to_numeric(cs_df["lrr_m_yr"],      errors="coerce")
        cs_df = cs_df.dropna(subset=["domain_number", "lrr_m_yr"])
        cs_df = cs_df[(cs_df["domain_number"] >= FIRST_FILE_NUMBER) &
                      (cs_df["domain_number"] <= LAST_FILE_NUMBER)]
        means = (cs_df.groupby("domain_number")["lrr_m_yr"]
                   .mean().reset_index())
        means["domain_number"] = means["domain_number"].astype(int)
        x    = np.array([_gis_to_pad(g) for g in means["domain_number"]])
        rate = means["lrr_m_yr"].values
        print(f"  ✓ CoastSat: {len(rate)} domain means  "
              f"(rate {rate.min():.2f}–{rate.max():.2f} m/yr)")
        return x, rate
    except Exception as e:
        print(f"  ⚠️  CoastSat load failed: {e}")
        return None, None


# =============================================================================
# ANNOTATION HELPERS
# =============================================================================

ANN_COMMUNITY_SPANS = [
    (_gis_to_pad(1),  _gis_to_pad(6),  "Cape Point"),
    (_gis_to_pad(7),  _gis_to_pad(8),  "Buxton"),
    (_gis_to_pad(21), _gis_to_pad(31), "Avon"),
    (_gis_to_pad(68), _gis_to_pad(83), "Tri-Village"),
]
ANN_VILLAGE_LINES = {"Salvo": _gis_to_pad(69), "Waves": _gis_to_pad(74),
                     "Rodanthe": _gis_to_pad(80)}
ANN_PIERS  = {"Avon Pier": _gis_to_pad(26), "Rodanthe Pier": _gis_to_pad(79)}
ANN_GROIN  = _gis_to_pad(6)
ANN_WIMBLE = (_gis_to_pad(60), _gis_to_pad(74))


def add_geographic_annotations(ax):
    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin

    for pad_lo, pad_hi, label in ANN_COMMUNITY_SPANS:
        x0 = FIRST_FILE_NUMBER + (pad_lo - START_REAL_INDEX)
        x1 = FIRST_FILE_NUMBER + (pad_hi - START_REAL_INDEX)
        ax.axvspan(x0, x1, alpha=0.12,
                   color=COMM_COLORS.get(label, "#eeeeee"), zorder=0)
        ax.text((x0 + x1) / 2, ymax - 0.04 * yrange, label,
                fontsize=7.5, color="0.35", ha="center", va="top", zorder=1)

    wx0 = FIRST_FILE_NUMBER + (ANN_WIMBLE[0] - START_REAL_INDEX)
    wx1 = FIRST_FILE_NUMBER + (ANN_WIMBLE[1] - START_REAL_INDEX)
    ax.axvspan(wx0, wx1, alpha=0.10, color="#f5deb3", zorder=0)

    for vname, pad_idx in ANN_VILLAGE_LINES.items():
        x = FIRST_FILE_NUMBER + (pad_idx - START_REAL_INDEX)
        ax.axvline(x, color="0.40", lw=0.9, ls="--", alpha=0.65, zorder=1)
        ax.text(x, ymax - 0.10 * yrange, vname, fontsize=7, color="0.30",
                ha="center", va="top", rotation=90, zorder=2)

    for _, pad_idx in ANN_PIERS.items():
        x = FIRST_FILE_NUMBER + (pad_idx - START_REAL_INDEX)
        ax.axvline(x, color="#2166ac", lw=0.9, ls="-.", alpha=0.70, zorder=1)

    gx = FIRST_FILE_NUMBER + (ANN_GROIN - START_REAL_INDEX)
    ax.axvline(gx, color="#b22222", lw=0.9, ls=":", alpha=0.70, zorder=1)

    xa = LAST_FILE_NUMBER + 0.5
    ax.annotate("Accretion ▲", xy=(xa, ymax - 0.18 * yrange),
                xytext=(xa, ymax - 0.32 * yrange), fontsize=8.5, color="0.30",
                ha="center", va="center",
                arrowprops=dict(arrowstyle="-|>", color="0.40", lw=1.0))
    ax.annotate("Erosion ▼", xy=(xa, ymin + 0.18 * yrange),
                xytext=(xa, ymin + 0.32 * yrange), fontsize=8.5, color="0.30",
                ha="center", va="center",
                arrowprops=dict(arrowstyle="-|>", color="0.40", lw=1.0))


def configure_gis_xaxis(ax):
    ax.set_xlim(FIRST_FILE_NUMBER - 0.5, LAST_FILE_NUMBER + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(DOMAIN_TICK_STEP))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="x", which="major", labelsize=9, direction="in", length=5)
    ax.tick_params(axis="x", which="minor", direction="in", length=3)
    ax.set_xlabel("GIS Domain (500 m alongshore, S→N)",
                  fontsize=11, fontweight="bold")


# =============================================================================
# PLOT 1: SHORELINE CHANGE RATES (end-of-period summary)
# =============================================================================

def plot_shoreline_change_rates(rate_profiles, coastsat_x, coastsat_rate):
    gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)

    fig, ax = plt.subplots(figsize=(20, 6), constrained_layout=True)

    for label, rate in rate_profiles.items():
        real_rate = rate[START_REAL_INDEX:END_REAL_INDEX]
        ax.plot(gis_ids, real_rate,
                color=SCENARIO_COLORS.get(label, "gray"),
                lw=SCENARIO_LW, label=label, zorder=5)

    if coastsat_x is not None and coastsat_rate is not None:
        cs_gis = coastsat_x - START_REAL_INDEX + FIRST_FILE_NUMBER
        ax.plot(cs_gis, coastsat_rate,
                color=COASTSAT_COLOR, lw=COASTSAT_LW,
                label=COASTSAT_LABEL, zorder=6)

    ax.axhline(0, color="#2c2c2c", lw=1.2, ls="--", alpha=0.65, zorder=3)

    configure_gis_xaxis(ax)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="y", labelsize=10, direction="in", length=5)
    ax.grid(True, which="major", ls=":", lw=0.6, alpha=0.4, color="gray")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Shoreline Change Rate (m/yr)", fontsize=11, fontweight="bold")
    ax.set_title(
        f"Hatteras Island — Management Scenarios | {START_YEAR}–{END_YEAR}",
        fontsize=13, fontweight="bold", pad=10)

    add_geographic_annotations(ax)
    ax.legend(loc="lower left", fontsize=9, frameon=True, framealpha=0.92)

    out = os.path.join(OUTPUT_DIR,
                       f"HAT_{START_YEAR}_{END_YEAR}_management_shoreline_rates.png")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"  ✓ Saved: {out}")
    plt.show()


# =============================================================================
# PLOT 2: YEARLY RELATIVE SHORELINE CHANGE  (+ optional GIF)
# =============================================================================

def plot_yearly_relative_shoreline(sc_by_label, bn_vol_by_year):
    yearly_dir = os.path.join(OUTPUT_DIR,
                              f"HAT_{START_YEAR}_{END_YEAR}_yearly_frames")
    os.makedirs(yearly_dir, exist_ok=True)

    gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)
    max_nt  = min(sc.shape[0] for sc in sc_by_label.values())

    # Fixed y-limits across all frames
    all_vals = np.concatenate([sc[:max_nt, START_REAL_INDEX:END_REAL_INDEX].ravel()
                               for sc in sc_by_label.values()])
    all_vals = all_vals[np.isfinite(all_vals)]
    y_pad = 0.15 * (all_vals.max() - all_vals.min()) if all_vals.size else 5.0
    y_min = (all_vals.min() - y_pad) if all_vals.size else -10.0
    y_max = (all_vals.max() + y_pad) if all_vals.size else  10.0

    # BN bar panel scale
    all_bn_m3 = [v for yr_arr in bn_vol_by_year.values()
                 for v in yr_arr if v > 0]
    max_bn = max(all_bn_m3) if all_bn_m3 else 1.0

    # Which label to use for BN dot placement
    bn_label = ("Roadway + Historical BN"
                if "Roadway + Historical BN" in sc_by_label
                else list(sc_by_label.keys())[-1])

    png_files = []

    for t in range(max_nt):
        model_year = START_YEAR + t
        bn_m3 = bn_vol_by_year.get(model_year, np.zeros(TOTAL_DOMAINS))
        bn_m3 = np.asarray(bn_m3, dtype=float)
        active_bn = np.where(bn_m3 > 0)[0]

        fig, (ax, ax_bn) = plt.subplots(
            2, 1, figsize=(22, 9), sharex=False,
            gridspec_kw={"height_ratios": [3.2, 1.15]},
            constrained_layout=True)

        for label, sc in sc_by_label.items():
            real_sc = sc[t, START_REAL_INDEX:END_REAL_INDEX]
            ax.plot(gis_ids, real_sc,
                    color=SCENARIO_COLORS.get(label, "gray"),
                    lw=SCENARIO_LW, label=label)

        ax.axhline(0, color="#2c2c2c", lw=1.2, ls="--", alpha=0.65)

        # Mark BN domains on the shoreline panel
        if len(active_bn) > 0:
            bn_gis = [FIRST_FILE_NUMBER + (p - START_REAL_INDEX)
                      for p in active_bn if START_REAL_INDEX <= p < END_REAL_INDEX]
            bn_y   = [sc_by_label[bn_label][t, p]
                      for p in active_bn if START_REAL_INDEX <= p < END_REAL_INDEX]
            ax.scatter(bn_gis, bn_y, s=40, marker="o",
                       color="red", edgecolor="black", lw=0.6, zorder=10,
                       label="BN domain")
            for gx in bn_gis:
                ax.axvline(gx, ls=":", lw=0.9, color="red", alpha=0.40, zorder=2)

        configure_gis_xaxis(ax)
        ax.set_ylim(y_min, y_max)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax.tick_params(axis="y", labelsize=10, direction="in", length=5)
        ax.grid(True, which="major", ls=":", lw=0.6, alpha=0.35, color="gray")
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylabel(f"Relative shoreline change from {START_YEAR} (m)\n"
                      "Accretion +, erosion −",
                      fontsize=10, fontweight="bold")
        ax.set_title(f"Hatteras Island — Management Comparison | {model_year}",
                     fontsize=13, fontweight="bold", pad=10)
        ax.legend(loc="lower left", fontsize=9, frameon=True, framealpha=0.92)
        add_geographic_annotations(ax)

        # BN bar panel
        bn_gis_all   = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)
        bn_vols_real = bn_m3[START_REAL_INDEX:END_REAL_INDEX]
        ax_bn.bar(bn_gis_all, bn_vols_real, width=0.85,
                  color="tab:blue", alpha=0.85)

        bn_title = (f"Historical BN in {model_year}: YES | "
                    f"Total = {bn_m3.sum():,.0f} m³"
                    if bn_m3.sum() > 0 else
                    f"Historical BN in {model_year}: None")
        ax_bn.set_title(bn_title, fontsize=11, fontweight="bold")
        ax_bn.set_ylabel("BN volume\n(m³/domain)", fontsize=9, fontweight="bold")
        ax_bn.set_ylim(0, max_bn * 1.25)
        ax_bn.grid(True, axis="y", ls=":", lw=0.6, alpha=0.35)
        configure_gis_xaxis(ax_bn)

        for p in active_bn:
            if START_REAL_INDEX <= p < END_REAL_INDEX:
                gx  = FIRST_FILE_NUMBER + (p - START_REAL_INDEX)
                vol = bn_m3[p]
                ax_bn.text(gx, vol + max_bn * 0.025,
                           f"GIS {gx}\n{vol/1000:.0f}k",
                           ha="center", va="bottom", fontsize=7.5, color="0.25")

        fig_out = os.path.join(yearly_dir,
                               f"HAT_{START_YEAR}_{END_YEAR}_{model_year}.png")
        fig.savefig(fig_out, dpi=200, bbox_inches="tight")
        plt.close(fig)
        png_files.append(fig_out)

    print(f"  ✓ Saved {len(png_files)} yearly frames to: {yearly_dir}")

    if MAKE_YEARLY_GIF and png_files:
        try:
            import imageio.v2 as imageio
            images  = [imageio.imread(f) for f in png_files]
            gif_out = os.path.join(OUTPUT_DIR,
                                   f"HAT_{START_YEAR}_{END_YEAR}_yearly.gif")
            imageio.mimsave(gif_out, images,
                            duration=GIF_DURATION_SECONDS, loop=0)
            print(f"  ✓ Saved GIF: {gif_out}")
        except ImportError:
            print("  imageio not installed — GIF skipped (pip install imageio)")

    return png_files


# =============================================================================
# MAIN
# =============================================================================

def main():
    # ── Load cascade objects ──────────────────────────────────────────────────
    print("Loading CASCADE runs...")
    cascades = {}
    for label, folder in RUN_PATHS.items():
        cascade = load_cascade(folder)
        if cascade is None:
            print(f"  ❌ Could not load '{label}' — skipping")
            continue
        cascades[label] = cascade

    if not cascades:
        print("❌ No runs loaded successfully. Check RUN_PATHS.")
        sys.exit(1)

    print(f"\n✓ Loaded {len(cascades)} run(s): {list(cascades.keys())}\n")

    # ── Load CoastSat ─────────────────────────────────────────────────────────
    print("Loading CoastSat data...")
    coastsat_x, coastsat_rate = load_coastsat()
    print()

    # ── Build BN volume arrays for bar panel ─────────────────────────────────
    bn_vol_by_year = build_bn_arrays()

    # ── Compute shoreline change rate profiles ────────────────────────────────
    print("Computing shoreline change rates...")
    rate_profiles   = {}
    sc_by_label     = {}
    time_span = float(END_YEAR - START_YEAR)

    for label, cascade in cascades.items():
        try:
            sm        = build_shoreline_matrix(cascade, to_meters=TO_METERS)
            total_chg = sm[-1, :] - sm[0, :]
            rate      = total_chg / time_span
            if FLIP_SIGN_MODEL:
                rate *= -1.0
            rate_profiles[label] = rate

            sc = build_relative_shoreline_change_matrix(
                cascade, to_meters=True, flip_sign=FLIP_SIGN_MODEL)
            sc_by_label[label] = sc

            real_rate = rate[START_REAL_INDEX:END_REAL_INDEX]
            print(f"  {label}: rate range "
                  f"{real_rate[np.isfinite(real_rate)].min():.2f}–"
                  f"{real_rate[np.isfinite(real_rate)].max():.2f} m/yr")
        except Exception as e:
            print(f"  ❌ Could not compute rates for '{label}': {e}")

    print()

    # ── Plot 1: Shoreline change rates ────────────────────────────────────────
    print("Generating shoreline rate plot...")
    plot_shoreline_change_rates(rate_profiles, coastsat_x, coastsat_rate)

    # ── Plot 2: Yearly relative shoreline change ──────────────────────────────
    print("Generating yearly frames...")
    plot_yearly_relative_shoreline(sc_by_label, bn_vol_by_year)

    print("\n✓ Management investigation plots complete.")
    print(f"  All outputs saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
