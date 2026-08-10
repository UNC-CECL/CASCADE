"""
HAT_plot_groin_runs.py
======================
Plots for the 1967-1997 groin-test runs. Loads saved shoreline matrices (the
*_shoreline_matrix.npy files written by the base / groin runners) and produces:

  FIG 1  Shoreline CHANGE (m, end - start) per domain          [position change]
  FIG 2  Shoreline change RATE (m/yr) per domain               [LRR-style]
  FIG 3  Shoreline POSITION over time, selected domains        [trajectories]
  FIG 4  (multi-run only) run-vs-run DIFFERENCE, e.g.
         groin_be minus no_groin -> the isolated groin signal

Matches your main hindcast's conventions: model orange (#FF8C00), groin red
(#B71C1C) at the D5/D6 boundary, GIS domain x-axis, real-domain focus with
buffers shaded. No CoastSat dependency -- this is for inspecting the model runs
themselves. Point RUNS at one or more saved run folders.

Usage
-----
Edit RUNS below to list the run folder name(s) under comparison/raw_runs/. One run
gives Figs 1-3; two or more also give Fig 4 (differences vs the first as baseline).
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================== CONFIG ==============================

PROJECT_BASE_DIR = r"/"
OUTPUT_BASE_DIR  = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs")

# Run folder name(s) under output/raw_runs/. First entry is the baseline for Fig 4.
RUNS = [
    "HAT_1967_1997_groinTest_no_groin",
    # "HAT_1967_1997_groinTest_groin",
]

# ── Geometry (must match the run) ──
NUM_REAL_DOMAINS   = 11
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 2
LAST_FILE_NUMBER   = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1     # 12
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS

START_YEAR = 1967
END_YEAR   = 1997

# ── Sign / rate conventions (match main hindcast exactly) ──
# Main script: change_rate = (x_s[-1] - x_s[0]) / (nt-1); then *= -1 if FLIP.
# Raw x_s increases LANDWARD (retreat). With FLIP_SIGN_MODEL=True, plotted
# rate is + = seaward/accretion, - = landward/erosion (same as your main run).
FLIP_SIGN_MODEL = True
SEA_LEVEL_RISE_RATE = 0.004   # for title only (matches the run)

# ── Styling (from your main hindcast) ──
MODEL_COLOR   = "#FF8C00"   # warm orange
GROIN_COLOR   = "#B71C1C"   # groin red
GROIN_BOUNDARY_GIS = 5.5    # D5/D6 interface (Buxton groin)
DOMAIN_SPACING_M   = 500
DOMAIN_TICK_STEP   = 5

# Domains whose position-vs-time trajectory to draw in Fig 3 (GIS ids).
TRAJECTORY_DOMAINS_GIS = [3, 5, 6, 7, 9, 12]

REAL_DOMAINS_ONLY = True    # focus x-axis on D2-D12; if False, show buffers too
SAVE_FIGS = True
SHOW_FIGS = True

# --- Observed reference (raw ArcGIS dune-line offset files) for the model-vs-
#     observed comparison figure. Used to overlay the real 1967->1997 change.
RAW_OFFSET_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "groin_module", "hindcast_groin_test",
    "groin_init", "2-brie-offset", "raw_offsets",
)
OBS_START_YEAR = 1967
OBS_END_YEAR   = 1997
OBS_RAW_START  = os.path.join(RAW_OFFSET_DIR, f"{OBS_START_YEAR}_duneline_offset_raw.csv")
OBS_RAW_END    = os.path.join(RAW_OFFSET_DIR, f"{OBS_END_YEAR}_duneline_offset_raw.csv")
OBS_DOMAIN_COL = "domain_id"
OBS_POS_COL    = "ORIG_LEN"

# --- Position-mode figures (match main hindcast convention) ---
# Real island planform with year-0 as the 0 reference, ocean-at-bottom layout.
OCEAN_AT_BOTTOM = True     # seaward plots downward (matches Hatteras cross-shore)


# ============================ HELPERS ============================

def _gis_to_pad(gis_id):
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


def _load_shoreline(run_name):
    """Load a run's shoreline matrix (nt, ndomain) in meters."""
    path = os.path.join(OUTPUT_BASE_DIR, run_name, f"{run_name}_shoreline_matrix.npy")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"shoreline matrix not found:\n  {path}")
    m = np.load(path)
    print(f"  Loaded {run_name}: shape {m.shape}")
    return m


def _flip(v):
    """Apply FLIP_SIGN_MODEL exactly as the main hindcast does."""
    return v * (-1.0 if FLIP_SIGN_MODEL else 1.0)


def _total_change(m):
    """(x_s[-1] - x_s[0]) per domain, flipped -> + = seaward/accretion. Real slice."""
    return _flip(_real_slice(m[-1]) - _real_slice(m[0]))


def _change_rate(m):
    """Main script convention: total_change / (nt - 1), flipped. Real slice."""
    nt = m.shape[0]
    denom = max(nt - 1, 1)
    return _flip(_real_slice(m[-1]) - _real_slice(m[0])) / float(denom)


def _real_slice(arr_1d):
    return arr_1d[START_REAL_INDEX:END_REAL_INDEX]


def _gis_axis():
    return np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)


def _mark_groin(ax, y_for_label=None):
    ax.axvline(GROIN_BOUNDARY_GIS, color=GROIN_COLOR, lw=1.5, ls="--",
               alpha=0.9, zorder=5)
    yl = ax.get_ylim()
    y = y_for_label if y_for_label is not None else yl[0] + 0.9 * (yl[1] - yl[0])
    ax.text(GROIN_BOUNDARY_GIS, y, " Buxton groin", color=GROIN_COLOR,
            fontsize=8, rotation=90, va="top", ha="left", alpha=0.9)


def _updrift_downdrift_shading(ax):
    """Light shading: updrift (D6+) vs downdrift (<=D5, not validated)."""
    ax.axvspan(FIRST_FILE_NUMBER - 0.5, GROIN_BOUNDARY_GIS,
               alpha=0.06, color="firebrick", zorder=0)   # downdrift
    ax.axvspan(GROIN_BOUNDARY_GIS, LAST_FILE_NUMBER + 0.5,
               alpha=0.06, color="seagreen", zorder=0)     # updrift


# ============================ FIGURES ============================

def fig_position_change(runs_data):
    """FIG 1: shoreline change (m, end - start) per domain."""
    gis = _gis_axis()
    fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
    for run_name, m in runs_data.items():
        ax.plot(gis, _total_change(m), marker="o", ms=4, lw=1.8, label=run_name)
    _updrift_downdrift_shading(ax)
    ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}\u2013{LAST_FILE_NUMBER})")
    ax.set_ylabel("Shoreline change (m)  [erosion \u25b2]")
    ax.set_title(f"Shoreline change {START_YEAR}\u2013{END_YEAR} (end \u2212 start)   |   "
                 f"erosion up / accretion down   |   updrift = D6+  downdrift = D5\u2212")
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom: erosion up
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return fig, "position_change"


def fig_change_rate(runs_data):
    """FIG 2: shoreline change RATE (m/yr) per domain -- LRR-style."""
    gis = _gis_axis()
    fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
    for run_name, m in runs_data.items():
        color = MODEL_COLOR if len(runs_data) == 1 else None
        ax.plot(gis, _change_rate(m), marker="o", ms=4, lw=2, color=color,
                label=run_name)
    _updrift_downdrift_shading(ax)
    ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}\u2013{LAST_FILE_NUMBER})")
    ax.set_ylabel("Shoreline change rate (m/yr)  [erosion \u25b2]")
    ax.set_title(f"Modeled Shoreline Change Rate \u2013 Hatteras Island | "
                 f"SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
                 f"{START_YEAR}\u2013{END_YEAR}")
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom: erosion up
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return fig, "change_rate"


def fig_trajectories(runs_data):
    """FIG 3: shoreline position over time for selected domains."""
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    years = np.arange(START_YEAR, START_YEAR + next(iter(runs_data.values())).shape[0])
    cmap = plt.cm.viridis(np.linspace(0, 0.9, len(TRAJECTORY_DOMAINS_GIS)))
    for run_idx, (run_name, m) in enumerate(runs_data.items()):
        ls = "-" if run_idx == 0 else "--"
        for c, gis_id in zip(cmap, TRAJECTORY_DOMAINS_GIS):
            pad = _gis_to_pad(gis_id)
            pos = _flip(m[:, pad])       # + = seaward when FLIP_SIGN_MODEL
            pos = pos - pos[0]  # relative to start, so trajectories share an origin
            lbl = f"D{gis_id}" if run_idx == 0 else None
            ax.plot(years[:len(pos)], pos, color=c, ls=ls, lw=1.8, label=lbl)
    ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    ax.set_xlabel("Year")
    ax.set_ylabel("Shoreline position change since start (m)   [landward \u25b2]")
    ax.set_title("Shoreline trajectories by domain"
                 + ("   (solid = 1st run, dashed = others)" if len(runs_data) > 1 else ""))
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom: landward/erosion up
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, title="Domain", ncol=2)
    return fig, "trajectories"


def fig_difference(runs_data):
    """FIG 4: run-vs-baseline difference -- isolates the groin signal."""
    names = list(runs_data.keys())
    baseline = names[0]
    base_m = runs_data[baseline]
    gis = _gis_axis()
    fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
    for run_name in names[1:]:
        m = runs_data[run_name]
        ax.plot(gis, _total_change(m) - _total_change(base_m), marker="o", ms=4, lw=2,
                label=f"{run_name}\n minus {baseline}")
    _updrift_downdrift_shading(ax)
    ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    _mark_groin(ax)
    ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}\u2013{LAST_FILE_NUMBER})")
    ax.set_ylabel("\u0394 shoreline change vs baseline (m)  [erosion \u25b2]")
    ax.set_title("Isolated groin signal   (run \u2212 baseline)   "
                 "|   validate UPDRIFT D6\u2013D12")
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom: erosion up
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return fig, "difference_vs_baseline"


# ============================== MAIN ==============================

def _load_observed_change():
    """Observed per-domain dune-line change 1967->1997 from raw offset files.
    Returns (gis_array, obs_change_m) where + = landward/erosion (RAW sign),
    or (None, None) if the raw files aren't found. Matches HAT_target extraction."""
    missing = [p for p in (OBS_RAW_START, OBS_RAW_END) if not os.path.isfile(p)]
    if missing:
        print("  [observed] MISSING file(s) -- model-vs-observed panel will be model-only:")
        for p in missing:
            print(f"    {p}")
        return None, None

    def per_domain(path):
        df = pd.read_csv(path, encoding="utf-8-sig")
        df = df[df[OBS_DOMAIN_COL].isin(range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1))]
        return df.groupby(OBS_DOMAIN_COL)[OBS_POS_COL].mean()

    p0 = per_domain(OBS_RAW_START)
    p1 = per_domain(OBS_RAW_END)
    gis = _gis_axis()
    obs_change = np.array([p1.get(d, np.nan) - p0.get(d, np.nan) for d in gis])

    n_ok = int(np.isfinite(obs_change).sum())
    print(f"  [observed] Loaded {OBS_START_YEAR}/{OBS_END_YEAR} raw offsets OK -- "
          f"{n_ok}/{len(gis)} domains (D{FIRST_FILE_NUMBER}-D{LAST_FILE_NUMBER}) have "
          f"both years present.")
    if n_ok < len(gis):
        missing_domains = [int(d) for d, v in zip(gis, obs_change) if not np.isfinite(v)]
        print(f"  [observed] WARNING: no data for domain(s) {missing_domains} in one "
              f"or both years -- gap(s) will show as a break in the observed line.")

    return gis, obs_change


def fig_model_vs_observed(runs_data):
    """Two-panel model-vs-observed comparison:
      TOP    absolute shoreline POSITION -- start (1967) and modeled end (1997)
      BOTTOM CHANGE since 1967 -- modeled vs observed (the validity check)
    The change panel is the apples-to-apples overlay: both are meters of change
    from 1967, + = seaward/accretion (flipped from raw)."""
    gis = _gis_axis()
    obs_gis, obs_change_raw = _load_observed_change()
    have_obs = obs_gis is not None

    fig, (axP, axC) = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)

    # ── TOP: absolute positions (start + modeled end) ──
    for run_name, m in runs_data.items():
        start_pos = _flip(_real_slice(m[0]))     # + = seaward
        end_pos   = _flip(_real_slice(m[-1]))
        if run_name == list(runs_data)[0]:
            axP.plot(gis, start_pos, marker="o", ms=4, lw=1.6, color="0.5",
                     ls="--", label=f"Start ({START_YEAR}) shoreline")
        axP.plot(gis, end_pos, marker="o", ms=4, lw=2, label=f"{run_name} end ({END_YEAR})")
    _updrift_downdrift_shading(axP)
    _mark_groin(axP)
    axP.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    axP.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}\u2013D{LAST_FILE_NUMBER})")
    axP.set_ylabel("Shoreline position (m)  [landward \u25b2]")
    axP.set_title(f"Absolute shoreline position: start vs modeled end")
    axP.set_ylim(axP.get_ylim()[::-1])   # invert: seaward down, landward/erosion up
    axP.grid(alpha=0.3); axP.legend(fontsize=8)

    # ── BOTTOM: change since 1967, model vs observed ──
    for run_name, m in runs_data.items():
        model_change = _total_change(m)          # + = seaward, real slice
        axC.plot(gis, model_change, marker="o", ms=4, lw=2,
                 label=f"{run_name} (modeled)")
    if have_obs:
        obs_change = -obs_change_raw             # flip: + = seaward, match model
        axC.plot(obs_gis, obs_change, marker="s", ms=6, lw=2.4, color="black",
                 ls="--", label=f"Observed {OBS_START_YEAR}\u2013{OBS_END_YEAR} (target)",
                 zorder=6)
    else:
        axC.text(0.5, 0.9, "observed raw files not found -- model only",
                 transform=axC.transAxes, ha="center", color="firebrick", fontsize=9)
    _updrift_downdrift_shading(axC)
    axC.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    _mark_groin(axC)
    axC.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    axC.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}\u2013D{LAST_FILE_NUMBER})")
    axC.set_ylabel("Shoreline change since 1967 (m)  [erosion \u25b2]")
    axC.set_title("Change vs observed target   |   validate UPDRIFT D6\u2013D12 "
                  "(downdrift D2\u2013D5 not a target: Cape dynamics)")
    axC.set_ylim(axC.get_ylim()[::-1])   # invert: erosion (negative) up
    axC.grid(alpha=0.3); axC.legend(fontsize=8)

    return fig, "model_vs_observed"


def fig_position_planform(runs_data):
    """Position-mode plot (main-hindcast convention): the REAL island planform,
    with the year-0 (1967) shoreline as the 0 reference, and how it changes.
    Cross-shore position relative to the year-0 alongshore mean, ocean at bottom.
    This is the 'real island orientation as position 0' view."""
    gis = _gis_axis()
    fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)

    for i, (run_name, m) in enumerate(runs_data.items()):
        pos = _flip(_real_slice(m[0]))              # year-0 position, + = seaward
        ref_mean = np.nanmean(pos)
        # year-0 planform (relative to its alongshore mean) -- the reference shape
        if i == 0:
            planform0 = pos - ref_mean
            ax.plot(gis, planform0, marker="o", ms=4, lw=1.8, color="0.5", ls="--",
                    label=f"{START_YEAR} shoreline (reference)", zorder=3)
        # end-of-run planform, same reference
        end = _flip(_real_slice(m[-1])) - ref_mean
        ax.plot(gis, end, marker="o", ms=4, lw=2.2, zorder=5,
                label=f"{run_name} ({END_YEAR})")

    _updrift_downdrift_shading(ax)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}\u2013D{LAST_FILE_NUMBER})")
    up_word = "landward" if OCEAN_AT_BOTTOM else "seaward"
    ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)\n{up_word} \u25b2")
    ax.set_title(f"Island planform: {START_YEAR} reference vs modeled end   "
                 f"(real orientation, position 0 = {START_YEAR} mean)")
    if OCEAN_AT_BOTTOM:
        ax.set_ylim(ax.get_ylim()[::-1])   # invert: seaward downward
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return fig, "position_planform"


def main():
    print("=" * 70)
    print("Plotting groin-test runs")
    print("=" * 70)

    runs_data = {}
    for run_name in RUNS:
        try:
            runs_data[run_name] = _load_shoreline(run_name)
        except FileNotFoundError as e:
            print(f"  [SKIP] {e}")
    if not runs_data:
        print("No runs loaded. Edit RUNS to point at saved run folders.")
        return

    figs = [fig_position_change(runs_data),
            fig_change_rate(runs_data),
            fig_trajectories(runs_data),
            fig_model_vs_observed(runs_data),
            fig_position_planform(runs_data)]
    if len(runs_data) > 1:
        figs.append(fig_difference(runs_data))

    if SAVE_FIGS:
        # Save into the first run's folder.
        out_dir = os.path.join(OUTPUT_BASE_DIR, RUNS[0])
        os.makedirs(out_dir, exist_ok=True)
        for fig, name in figs:
            out = os.path.join(out_dir, f"PLOT_{name}.png")
            fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
            print(f"  Saved: {out}")

    if SHOW_FIGS:
        plt.show()


if __name__ == "__main__":
    main()
