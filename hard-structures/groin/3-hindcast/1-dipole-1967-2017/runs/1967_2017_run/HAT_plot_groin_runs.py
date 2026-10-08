"""
Figures of saved 1967-2017 rig runs: change, rate, trajectories, model vs observed wet/dry change, planform.

    python HAT_plot_groin_runs.py

Reads each run's *_shoreline_matrix.npy under output/raw_runs/<run>/ (RUNS)
and the wet/dry change table; saves PLOT_<name>.png into the first run's
folder. The rig runners import its fig_ functions for per-run figures.
Final years come from each matrix's length, never END_YEAR.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

import os
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# --- CONFIG ------------------------------------------------------------------
# Repo root, found by searching upward (ORGANIZATION.md rule 5)
PROJECT_BASE_DIR = str(next(
    p for p in pathlib.Path(__file__).resolve().parents
    if (p / "pyproject.toml").exists()))
OUTPUT_BASE_DIR  = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs")

# Run folder name(s) under output/raw_runs/; the first is the baseline for the difference figure
RUNS = [
    # The sweep's best cell, M = 60, f = 0.6, on 1984-start/v1 (set 2026-08-30)
    "HAT_1967_2018_edge_calibrated_groin",
    # "HAT_1967_2018_no_BE_no_groin",
]

# Geometry, must match the run
NUM_REAL_DOMAINS   = 11
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 2
LAST_FILE_NUMBER   = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1     # 12
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS

START_YEAR = 1967
END_YEAR   = 2018   # EXCLUSIVE convention -- see module docstring; NOT a label year

# Sign: raw x_s grows landward; True plots + = seaward, as the main hindcast (README)
FLIP_SIGN_MODEL = True
SEA_LEVEL_RISE_RATE = 0.004   # for title only (matches the run)

# Styling, as the main hindcast
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

# Observed checkpoint years, ~10 apart with full D2-D12 cover; 2018 stands in for 2017 (README)
OBSERVED_YEARS = [1978, 1987, 1997, 2008, 2018]
# The validated wet/dry change table, not the dune-line CSVs (README)
WETDRY_CHANGE_TABLE = os.path.join(
    PROJECT_BASE_DIR, "hard-structures", "groin", "1-observations",
    "wetdry_photo_positions",
    "Change_from_wetdry_1967_D2_D12.csv",
)

WETDRY_DOMAIN_COL = "Domain_ID"

# Position-mode figures: year-0 planform as the 0 reference, seaward plotted downward
OCEAN_AT_BOTTOM = True     # seaward plots downward (matches Hatteras cross-shore)
# -----------------------------------------------------------------------------


# Refuse to draw model-only figures that look like comparisons
if not os.path.isfile(WETDRY_CHANGE_TABLE):
    raise SystemExit(
        "observed wet/dry table not found at " + WETDRY_CHANGE_TABLE
        + " -- refusing to draw model-only figures that look like comparisons.")


# GIS domain id -> padded index
def _gis_to_pad(gis_id):
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


# A run's shoreline matrix (nt, ndomain), metres
def _load_shoreline(run_name):
    path = os.path.join(OUTPUT_BASE_DIR, run_name, f"{run_name}_shoreline_matrix.npy")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"shoreline matrix not found:\n  {path}")
    m = np.load(path)
    print(f"  Loaded {run_name}: shape {m.shape}")
    return m


# True final modelled year from the data's length (END_YEAR is exclusive); warns if runs differ
def _final_model_year(runs_data):
    lengths = {name: m.shape[0] for name, m in runs_data.items()}
    if len(set(lengths.values())) > 1:
        print(f"  WARNING: runs have different lengths {lengths} -- "
              f"using the first run's length for figure labels.")
    nt = next(iter(lengths.values()))
    return START_YEAR + nt - 1


# Apply FLIP_SIGN_MODEL exactly as the main hindcast does
def _flip(v):
    return v * (-1.0 if FLIP_SIGN_MODEL else 1.0)


# End minus start per real domain, flipped: + = seaward/accretion
def _total_change(m):
    return _flip(_real_slice(m[-1]) - _real_slice(m[0]))


# Total change / (nt - 1), flipped, per real domain (the main script's convention)
def _change_rate(m):
    nt = m.shape[0]
    denom = max(nt - 1, 1)
    return _flip(_real_slice(m[-1]) - _real_slice(m[0])) / float(denom)


# The real-domain slice of a padded row
def _real_slice(arr_1d):
    return arr_1d[START_REAL_INDEX:END_REAL_INDEX]


# GIS ids of the real domains
def _gis_axis():
    return np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)


# The groin line and its label
def _mark_groin(ax, y_for_label=None):
    ax.axvline(GROIN_BOUNDARY_GIS, color=GROIN_COLOR, lw=1.5, ls="--",
               alpha=0.9, zorder=5)
    yl = ax.get_ylim()
    y = y_for_label if y_for_label is not None else yl[0] + 0.9 * (yl[1] - yl[0])
    ax.text(GROIN_BOUNDARY_GIS, y, " Buxton groin", color=GROIN_COLOR,
            fontsize=8, rotation=90, va="top", ha="left", alpha=0.9)


# Light shading: updrift (D6+) vs downdrift (D5 and south, not validated)
def _updrift_downdrift_shading(ax):
    ax.axvspan(FIRST_FILE_NUMBER - 0.5, GROIN_BOUNDARY_GIS,
               alpha=0.06, color="firebrick", zorder=0)   # downdrift
    ax.axvspan(GROIN_BOUNDARY_GIS, LAST_FILE_NUMBER + 0.5,
               alpha=0.06, color="seagreen", zorder=0)     # updrift


# Fig 1: shoreline change (m, end - start) per domain
def fig_position_change(runs_data):
    gis = _gis_axis()
    final_year = _final_model_year(runs_data)
    fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
    for run_name, m in runs_data.items():
        ax.plot(gis, _total_change(m), marker="o", ms=4, lw=1.8, label=run_name)
    _updrift_downdrift_shading(ax)
    ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}–{LAST_FILE_NUMBER})")
    ax.set_ylabel("Shoreline change (m)  [erosion ▲]")
    ax.set_title(f"Shoreline change {START_YEAR}–{final_year} (end − start)   |   "
                 f"erosion up / accretion down   |   updrift = D6+  downdrift = D5−")
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom: erosion up
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return fig, "position_change"


# Fig 2: shoreline change rate (m/yr) per domain
def fig_change_rate(runs_data):
    gis = _gis_axis()
    final_year = _final_model_year(runs_data)
    fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
    for run_name, m in runs_data.items():
        color = MODEL_COLOR if len(runs_data) == 1 else None
        ax.plot(gis, _change_rate(m), marker="o", ms=4, lw=2, color=color,
                label=run_name)
    _updrift_downdrift_shading(ax)
    ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}–{LAST_FILE_NUMBER})")
    ax.set_ylabel("Shoreline change rate (m/yr)  [erosion ▲]")
    ax.set_title(f"Modeled Shoreline Change Rate – Hatteras Island | "
                 f"SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
                 f"{START_YEAR}–{final_year}")
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom: erosion up
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return fig, "change_rate"


# Fig 3: shoreline position over time for TRAJECTORY_DOMAINS_GIS
def fig_trajectories(runs_data):
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
    ax.set_ylabel("Shoreline position change since start (m)   [landward ▲]")
    ax.set_title("Shoreline trajectories by domain"
                 + ("   (solid = 1st run, dashed = others)" if len(runs_data) > 1 else ""))
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom: landward/erosion up
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, title="Domain", ncol=2)
    return fig, "trajectories"


# Observed change START_YEAR -> year per OBSERVED_YEARS, raw sign (+ = landward); None where missing
def _load_observed_changes():
    gis = _gis_axis()

    if not os.path.isfile(WETDRY_CHANGE_TABLE):
        print(f"  [observed] MISSING wet/dry change table -- model-vs-observed "
              f"panel will be model-only:\n    {WETDRY_CHANGE_TABLE}")
        return {y: None for y in OBSERVED_YEARS}

    df = pd.read_csv(WETDRY_CHANGE_TABLE).set_index(WETDRY_DOMAIN_COL)

    changes = {}
    for year in OBSERVED_YEARS:
        col = f"change_from_wetdry_1967_wetdry_{year}_m"
        if col not in df.columns:
            print(f"  [observed] MISSING column '{col}' -- {year} omitted.")
            changes[year] = None
            continue
        obs_change = np.array([df[col].get(d, np.nan) for d in gis])
        n_ok = int(np.isfinite(obs_change).sum())
        print(f"  [observed] Loaded {START_YEAR}->{year} wet/dry change OK -- "
              f"{n_ok}/{len(gis)} domains.")
        if n_ok < len(gis):
            missing_domains = [int(d) for d, v in zip(gis, obs_change) if not np.isfinite(v)]
            print(f"    WARNING: no data for domain(s) {missing_domains} -- "
                  f"gap(s) will show as a break in that year's line.")
        changes[year] = obs_change

    return changes


# Fig 4: start and modelled-end positions above; change since 1967, model vs every observed year, below
def fig_model_vs_observed(runs_data):
    gis = _gis_axis()
    final_year = _final_model_year(runs_data)
    observed_changes = _load_observed_changes()
    years_with_data = [y for y, v in observed_changes.items() if v is not None]

    fig, (axP, axC) = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)

    # Top: absolute positions, start and modelled end
    for run_name, m in runs_data.items():
        start_pos = _flip(_real_slice(m[0]))     # + = seaward
        end_pos   = _flip(_real_slice(m[-1]))
        if run_name == list(runs_data)[0]:
            axP.plot(gis, start_pos, marker="o", ms=4, lw=1.6, color="0.5",
                     ls="--", label=f"Start ({START_YEAR}) shoreline")
        axP.plot(gis, end_pos, marker="o", ms=4, lw=2, label=f"{run_name} end ({final_year})")
    _updrift_downdrift_shading(axP)
    _mark_groin(axP)
    axP.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    axP.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}–D{LAST_FILE_NUMBER})")
    axP.set_ylabel("Shoreline position (m)  [landward ▲]")
    axP.set_title(f"Absolute shoreline position: start vs modeled end")
    axP.set_ylim(axP.get_ylim()[::-1])   # invert: seaward down, landward/erosion up
    axP.grid(alpha=0.3); axP.legend(fontsize=8)

    # Bottom: change since 1967, model vs observed at every available year
    for run_name, m in runs_data.items():
        model_change = _total_change(m)          # + = seaward, real slice
        axC.plot(gis, model_change, marker="o", ms=4, lw=2.5, color=MODEL_COLOR,
                 zorder=6, label=f"{run_name} (modeled)")

    if years_with_data:
        vmin, vmax = min(years_with_data), max(years_with_data)
        norm = plt.Normalize(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1)
        cmap = plt.cm.coolwarm
        for year in years_with_data:
            obs_change = -observed_changes[year]     # flip: + = seaward, match model
            axC.plot(gis, obs_change, marker="s", ms=5, lw=1.8, color=cmap(norm(year)),
                      ls="--", alpha=0.85, label=f"Observed {START_YEAR}–{year}", zorder=5)
    else:
        axC.text(0.5, 0.9, "observed wet/dry change table not found -- model only",
                 transform=axC.transAxes, ha="center", color="firebrick", fontsize=9)
    _updrift_downdrift_shading(axC)
    axC.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    _mark_groin(axC)
    axC.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    axC.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}–D{LAST_FILE_NUMBER})")
    axC.set_ylabel("Shoreline change since 1967 (m)  [erosion ▲]")
    axC.set_title("Change vs observed targets (~10-yr increments)   |   validate "
                  "UPDRIFT D6–D12 (downdrift D2–D5 not a target: Cape dynamics)")
    axC.set_ylim(axC.get_ylim()[::-1])   # invert: erosion (negative) up
    axC.grid(alpha=0.3); axC.legend(fontsize=7.5, ncol=2)

    return fig, "model_vs_observed"


# Fig 5: the real planform against the 1967 alongshore mean, ocean at bottom
def fig_position_planform(runs_data):
    gis = _gis_axis()
    final_year = _final_model_year(runs_data)
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
                label=f"{run_name} ({final_year})")

    _updrift_downdrift_shading(ax)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}–D{LAST_FILE_NUMBER})")
    up_word = "landward" if OCEAN_AT_BOTTOM else "seaward"
    ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)\n{up_word} ▲")
    ax.set_title(f"Island planform: {START_YEAR} reference vs modeled end   "
                 f"(real orientation, position 0 = {START_YEAR} mean)")
    if OCEAN_AT_BOTTOM:
        ax.set_ylim(ax.get_ylim()[::-1])   # invert: seaward downward
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return fig, "position_planform"


# Fig 6 (two or more runs): each run minus the first, the isolated groin signal
def fig_difference(runs_data):
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
    ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}–{LAST_FILE_NUMBER})")
    ax.set_ylabel("Δ shoreline change vs baseline (m)  [erosion ▲]")
    ax.set_title("Isolated groin signal   (run − baseline)   "
                 "|   validate UPDRIFT D6–D12")
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom: erosion up
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return fig, "difference_vs_baseline"


# Run: load RUNS, draw every figure, save into the first run's folder
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
