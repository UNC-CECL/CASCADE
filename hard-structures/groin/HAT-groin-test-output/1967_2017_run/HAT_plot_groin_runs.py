"""
HAT_plot_groin_runs.py
======================
Plots for the 1967-2017 groin-DETERIORATION-test runs. Loads saved shoreline
matrices (the *_shoreline_matrix.npy files written by the base / groin
runners) and produces:

  FIG 1  Shoreline CHANGE (m, end - start) per domain          [position change]
  FIG 2  Shoreline change RATE (m/yr) per domain               [LRR-style]
  FIG 3  Shoreline POSITION over time, selected domains        [trajectories]
  FIG 4  (multi-run only) run-vs-run DIFFERENCE, e.g.
         groin_be minus no_groin -> the isolated groin signal

Matches your main hindcast's conventions: model orange (#FF8C00), groin red
(#B71C1C) at the D5/D6 boundary, GIS domain x-axis, real-domain focus with
buffers shaded. No CoastSat dependency -- this is for inspecting the model runs
themselves. Point RUNS at one or more saved run folders.

IMPORTANT -- END_YEAR is EXCLUSIVE (RUN_YEARS = END_YEAR - START_YEAR), so
END_YEAR=2018 here means the run's actual final modeled year is 2017, not
2018. Every figure label below derives the true final year from the loaded
data's own shape (_final_model_year()) rather than trusting END_YEAR
directly -- END_YEAR is used only for RUN_NAME_SUFFIX-adjacent bookkeeping,
never for a figure title.

OBSERVED COMPARISON: fig_model_vs_observed() now reads from the validated
wet/dry change table (Change_from_wetdry_1967_D2_D12.csv, built by
HAT_geometric_distance_sanity_check.py) instead of the raw dune-line CSVs --
see HAT_groin_effect_comparison.py's module docstring for the full reasoning
(x_s is conceptually closer to a water-line proxy than to the dune line, and
the wet/dry line was confirmed seaward of the dune line at every domain).

Usage
-----
Edit RUNS below to list the run folder name(s) under output/raw_runs/. One run
gives Figs 1-3; two or more also give Fig 4 (differences vs the first as baseline).
"""

import os
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================== CONFIG ==============================

# Derived, not hardcoded. This was `r"/"`, which silently resolved every data
# path to the filesystem root -- the scripts imported fine and then reported
# every input as missing. Walking up to pyproject.toml survives the repo
# reorganisation that moved this tree from scripts/groin/ to hard-structures/.
PROJECT_BASE_DIR = str(next(
    p for p in pathlib.Path(__file__).resolve().parents
    if (p / "pyproject.toml").exists()))
OUTPUT_BASE_DIR  = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs")

# Run folder name(s) under output/raw_runs/. First entry is the baseline for Fig 4.
RUNS = [
    "HAT_1967_2018_no_BE_no_groin",
    # "HAT_1967_2018_M60_deterioration_groin",
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
END_YEAR   = 2018   # EXCLUSIVE convention -- see module docstring; NOT a label year

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

# --- Observed reference for the model-vs-observed comparison (fig_model_vs_observed):
# the validated wet/dry change table, NOT the raw dune-line CSVs -- see module
# docstring. OBSERVED_YEARS are the checkpoint years to overlay, roughly every
# ~10 years, chosen from years with full (or near-full) D2-D12 coverage in the
# change table -- 2017 itself only has 6/11 domains (D8-D12 missing), so 2018
# (11/11, one year later than the model's actual endpoint) is used instead.
OBSERVED_YEARS = [1978, 1987, 1997, 2008, 2018]
WETDRY_CHANGE_TABLE = os.path.join(
    PROJECT_BASE_DIR, "hard-structures", "groin", "HAT-hindcast-groin-test",
    "input_prep", "shoreline_position", "output",
    "Change_from_wetdry_1967_D2_D12.csv",
)
WETDRY_DOMAIN_COL = "Domain_ID"

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


def _final_model_year(runs_data):
    """True final modeled year, derived from the loaded data's own shape
    (START_YEAR + nt - 1) rather than trusting END_YEAR directly -- END_YEAR
    is EXCLUSIVE (see module docstring), so using it as a label would be off
    by one. Warns (doesn't crash) if runs disagree in length."""
    lengths = {name: m.shape[0] for name, m in runs_data.items()}
    if len(set(lengths.values())) > 1:
        print(f"  WARNING: runs have different lengths {lengths} -- "
              f"using the first run's length for figure labels.")
    nt = next(iter(lengths.values()))
    return START_YEAR + nt - 1


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
    final_year = _final_model_year(runs_data)
    fig, ax = plt.subplots(figsize=(12, 5), constrained_layout=True)
    for run_name, m in runs_data.items():
        ax.plot(gis, _total_change(m), marker="o", ms=4, lw=1.8, label=run_name)
    _updrift_downdrift_shading(ax)
    ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}\u2013{LAST_FILE_NUMBER})")
    ax.set_ylabel("Shoreline change (m)  [erosion \u25b2]")
    ax.set_title(f"Shoreline change {START_YEAR}\u2013{final_year} (end \u2212 start)   |   "
                 f"erosion up / accretion down   |   updrift = D6+  downdrift = D5\u2212")
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom: erosion up
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    return fig, "position_change"


def fig_change_rate(runs_data):
    """FIG 2: shoreline change RATE (m/yr) per domain -- LRR-style."""
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
    ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}\u2013{LAST_FILE_NUMBER})")
    ax.set_ylabel("Shoreline change rate (m/yr)  [erosion \u25b2]")
    ax.set_title(f"Modeled Shoreline Change Rate \u2013 Hatteras Island | "
                 f"SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
                 f"{START_YEAR}\u2013{final_year}")
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

def _load_observed_changes():
    """Observed per-domain shoreline change START_YEAR->year, for every year in
    OBSERVED_YEARS, from the validated wet/dry change table. Returns
    dict {year: np.ndarray or None} -- + = landward/erosion (RAW sign)."""
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


def fig_model_vs_observed(runs_data):
    """Two-panel model-vs-observed comparison:
      TOP    absolute shoreline POSITION -- start (1967) and modeled end
      BOTTOM CHANGE since 1967 -- modeled vs observed AT EVERY YEAR IN
             OBSERVED_YEARS (the validity check), color-coded chronologically
             (blue=older, red=newer) so multiple decades overlay readably.
    The change panel is the apples-to-apples overlay: both are meters of change
    from 1967, + = seaward/accretion (flipped from raw)."""
    gis = _gis_axis()
    final_year = _final_model_year(runs_data)
    observed_changes = _load_observed_changes()
    years_with_data = [y for y, v in observed_changes.items() if v is not None]

    fig, (axP, axC) = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)

    # ── TOP: absolute positions (start + modeled end) ──
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
    axP.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}\u2013D{LAST_FILE_NUMBER})")
    axP.set_ylabel("Shoreline position (m)  [landward \u25b2]")
    axP.set_title(f"Absolute shoreline position: start vs modeled end")
    axP.set_ylim(axP.get_ylim()[::-1])   # invert: seaward down, landward/erosion up
    axP.grid(alpha=0.3); axP.legend(fontsize=8)

    # ── BOTTOM: change since 1967, model vs observed at every available year ──
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
                      ls="--", alpha=0.85, label=f"Observed {START_YEAR}\u2013{year}", zorder=5)
    else:
        axC.text(0.5, 0.9, "observed wet/dry change table not found -- model only",
                 transform=axC.transAxes, ha="center", color="firebrick", fontsize=9)
    _updrift_downdrift_shading(axC)
    axC.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    _mark_groin(axC)
    axC.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    axC.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}\u2013D{LAST_FILE_NUMBER})")
    axC.set_ylabel("Shoreline change since 1967 (m)  [erosion \u25b2]")
    axC.set_title("Change vs observed targets (~10-yr increments)   |   validate "
                  "UPDRIFT D6\u2013D12 (downdrift D2\u2013D5 not a target: Cape dynamics)")
    axC.set_ylim(axC.get_ylim()[::-1])   # invert: erosion (negative) up
    axC.grid(alpha=0.3); axC.legend(fontsize=7.5, ncol=2)

    return fig, "model_vs_observed"


def fig_position_planform(runs_data):
    """Position-mode plot (main-hindcast convention): the REAL island planform,
    with the year-0 (1967) shoreline as the 0 reference, and how it changes.
    Cross-shore position relative to the year-0 alongshore mean, ocean at bottom.
    This is the 'real island orientation as position 0' view."""
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
