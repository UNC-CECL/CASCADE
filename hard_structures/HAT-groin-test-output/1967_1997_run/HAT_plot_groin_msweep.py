"""
HAT_plot_groin_msweep.py
========================
Compare multiple groin runs (an M-sweep, or any set of runs) against the observed
1967->1997 target on a SINGLE axis, so you can see convergence at a glance instead
of flipping between per-run figures.

Produces:
  FIG A  change-vs-observed overlay -- every run's D2-D12 change curve, color-
         graded, with the observed target in bold black. Ocean-at-bottom.
  FIG B  updrift-only skill vs M -- RMSE (and mean bias) between modeled and
         observed over the VALIDATION zone D6-D12, one point per run. The
         minimum-RMSE run is the best-fit M. Downdrift is excluded on purpose.

Reads the same saved shoreline matrices the runner writes, and the same raw
offset files the target script uses (self-contained; no dependency on other
scripts having been run).

Convention matches HAT_plot_groin_runs.py: FLIP_SIGN_MODEL applied, + = seaward,
erosion plotted upward (ocean at bottom).
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================== CONFIG ==============================

PROJECT_BASE_DIR = r"/"
# WHERE THE RUNS ARE (must match the hindcast runner's output folder).
OUTPUT_BASE_DIR  = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs")

# Runs to compare. List the folder names under output/raw_runs/. Order sets the
# color ramp (low -> high M). Include the no_groin baseline if you want it shown.
RUNS = [
    "HAT_1967_1997_noBE_no_groin",
    "HAT_1967_1997_noBE_40M_groin",
    "HAT_1967_1997_noBE_50M_groin",
    "HAT_1967_1997_noBE_60M_groin",
    "HAT_1967_1997_noBE_70M_groin",
]

# Optional: pull the M value out of each run name for the skill-vs-M axis.
# Matches "..._<number>M_groin". no_groin -> M = 0.
M_FROM_NAME = re.compile(r"_(\d+)M_groin")

# ── Geometry (must match the runs) ──
NUM_REAL_DOMAINS   = 11
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 2
LAST_FILE_NUMBER   = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1     # 12
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS

START_YEAR = 1967
END_YEAR   = 1997
FLIP_SIGN_MODEL = True

# Validation zone (updrift only) for the skill metric.
VALIDATE_LO_GIS = 6
VALIDATE_HI_GIS = 12

GROIN_BOUNDARY_GIS = 5.5
GROIN_COLOR = "#B71C1C"
DOMAIN_TICK_STEP = 5

# Observed reference-line styling (position figure). Clean solid line, distinct
# accent color, white-edged markers -- reads as "the real world" without a halo.
OBS_TARGET_COLOR = "#8B0000"   # dark red, ties to the groin-red palette
OBS_TARGET_LW    = 2.6
OBS_TARGET_MS    = 8
OBS_START_COLOR  = "0.55"      # grey for the 1967 start reference

# Observed target (raw offset files), same extraction as HAT_target.
RAW_OFFSET_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "groin_module_noBE", "hindcast_groin_test",
    "groin_init", "2-brie-offset", "raw_offsets",
)
OBS_RAW_START = os.path.join(RAW_OFFSET_DIR, f"{START_YEAR}_duneline_offset_raw.csv")
OBS_RAW_END   = os.path.join(RAW_OFFSET_DIR, f"{END_YEAR}_duneline_offset_raw.csv")
OBS_DOMAIN_COL = "domain_id"
OBS_POS_COL    = "ORIG_LEN"

# --- Observed -> model coordinate conversion ---------------------------------
# The offset files share a fixed offshore baseline, and CASCADE initializes x_s
# from the 1967 file. Empirically model x_s = A*(ORIG_LEN - MIN1967) + B, where
# A ~= 0.1 (the m->decameter factor CASCADE stores x_s in) and B, MIN1967 are the
# 1967 reference constants. This maps ANY same-baseline year into the model x_s
# frame, so observed positions plot directly (verified: converted 1967 lands on
# model year-0 to <0.5 m). Recompute with calibrate_obs_to_model() if the init
# pipeline changes.
OBS2MODEL_A   = 0.10017      # ORIG_LEN(m) -> x_s(model): the m/dam factor
OBS2MODEL_MIN = 7096.00      # min ORIG_LEN over D2-D12 in the 1967 file
OBS2MODEL_B   = 1590.228     # x_s intercept (model frame, meters)

# Observed years to overlay as reference positions (need <year>_duneline_offset_
# raw.csv in RAW_OFFSET_DIR). All share the baseline, so all are convertible.
OBS_REFERENCE_YEARS = [1967, 1997]     # e.g. [1967, 1978, 1984, 1997] for all

# WHERE TO SAVE the comparison figures (independent of where runs are read from).
OUT_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "groin_module_noBE", "hindcast_groin_test", "comparison"
)
OUT_BASENAME = "HAT_1967_1997_noBE_Msweep"


# ============================ HELPERS ============================

def _gis_axis():
    return np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)


def _flip(v):
    return v * (-1.0 if FLIP_SIGN_MODEL else 1.0)


def _real_slice(a):
    return a[START_REAL_INDEX:END_REAL_INDEX]


def _load_matrix(run_name):
    p = os.path.join(OUTPUT_BASE_DIR, run_name, f"{run_name}_shoreline_matrix.npy")
    if not os.path.isfile(p):
        raise FileNotFoundError(p)
    return np.load(p)


def _model_change(run_name):
    """Per-domain modeled change since start, + = seaward (flipped). Real slice."""
    m = _load_matrix(run_name)
    return _flip(_real_slice(m[-1]) - _real_slice(m[0]))


def _n_years(run_name):
    """Number of years spanned by a run's matrix (rows - 1 = annual steps)."""
    m = _load_matrix(run_name)
    return max(1, m.shape[0] - 1)


def _obs_orig_len(year):
    """Mean ORIG_LEN per domain (D2-D12) for a given observed year, or None."""
    path = os.path.join(RAW_OFFSET_DIR, f"{year}_duneline_offset_raw.csv")
    if not os.path.isfile(path):
        return None
    df = pd.read_csv(path, encoding="utf-8-sig")
    df = df[df[OBS_DOMAIN_COL].isin(range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1))]
    s = df.groupby(OBS_DOMAIN_COL)[OBS_POS_COL].mean()
    gis = _gis_axis()
    return np.array([s.get(d, np.nan) for d in gis])


def _obs_to_model_xs(orig_len):
    """Convert observed ORIG_LEN (m, shared baseline) into model x_s coordinates
    (m), so observed positions can be plotted directly against model output."""
    return OBS2MODEL_A * (orig_len - OBS2MODEL_MIN) + OBS2MODEL_B


def calibrate_obs_to_model(run_name, year=1967):
    """One-off helper: fit model_x_s = A*(ORIG_LEN - MIN) + B from a run's year-0
    shoreline vs the observed <year> file, and print the constants to paste into
    the config. Run this if the init pipeline changes. Not called automatically."""
    m = _load_matrix(run_name)
    model_yr0 = _real_slice(m[0])          # model x_s, year 0, meters
    orig = _obs_orig_len(year)
    if orig is None:
        print(f"  [calibrate] no observed file for {year}")
        return
    mn = np.nanmin(orig)
    A, B = np.polyfit(orig - mn, model_yr0, 1)
    resid = model_yr0 - (A * (orig - mn) + B)
    print("calibrate_obs_to_model:")
    print(f"  OBS2MODEL_A   = {A:.5f}")
    print(f"  OBS2MODEL_MIN = {mn:.2f}")
    print(f"  OBS2MODEL_B   = {B:.3f}")
    print(f"  (residual std {np.nanstd(resid):.3f} m -- should be ~0)")


def _model_positions(run_name):
    """Per-domain absolute start and end positions, + = seaward (flipped). Real slice."""
    m = _load_matrix(run_name)
    start = _flip(_real_slice(m[0]))
    end   = _flip(_real_slice(m[-1]))
    return start, end


def _observed_change():
    """Observed 1967->1997 change, + = seaward (flipped to match model), or None."""
    if not (os.path.isfile(OBS_RAW_START) and os.path.isfile(OBS_RAW_END)):
        return None

    def per_domain(path):
        df = pd.read_csv(path, encoding="utf-8-sig")
        df = df[df[OBS_DOMAIN_COL].isin(range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1))]
        return df.groupby(OBS_DOMAIN_COL)[OBS_POS_COL].mean()

    p0 = per_domain(OBS_RAW_START)
    p1 = per_domain(OBS_RAW_END)
    gis = _gis_axis()
    raw = np.array([p1.get(d, np.nan) - p0.get(d, np.nan) for d in gis])  # + = landward
    return -raw   # flip: + = seaward, match model


def _m_value(run_name):
    m = M_FROM_NAME.search(run_name)
    if m:
        return float(m.group(1))
    if "no_groin" in run_name:
        return 0.0
    return np.nan


def _updrift_mask():
    gis = _gis_axis()
    return (gis >= VALIDATE_LO_GIS) & (gis <= VALIDATE_HI_GIS)


# ============================ FIGURES ============================

def fig_overlay(runs_rate, obs_rate, period_years):
    """All runs' change-RATE curves (m/yr) on one axis, vs the observed rate.
    Ocean at bottom.

    Rate rather than cumulative change: this is the LRR-native quantity you
    calibrate against and it's directly comparable to CoastSat rates. Model runs
    use a single-hue Blues ramp (light=low M, dark=high M); the observed target is
    a clean solid accent line. Set FADE_DOWNDRIFT=True to mute model lines left of
    the groin (not a target).
    """
    FADE_DOWNDRIFT = False   # True -> de-emphasize model lines in the D2-D5 zone

    gis = _gis_axis()
    fig, ax = plt.subplots(figsize=(13, 6), constrained_layout=True)

    n = len(runs_rate)
    cmap = plt.cm.Blues(np.linspace(0.40, 0.95, n))

    updrift = gis >= GROIN_BOUNDARY_GIS
    downdrift = ~updrift

    for c, (run_name, rate) in zip(cmap, runs_rate.items()):
        mval = _m_value(run_name)
        lbl = f"M={mval:.0f}" if not np.isnan(mval) else run_name
        if FADE_DOWNDRIFT:
            ax.plot(gis[downdrift], rate[downdrift], marker="o", ms=3, lw=1.2,
                    color=c, alpha=0.25)
            ax.plot(gis[updrift], rate[updrift], marker="o", ms=4, lw=2.2,
                    color=c, label=lbl)
        else:
            ax.plot(gis, rate, marker="o", ms=4, lw=2.2, color=c, label=lbl)

    if obs_rate is not None:
        ax.plot(gis, obs_rate, marker="s", ms=OBS_TARGET_MS,
                lw=OBS_TARGET_LW, color=OBS_TARGET_COLOR, ls="-",
                markeredgecolor="white", markeredgewidth=1.1,
                label=f"Observed {START_YEAR}\u2013{END_YEAR}", zorder=10)

    # shading + groin
    ax.axvspan(FIRST_FILE_NUMBER - 0.5, GROIN_BOUNDARY_GIS, alpha=0.06, color="firebrick")
    ax.axvspan(GROIN_BOUNDARY_GIS, LAST_FILE_NUMBER + 0.5, alpha=0.06, color="seagreen")
    ax.axvline(GROIN_BOUNDARY_GIS, color=GROIN_COLOR, ls="--", lw=1.6)
    ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}\u2013D{LAST_FILE_NUMBER})")
    ax.set_ylabel("Shoreline change rate (m/yr)  [erosion \u25b2]")
    ax.set_title(f"Groin M-sweep change RATE vs observed   |   validate UPDRIFT "
                 f"D6\u2013D12   ({START_YEAR}\u2013{END_YEAR}, {period_years} yr)")
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, title="Trapping rate", loc="upper right",
              framealpha=0.95, ncol=1)
    return fig, "rate_overlay"


def fig_skill_vs_m(runs_change, obs_change):
    """Updrift-only RMSE and mean bias vs M. Minimum RMSE = best-fit M."""
    if obs_change is None:
        return None, None
    mask = _updrift_mask()
    obs_up = obs_change[mask]

    rows = []
    for run_name, change in runs_change.items():
        mval = _m_value(run_name)
        resid = change[mask] - obs_up
        rmse = float(np.sqrt(np.nanmean(resid ** 2)))
        bias = float(np.nanmean(resid))   # + = model too seaward (under-eroded)
        rows.append((mval, rmse, bias, run_name))
    rows.sort(key=lambda r: (np.isnan(r[0]), r[0]))
    M = np.array([r[0] for r in rows])
    RMSE = np.array([r[1] for r in rows])
    BIAS = np.array([r[2] for r in rows])

    best_i = int(np.nanargmin(RMSE))
    best_M = M[best_i]

    fig, ax1 = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    ax1.plot(M, RMSE, marker="o", ms=7, lw=2, color="#08519C", label="RMSE (updrift)")
    ax1.axvline(best_M, color="#08519C", ls=":", lw=1.5, alpha=0.7)
    ax1.annotate(f"best M \u2248 {best_M:.0f}\nRMSE={RMSE[best_i]:.1f} m",
                 (best_M, RMSE[best_i]), textcoords="offset points",
                 xytext=(10, 10), fontsize=9, color="#08519C", fontweight="bold")
    ax1.set_xlabel("Groin trapping rate M (m/yr)")
    ax1.set_ylabel("Updrift RMSE vs observed (m)", color="#08519C")
    ax1.tick_params(axis="y", labelcolor="#08519C")
    ax1.grid(alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(M, BIAS, marker="s", ms=6, lw=1.8, color="#B71C1C", ls="--",
             label="mean bias")
    ax2.axhline(0, color="#B71C1C", lw=1, ls=":", alpha=0.6)
    ax2.set_ylabel("Mean bias (m)  [+ = model under-eroded]", color="#B71C1C")
    ax2.tick_params(axis="y", labelcolor="#B71C1C")

    ax1.set_title(f"Updrift skill (D{VALIDATE_LO_GIS}\u2013D{VALIDATE_HI_GIS}) vs M   "
                  f"|   best-fit M \u2248 {best_M:.0f}")
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc="upper center")
    return fig, "skill_vs_M"


def fig_positions(runs_positions, obs_change=None):
    """Final-POSITION comparison: each run's 1997 planform plus OBSERVED reference
    positions, all in the model frame. Positions shown relative to the 1967
    alongshore mean. Ocean at bottom.

    Observed positions are ANCHORED, not raw-converted: the offset pipeline
    compresses the 1967 planform SHAPE by ~10x when it initializes the model, but
    the model's shoreline CHANGE is in true meters. So a single linear ORIG_LEN->
    x_s map can't honor both (it would shrink the observed change 10x). Instead we
    anchor: observed_position(year) = (model 1967 line) + (observed change
    1967->year, in true meters). The model 1967 line and the observed change are
    both in real meters, so the result is correct. Verified: the observed change
    here matches the change-vs-observed figure exactly."""
    gis = _gis_axis()
    fig, ax = plt.subplots(figsize=(13, 6), constrained_layout=True)

    first_run = next(iter(runs_positions))
    start0, _ = runs_positions[first_run]      # model 1967, flipped (seaward+), m
    ref_mean = np.nanmean(start0)
    ref_planform = start0 - ref_mean           # model 1967 in the plot frame

    n = len(runs_positions)
    cmap = plt.cm.Blues(np.linspace(0.40, 0.95, n))
    for c, (run_name, (start, end)) in zip(cmap, runs_positions.items()):
        mval = _m_value(run_name)
        lbl = f"M={mval:.0f}" if not np.isnan(mval) else run_name
        ax.plot(gis, end - ref_mean, marker="o", ms=4, lw=2.2, color=c, label=lbl)

    # Observed reference positions = model 1967 + observed change (true meters).
    obs_67 = _obs_orig_len(START_YEAR)
    # solid, then dashed, then dash-dot for successive target years (if >1)
    obs_target_styles = ["-", (0, (6, 2)), (0, (3, 1, 1, 1))]
    t = 0
    for year in OBS_REFERENCE_YEARS:
        if year == START_YEAR:
            # start reference: quiet grey dashed line
            ax.plot(gis, ref_planform, marker="o", ms=4, lw=1.6,
                    color=OBS_START_COLOR, ls="--",
                    label=f"Observed {year} (start ref)", zorder=3)
            continue
        obs_yr = _obs_orig_len(year)
        if obs_yr is None or obs_67 is None:
            print(f"  [positions] observed file for {year} not found; skipping.")
            continue
        # observed change in TRUE meters, + = seaward (flip landward-positive raw)
        obs_change_m = -(obs_yr - obs_67)
        obs_position = ref_planform + obs_change_m
        ax.plot(gis, obs_position, marker="s", ms=OBS_TARGET_MS,
                lw=OBS_TARGET_LW, color=OBS_TARGET_COLOR,
                ls=obs_target_styles[t % len(obs_target_styles)],
                markeredgecolor="white", markeredgewidth=1.1,
                label=f"Observed {year} (target)", zorder=10)
        t += 1

    ax.axvspan(FIRST_FILE_NUMBER - 0.5, GROIN_BOUNDARY_GIS, alpha=0.06, color="firebrick")
    ax.axvspan(GROIN_BOUNDARY_GIS, LAST_FILE_NUMBER + 0.5, alpha=0.06, color="seagreen")
    ax.axvline(GROIN_BOUNDARY_GIS, color=GROIN_COLOR, ls="--", lw=1.6)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}\u2013D{LAST_FILE_NUMBER})")
    ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)  [landward \u25b2]")
    ax.set_title(f"Final ({END_YEAR}) shoreline POSITION by M   "
                 f"vs observed reference positions")
    ax.set_ylim(ax.get_ylim()[::-1])   # ocean at bottom
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, title="Trapping rate", loc="upper right", framealpha=0.95)
    return fig, "positions"


# ============================== MAIN ==============================

def main():
    print("=" * 70)
    print("Groin M-sweep comparison")
    print("=" * 70)

    runs_change = {}
    for run_name in RUNS:
        try:
            runs_change[run_name] = _model_change(run_name)
            print(f"  loaded {run_name}  (M={_m_value(run_name):.0f})")
        except FileNotFoundError as e:
            print(f"  [SKIP] matrix not found: {e}")
    if not runs_change:
        print("No runs loaded; check RUNS and OUTPUT_BASE_DIR.")
        return

    obs = _observed_change()
    if obs is None:
        print("  [warn] observed raw files not found; overlay will omit target.")

    # Rates (m/yr): model change / run years; observed change / observed period.
    period_years = END_YEAR - START_YEAR
    runs_rate = {rn: ch / _n_years(rn) for rn, ch in runs_change.items()}
    obs_rate = (obs / period_years) if obs is not None else None

    os.makedirs(OUT_DIR, exist_ok=True)
    figA, tagA = fig_overlay(runs_rate, obs_rate, period_years)
    outA = os.path.join(OUT_DIR, f"{OUT_BASENAME}_{tagA}.png")
    figA.savefig(outA, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"\nSaved: {outA}")

    # position comparison (final planform by M)
    runs_positions = {}
    for run_name in runs_change:      # same runs that loaded successfully
        runs_positions[run_name] = _model_positions(run_name)
    figP, tagP = fig_positions(runs_positions, obs_change=obs)
    outP = os.path.join(OUT_DIR, f"{OUT_BASENAME}_{tagP}.png")
    figP.savefig(outP, dpi=200, bbox_inches="tight", facecolor="white")
    print(f"Saved: {outP}")

    figB, tagB = fig_skill_vs_m(runs_change, obs)
    if figB is not None:
        outB = os.path.join(OUT_DIR, f"{OUT_BASENAME}_{tagB}.png")
        figB.savefig(outB, dpi=200, bbox_inches="tight", facecolor="white")
        print(f"Saved: {outB}")

        # also print the skill table
        mask = _updrift_mask()
        print(f"\nUpdrift skill (D{VALIDATE_LO_GIS}-D{VALIDATE_HI_GIS}):")
        print(f"  {'M':>5} {'RMSE':>8} {'bias':>8}   run")
        rows = []
        for run_name, change in runs_change.items():
            resid = change[mask] - obs[mask]
            rows.append((_m_value(run_name),
                         float(np.sqrt(np.nanmean(resid**2))),
                         float(np.nanmean(resid)), run_name))
        for mval, rmse, bias, rn in sorted(rows, key=lambda r: (np.isnan(r[0]), r[0])):
            print(f"  {mval:>5.0f} {rmse:>8.1f} {bias:>+8.1f}   {rn}")


if __name__ == "__main__":
    main()
