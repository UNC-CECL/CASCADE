"""
Grid search over the groin's M and deterioration fraction on the 1967-2017 rig, scored against the 2018 wet/dry shoreline.

    python HAT_groin_sensitivity_sweep.py

A coarse grid, then a 3 x 3 fine grid around its best cell; every cell runs
in its own process (HAT_groin_sweep_single_combo.py) and is appended to the
results CSV at once, so a rerun resumes. Writes the CSV, a heatmap and two
profile figures to HAT-buxton-hindcast-groin-test/sensitivity_sweep/.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

import os
import sys
import re
import subprocess
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import HAT_groin_hindcast_1967_2017 as hc


# --- CONFIG ------------------------------------------------------------------
# Stage 1: coarse grid, wide and sparse (illustrative ranges, not validated)
COARSE_M_VALUES        = [20, 40, 60, 80, 100, 120]
COARSE_FRACTION_VALUES = [0.1, 0.3, 0.5, 0.7, 0.9]

# Stage 2: fine grid around the coarse best; half-width == step gives 3 x 3 (README)
FINE_M_HALF_WIDTH        = 10     # +/- around coarse best M
FINE_M_STEP              = 10
FINE_FRACTION_HALF_WIDTH = 0.1
FINE_FRACTION_STEP       = 0.1

# Validation target; must match HAT_groin_sweep_single_combo.py
MODEL_FIT_YEAR    = 2017   # model's TRUE final simulated year
OBSERVED_FIT_YEAR = 2018   # wet/dry column used for the observed target
FIT_DOMAINS_GIS = list(range(2, 13))   # D2-D12, full range

WETDRY_CHANGE_TABLE = os.path.join(
    hc.PROJECT_BASE_DIR, "hard-structures", "groin", "HAT-groin-buxton-output",
    "shoreline_position_output",
    "Change_from_wetdry_1967_D2_D12.csv",
)

OUTPUT_DIR = os.path.join(
    hc.PROJECT_BASE_DIR, "hard-structures", "groin", "HAT-buxton-hindcast-groin-test",
    "sensitivity_sweep",
)
RESULTS_CSV = os.path.join(OUTPUT_DIR, "HAT_groin_sweep_results.csv")
PROFILE_DIR = os.path.join(OUTPUT_DIR, "profiles")   # must match the worker script
WORKER_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "HAT_groin_sweep_single_combo.py")

# Styling, as HAT_groin_effect_comparison.py
MODEL_COLOR = "#FF8C00"
GROIN_COLOR = "#B71C1C"
GROIN_BOUNDARY_GIS = 5.5
DOMAIN_TICK_STEP = 2
# -----------------------------------------------------------------------------


# Observed 2018 change per fit domain, for printing only (each worker loads its own)
def load_observed_target():
    if not os.path.isfile(WETDRY_CHANGE_TABLE):
        raise FileNotFoundError(f"Wet/dry change table not found:\n  {WETDRY_CHANGE_TABLE}")
    df = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
    col = f"change_from_wetdry_1967_wetdry_{OBSERVED_FIT_YEAR}_m"
    if col not in df.columns:
        raise KeyError(f"'{col}' not found in change table. Available: "
                        f"{[c for c in df.columns if 'wetdry' in c]}")
    return np.array([df[col].get(d, np.nan) for d in FIT_DOMAINS_GIS])


# Cells already done in an earlier (possibly crashed) sweep; typed empty frame otherwise (README)
def load_existing_results():
    if os.path.isfile(RESULTS_CSV):
        return pd.read_csv(RESULTS_CSV)
    return pd.DataFrame({
        "M": pd.Series(dtype="float64"),
        "fraction": pd.Series(dtype="float64"),
        "rmse": pd.Series(dtype="float64"),
        "stage": pd.Series(dtype="object"),
    })


# Append one cell's result and rewrite the CSV at once
def append_result(results_df, M, fraction, rmse_val, stage):
    new_row = pd.DataFrame([dict(M=M, fraction=fraction, rmse=rmse_val, stage=stage)])
    updated = pd.concat([results_df, new_row], ignore_index=True)
    updated.to_csv(RESULTS_CSV, index=False)
    return updated


# Run one (M, fraction) cell in a fresh process; (rmse, ok), ok False on any failure
def run_one_combo_subprocess(M, fraction):
    try:
        result = subprocess.run(
            [sys.executable, WORKER_SCRIPT, str(M), str(fraction)],
            capture_output=True, text=True,
        )
    except Exception as e:
        print(f"    [FAILED] could not launch subprocess: {e}")
        return np.nan, False

    if result.returncode != 0:
        print(f"    [FAILED] worker exited with code {result.returncode} "
              f"(this is how a crash shows up here -- the orchestrator "
              f"survives it and moves on to the next combination)")
        tail = result.stderr.strip().splitlines()[-15:]
        if tail:
            print("    ---- worker stderr (last 15 lines) ----")
            for line in tail:
                print(f"    {line}")
        return np.nan, False

    match = re.search(r"RESULT_RMSE=([\w.+\-]+)", result.stdout)
    if not match:
        print(f"    [FAILED] worker exited cleanly but printed no RESULT_RMSE "
              f"line -- check its stdout/stderr manually if this recurs.")
        return np.nan, False

    try:
        return float(match.group(1)), True
    except ValueError:
        return np.nan, False


# Every cell of one grid, skipping cells already done
def run_grid(m_values, fraction_values, stage_label, results_df):
    combos = list(itertools.product(m_values, fraction_values))
    print(f"\n{'=' * 70}\nSTAGE: {stage_label} ({len(m_values)} x "
          f"{len(fraction_values)} = {len(combos)} runs)\n{'=' * 70}")

    for i, (M, frac) in enumerate(combos):
        already_ok = results_df[
            (results_df["M"] == M)
            & np.isclose(results_df["fraction"], frac)
            & (results_df["stage"] == stage_label)
            & results_df["rmse"].notna()
        ]
        if len(already_ok) > 0:
            print(f"  [{i + 1}/{len(combos)}] M={M:>4}, fraction={frac:.3f} "
                  f"-> already done (RMSE={already_ok.iloc[0]['rmse']:.2f} m), skipping")
            continue

        print(f"  [{i + 1}/{len(combos)}] M={M:>4}, fraction={frac:.3f} -> running...")
        err, ok = run_one_combo_subprocess(M, frac)
        print(f"    -> RMSE={err} ({'OK' if ok else 'FAILED'})")
        results_df = append_result(results_df, M, frac, err, stage_label)

    return results_df


# The groin line and its label
def _mark_groin(ax):
    ax.axvline(GROIN_BOUNDARY_GIS, color=GROIN_COLOR, lw=1.5, ls="--", alpha=0.9, zorder=5)
    yl = ax.get_ylim()
    ax.text(GROIN_BOUNDARY_GIS, yl[0] + 0.9 * (yl[1] - yl[0]), " Buxton groin",
            color=GROIN_COLOR, fontsize=8, rotation=90, va="top", ha="left", alpha=0.9)


# Light shading: downdrift (D5 and south) vs updrift (D6 and north)
def _updrift_downdrift_shading(ax):
    ax.axvspan(FIT_DOMAINS_GIS[0] - 0.5, GROIN_BOUNDARY_GIS,
               alpha=0.06, color="firebrick", zorder=0)   # downdrift
    ax.axvspan(GROIN_BOUNDARY_GIS, FIT_DOMAINS_GIS[-1] + 0.5,
               alpha=0.06, color="seagreen", zorder=0)     # updrift


# A cell's saved modelled profile, or None if the worker never wrote one
def load_profile(M, fraction):
    path = os.path.join(PROFILE_DIR, f"M{M:g}_frac{fraction:g}.npy")
    if not os.path.isfile(path):
        return None
    return np.load(path)


# Best cell's modelled change against the observed, D2-D12
def fig_best_fit_profile(best, observed):
    modeled = load_profile(best.M, best.fraction)
    if modeled is None:
        print(f"  [figure] No saved profile for best combo M={best.M}, "
              f"fraction={best.fraction} -- skipping best-fit profile figure.")
        return

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    ax.plot(FIT_DOMAINS_GIS, observed, "s--", color="black", lw=2.2, ms=6,
            label=f"Observed {OBSERVED_FIT_YEAR}", zorder=6)
    ax.plot(FIT_DOMAINS_GIS, modeled, "D-", color=MODEL_COLOR, lw=2.2, ms=6,
            label=f"Best fit (M={best.M:g}, frac={best.fraction:.2f}, "
                  f"RMSE={best.rmse:.1f} m)", zorder=5)
    _updrift_downdrift_shading(ax)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIT_DOMAINS_GIS[0], FIT_DOMAINS_GIS[-1] + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIT_DOMAINS_GIS[0]}-D{FIT_DOMAINS_GIS[-1]})")
    ax.set_ylabel(f"Shoreline change {hc.START_YEAR}->{MODEL_FIT_YEAR}/"
                  f"{OBSERVED_FIT_YEAR} (m)  [+ = landward]")
    ax.set_title("Best-fit sweep result vs observed shoreline position")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig_out = os.path.join(OUTPUT_DIR, "HAT_groin_sweep_best_fit_profile.png")
    fig.savefig(fig_out, dpi=200, facecolor="white")
    plt.close(fig)
    print(f"Saved: {fig_out}")


# The top-N cells' profiles over the observed: one winner, or a ridge of equals
def fig_top_n_profiles(all_ok, observed, n=5):
    top = all_ok.sort_values("rmse").head(n)
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)
    ax.plot(FIT_DOMAINS_GIS, observed, "s--", color="black", lw=2.5, ms=7,
            label=f"Observed {OBSERVED_FIT_YEAR}", zorder=10)

    cmap = plt.cm.viridis(np.linspace(0, 0.85, len(top)))
    plotted_any = False
    for color, (_, row) in zip(cmap, top.iterrows()):
        modeled = load_profile(row.M, row.fraction)
        if modeled is None:
            continue
        plotted_any = True
        ax.plot(FIT_DOMAINS_GIS, modeled, "-", color=color, lw=1.8, alpha=0.85,
                label=f"M={row.M:g}, frac={row.fraction:.2f} (RMSE={row.rmse:.1f})")

    if not plotted_any:
        plt.close(fig)
        print("  [figure] No saved profiles found for top-N combos -- skipping.")
        return

    _updrift_downdrift_shading(ax)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIT_DOMAINS_GIS[0], FIT_DOMAINS_GIS[-1] + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIT_DOMAINS_GIS[0]}-D{FIT_DOMAINS_GIS[-1]})")
    ax.set_ylabel(f"Shoreline change {hc.START_YEAR}->{MODEL_FIT_YEAR}/"
                  f"{OBSERVED_FIT_YEAR} (m)  [+ = landward]")
    ax.set_title(f"Top {len(top)} sweep results vs observed shoreline position")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig_out = os.path.join(OUTPUT_DIR, "HAT_groin_sweep_top_n_profiles.png")
    fig.savefig(fig_out, dpi=200, facecolor="white")
    plt.close(fig)
    print(f"Saved: {fig_out}")


# Run: observed target, coarse grid, fine grid, best cell, heatmap and profile figures
def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    observed = load_observed_target()
    n_ok = int(np.isfinite(observed).sum())
    print(f"Observed target ({OBSERVED_FIT_YEAR}, D{FIT_DOMAINS_GIS[0]}-"
          f"D{FIT_DOMAINS_GIS[-1]}): {n_ok}/{len(FIT_DOMAINS_GIS)} domains have data.")
    print(f"Values: {np.round(observed, 1)}")

    results_df = load_existing_results()
    if len(results_df) > 0:
        n_prev_ok = results_df["rmse"].notna().sum()
        print(f"\nResuming from a previous run: {n_prev_ok} combination(s) "
              f"already completed successfully -- these will be skipped.")

    # Stage 1: coarse grid
    results_df = run_grid(COARSE_M_VALUES, COARSE_FRACTION_VALUES, "coarse", results_df)

    coarse_ok = results_df[(results_df["stage"] == "coarse") & results_df["rmse"].notna()]
    if len(coarse_ok) == 0:
        sys.exit("\nNo coarse-grid combination completed successfully -- "
                  "cannot proceed to the fine grid. Check the [FAILED] messages above.")
    best_coarse = coarse_ok.loc[coarse_ok["rmse"].idxmin()]
    print(f"\nBest coarse combo: M={best_coarse.M}, fraction={best_coarse.fraction:.3f}, "
          f"RMSE={best_coarse.rmse:.2f} m")

    # Stage 2: fine grid around the coarse best
    fine_M_values = list(range(
        max(0, int(best_coarse.M - FINE_M_HALF_WIDTH)),
        int(best_coarse.M + FINE_M_HALF_WIDTH) + 1,
        FINE_M_STEP,
    ))
    fine_fraction_values = sorted(set(
        round(f, 3) for f in np.arange(
            max(0.0, best_coarse.fraction - FINE_FRACTION_HALF_WIDTH),
            min(1.0, best_coarse.fraction + FINE_FRACTION_HALF_WIDTH) + 1e-9,
            FINE_FRACTION_STEP,
        )
    ))
    results_df = run_grid(fine_M_values, fine_fraction_values, "fine", results_df)

    fine_ok = results_df[(results_df["stage"] == "fine") & results_df["rmse"].notna()]
    all_ok = results_df[results_df["rmse"].notna()].copy()
    if len(all_ok) == 0:
        sys.exit("\nNo combination completed successfully across either stage.")
    all_ok["M"] = pd.to_numeric(all_ok["M"])
    all_ok["fraction"] = pd.to_numeric(all_ok["fraction"])
    all_ok["rmse"] = pd.to_numeric(all_ok["rmse"])

    best = all_ok.loc[all_ok["rmse"].idxmin()]
    print(f"\n{'=' * 70}\nBEST OVERALL: M={best.M}, fraction={best.fraction:.3f}, "
          f"RMSE={best.rmse:.2f} m ({best.stage} grid)\n{'=' * 70}")
    print(f"Full results (including any failed attempts): {RESULTS_CSV}")

    n_failed = results_df["rmse"].isna().sum()
    if n_failed > 0:
        print(f"\n{n_failed} attempt(s) failed (crashed or errored). Re-run this "
              f"script to retry just those -- already-successful combos are "
              f"skipped automatically.")

    # Heatmap of the fine grid
    if len(fine_ok) > 0:
        fine_ok = fine_ok.copy()
        fine_ok["M"] = pd.to_numeric(fine_ok["M"])
        fine_ok["fraction"] = pd.to_numeric(fine_ok["fraction"])
        fine_ok["rmse"] = pd.to_numeric(fine_ok["rmse"])
        pivot = fine_ok.pivot(index="fraction", columns="M", values="rmse")
        fig, ax = plt.subplots(figsize=(9, 6))
        im = ax.imshow(pivot.values, aspect="auto", origin="lower", cmap="viridis_r",
                        extent=[min(fine_M_values), max(fine_M_values),
                                min(fine_fraction_values), max(fine_fraction_values)])
        ax.set_xlabel("GROIN_TRAPPING_RATE_M_YR (M, m/yr)")
        ax.set_ylabel("GROIN_DETERIORATION_FRACTION")
        ax.set_title(f"RMSE vs observed {OBSERVED_FIT_YEAR} shoreline "
                     f"(D{FIT_DOMAINS_GIS[0]}-D{FIT_DOMAINS_GIS[-1]}, fine grid)")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("RMSE (m)")
        ax.plot(best.M, best.fraction, marker="*", ms=22, color="red",
                markeredgecolor="white", label=f"Best: M={best.M}, frac={best.fraction:.2f}")
        ax.legend()
        fig.tight_layout()
        fig_out = os.path.join(OUTPUT_DIR, "HAT_groin_sweep_heatmap.png")
        fig.savefig(fig_out, dpi=200, facecolor="white")
        print(f"Saved: {fig_out}")
    else:
        print("\nNo fine-grid combination succeeded -- skipping heatmap.")

    # Modelled vs observed profile figures
    fig_best_fit_profile(best, observed)
    fig_top_n_profiles(all_ok, observed, n=5)


if __name__ == "__main__":
    main()
