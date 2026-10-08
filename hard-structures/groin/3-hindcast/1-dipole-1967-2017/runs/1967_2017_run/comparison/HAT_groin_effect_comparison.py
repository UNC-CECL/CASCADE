"""
The groin's effect at each checkpoint year: a no-groin and a groin 1967-2017 run against the observed wet/dry shoreline.

    python HAT_groin_effect_comparison.py

Reads both runs' *_shoreline_matrix.npy (RUN_NO_GROIN, RUN_GROIN) and the
wet/dry change table; writes an overview, one figure per CHECKPOINT_YEARS
(each also as a _v2 without title and footer) and a shoreline GIF to
COMPARISON_OUTPUT_DIR. Every comparison is change from each series' own 1967.
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
import matplotlib.ticker as mticker
import matplotlib.animation as animation


# --- CONFIG ------------------------------------------------------------------
# Repo root, found by searching upward (ORGANIZATION.md rule 5)
PROJECT_BASE_DIR = str(next(p for p in Path(__file__).resolve().parents
                            if (p / "pyproject.toml").exists()))

# Where the two runs' shoreline matrices live (the runner's OUTPUT_BASE_DIR)
RUN_DATA_DIR = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs")

# The two runs compared
RUN_NO_GROIN = "HAT_1967_2018_M60_deterioration_no_groin"
RUN_GROIN    = "HAT_1967_2018_M60_deterioration_groin"

# Name each comparison's subfolder so comparisons never overwrite each other
COMPARISON_SUBFOLDER = "deterioration_1995_2003_Mover3"   # <- edit per comparison
COMPARISON_OUTPUT_DIR = os.path.join(
    PROJECT_BASE_DIR, "hard-structures", "groin", "3-hindcast", "1-dipole-1967-2017", "results",
    "comparison", COMPARISON_SUBFOLDER,
)

# Geometry, must match the run
NUM_REAL_DOMAINS   = 11
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 2
LAST_FILE_NUMBER   = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1     # 12
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS

# Row 0 of every matrix; checkpoint years need a matching wetdry_<year> column
START_YEAR = 1967
CHECKPOINT_YEARS = [1997, 2017]

# Per-checkpoint style; colour stays by series type (observed black, no groin orange, groin red)
CHECKPOINT_STYLE = {
    1997: dict(alpha=0.55, marker="s", ls="--"),
    2017: dict(alpha=1.00, marker="D", ls="-"),
}
DEFAULT_CHECKPOINT_STYLE = dict(alpha=0.8, marker="o", ls="-.")

# Sign convention, as the main hindcast and HAT_plot_groin_runs.py
FLIP_SIGN_MODEL = True

# Styling
MODEL_COLOR        = "#FF8C00"   # warm orange
GROIN_COLOR        = "#B71C1C"   # groin red
GROIN_BOUNDARY_GIS = 5.5         # D5/D6 interface (Buxton groin)
DOMAIN_TICK_STEP   = 2
DOMAIN_SPACING_M   = 500.0       # alongshore domain width (m) -- project convention
OCEAN_AT_BOTTOM    = True        # seaward plots downward
AXIS_LABEL_FONTSIZE = 12         # shared across every figure variant

# Observed: the wet/dry-referenced change table, never a dune-line column (README)
WETDRY_CHANGE_TABLE = os.path.join(
    PROJECT_BASE_DIR, "hard-structures", "groin", "1-observations",
    "wetdry_photo_positions",
    "Change_from_wetdry_1967_D2_D12.csv",
)
WETDRY_DOMAIN_COL = "Domain_ID"

FIGURE_DPI   = 300     # publication-quality raster
SHOW_FIGURE  = True

GIF_FPS = 4       # frames per second in the shoreline-evolution GIF
GIF_DPI = 100     # GIF resolution -- kept lower than FIGURE_DPI to control file size
# -----------------------------------------------------------------------------


# GIS ids of the real domains
def _gis_axis():
    return np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)


# The real-domain slice of a padded row
def _real_slice(arr_1d):
    return arr_1d[START_REAL_INDEX:END_REAL_INDEX]


# Apply FLIP_SIGN_MODEL exactly as the main hindcast and plotter do
def _flip(v):
    return v * (-1.0 if FLIP_SIGN_MODEL else 1.0)


# Calendar year -> matrix row; None, with a warning, if the run does not reach it
def _year_to_row(year, nt):
    row = year - START_YEAR
    if not (0 <= row < nt):
        print(f"  WARNING: year {year} (row {row}) is outside this run's "
              f"{nt} modeled years (rows 0-{nt - 1}) -- checkpoint skipped.")
        return None
    return row


# A run's shoreline matrix (nt, ndomain), metres
def _load_shoreline(run_name):
    path = os.path.join(RUN_DATA_DIR, run_name, f"{run_name}_shoreline_matrix.npy")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"shoreline matrix not found:\n  {path}")
    m = np.load(path)
    print(f"  Loaded {run_name}: shape {m.shape}")
    return m


# Observed change 1967 -> each checkpoint year, raw sign (+ = landward); None where missing
def _load_observed_changes():
    gis = _gis_axis()

    if not os.path.isfile(WETDRY_CHANGE_TABLE):
        print(f"  [observed] MISSING wet/dry change table -- ALL observed "
              f"checkpoints will be omitted:\n    {WETDRY_CHANGE_TABLE}")
        return {y: None for y in CHECKPOINT_YEARS}

    df = pd.read_csv(WETDRY_CHANGE_TABLE)
    df = df.set_index(WETDRY_DOMAIN_COL)

    changes = {}
    for year in CHECKPOINT_YEARS:
        col = f"change_from_wetdry_1967_wetdry_{year}_m"
        if col not in df.columns:
            print(f"  [observed] MISSING column '{col}' -- checkpoint omitted. "
                  f"(Available wetdry columns: "
                  f"{[c for c in df.columns if 'wetdry' in c]})")
            changes[year] = None
            continue

        change = np.array([df[col].get(d, np.nan) for d in gis])
        n_ok = int(np.isfinite(change).sum())
        print(f"  [observed] {START_YEAR}->{year}: {n_ok}/{len(gis)} domains "
              f"(D{FIRST_FILE_NUMBER}-D{LAST_FILE_NUMBER}) have data.")
        if n_ok < len(gis):
            missing = [int(d) for d, v in zip(gis, change) if not np.isfinite(v)]
            print(f"    WARNING: no data for domain(s) {missing} -- "
                  f"gap(s) will show as a break in the observed line.")
        changes[year] = change

    return changes


# The groin line, labelled near the top of the final (inverted) axis (README)
def _mark_groin(ax, near_top_frac=0.08):
    ax.axvline(GROIN_BOUNDARY_GIS, color=GROIN_COLOR, lw=1.5, ls="--",
               alpha=0.9, zorder=5)
    yl = ax.get_ylim()   # pre-inversion auto limits at call time
    # With OCEAN_AT_BOTTOM the final top is the current bottom, so a small fraction
    frac = near_top_frac if OCEAN_AT_BOTTOM else (1.0 - near_top_frac)
    y = yl[0] + frac * (yl[1] - yl[0])
    ax.text(GROIN_BOUNDARY_GIS + 0.15, y, "Buxton groin", color=GROIN_COLOR,
            fontsize=8, rotation=90, va="top", ha="left", alpha=0.9, zorder=6)


# Light shading: updrift (D6+) vs downdrift (D5 and south, not validated)
def _updrift_downdrift_shading(ax):
    ax.axvspan(FIRST_FILE_NUMBER - 0.5, GROIN_BOUNDARY_GIS,
               alpha=0.06, color="firebrick", zorder=0)   # downdrift
    ax.axvspan(GROIN_BOUNDARY_GIS, LAST_FILE_NUMBER + 0.5,
               alpha=0.06, color="seagreen", zorder=0)     # updrift


# GIS domain -> km alongshore from the window's left edge, for the top axis
def _dom_to_km(x):
    return (np.asarray(x, dtype=float) - FIRST_FILE_NUMBER) * (DOMAIN_SPACING_M / 1000.0)


# km -> GIS domain, the inverse secondary_xaxis needs
def _km_to_dom(x_km):
    return np.asarray(x_km, dtype=float) / (DOMAIN_SPACING_M / 1000.0) + FIRST_FILE_NUMBER


# All checkpoints on one frame: 1967, each year's observed target and both modelled shorelines
def fig_combined_overview(no_groin_m, groin_m, no_groin_name, groin_name,
                           observed_changes, show_title=True, figsize=(12, 7)):
    gis = _gis_axis()

    pos0 = _flip(_real_slice(no_groin_m[0]))
    ref_mean = np.nanmean(pos0)
    planform_1967 = pos0 - ref_mean
    nt = no_groin_m.shape[0]

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    ax.plot(gis, planform_1967, color="0.45", ls="--", lw=2.0, marker="o", ms=5,
            label=f"{START_YEAR} shoreline (observed, model-initialized)", zorder=3)

    for year in CHECKPOINT_YEARS:
        style = CHECKPOINT_STYLE.get(year, DEFAULT_CHECKPOINT_STYLE)
        row = _year_to_row(year, nt)
        observed_change = observed_changes.get(year)

        if observed_change is not None:
            observed_target = planform_1967 - observed_change
            ax.plot(gis, observed_target, color="black", ls=style["ls"], lw=2.2,
                    marker=style["marker"], ms=6, alpha=style["alpha"],
                    label=f"{year} shoreline (observed, target)", zorder=6)
        else:
            print(f"  [overview] no observed data for {year} -- omitted from overview.")

        if row is not None:
            end_no_groin = _flip(_real_slice(no_groin_m[row])) - ref_mean
            end_groin    = _flip(_real_slice(groin_m[row]))    - ref_mean
            ax.plot(gis, end_no_groin, color=MODEL_COLOR, ls=style["ls"], lw=2.0,
                    marker=style["marker"], ms=5, alpha=style["alpha"] * 0.9,
                    label=f"{year} modeled — no groin", zorder=4)
            ax.plot(gis, end_groin, color=GROIN_COLOR, ls=style["ls"], lw=2.2,
                    marker=style["marker"], ms=5, alpha=style["alpha"],
                    label=f"{year} modeled — with groin", zorder=5)
        else:
            print(f"  [overview] no modeled data for {year} -- omitted from overview.")

    _updrift_downdrift_shading(ax)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}–D{LAST_FILE_NUMBER})",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

    secax = ax.secondary_xaxis("top", functions=(_dom_to_km, _km_to_dom))
    secax.set_xlabel(f"Alongshore distance from D{FIRST_FILE_NUMBER} (km)",
                      fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    secax.xaxis.set_major_locator(mticker.MultipleLocator(1.0))

    ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    if show_title:
        year_span = f"{START_YEAR}–{CHECKPOINT_YEARS[-1]}" if CHECKPOINT_YEARS else f"{START_YEAR}"
        ax.set_title(
            f"Effect of Groin Deterioration on Modeled Shoreline Position — All Checkpoints\n"
            f"Buxton Groin Field, Hatteras Island, NC  |  {year_span}",
            fontsize=13, fontweight="bold"
        )
    if OCEAN_AT_BOTTOM:
        ax.set_ylim(ax.get_ylim()[::-1])
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7.5, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              framealpha=0.95, ncol=1)
    ax.spines[["top", "right"]].set_visible(False)
    if show_title:
        ax.annotate(
            f"Model: CASCADE (Barrier3D + BRIE)  |  Obs: digitized dune-line offsets "
            f"({START_YEAR}, {', '.join(str(y) for y in CHECKPOINT_YEARS)})  |  "
            f"no-groin: {no_groin_name}  |  groin: {groin_name}",
            xy=(0, 0), xycoords="axes fraction", xytext=(0, -0.16),
            textcoords="axes fraction", fontsize=7.5, color="#666666",
            ha="left", va="top", style="italic", annotation_clip=False,
        )

    tag = "groin_effect_overview" if show_title else "groin_effect_overview_v2"
    return fig, tag


# One checkpoint year: 1967, the observed target, and both modelled shorelines; show_title=False is the _v2 (README)
def fig_checkpoint(no_groin_m, groin_m, no_groin_name, groin_name,
                    checkpoint_year, observed_change,
                    show_title=True, figsize=(11, 6)):
    gis = _gis_axis()

    # Shared 1967 reference: both runs share the same initial planform
    pos0 = _flip(_real_slice(no_groin_m[0]))
    ref_mean = np.nanmean(pos0)
    planform_1967 = pos0 - ref_mean

    nt = no_groin_m.shape[0]
    row = _year_to_row(checkpoint_year, nt)
    have_model_row = row is not None

    end_no_groin = _flip(_real_slice(no_groin_m[row])) - ref_mean if have_model_row else None
    end_groin    = _flip(_real_slice(groin_m[row]))    - ref_mean if have_model_row else None

    # Observed target: 1967 reference minus the landward-positive observed change
    have_obs = observed_change is not None
    observed_target = (planform_1967 - observed_change) if have_obs else None

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    ax.plot(gis, planform_1967, color="0.45", ls="--", lw=2.0, marker="o", ms=5,
            label=f"{START_YEAR} shoreline (observed, model-initialized)", zorder=3)

    if have_obs:
        ax.plot(gis, observed_target, color="black", ls="-", lw=2.4, marker="o", ms=6,
                label=f"{checkpoint_year} shoreline (observed, target)", zorder=6)
    else:
        ax.text(0.5, 0.92, f"observed {checkpoint_year} raw file not found -- omitted",
                transform=ax.transAxes, ha="center", color="firebrick", fontsize=9)

    if have_model_row:
        ax.plot(gis, end_no_groin, color=MODEL_COLOR, ls="-", lw=2.2, marker="D", ms=5,
                alpha=0.85, label=f"{checkpoint_year} modeled — no groin", zorder=4)
        ax.plot(gis, end_groin, color=GROIN_COLOR, ls="-", lw=2.4, marker="D", ms=5,
                label=f"{checkpoint_year} modeled — with groin", zorder=5)
    else:
        ax.text(0.5, 0.85, f"modeled {checkpoint_year} outside this run's simulated years",
                transform=ax.transAxes, ha="center", color="firebrick", fontsize=9)

    _updrift_downdrift_shading(ax)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}–D{LAST_FILE_NUMBER})",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

    # Top axis: km alongshore from the window's left edge
    secax = ax.secondary_xaxis("top", functions=(_dom_to_km, _km_to_dom))
    secax.set_xlabel(f"Alongshore distance from D{FIRST_FILE_NUMBER} (km)",
                      fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    secax.xaxis.set_major_locator(mticker.MultipleLocator(1.0))

    ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    if show_title:
        ax.set_title(
            f"Effect of Simplified Groin Parameterization on Modeled Shoreline Position\n"
            f"Buxton Groin Field, Hatteras Island, NC  |  {START_YEAR}–{checkpoint_year}",
            fontsize=13, fontweight="bold"
        )
    if OCEAN_AT_BOTTOM:
        ax.set_ylim(ax.get_ylim()[::-1])
    ax.grid(alpha=0.25)
    # Legend's left edge at the 3 km mark, clear of the lines further right
    legend_start_domain = _km_to_dom(3.0)
    xlim = ax.get_xlim()
    legend_x_frac = (legend_start_domain - xlim[0]) / (xlim[1] - xlim[0])
    ax.legend(fontsize=9, loc="upper left", bbox_to_anchor=(legend_x_frac, 0.98),
              framealpha=0.95)
    ax.spines[["top", "right"]].set_visible(False)
    if show_title:
        ax.annotate(
            f"Model: CASCADE (Barrier3D + BRIE)  |  Obs: digitized dune-line offsets "
            f"({START_YEAR}, {checkpoint_year})  |  no-groin: {no_groin_name}  |  "
            f"groin: {groin_name}",
            xy=(0, 0), xycoords="axes fraction", xytext=(0, -0.16),
            textcoords="axes fraction", fontsize=7.5, color="#666666",
            ha="left", va="top", style="italic", annotation_clip=False,
        )

    tag = f"groin_effect_{checkpoint_year}" if show_title else f"groin_effect_{checkpoint_year}_v2"
    return fig, tag


# GIF: both modelled shorelines year by year over static 1967 and observed checkpoint lines
def make_shoreline_evolution_gif(no_groin_m, groin_m, no_groin_name, groin_name,
                                  observed_changes, out_path):
    gis = _gis_axis()

    # Shared 1967 reference (identical convention to the static figures).
    pos0 = _flip(_real_slice(no_groin_m[0]))
    ref_mean = np.nanmean(pos0)
    planform_1967 = pos0 - ref_mean

    nt = no_groin_m.shape[0]
    no_groin_traj = np.array([_flip(_real_slice(no_groin_m[t])) - ref_mean
                               for t in range(nt)])
    groin_traj    = np.array([_flip(_real_slice(groin_m[t])) - ref_mean
                               for t in range(nt)])

    # Year per frame: row t = START_YEAR + t, exact integers
    years = START_YEAR + np.arange(nt)

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)

    # Static references, the same in every frame
    ax.plot(gis, planform_1967, color="0.6", ls="--", lw=1.2, marker="o", ms=3,
            label=f"{START_YEAR} shoreline (reference)", zorder=2)

    all_static_vals = [planform_1967]
    for year in CHECKPOINT_YEARS:
        style = CHECKPOINT_STYLE.get(year, DEFAULT_CHECKPOINT_STYLE)
        observed_change = observed_changes.get(year)
        if observed_change is None:
            print(f"  [GIF] no observed data for {year} -- reference line omitted.")
            continue
        observed_target = planform_1967 - observed_change
        ax.plot(gis, observed_target, color="black", ls="--", lw=1.2,
                marker=style["marker"], ms=3, alpha=0.65,
                label=f"{year} shoreline (observed, target)", zorder=2)
        all_static_vals.append(observed_target)

    # Animated lines, ydata rewritten each frame
    line_no_groin, = ax.plot(gis, no_groin_traj[0], color=MODEL_COLOR, ls="-", lw=2.2,
                              marker="D", ms=5, alpha=0.9,
                              label="Modeled — no groin", zorder=4)
    line_groin, = ax.plot(gis, groin_traj[0], color=GROIN_COLOR, ls="-", lw=2.4,
                           marker="D", ms=5, label="Modeled — with groin", zorder=5)

    _updrift_downdrift_shading(ax)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}–D{LAST_FILE_NUMBER})",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

    secax = ax.secondary_xaxis("top", functions=(_dom_to_km, _km_to_dom))
    secax.set_xlabel(f"Alongshore distance from D{FIRST_FILE_NUMBER} (km)",
                      fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    secax.xaxis.set_major_locator(mticker.MultipleLocator(1.0))

    ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

    # y-limits fixed for the whole animation from every series, so the axis never rescales
    all_vals = np.concatenate(
        all_static_vals + [no_groin_traj.ravel(), groin_traj.ravel()]
    )
    pad = 0.08 * (np.nanmax(all_vals) - np.nanmin(all_vals) + 1e-9)
    ax.set_ylim(np.nanmin(all_vals) - pad, np.nanmax(all_vals) + pad)
    if OCEAN_AT_BOTTOM:
        ax.set_ylim(ax.get_ylim()[::-1])

    ax.grid(alpha=0.25)
    legend_start_domain = _km_to_dom(3.0)
    xlim = ax.get_xlim()
    legend_x_frac = (legend_start_domain - xlim[0]) / (xlim[1] - xlim[0])
    ax.legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(legend_x_frac, 0.98),
              framealpha=0.95)
    ax.spines[["top", "right"]].set_visible(False)

    year_text = ax.text(0.02, 0.03, "", transform=ax.transAxes, fontsize=12,
                         fontweight="bold", va="bottom", ha="left", zorder=10)

    # One frame: move both modelled lines and the year label
    def update(frame):
        line_no_groin.set_ydata(no_groin_traj[frame])
        line_groin.set_ydata(groin_traj[frame])
        year_text.set_text(f"Year: {years[frame]:.0f}")
        return line_no_groin, line_groin, year_text

    anim = animation.FuncAnimation(fig, update, frames=nt, blit=False)
    anim.save(out_path, writer=animation.PillowWriter(fps=GIF_FPS), dpi=GIF_DPI)
    plt.close(fig)
    print(f"  Saved GIF: {out_path}")


# Run: load both runs and the observed checkpoints, write the overview, per-year figures and GIF
def main():
    print("=" * 78)
    print("GROIN-EFFECT COMPARISON FIGURE")
    print(f"  no_groin run: {RUN_NO_GROIN}")
    print(f"  groin run:    {RUN_GROIN}")
    print(f"  Saving to:    {COMPARISON_OUTPUT_DIR}")
    print("=" * 78)

    no_groin_m = _load_shoreline(RUN_NO_GROIN)
    groin_m    = _load_shoreline(RUN_GROIN)

    os.makedirs(COMPARISON_OUTPUT_DIR, exist_ok=True)

    print(f"\nLoading observed checkpoints ({START_YEAR} -> {CHECKPOINT_YEARS})...")
    observed_changes = _load_observed_changes()

    # Overview: every checkpoint on one figure
    fig_ov1, tag_ov1 = fig_combined_overview(no_groin_m, groin_m, RUN_NO_GROIN, RUN_GROIN,
                                              observed_changes, show_title=True)
    out_ov1 = os.path.join(COMPARISON_OUTPUT_DIR, f"{tag_ov1}.png")
    fig_ov1.savefig(out_ov1, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {out_ov1}")

    fig_ov2, tag_ov2 = fig_combined_overview(no_groin_m, groin_m, RUN_NO_GROIN, RUN_GROIN,
                                              observed_changes, show_title=False,
                                              figsize=(12, 6))
    out_ov2 = os.path.join(COMPARISON_OUTPUT_DIR, f"{tag_ov2}.png")
    fig_ov2.savefig(out_ov2, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {out_ov2}")

    # One figure per checkpoint year
    for year in CHECKPOINT_YEARS:
        observed_change = observed_changes.get(year)

        fig1, tag1 = fig_checkpoint(no_groin_m, groin_m, RUN_NO_GROIN, RUN_GROIN,
                                     year, observed_change, show_title=True)
        out1 = os.path.join(COMPARISON_OUTPUT_DIR, f"{tag1}.png")
        fig1.savefig(out1, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out1}")

        # v2: no title, subtitle or footer, for a captioned figure
        fig2, tag2 = fig_checkpoint(no_groin_m, groin_m, RUN_NO_GROIN, RUN_GROIN,
                                     year, observed_change, show_title=False,
                                     figsize=(11, 5))
        out2 = os.path.join(COMPARISON_OUTPUT_DIR, f"{tag2}.png")
        fig2.savefig(out2, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out2}")

    # GIF: shoreline evolution, both runs, observed checkpoints as static references
    gif_out = os.path.join(COMPARISON_OUTPUT_DIR, "shoreline_evolution.gif")
    make_shoreline_evolution_gif(no_groin_m, groin_m, RUN_NO_GROIN, RUN_GROIN,
                                  observed_changes, gif_out)

    if SHOW_FIGURE:
        plt.show()


if __name__ == "__main__":
    main()
