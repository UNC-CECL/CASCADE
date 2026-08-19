"""
HAT_groin_effect_comparison.py
===============================
Standalone publication-figure script. Compares a 'no_groin' and a 'groin'
1967-1997 hindcast run against the real observed 1967/1997 shoreline,
isolating the groin's effect on modeled shoreline position. Single dedicated
script -- no dependency on HAT_plot_groin_runs.py's full multi-figure suite.

ONE FIGURE, FOUR SERIES, ONE SHARED REFERENCE FRAME
----------------------------------------------------
  - 1967 shoreline    (model's own initial planform -- itself derived from the
                        real 1967 offset file at CASCADE init, so it stands in
                        directly for the observed 1967 shape)
  - 1997 OBSERVED      (real target: 1967->1997 raw-offset CHANGE added onto
                        the shared 1967 reference -- keeps the raw offset
                        file's own coordinate system from ever touching
                        CASCADE's x_s directly; only the delta is borrowed)
  - 1997 modeled, no groin
  - 1997 modeled, with groin

Usage
-----
Set RUN_NO_GROIN / RUN_GROIN below to your two run folder names under
output/raw_runs/. Set COMPARISON_SUBFOLDER to whatever label you want for
THIS comparison (e.g. "M70_install1970") -- different comparisons (different
M, different install year, etc.) each get their own subfolder and never
overwrite each other. Then run.

Saves to:
    scripts/groin/HAT-hindcast-groin-test/comparison/{COMPARISON_SUBFOLDER}/
        groin_effect_final_position.png      (with title/subtitle)
        groin_effect_final_position_v2.png   (same plot, no title/subtitle --
                                               for a captioned dissertation figure)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.animation as animation

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_BASE_DIR = r"/"

# Where the two runs' shoreline matrices live (matches the hindcast runner's
# own OUTPUT_BASE_DIR -- do not edit unless that convention changes).
RUN_DATA_DIR = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs")

# --- The two runs being compared -- edit these to your actual run folder names
RUN_NO_GROIN = "HAT_1967_1997_base_no_groin"
RUN_GROIN    = "HAT_1967_1997_60M_groin"

# --- Where THIS comparison's figure is saved. Name the subfolder yourself so
# different comparisons never overwrite each other. ---
COMPARISON_SUBFOLDER = "base_vs_60M"          # <- edit this per comparison
COMPARISON_OUTPUT_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "groin", "HAT-hindcast-groin-test",
    "comparison", COMPARISON_SUBFOLDER,
)

# --- Geometry (must match the run) ---
NUM_REAL_DOMAINS   = 11
NUM_BUFFER_DOMAINS = 15
FIRST_FILE_NUMBER  = 2
LAST_FILE_NUMBER   = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1     # 12
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS

START_YEAR = 1967
END_YEAR   = 1997

# ── Sign convention (matches main hindcast + HAT_plot_groin_runs.py exactly) ──
FLIP_SIGN_MODEL = True

# --- Styling (from your established figure conventions) ---
MODEL_COLOR        = "#FF8C00"   # warm orange
GROIN_COLOR        = "#B71C1C"   # groin red
GROIN_BOUNDARY_GIS = 5.5         # D5/D6 interface (Buxton groin)
DOMAIN_TICK_STEP   = 2
DOMAIN_SPACING_M   = 500.0       # alongshore domain width (m) -- project convention
OCEAN_AT_BOTTOM    = True        # seaward plots downward
AXIS_LABEL_FONTSIZE = 12         # shared by both the titled and simplified (v2) figure

# --- Observed reference (raw ArcGIS dune-line offset files) ---
RAW_OFFSET_DIR = os.path.join(
    PROJECT_BASE_DIR, "scripts", "groin", "HAT-hindcast-groin-test",
    "groin_init", "2-brie-offset", "raw_offsets",
)
OBS_RAW_START  = os.path.join(RAW_OFFSET_DIR, f"{START_YEAR}_duneline_offset_raw.csv")
OBS_RAW_END    = os.path.join(RAW_OFFSET_DIR, f"{END_YEAR}_duneline_offset_raw.csv")
OBS_DOMAIN_COL = "domain_id"
OBS_POS_COL    = "ORIG_LEN"

FIGURE_DPI   = 300     # publication-quality raster
SHOW_FIGURE  = True

GIF_FPS = 4       # frames per second in the shoreline-evolution GIF
GIF_DPI = 100     # GIF resolution -- kept lower than FIGURE_DPI to control file size

# =============================================================================
# HELPERS
# =============================================================================

def _gis_axis():
    return np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)


def _real_slice(arr_1d):
    return arr_1d[START_REAL_INDEX:END_REAL_INDEX]


def _flip(v):
    """Apply FLIP_SIGN_MODEL exactly as the main hindcast / plotter do."""
    return v * (-1.0 if FLIP_SIGN_MODEL else 1.0)


def _load_shoreline(run_name):
    """Load a run's shoreline matrix (nt, ndomain) in meters."""
    path = os.path.join(RUN_DATA_DIR, run_name, f"{run_name}_shoreline_matrix.npy")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"shoreline matrix not found:\n  {path}")
    m = np.load(path)
    print(f"  Loaded {run_name}: shape {m.shape}")
    return m


def _load_observed_change():
    """Observed per-domain dune-line change 1967->1997 from raw offset files.
    Returns (gis_array, obs_change_raw) where + = landward/erosion (RAW sign),
    or (None, None) if the raw files aren't found."""
    missing = [p for p in (OBS_RAW_START, OBS_RAW_END) if not os.path.isfile(p)]
    if missing:
        print("  [observed] MISSING file(s) -- observed target line will be omitted:")
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
    obs_change_raw = np.array([p1.get(d, np.nan) - p0.get(d, np.nan) for d in gis])

    n_ok = int(np.isfinite(obs_change_raw).sum())
    print(f"  [observed] Loaded {START_YEAR}/{END_YEAR} raw offsets OK -- "
          f"{n_ok}/{len(gis)} domains (D{FIRST_FILE_NUMBER}-D{LAST_FILE_NUMBER}) have "
          f"both years present.")
    if n_ok < len(gis):
        missing_domains = [int(d) for d, v in zip(gis, obs_change_raw) if not np.isfinite(v)]
        print(f"  [observed] WARNING: no data for domain(s) {missing_domains} in one "
              f"or both years -- gap(s) will show as a break in the observed line.")

    return gis, obs_change_raw


def _mark_groin(ax, near_top_frac=0.08):
    """Draw the Buxton groin boundary line + label. The label is placed near
    the top of the FINAL displayed axis (after fig_groin_effect_final_position's
    later ax.set_ylim(ax.get_ylim()[::-1]) inversion when OCEAN_AT_BOTTOM), and
    offset to the right of the dashed line so it never overlaps either the
    line itself or the x-axis."""
    ax.axvline(GROIN_BOUNDARY_GIS, color=GROIN_COLOR, lw=1.5, ls="--",
               alpha=0.9, zorder=5)
    yl = ax.get_ylim()   # pre-inversion auto limits at call time
    # When OCEAN_AT_BOTTOM, the eventual TOP of the chart corresponds to the
    # current (pre-inversion) BOTTOM of the y-range -- so "near top of final"
    # means a SMALL fraction here; otherwise it means a large one.
    frac = near_top_frac if OCEAN_AT_BOTTOM else (1.0 - near_top_frac)
    y = yl[0] + frac * (yl[1] - yl[0])
    ax.text(GROIN_BOUNDARY_GIS + 0.15, y, "Buxton groin", color=GROIN_COLOR,
            fontsize=8, rotation=90, va="top", ha="left", alpha=0.9, zorder=6)


def _updrift_downdrift_shading(ax):
    """Light shading: updrift (D6+) vs downdrift (<=D5, not validated)."""
    ax.axvspan(FIRST_FILE_NUMBER - 0.5, GROIN_BOUNDARY_GIS,
               alpha=0.06, color="firebrick", zorder=0)   # downdrift
    ax.axvspan(GROIN_BOUNDARY_GIS, LAST_FILE_NUMBER + 0.5,
               alpha=0.06, color="seagreen", zorder=0)     # updrift


def _dom_to_km(x):
    """GIS domain ID -> alongshore distance (km) from the left edge of this
    test window (D{FIRST_FILE_NUMBER} = 0 km), using the project's 500 m
    domain-width convention. Used for the secondary top x-axis."""
    return (np.asarray(x, dtype=float) - FIRST_FILE_NUMBER) * (DOMAIN_SPACING_M / 1000.0)


def _km_to_dom(x_km):
    """Inverse of _dom_to_km, required by matplotlib's secondary_xaxis."""
    return np.asarray(x_km, dtype=float) / (DOMAIN_SPACING_M / 1000.0) + FIRST_FILE_NUMBER


# =============================================================================
# FIGURE
# =============================================================================

def fig_groin_effect_final_position(no_groin_m, groin_m, no_groin_name, groin_name,
                                     show_title=True, figsize=(11, 6)):
    """PUBLICATION FIGURE: isolates the groin's effect on final (1997) shoreline
    position against a real observed target, all on one shared reference frame.
    See module docstring for the four series and why they're safe to overlay
    despite the model and raw GIS offset files living in different absolute
    coordinate systems.

    show_title=True  -> full figure with title/subtitle AND the gray italic
                         source-citation footer (for standalone use, e.g. slides).
    show_title=False -> "v2": no title/subtitle and no footer -- just the plot
                         -- for figures that will carry their caption in
                         surrounding text (e.g. a dissertation figure with a
                         numbered caption).
    figsize           -> pass a shorter height for v2, e.g. (11, 5), since it
                         no longer needs vertical room for the title block.
    Both versions share identical axis label font sizes (AXIS_LABEL_FONTSIZE)
    so they drop into the same document/deck without looking mismatched."""
    gis = _gis_axis()

    # Shared 1967 reference: the no-groin run's own initial planform (both runs
    # share identical 1967 init, so either run's [0] gives the same result).
    pos0 = _flip(_real_slice(no_groin_m[0]))
    ref_mean = np.nanmean(pos0)
    planform_1967 = pos0 - ref_mean

    # Modeled 1997 end positions, same shared reference.
    end_no_groin = _flip(_real_slice(no_groin_m[-1])) - ref_mean
    end_groin    = _flip(_real_slice(groin_m[-1]))    - ref_mean

    # Observed 1997 target: real 1967->1997 CHANGE (from raw offset files)
    # added onto the shared 1967 reference -- see module docstring.
    obs_gis, obs_change_raw = _load_observed_change()
    have_obs = obs_gis is not None
    observed_1997 = (planform_1967 - obs_change_raw) if have_obs else None
    # (obs_change_raw is "+ = landward"; subtracting it flips to "+ = seaward",
    # matching the model's _flip() convention, before adding onto the 1967 ref.)

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    ax.plot(gis, planform_1967, color="0.45", ls="--", lw=2.0, marker="o", ms=5,
            label=f"{START_YEAR} shoreline (observed, model-initialized)", zorder=3)

    if have_obs:
        ax.plot(gis, observed_1997, color="black", ls="-", lw=2.4, marker="o", ms=6,
                label=f"{END_YEAR} shoreline (observed, target)", zorder=6)
    else:
        ax.text(0.5, 0.92, "observed 1997 raw file not found -- omitted from figure",
                transform=ax.transAxes, ha="center", color="firebrick", fontsize=9)

    ax.plot(gis, end_no_groin, color=MODEL_COLOR, ls="-", lw=2.2, marker="D", ms=5,
            alpha=0.85, label=f"{END_YEAR} modeled \u2014 no groin", zorder=4)
    ax.plot(gis, end_groin, color=GROIN_COLOR, ls="-", lw=2.4, marker="D", ms=5,
            label=f"{END_YEAR} modeled \u2014 with groin", zorder=5)

    _updrift_downdrift_shading(ax)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}\u2013D{LAST_FILE_NUMBER})",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

    # Secondary top x-axis: alongshore distance in km, sharing the same domain
    # axis (500 m/domain), zeroed at D{FIRST_FILE_NUMBER} (left edge of window).
    secax = ax.secondary_xaxis("top", functions=(_dom_to_km, _km_to_dom))
    secax.set_xlabel(f"Alongshore distance from D{FIRST_FILE_NUMBER} (km)",
                      fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    secax.xaxis.set_major_locator(mticker.MultipleLocator(1.0))

    ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    if show_title:
        ax.set_title(
            f"Effect of Simplified Groin Parameterization on Modeled Shoreline Position\n"
            f"Buxton Groin Field, Hatteras Island, NC  |  {START_YEAR}\u2013{END_YEAR}",
            fontsize=13, fontweight="bold"
        )
    if OCEAN_AT_BOTTOM:
        ax.set_ylim(ax.get_ylim()[::-1])
    ax.grid(alpha=0.25)
    # Anchor the legend's LEFT edge at the 3 km mark on the top axis (i.e. GIS
    # domain _km_to_dom(3.0)) rather than flush against the right border --
    # keeps it clear of the with-groin/no-groin lines further right.
    legend_start_domain = _km_to_dom(3.0)
    xlim = ax.get_xlim()
    legend_x_frac = (legend_start_domain - xlim[0]) / (xlim[1] - xlim[0])
    ax.legend(fontsize=9, loc="upper left", bbox_to_anchor=(legend_x_frac, 0.98),
              framealpha=0.95)
    ax.spines[["top", "right"]].set_visible(False)
    if show_title:
        ax.annotate(
            f"Model: CASCADE (Barrier3D + BRIE)  |  Obs: digitized dune-line offsets "
            f"({START_YEAR}, {END_YEAR})  |  no-groin: {no_groin_name}  |  "
            f"groin: {groin_name}",
            xy=(0, 0), xycoords="axes fraction", xytext=(0, -0.16),
            textcoords="axes fraction", fontsize=7.5, color="#666666",
            ha="left", va="top", style="italic", annotation_clip=False,
        )

    tag = "groin_effect_final_position" if show_title else "groin_effect_final_position_v2"
    return fig, tag


# =============================================================================
# ANIMATED GIF
# =============================================================================

def make_shoreline_evolution_gif(no_groin_m, groin_m, no_groin_name, groin_name, out_path):
    """Animated GIF: shoreline position evolving year-by-year for both the
    no-groin (orange) and with-groin (red) runs, on the same shared
    1967-anchored reference frame as fig_groin_effect_final_position().

    The 1967 initial shoreline and the 1997 OBSERVED target are drawn as thin,
    static dashed reference lines that stay put for the whole animation; only
    the two modeled trajectories move.
    """
    gis = _gis_axis()

    # Shared 1967 reference (identical convention to the static figure).
    pos0 = _flip(_real_slice(no_groin_m[0]))
    ref_mean = np.nanmean(pos0)
    planform_1967 = pos0 - ref_mean

    nt = no_groin_m.shape[0]
    no_groin_traj = np.array([_flip(_real_slice(no_groin_m[t])) - ref_mean
                               for t in range(nt)])
    groin_traj    = np.array([_flip(_real_slice(groin_m[t])) - ref_mean
                               for t in range(nt)])

    # Year label per frame, anchored so frame 0 = START_YEAR and the last
    # frame = END_YEAR regardless of exactly how many rows nt is.
    years = np.linspace(START_YEAR, END_YEAR, nt)

    obs_gis, obs_change_raw = _load_observed_change()
    have_obs = obs_gis is not None
    observed_1997 = (planform_1967 - obs_change_raw) if have_obs else None

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)

    # --- Static thin dashed references (do not move across frames) ---
    ax.plot(gis, planform_1967, color="0.6", ls="--", lw=1.2, marker="o", ms=3,
            label=f"{START_YEAR} shoreline (reference)", zorder=2)
    if have_obs:
        ax.plot(gis, observed_1997, color="black", ls="--", lw=1.2, marker="o", ms=3,
                alpha=0.75, label=f"{END_YEAR} shoreline (observed, target)", zorder=2)

    # --- Animated lines (ydata rewritten each frame in update()) ---
    line_no_groin, = ax.plot(gis, no_groin_traj[0], color=MODEL_COLOR, ls="-", lw=2.2,
                              marker="D", ms=5, alpha=0.9,
                              label="Modeled \u2014 no groin", zorder=4)
    line_groin, = ax.plot(gis, groin_traj[0], color=GROIN_COLOR, ls="-", lw=2.4,
                           marker="D", ms=5, label="Modeled \u2014 with groin", zorder=5)

    _updrift_downdrift_shading(ax)
    _mark_groin(ax)
    ax.set_xticks(np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP))
    ax.set_xlabel(f"GIS Domain ID (D{FIRST_FILE_NUMBER}\u2013D{LAST_FILE_NUMBER})",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

    secax = ax.secondary_xaxis("top", functions=(_dom_to_km, _km_to_dom))
    secax.set_xlabel(f"Alongshore distance from D{FIRST_FILE_NUMBER} (km)",
                      fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    secax.xaxis.set_major_locator(mticker.MultipleLocator(1.0))

    ax.set_ylabel(f"Cross-shore position (m, rel. {START_YEAR} mean)",
                  fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

    # Fix y-limits across the WHOLE animation up front (using the full range
    # of every series, static and animated) so the axis never rescales
    # frame-to-frame -- a rescaling axis is the #1 way animated plots end up
    # looking jittery/misleading.
    all_vals = np.concatenate(
        [planform_1967, no_groin_traj.ravel(), groin_traj.ravel()]
        + ([observed_1997] if have_obs else [])
    )
    pad = 0.08 * (np.nanmax(all_vals) - np.nanmin(all_vals) + 1e-9)
    ax.set_ylim(np.nanmin(all_vals) - pad, np.nanmax(all_vals) + pad)
    if OCEAN_AT_BOTTOM:
        ax.set_ylim(ax.get_ylim()[::-1])

    ax.grid(alpha=0.25)
    legend_start_domain = _km_to_dom(3.0)
    xlim = ax.get_xlim()
    legend_x_frac = (legend_start_domain - xlim[0]) / (xlim[1] - xlim[0])
    ax.legend(fontsize=9, loc="upper left", bbox_to_anchor=(legend_x_frac, 0.98),
              framealpha=0.95)
    ax.spines[["top", "right"]].set_visible(False)

    year_text = ax.text(0.02, 0.03, "", transform=ax.transAxes, fontsize=12,
                         fontweight="bold", va="bottom", ha="left", zorder=10)

    def update(frame):
        line_no_groin.set_ydata(no_groin_traj[frame])
        line_groin.set_ydata(groin_traj[frame])
        year_text.set_text(f"Year: {years[frame]:.0f}")
        return line_no_groin, line_groin, year_text

    anim = animation.FuncAnimation(fig, update, frames=nt, blit=False)
    anim.save(out_path, writer=animation.PillowWriter(fps=GIF_FPS), dpi=GIF_DPI)
    plt.close(fig)
    print(f"  Saved GIF: {out_path}")


# =============================================================================
# MAIN
# =============================================================================

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

    # v1: full figure with title/subtitle (e.g. for slides / standalone use)
    fig1, tag1 = fig_groin_effect_final_position(no_groin_m, groin_m,
                                                  RUN_NO_GROIN, RUN_GROIN,
                                                  show_title=True)
    out1 = os.path.join(COMPARISON_OUTPUT_DIR, f"{tag1}.png")
    fig1.savefig(out1, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {out1}")

    # v2: same plot, no title/subtitle/footer (e.g. for a numbered dissertation
    # figure with its caption set in the surrounding document text) -- a bit
    # shorter since it no longer needs vertical room for the title block.
    fig2, tag2 = fig_groin_effect_final_position(no_groin_m, groin_m,
                                                  RUN_NO_GROIN, RUN_GROIN,
                                                  show_title=False,
                                                  figsize=(11, 5))
    out2 = os.path.join(COMPARISON_OUTPUT_DIR, f"{tag2}.png")
    fig2.savefig(out2, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
    print(f"  Saved: {out2}")

    # Animated GIF: shoreline evolution over time, no-groin vs with-groin
    gif_out = os.path.join(COMPARISON_OUTPUT_DIR, "shoreline_evolution.gif")
    make_shoreline_evolution_gif(no_groin_m, groin_m, RUN_NO_GROIN, RUN_GROIN, gif_out)

    if SHOW_FIGURE:
        plt.show()


if __name__ == "__main__":
    main()
