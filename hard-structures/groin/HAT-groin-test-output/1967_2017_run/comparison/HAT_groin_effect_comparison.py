"""
HAT_groin_effect_comparison.py
===============================
Standalone publication-figure script. Compares a 'no_groin' and a 'groin'
1967-2017 hindcast run against the real observed shoreline at MULTIPLE
checkpoint years (1997 mid-run, 2017 endpoint), isolating the groin's
(and its deterioration's) effect on modeled shoreline position over time.

Generalized from the original single-endpoint (1967->1997) version to an
arbitrary list of checkpoint years, so the 1997 mid-run point is a genuine
validation target -- not just the final year -- which matters specifically
for the 1995->2003 deterioration ramp this extended run is testing.

OBSERVED CHANGE METHODOLOGY -- IMPORTANT
------------------------------------------
Observed change comes from the validated WET/DRY change table
(Change_from_wetdry_1967_D2_D12.csv, built by
HAT_geometric_distance_sanity_check.py) -- NOT the dune-line raw CSVs used
in earlier versions of this script. Two decisions behind that switch, both
already validated empirically, not just asserted:
  1. CASCADE's own x_s is defined (barrier3d.py) as x_t + LShoreface -- a
     shoreface-toe-based "shoreline" conceptually much closer to a water-
     line/MHW-type contour than to the dune/vegetation line. Confirmed: the
     1967 wet/dry line sits seaward of the 1967 dune line at every domain,
     exactly as physically expected.
  2. The change table is internally self-consistent -- wetdry_1967 is the
     baseline, every other wetdry_YYYY is differenced against that SAME
     fixed reference throughout, never re-baselined per year (re-baselining
     to each year's own minimum domain was shown to silently erase real
     signal -- D12's own raw distance moved ~33 m between 1967 and 2017,
     which a per-year re-zeroing would hide).
Only the wetdry_* columns of that table are used here -- the table also
carries duneline_* columns (referenced to the same wetdry_1967 baseline)
for diagnostic purposes, but mixing one into the observed series would
reintroduce the exact feature-mismatch this switch fixed.

ONE SHARED REFERENCE FRAME, MULTIPLE CHECKPOINTS
--------------------------------------------------
  - 1967 shoreline    (model's own dune-line-based initial planform -- a
                        DIFFERENT feature than the wet/dry observed series,
                        but that's fine: every comparison here is CHANGE
                        from each series' own 1967 reference, never
                        absolute position across the two)
  - For each year in CHECKPOINT_YEARS (e.g. 1997, 2017):
      - OBSERVED     (wetdry_1967 baseline + observed wetdry 1967->year change)
      - modeled, no groin  (model's own row at that calendar year)
      - modeled, with groin

Usage
-----
Set RUN_NO_GROIN / RUN_GROIN below to your two run folder names under
output/raw_runs/. Set COMPARISON_SUBFOLDER to whatever label you want for
THIS comparison. Then run. Produces a combined overview figure (all
checkpoints together) AND one individual figure per checkpoint year.

Saves to:
    scripts/groin/HAT-hindcast-groin-test/comparison/{COMPARISON_SUBFOLDER}/
        groin_effect_overview.png             (all checkpoints, one figure)
        groin_effect_overview_v2.png           (same, no title/footer)
        groin_effect_{year}.png                (per checkpoint year)
        groin_effect_{year}_v2.png              (same, no title/footer)
        shoreline_evolution.gif
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
RUN_NO_GROIN = "HAT_1967_2018_M60_deterioration_no_groin"
RUN_GROIN    = "HAT_1967_2018_M60_deterioration_groin"

# --- Where THIS comparison's figures are saved. Name the subfolder yourself so
# different comparisons never overwrite each other. ---
COMPARISON_SUBFOLDER = "deterioration_1995_2003_Mover3"   # <- edit per comparison
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

# START_YEAR is the shared reference year (row 0 of every shoreline matrix).
# CHECKPOINT_YEARS are the observed validation targets -- add/remove years
# here; everything else (figures, GIF markers, observed-change loading)
# follows automatically, AS LONG AS a matching "wetdry_<year>" column exists
# in WETDRY_CHANGE_TABLE (below).
START_YEAR = 1967
CHECKPOINT_YEARS = [1997, 2017]

# Per-checkpoint styling so 1997 vs 2017 are visually distinguishable while
# keeping color-by-series-type (observed=black, no_groin=orange, groin=red)
# consistent across years in the combined overview figure.
CHECKPOINT_STYLE = {
    1997: dict(alpha=0.55, marker="s", ls="--"),
    2017: dict(alpha=1.00, marker="D", ls="-"),
}
DEFAULT_CHECKPOINT_STYLE = dict(alpha=0.8, marker="o", ls="-.")

# ── Sign convention (matches main hindcast + HAT_plot_groin_runs.py exactly) ──
FLIP_SIGN_MODEL = True

# --- Styling (from your established figure conventions) ---
MODEL_COLOR        = "#FF8C00"   # warm orange
GROIN_COLOR        = "#B71C1C"   # groin red
GROIN_BOUNDARY_GIS = 5.5         # D5/D6 interface (Buxton groin)
DOMAIN_TICK_STEP   = 2
DOMAIN_SPACING_M   = 500.0       # alongshore domain width (m) -- project convention
OCEAN_AT_BOTTOM    = True        # seaward plots downward
AXIS_LABEL_FONTSIZE = 12         # shared across every figure variant

# --- Observed reference: the validated wet/dry change table from
# HAT_geometric_distance_sanity_check.py, NOT the raw duneline CSVs -- see
# module docstring for why wet/dry is the better-matched proxy for CASCADE's
# own x_s ("shoreline" = x_t + LShoreface, conceptually closer to a water-
# line than to the dune/vegetation line), and why it must be the WETDRY-
# referenced table specifically (Change_from_wetdry_1967_D2_D12.csv), not
# the dune-line-referenced one -- mixing a dune-line point into an otherwise
# all-wet/dry observed series reintroduces the exact feature-mismatch this
# switch was meant to fix.
WETDRY_CHANGE_TABLE = os.path.join(
    PROJECT_BASE_DIR, "scripts", "groin", "HAT-hindcast-groin-test",
    "input_prep", "shoreline_position", "output",
    "Change_from_wetdry_1967_D2_D12.csv",
)
WETDRY_DOMAIN_COL = "Domain_ID"

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


def _year_to_row(year, nt):
    """Calendar year -> row index in a run's shoreline matrix (row 0 =
    START_YEAR, row t = START_YEAR + t). Returns None (with a printed
    warning) if the run doesn't reach that year -- guards against the
    END_YEAR-is-exclusive convention silently pulling the wrong row."""
    row = year - START_YEAR
    if not (0 <= row < nt):
        print(f"  WARNING: year {year} (row {row}) is outside this run's "
              f"{nt} modeled years (rows 0-{nt - 1}) -- checkpoint skipped.")
        return None
    return row


def _load_shoreline(run_name):
    """Load a run's shoreline matrix (nt, ndomain) in meters."""
    path = os.path.join(RUN_DATA_DIR, run_name, f"{run_name}_shoreline_matrix.npy")
    if not os.path.isfile(path):
        raise FileNotFoundError(f"shoreline matrix not found:\n  {path}")
    m = np.load(path)
    print(f"  Loaded {run_name}: shape {m.shape}")
    return m


def _load_observed_changes():
    """Load observed shoreline CHANGE for every year in CHECKPOINT_YEARS from
    the validated wet/dry change table (Change_from_wetdry_1967_D2_D12.csv --
    already computed as raw_distance[year] - raw_distance[1967], sharing ONE
    fixed absolute baseline throughout; see that script's module docstring
    for why re-baselining per year would silently erase real signal).

    Returns dict {year: np.ndarray or None} -- + = landward/erosion (RAW
    sign, matching what _flip() expects downstream). Missing file/column
    gives None for that year (omitted from plots, run continues)."""
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

def fig_checkpoint(no_groin_m, groin_m, no_groin_name, groin_name,
                    checkpoint_year, observed_change,
                    show_title=True, figsize=(11, 6)):
    """PUBLICATION FIGURE: isolates the groin's effect on shoreline position
    AT ONE CHECKPOINT YEAR against a real observed target, all on one shared
    reference frame. Generalizes the original single-endpoint version to any
    checkpoint_year in CHECKPOINT_YEARS (e.g. 1997 mid-run, 2017 endpoint).

    observed_change : the 1967->checkpoint_year array from
        _load_observed_changes(), or None if that year's raw file is missing.

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

    nt = no_groin_m.shape[0]
    row = _year_to_row(checkpoint_year, nt)
    have_model_row = row is not None

    end_no_groin = _flip(_real_slice(no_groin_m[row])) - ref_mean if have_model_row else None
    end_groin    = _flip(_real_slice(groin_m[row]))    - ref_mean if have_model_row else None

    # Observed checkpoint target: real 1967->checkpoint_year CHANGE (from raw
    # offset files) added onto the shared 1967 reference -- see module docstring.
    have_obs = observed_change is not None
    observed_target = (planform_1967 - observed_change) if have_obs else None
    # (observed_change is "+ = landward"; subtracting it flips to "+ = seaward",
    # matching the model's _flip() convention, before adding onto the 1967 ref.)

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
                alpha=0.85, label=f"{checkpoint_year} modeled \u2014 no groin", zorder=4)
        ax.plot(gis, end_groin, color=GROIN_COLOR, ls="-", lw=2.4, marker="D", ms=5,
                label=f"{checkpoint_year} modeled \u2014 with groin", zorder=5)
    else:
        ax.text(0.5, 0.85, f"modeled {checkpoint_year} outside this run's simulated years",
                transform=ax.transAxes, ha="center", color="firebrick", fontsize=9)

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
            f"Buxton Groin Field, Hatteras Island, NC  |  {START_YEAR}\u2013{checkpoint_year}",
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
            f"({START_YEAR}, {checkpoint_year})  |  no-groin: {no_groin_name}  |  "
            f"groin: {groin_name}",
            xy=(0, 0), xycoords="axes fraction", xytext=(0, -0.16),
            textcoords="axes fraction", fontsize=7.5, color="#666666",
            ha="left", va="top", style="italic", annotation_clip=False,
        )

    tag = f"groin_effect_{checkpoint_year}" if show_title else f"groin_effect_{checkpoint_year}_v2"
    return fig, tag


def fig_combined_overview(no_groin_m, groin_m, no_groin_name, groin_name,
                           observed_changes, show_title=True, figsize=(12, 7)):
    """PUBLICATION FIGURE: ALL checkpoint years on one shared reference frame
    -- 1967 reference plus, for every year in CHECKPOINT_YEARS, the observed
    target and both modeled trajectories. Busier than fig_checkpoint() by
    design (up to 3 + 3*len(CHECKPOINT_YEARS) series); each checkpoint year
    gets a distinct alpha/marker/linestyle combo (CHECKPOINT_STYLE) so 1997
    vs 2017 stay visually distinguishable while color still encodes series
    TYPE (observed=black, no_groin=orange, groin=red) consistently."""
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
                    label=f"{year} modeled \u2014 no groin", zorder=4)
            ax.plot(gis, end_groin, color=GROIN_COLOR, ls=style["ls"], lw=2.2,
                    marker=style["marker"], ms=5, alpha=style["alpha"],
                    label=f"{year} modeled \u2014 with groin", zorder=5)
        else:
            print(f"  [overview] no modeled data for {year} -- omitted from overview.")

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
    if show_title:
        year_span = f"{START_YEAR}\u2013{CHECKPOINT_YEARS[-1]}" if CHECKPOINT_YEARS else f"{START_YEAR}"
        ax.set_title(
            f"Effect of Groin Deterioration on Modeled Shoreline Position \u2014 All Checkpoints\n"
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


# =============================================================================
# ANIMATED GIF
# =============================================================================

def make_shoreline_evolution_gif(no_groin_m, groin_m, no_groin_name, groin_name,
                                  observed_changes, out_path):
    """Animated GIF: shoreline position evolving year-by-year for both the
    no-groin (orange) and with-groin (red) runs, on the same shared
    1967-anchored reference frame as fig_checkpoint()/fig_combined_overview().

    The 1967 initial shoreline and EVERY checkpoint year's OBSERVED target
    (1997, 2017, ...) are drawn as thin, static dashed reference lines that
    stay put for the whole animation; only the two modeled trajectories move.
    """
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

    # Year label per frame: row 0 = START_YEAR, row t = START_YEAR + t (exact
    # integers -- the previous np.linspace(START_YEAR, END_YEAR, nt) approach
    # didn't land on clean years whenever nt-1 didn't evenly divide the span).
    years = START_YEAR + np.arange(nt)

    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)

    # --- Static thin dashed references (do not move across frames) ---
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

    print(f"\nLoading observed checkpoints ({START_YEAR} -> {CHECKPOINT_YEARS})...")
    observed_changes = _load_observed_changes()

    # --- Combined overview: all checkpoints on one figure ---
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

    # --- Individual figure per checkpoint year ---
    for year in CHECKPOINT_YEARS:
        observed_change = observed_changes.get(year)

        fig1, tag1 = fig_checkpoint(no_groin_m, groin_m, RUN_NO_GROIN, RUN_GROIN,
                                     year, observed_change, show_title=True)
        out1 = os.path.join(COMPARISON_OUTPUT_DIR, f"{tag1}.png")
        fig1.savefig(out1, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out1}")

        # v2: same plot, no title/subtitle/footer (e.g. for a numbered
        # dissertation figure with its caption set in the surrounding text).
        fig2, tag2 = fig_checkpoint(no_groin_m, groin_m, RUN_NO_GROIN, RUN_GROIN,
                                     year, observed_change, show_title=False,
                                     figsize=(11, 5))
        out2 = os.path.join(COMPARISON_OUTPUT_DIR, f"{tag2}.png")
        fig2.savefig(out2, dpi=FIGURE_DPI, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out2}")

    # Animated GIF: shoreline evolution over time, no-groin vs with-groin,
    # with every checkpoint year's observed target as a static reference.
    gif_out = os.path.join(COMPARISON_OUTPUT_DIR, "shoreline_evolution.gif")
    make_shoreline_evolution_gif(no_groin_m, groin_m, RUN_NO_GROIN, RUN_GROIN,
                                  observed_changes, gif_out)

    if SHOW_FIGURE:
        plt.show()


if __name__ == "__main__":
    main()
