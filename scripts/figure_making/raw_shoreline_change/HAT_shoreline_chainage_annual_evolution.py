"""
HAT_shoreline_chainage_alldata_evolution.py
====================================
Visualises CoastSat shoreline chainage across Hatteras Island over 40 years.

Outputs
-------
1.  Four 10-year panel figures (one PNG each):
        1984–1994, 1994–2004, 2004–2014, 2014–2024
    Each shows annual median chainage deviation from the period-start baseline
    for every individual transect.  Lines coloured light→dark as years pass.

2.  A GIF animating spring (Apr–May) and fall (Oct–Nov) seasonal median
    shoreline profiles year by year, 1984–2024.

Usage
-----
    Edit the CONFIG section, then run:
        python HAT_shoreline_chainage_alldata_evolution.py

Dependencies
------------
    pip install pandas numpy matplotlib scipy imageio tqdm
"""

import os
import glob

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import imageio.v2 as imageio
from tqdm import tqdm

# ============================================================
# CONFIG
# ============================================================

ROOT_DATA_DIR = r"/scripts/input_prep/5-scr/CoastSat\coastsat_timeseries"
LOOKUP_CSV    = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\transect_domain_lookup.csv"
OUTPUT_DIR    = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\figure_making\raw_shoreline_change\annual_output"
SITE_FILTER   = "usa_NC"

FULL_START = "1984-01-01"
FULL_END   = "2024-12-31"

PANELS = [
    ("1984–1994", "1984-01-01", "1994-12-31"),
    ("1994–2004", "1994-01-01", "2004-12-31"),
    ("2004–2014", "2004-01-01", "2014-12-31"),
    ("2014–2024", "2014-01-01", "2024-12-31"),
]

MIN_OBS_PER_YEAR      = 2     # min observations per transect per year/window to use
MIN_TRANSECT_FRACTION = 0.20  # min fraction of transects valid to draw a line/frame
SMOOTH_SIGMA          = 2     # Gaussian sigma in transect-index units (0 = off)

# GIF frame durations — set in seconds, both values are straightforward
GIF_SEASONAL_DURATION  = 0.65   # seconds per frame for the spring/fall 40-year GIF
GIF_DECADAL_DURATION   = 3   # seconds per frame for the 4-decade panel GIF
GIF_FIG_WIDTH      = 14    # inches
PANEL_CMAP         = "Blues"
BG_COLOR           = "white"

NUM_REAL_DOMAINS = 90

# ── Geographic annotation (matches publication figure style) ──────────────────

ANN_TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
ANN_VILLAGE_LINES = {"Salvo": 69, "Waves": 74, "Rodanthe": 80}
ANN_PIERS         = {"Avon Pier": 26, "Rodanthe Pier": 79}
ANN_GROINS        = {"Buxton Groin": 5.5}
ANN_WIMBLE_SHOALS = (60, 74)

ANN_PIER_LABEL_Y  = 0.76
ANN_GROIN_LABEL_Y = 0.68

ANN_C_TOWN_SPAN    = "#90AFC5"   # steel blue — community span fill
ANN_C_WIMBLE       = "#E0A800"   # amber     — Wimble Shoals fill
ANN_C_VILLAGE_LINE = "0.40"      # mid-grey  — Salvo/Waves/Rodanthe lines
ANN_C_PIER         = "#1565C0"   # dark blue — pier lines
ANN_C_GROIN        = "#B71C1C"   # dark red  — groin lines

# ============================================================
# HELPERS
# ============================================================

def find_csvs(root_dir, site_filter):
    csv_map = {}
    for folder in os.listdir(root_dir):
        if not folder.startswith(site_filter):
            continue
        folder_path = os.path.join(root_dir, folder)
        if not os.path.isdir(folder_path):
            continue
        for csv_file in glob.glob(os.path.join(folder_path, "*.csv")):
            tid = os.path.splitext(os.path.basename(csv_file))[0]
            csv_map[tid] = csv_file
    return csv_map


def load_transect(csv_path):
    try:
        df = pd.read_csv(csv_path)
        df.columns = [c.strip() for c in df.columns]
        date_col  = next((c for c in df.columns if c.lower().startswith("date")), None)
        chain_col = next((c for c in df.columns
                          if "chainage" in c.lower() or "cross_dist" in c.lower()), None)
        if date_col is None or chain_col is None:
            return None
        df[date_col] = pd.to_datetime(df[date_col], utc=True).dt.tz_localize(None)
        df = df.dropna(subset=[date_col, chain_col]).sort_values(date_col)
        return df.set_index(date_col)[chain_col]
    except Exception:
        return None


def gaussian_smooth_1d(arr, sigma):
    if sigma <= 0:
        return arr
    from scipy.ndimage import gaussian_filter1d
    mask   = np.isnan(arr)
    filled = arr.copy(); filled[mask] = 0.0
    sm     = gaussian_filter1d(filled, sigma=sigma)
    wt     = gaussian_filter1d((~mask).astype(float), sigma=sigma)
    with np.errstate(invalid="ignore"):
        return np.where(wt > 0.05, sm / wt, np.nan)


def build_annual_matrix(ts_dict, transect_order, start, end, min_obs=2):
    p_start = pd.Timestamp(start)
    p_end   = pd.Timestamp(end)
    years   = list(range(p_start.year, p_end.year + 1))
    records = {}
    for tid in transect_order:
        ts = ts_dict.get(tid)
        if ts is None:
            records[tid] = [np.nan] * len(years); continue
        ts_p = ts[(ts.index >= p_start) & (ts.index <= p_end)]
        row = []
        for yr in years:
            d = ts_p[ts_p.index.year == yr]
            row.append(float(d.median()) if len(d) >= min_obs else np.nan)
        records[tid] = row
    df = pd.DataFrame(records, index=years).T
    return df, np.array(years)


def transect_index(domain_per_transect, domain):
    """Return the transect array index nearest to a given domain number."""
    return int(np.searchsorted(domain_per_transect, domain))


def annotate_axes_publication(ax, domain_per_transect, ylim, domain_tick_every=5):
    """
    Apply the full publication annotation suite to ax.
    ylim must be the (ymin, ymax) already set on ax so labels are placed correctly.
    """
    ymin, ymax = ylim
    yspan = ymax - ymin

    # ── Wimble Shoals zone ────────────────────────────────────────────────────
    ws0 = transect_index(domain_per_transect, ANN_WIMBLE_SHOALS[0])
    ws1 = transect_index(domain_per_transect, ANN_WIMBLE_SHOALS[1])
    ax.axvspan(ws0, ws1, color=ANN_C_WIMBLE, alpha=0.18, zorder=0)
    ax.text((ws0 + ws1) / 2, ymax - 0.01 * yspan,
            "Wimble Shoals", ha="center", va="top",
            fontsize=6.5, color=ANN_C_WIMBLE, style="italic", zorder=3)

    # ── Community spans ───────────────────────────────────────────────────────
    for name, (d0, d1) in ANN_TOWN_SPANS.items():
        i0 = transect_index(domain_per_transect, d0)
        i1 = transect_index(domain_per_transect, d1)
        ax.axvspan(i0, i1, color=ANN_C_TOWN_SPAN, alpha=0.25, zorder=0)
        ax.text((i0 + i1) / 2, ymax - 0.01 * yspan,
                name, ha="center", va="top",
                fontsize=7, color=ANN_C_TOWN_SPAN, fontweight="bold", zorder=3)

    # ── Village centre lines (Salvo, Waves, Rodanthe) ─────────────────────────
    for vname, vdom in ANN_VILLAGE_LINES.items():
        ix = transect_index(domain_per_transect, vdom)
        ax.axvline(ix, color=ANN_C_VILLAGE_LINE, lw=0.7, ls=(0, (3, 4)), zorder=2)

    # ── Piers ─────────────────────────────────────────────────────────────────
    for pname, pdom in ANN_PIERS.items():
        ix = transect_index(domain_per_transect, pdom)
        ax.axvline(ix, color=ANN_C_PIER, lw=1.0, ls="--", zorder=2)
        label_y = ymin + ANN_PIER_LABEL_Y * yspan
        ax.text(ix + 1, label_y, pname,
                rotation=90, va="bottom", ha="left",
                fontsize=6.5, color=ANN_C_PIER, zorder=3)

    # ── Groins ────────────────────────────────────────────────────────────────
    for gname, gdom in ANN_GROINS.items():
        ix = transect_index(domain_per_transect, gdom)
        ax.axvline(ix, color=ANN_C_GROIN, lw=1.0, ls="--", zorder=2)
        label_y = ymin + ANN_GROIN_LABEL_Y * yspan
        ax.text(ix + 1, label_y, gname,
                rotation=90, va="bottom", ha="left",
                fontsize=6.5, color=ANN_C_GROIN, zorder=3)

    # ── Top x-axis: transect index numbers ────────────────────────────────────
    transect_ticks = np.arange(0, len(domain_per_transect), 100)
    ax2 = ax.secondary_xaxis("top")
    ax2.set_xticks(transect_ticks)
    ax2.set_xticklabels([str(t) for t in transect_ticks], fontsize=6)
    ax2.set_xlabel("Transect index", fontsize=7)

    # ── Bottom x-axis: replace generic label with domain tick marks ───────────
    # Place domain number ticks on the primary axis
    domain_ticks = np.arange(1, NUM_REAL_DOMAINS + 1, domain_tick_every)
    tick_pos     = [transect_index(domain_per_transect, d) for d in domain_ticks]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels([str(d) for d in domain_ticks], fontsize=6)
    ax.set_xlabel("CASCADE domain  (1 = Buxton  →  90 = Rodanthe)", fontsize=7)


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    gif_dir = os.path.join(OUTPUT_DIR, "gif_frames")
    os.makedirs(gif_dir, exist_ok=True)

    # ── 1. Lookup table ──────────────────────────────────────────────────────
    print("Loading lookup table …")
    lookup = pd.read_csv(LOOKUP_CSV)
    lookup["transect_id"]   = lookup["transect_id"].astype(str)
    lookup["domain_number"] = lookup["domain_number"].astype(int)
    lookup = lookup[
        (lookup["domain_number"] >= 1) &
        (lookup["domain_number"] <= NUM_REAL_DOMAINS)
    ].sort_values(["domain_number", "transect_id"]).reset_index(drop=True)
    transect_order      = lookup["transect_id"].tolist()
    domain_per_transect = lookup["domain_number"].values
    n_transects         = len(transect_order)
    print(f"  {n_transects} transects across {lookup['domain_number'].nunique()} domains.\n")

    # ── 2. Load CSVs ─────────────────────────────────────────────────────────
    print("Discovering CSV files …")
    csv_map = find_csvs(ROOT_DATA_DIR, SITE_FILTER)
    print(f"  Found {len(csv_map)} CSVs on disk.\n")

    print("Loading transect time-series …")
    ts_dict = {}
    for tid in tqdm(transect_order, unit="transect"):
        if tid in csv_map:
            ts = load_transect(csv_map[tid])
            if ts is not None and len(ts) > 0:
                ts_dict[tid] = ts
    print(f"  Successfully loaded {len(ts_dict)} / {n_transects} transects.\n")

    x = np.arange(n_transects)

    # ============================================================
    # PART A — Four 10-year panel figures (shared 1984 baseline)
    # ============================================================
    print("=" * 60)
    print("Building 10-year panel figures …")
    print("=" * 60)

    # ── Compute 1984 baseline ONCE — reused by all panels and GIF ────────────
    # Use full-year 1984 median (not spring-only) for panels so early Landsat
    # years with sparse seasonal coverage still get a valid baseline.
    print("  Computing 1984 annual baseline for panels …")
    panel_baseline_1984 = {}
    for tid in transect_order:
        ts = ts_dict.get(tid)
        if ts is not None:
            fy84 = ts[ts.index.year == 1984]
            panel_baseline_1984[tid] = float(fy84.median()) if len(fy84) >= 1 else np.nan
        else:
            panel_baseline_1984[tid] = np.nan
    panel_baseline_ser = pd.Series(panel_baseline_1984, name=1984)

    # ── Pre-compute all annual deviations from 1984 baseline across full record
    # so we can set a shared y-axis range before drawing any panel.
    print("  Pre-computing full-record deviations for shared y-axis …")
    full_ann_df, full_years = build_annual_matrix(
        ts_dict, transect_order, FULL_START, FULL_END, MIN_OBS_PER_YEAR
    )
    full_dev_df = full_ann_df.subtract(panel_baseline_ser, axis=0)

    all_panel_vals = full_dev_df.values.flatten()
    all_panel_vals = all_panel_vals[~np.isnan(all_panel_vals)]
    if len(all_panel_vals):
        plo = all_panel_vals.min()
        phi = all_panel_vals.max()
        pad = (phi - plo) * 0.08   # 8% breathing room
        shared_ylim = (plo - pad, phi + pad)
    else:
        shared_ylim = (-150, 150)
    print(f"  Shared y-axis: {shared_ylim[0]:.0f} to {shared_ylim[1]:.0f} m\n")

    cmap_panel = matplotlib.colormaps.get_cmap(PANEL_CMAP)
    panel_frame_paths = []   # for decadal GIF

    for panel_label, p_start, p_end in PANELS:
        print(f"  Panel: {panel_label}")
        p_start_yr = int(p_start[:4])
        p_end_yr   = int(p_end[:4])

        # Slice full deviation matrix to this panel's years
        panel_years = [y for y in full_years if p_start_yr <= y <= p_end_yr]
        if not panel_years:
            print(f"    No data — skipping")
            continue

        norm_panel = mcolors.Normalize(vmin=0, vmax=max(len(panel_years) - 1, 1))

        fig, ax = plt.subplots(figsize=(16, 5))
        ax.set_facecolor(BG_COLOR)
        fig.patch.set_facecolor(BG_COLOR)

        for i, yr in enumerate(panel_years):
            if yr not in full_dev_df.columns:
                continue
            profile = full_dev_df[yr].values.astype(float)
            n_valid = np.sum(~np.isnan(profile))
            if n_valid < MIN_TRANSECT_FRACTION * n_transects:
                print(f"    Skipping {yr}: only {n_valid}/{n_transects} valid transects")
                continue
            if SMOOTH_SIGMA > 0:
                profile = gaussian_smooth_1d(profile, SMOOTH_SIGMA)
            color = cmap_panel(norm_panel(i))
            alpha = 0.65 + 0.35 * (i / max(len(panel_years) - 1, 1))
            ax.plot(x, profile, color=color, lw=0.9, alpha=alpha, zorder=2 + i)

        ax.axhline(0, color="black", lw=0.9, ls="--", alpha=0.6, zorder=1,
                   label="1984 baseline")

        ax.set_xlim(0, n_transects - 1)
        ax.set_ylim(shared_ylim)

        ax.set_xlabel(
            "CASCADE domain  (1 = Buxton  →  90 = Rodanthe)",
            fontsize=9)
        ax.set_ylabel(
            "Shoreline deviation from 1984 baseline (m)\n▲ accretion (seaward)   ▼ erosion (landward)",
            fontsize=9)
        ax.set_title(
            f"Hatteras Island  |  Shoreline chainage evolution  |  {panel_label}\n"
            "Individual CoastSat transects · annual median · deviation from 1984 baseline",
            fontsize=11, fontweight="bold")

        # Colorbar as year legend
        sm = cm.ScalarMappable(
            cmap=PANEL_CMAP,
            norm=mcolors.Normalize(vmin=p_start_yr, vmax=p_end_yr))
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, orientation="vertical", pad=0.01, shrink=0.75)
        cb.set_label("Year", fontsize=8)
        cb.set_ticks(np.arange(p_start_yr, p_end_yr + 1, 2))

        annotate_axes_publication(ax, domain_per_transect, shared_ylim, domain_tick_every=5)

        ax.tick_params(axis="both", labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)

        fig.tight_layout()
        safe_label = panel_label.replace("–", "_")
        out_path = os.path.join(OUTPUT_DIR, f"shoreline_evolution_{safe_label}.png")
        fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor=BG_COLOR)
        panel_frame_paths.append(out_path)
        plt.close(fig)
        print(f"    Saved → {out_path}")

    # ── Decadal GIF: cycle through the four panel PNGs slowly ─────────────────
    if panel_frame_paths:
        print("\n  Assembling decadal GIF …")
        decadal_gif_path = os.path.join(OUTPUT_DIR, "HAT_shoreline_decadal.gif")
        # Use PIL directly — imageio duration is unreliable across versions/viewers.
        # PIL duration is always milliseconds, consistent across all GIF viewers.
        from PIL import Image
        pil_frames = [Image.open(fp).convert("RGBA") for fp in panel_frame_paths]
        pil_frames[0].save(
            decadal_gif_path,
            save_all=True,
            append_images=pil_frames[1:],
            duration=int(GIF_DECADAL_DURATION * 500),  # *500 corrects Pillow doubling bug
            loop=0,
        )
        print(f"  Decadal GIF saved → {decadal_gif_path}")

    # ============================================================
    # PART B — GIF: spring (Apr–May) and fall (Oct–Nov) per year
    # ============================================================
    print("\n" + "=" * 60)
    print("Building GIF frames (spring + fall per year, 1984–2024) …")
    print("=" * 60)

    SEASONS = [
        ("Spring (Apr–May)", 4, 5,  "▲"),
        ("Fall  (Oct–Nov)",  10, 11, "▼"),
    ]

    gif_years = list(range(int(FULL_START[:4]) + 1, int(FULL_END[:4]) + 1))  # skip 1984 (baseline year)

    # 1984 baseline: spring median, fall back to full-year median
    print("  Computing 1984 spring baseline …")
    baseline_1984 = {}
    for tid in transect_order:
        ts = ts_dict.get(tid)
        if ts is not None:
            sp84 = ts[(ts.index.year == 1984) & ts.index.month.isin([4, 5])]
            fy84 = ts[ts.index.year == 1984]
            if len(sp84) >= 1:
                baseline_1984[tid] = float(sp84.median())
            elif len(fy84) >= 1:
                baseline_1984[tid] = float(fy84.median())
            else:
                baseline_1984[tid] = np.nan
        else:
            baseline_1984[tid] = np.nan
    baseline_arr = np.array([baseline_1984[tid] for tid in transect_order], dtype=float)

    # Build seasonal profiles
    print("  Computing seasonal window medians …")
    seasonal_profiles = {}
    frame_labels      = []
    for yr in gif_years:
        for season_label, m0, m1, sym in SEASONS:
            row = []
            for tid in transect_order:
                ts = ts_dict.get(tid)
                if ts is not None:
                    w = ts[(ts.index.year == yr) & ts.index.month.isin(range(m0, m1 + 1))]
                    row.append(float(w.median()) if len(w) >= 1 else np.nan)
                else:
                    row.append(np.nan)
            arr = np.array(row, dtype=float) - baseline_arr
            key = (yr, season_label, sym)
            seasonal_profiles[key] = arr
            frame_labels.append(key)

    # Global colour scale
    all_devs   = np.concatenate(list(seasonal_profiles.values()))
    valid_devs = all_devs[~np.isnan(all_devs)]
    if len(valid_devs) == 0:
        clim = 50.0
    else:
        vabs = np.ceil(max(abs(valid_devs.min()), abs(valid_devs.max())) / 5) * 5
        clim = vabs
    print(f"  Colour scale: ±{clim:.0f} m")
    print(f"  Rendering up to {len(frame_labels)} frames …\n")

    # Custom diverging colormap: saturated red (erosion/negative) → mid-grey
    # → saturated blue (accretion/positive). Stays visible near zero unlike
    # RdBu which passes through white. Negative = erosion = red, Positive = accretion = blue.
    from matplotlib.colors import LinearSegmentedColormap
    gif_cmap = LinearSegmentedColormap.from_list(
        "erosion_accretion",
        [
            (0.00, "#4a0010"),   # near-black deep crimson — extreme erosion
            (0.15, "#b2182b"),   # strong red              — heavy erosion
            (0.35, "#ef8a62"),   # light red-orange        — mild erosion
            (0.50, "#878787"),   # mid grey                — no change
            (0.65, "#67a9cf"),   # light blue              — mild accretion
            (0.85, "#2166ac"),   # strong blue             — heavy accretion
            (1.00, "#08306b"),   # near-black deep navy    — extreme accretion
        ]
    )
    norm_gif  = mcolors.Normalize(vmin=-clim, vmax=clim)

    frame_paths = []
    for di, (yr, season_label, sym) in enumerate(tqdm(frame_labels, desc="  GIF frames")):
        profile = seasonal_profiles[(yr, season_label, sym)].copy()

        if SMOOTH_SIGMA > 0:
            profile = gaussian_smooth_1d(profile, SMOOTH_SIGMA)

        n_valid = np.sum(~np.isnan(profile))
        if n_valid < MIN_TRANSECT_FRACTION * n_transects:
            continue

        fig, ax = plt.subplots(figsize=(GIF_FIG_WIDTH, 4.2))
        ax.set_facecolor(BG_COLOR)
        fig.patch.set_facecolor(BG_COLOR)

        # Draw the chainage profile coloured by deviation
        for i in range(n_transects - 1):
            if np.isnan(profile[i]) and np.isnan(profile[i + 1]):
                continue
            mid_val = float(np.nanmean([profile[i], profile[i + 1]]))
            ax.plot([i, i + 1], [profile[i], profile[i + 1]],
                    color=gif_cmap(norm_gif(mid_val)), lw=1.4, solid_capstyle="round")

        ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.55)
        ax.set_xlim(0, n_transects - 1)
        ax.set_ylim(-clim * 1.05, clim * 1.05)
        ylim = ax.get_ylim()

        ax.set_xlabel(
            "Along-island position  →  (south Cape Hatteras → north Rodanthe)",
            fontsize=8)
        ax.set_ylabel("Deviation from 1984 baseline (m)\n▲ accretion (seaward)  ▼ erosion (landward)", fontsize=8)

        # Title with season badge as coloured text prefix
        ax.set_title(
            f"Hatteras Island  —  CoastSat shoreline position  |  {yr}",
            fontsize=11, fontweight="bold", color="black")
        # Season shown as a quiet subtitle below the main title
        ax.text(0.50, 0.11, f"{sym} {season_label}",
                transform=ax.transAxes, ha="center", va="bottom",
                fontsize=8, color="#555555", style="italic")

        # Colorbar
        sm = cm.ScalarMappable(cmap=gif_cmap, norm=norm_gif)
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, orientation="vertical", pad=0.01, shrink=0.82)
        cb.set_label("Shoreline deviation (m)\nblue = accretion (seaward) · red = erosion (landward)", fontsize=7)
        cb.ax.tick_params(labelsize=7)

        annotate_axes_publication(ax, domain_per_transect, ylim, domain_tick_every=5)

        # Year + elapsed box — bottom centre
        ax.text(0.50, 0.03,
                f"{yr}  ({yr - 1984} yr{'s' if yr > 1984 else ''} from baseline)",
                transform=ax.transAxes, ha="center", va="bottom",
                fontsize=9, color="#222222",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc",
                          alpha=0.90, lw=0.8))

        ax.tick_params(axis="both", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()

        frame_path = os.path.join(gif_dir, f"frame_{di:04d}.png")
        fig.savefig(frame_path, dpi=110, bbox_inches="tight", facecolor=BG_COLOR)
        plt.close(fig)
        frame_paths.append(frame_path)

    # Assemble GIF
    if frame_paths:
        print(f"\n  Assembling {len(frame_paths)}-frame GIF …")
        gif_path = os.path.join(OUTPUT_DIR, "HAT_shoreline_evolution_40yr.gif")
        from PIL import Image
        pil_seasonal = [Image.open(fp).convert("RGBA") for fp in frame_paths]
        pil_seasonal[0].save(
            gif_path,
            save_all=True,
            append_images=pil_seasonal[1:],
            duration=int(GIF_SEASONAL_DURATION * 500),  # *500 corrects Pillow doubling bug
            loop=0,
        )
        print(f"  GIF saved → {gif_path}")
    else:
        print("  No valid frames — GIF not created.")

    print("\nDone.")


if __name__ == "__main__":
    main()
