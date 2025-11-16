# =============================================================================
# CASCADE Topography Figure Maker (Benton-parity outputs, presentation-ready)
# ---------------------------------------------------------------------------
# Supports 200-row topography. Filters stitched view to domains 30–134.
# Saves into <repo>/output/figures/topography_<YEAR>/
# Toggle views with VIEW_PRESET = "index" or "metric"
# =============================================================================

import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter, FormatStrFormatter
from pathlib import Path

# --------------------------- PATHS -------------------------------------------
# Repo root inferred from this file's path:
#   .../CASCADE/scripts/input_preperation/topography/topo_figures.py
PROJECT_ROOT = Path(__file__).resolve().parents[3]

# Data locations (adjust if your year subfolders differ)
BASE_PATH = PROJECT_ROOT / "data" / "hatteras_init"
TOPO_DIR  = BASE_PATH / "topography" / "2009"
DUNE_DIR  = BASE_PATH / "dunes"      / "2009"
YEAR_SUFFIX = "2009"                 # *_topography_2009.npy, *_dune_2009.npy

# Save directory INSIDE the repo
SAVE_DIR = PROJECT_ROOT / "output" / "figures" / f"topography_{YEAR_SUFFIX}"
SAVE_PDF = True
SAVE_DIR.mkdir(parents=True, exist_ok=True)

VERBOSE = True

# --------------------------- USER INPUTS -------------------------------------
# Choose ONE (None stitches island with filter below)
DOMAIN_ID   = None                  # e.g., "domain_84_resampled" or None

# Domain range filter for stitched plots
DOMAIN_MIN_ID = 30
DOMAIN_MAX_ID = 134                 # inclusive

# View preset
VIEW_PRESET = "index"               # "index" or "metric"

# Figure + style
FIGSIZE     = (12, 7)
DPI         = 300
CMAP        = "terrain"
SHOW_IN_METERS = True
TITLE_YEAR  = "2009"
MHW_VALUE_M = 0.26
MHW_NOTE    = f"MHW = {MHW_VALUE_M:.2f} m (NAVD88)"
SENTINEL_DAM = -0.3

# ----- orientation controls -----
FLIP_DOMAIN_ALONGSHORE = True   # prevent per-domain mirroring
DOMAIN_SORT_DESC       = True    # sort 134..30 so 134 renders at BOTTOM (origin="lower")
STACK_FLIP_VERTICAL    = False   # set True to invert whole stack (top↔bottom) after sorting

# --- Base pixel/metric controls (set by preset below) ------------------------
SQUARE_IN_INDEX  = True              # if VIEW_PRESET == "index"
SQUARE_IN_METERS = False             # if VIEW_PRESET == "metric"
CROSS_SHORE_DX_M = None              # set by preset (metric)
ALONGSHORE_DY_M  = None              # set by preset (metric)
SHOW_METRIC_AXES = True

# --- figure sizing for long, skinny islands ----------------------------------
AUTO_SIZE            = True
BASE_WIDTH_IN        = 7.0           # index-view base width
METRIC_BASE_WIDTH_IN = 9.5           # metric view base width
MIN_HEIGHT_IN        = 8.0
MAX_HEIGHT_IN        = 28.0

# --- domain separators / labels in stitched mode -----------------------------
DRAW_SEPARATORS      = True
LABEL_EVERY_N_DOMAINS= 10
EXCLUDE_DOMAIN_LABELS= {41, 51}

# --- metric x-axis readability controls --------------------------------------
METRIC_X_NTICKS      = 3
METRIC_X_LABEL_ROT   = 0

# -----------------------------------------------------------------------------


N_ROWS_EXPECT  = 200   # 200-row extractor
ALONG_EXPECT   = 50

def apply_view_preset():
    """Set flags for 'index' vs 'metric' views."""
    global SQUARE_IN_INDEX, SQUARE_IN_METERS, CROSS_SHORE_DX_M, ALONGSHORE_DY_M, SHOW_METRIC_AXES
    if VIEW_PRESET.lower() == "metric":
        SQUARE_IN_INDEX  = False
        SQUARE_IN_METERS = True
        CROSS_SHORE_DX_M = 10.0     # meters
        ALONGSHORE_DY_M  = 10.0     # meters
        SHOW_METRIC_AXES = True
    elif VIEW_PRESET.lower() == "index":
        SQUARE_IN_INDEX  = False
        SQUARE_IN_METERS = False
        SHOW_METRIC_AXES = True
    else:
        raise ValueError("VIEW_PRESET must be 'index' or 'metric'.")

def domain_num_from_name(stem: str) -> int:
    m = re.search(r"domain_(\d+)", stem)
    return int(m.group(1)) if m else 10**9

def nice_domain_label(stem: str, fallback_index: int) -> str:
    m = re.search(r"domain_(\d+)", stem)
    num = m.group(1) if m else str(fallback_index)
    return f"Domain {num}"

def mask_sentinel(a, sentinel=SENTINEL_DAM):
    out = a.astype(float).copy()
    out[out <= sentinel + 1e-6] = np.nan
    return out

def per_domain_flip(topo_dm: np.ndarray, dune_dm: np.ndarray):
    if not FLIP_DOMAIN_ALONGSHORE:
        return topo_dm, dune_dm
    return topo_dm[:, ::-1], dune_dm[::-1]

def _shape_warn(stem, top, dun):
    if top.ndim != 2 or dun.ndim != 1:
        print(f"[warn] {stem}: unexpected dims top.ndim={top.ndim}, dune.ndim={dun.ndim}")
    if top.shape[1] != dun.shape[0]:
        print(f"[warn] {stem}: alongshore mismatch top={top.shape[1]} vs dune={dun.shape[0]}")
    if top.shape[0] != N_ROWS_EXPECT:
        print(f"[warn] {stem}: N_ROWS={top.shape[0]} (expected {N_ROWS_EXPECT}).")

def load_single(stem: str):
    dnum = domain_num_from_name(stem)
    if dnum < DOMAIN_MIN_ID or dnum > DOMAIN_MAX_ID:
        print(f"[note] {stem} (#{dnum}) outside requested range {DOMAIN_MIN_ID}-{DOMAIN_MAX_ID}, plotting anyway.")
    topo_p = TOPO_DIR / f"{stem}_topography_{YEAR_SUFFIX}.npy"
    dune_p = DUNE_DIR  / f"{stem}_dune_{YEAR_SUFFIX}.npy"
    if VERBOSE:
        print(f"[load single] topo: {topo_p}")
        print(f"[load single] dune: {dune_p}")
    if not topo_p.exists():
        raise FileNotFoundError(f"Topography file not found: {topo_p}")
    if not dune_p.exists():
        raise FileNotFoundError(f"Dune file not found: {dune_p}")

    topo = np.load(topo_p)
    dune = np.load(dune_p)
    _shape_warn(stem, topo, dune)
    topo, dune = per_domain_flip(topo, dune)
    return topo, dune, [stem], [topo.shape[1]]

def load_all_stitched():
    tops, dunes, stems, widths = [], [], [], []
    cand = sorted(
        TOPO_DIR.glob(f"domain_*_topography_{YEAR_SUFFIX}.npy"),
        key=lambda x: domain_num_from_name(x.stem),
        reverse=DOMAIN_SORT_DESC,   # 134..30 so 134 stacks at bottom
    )
    if VERBOSE:
        print(f"[scan] found {len(cand)} topo files for YEAR_SUFFIX={YEAR_SUFFIX}")

    for p in cand:
        stem = p.stem.replace(f"_topography_{YEAR_SUFFIX}", "")
        dnum = domain_num_from_name(stem)
        if not (DOMAIN_MIN_ID <= dnum <= DOMAIN_MAX_ID):
            if VERBOSE:
                print(f"[skip] {stem} (#{dnum}) outside range {DOMAIN_MIN_ID}-{DOMAIN_MAX_ID}")
            continue

        topo_p = TOPO_DIR / f"{stem}_topography_{YEAR_SUFFIX}.npy"
        dune_p = DUNE_DIR  / f"{stem}_dune_{YEAR_SUFFIX}.npy"
        if not dune_p.exists():
            print(f"[skip] Missing dune for {stem}")
            continue

        if VERBOSE:
            print(f"[load] {stem} → topo={topo_p.name}, dune={dune_p.name}")

        top = np.load(topo_p)
        dun = np.load(dune_p)
        _shape_warn(stem, top, dun)

        top, dun = per_domain_flip(top, dun)

        tops.append(top); dunes.append(dun)
        stems.append(stem); widths.append(top.shape[1])

    if not tops:
        raise RuntimeError(f"No domain outputs in range {DOMAIN_MIN_ID}-{DOMAIN_MAX_ID}. "
                           f"Looked in {TOPO_DIR} and {DUNE_DIR} for '*_{YEAR_SUFFIX}.npy'.")

    if VERBOSE:
        kept = ", ".join(stems)
        print(f"[stitch] domains used ({len(stems)} kept in {DOMAIN_MIN_ID}-{DOMAIN_MAX_ID}): {kept}")

    topo_all = np.hstack(tops)      # stitch along alongshore
    dune_all = np.concatenate(dunes)
    return topo_all, dune_all, stems, widths

def reorient_map_like(topo_dm: np.ndarray, dune_dm: np.ndarray, scale: float):
    topo_dm_masked = mask_sentinel(topo_dm)
    # (N_ROWS x alongshore) -> (alongshore x cross_shore)
    Z_map = (topo_dm_masked * scale).T
    # Put dune at RIGHT (ocean on LEFT) to match your preferred orientation
    Z_map = np.flip(Z_map, axis=1)
    D_map = dune_dm * scale

    # Optional whole-stack vertical flip (top↔bottom) after sorting
    if STACK_FLIP_VERTICAL:
        Z_map = np.flip(Z_map, axis=0)
        D_map = np.flip(D_map, axis=0)

    return Z_map, D_map

def build_title(stitched: bool, domain_stems):
    if stitched:
        return f"Stitched island (domains {DOMAIN_MIN_ID}–{DOMAIN_MAX_ID}) — back-barrier topography, {TITLE_YEAR}"
    return f"{nice_domain_label(domain_stems[0], 1)} — back-barrier topography, {TITLE_YEAR}"

def draw_figure(Z_map, D_map, domain_stems, domain_widths, stitched: bool, out_name: str):
    unit_label = "m" if SHOW_IN_METERS else "dam"
    metric_mode = (SQUARE_IN_METERS and (CROSS_SHORE_DX_M and ALONGSHORE_DY_M))

    # dynamic figure size
    if AUTO_SIZE:
        along, xshore = Z_map.shape
        base_w = METRIC_BASE_WIDTH_IN if metric_mode else BASE_WIDTH_IN
        height = base_w * (along / max(xshore, 1))
        height = float(np.clip(height, MIN_HEIGHT_IN, MAX_HEIGHT_IN))
        figsize = (base_w, height)
    else:
        figsize = FIGSIZE

    fig, ax = plt.subplots(figsize=figsize, dpi=DPI, layout="constrained")

    # robust color limits
    if np.isfinite(Z_map).any():
        vmin = np.nanpercentile(Z_map, 1)
        vmax = np.nanpercentile(Z_map, 99)
    else:
        vmin, vmax = 0, 1

    # ----- draw raster: METRIC vs INDEX -----
    if metric_mode:
        extent = (0, Z_map.shape[1] * CROSS_SHORE_DX_M, 0, Z_map.shape[0] * ALONGSHORE_DY_M)
        im = ax.imshow(Z_map, origin="lower", cmap=CMAP, vmin=vmin, vmax=vmax,
                       interpolation="nearest", resample=False, rasterized=True, extent=extent)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("Cross-shore distance (m)", labelpad=10)
        ax.set_ylabel("Alongshore distance (m)")
        ax.axvline(extent[1], color="k", lw=0.8, alpha=0.35)
        draw_separators_here = False

        n = max(2, int(3))
        ticks = np.linspace(extent[0], extent[1], n)
        ax.set_xticks(ticks); ax.set_xticklabels([f"{t:.0f}" for t in ticks])
        ax.tick_params(axis="x", labelsize=9, rotation=0, pad=2)

        if SHOW_METRIC_AXES and ALONGSHORE_DY_M is not None:
            ax_right = ax.secondary_yaxis("right",
                functions=(lambda m: m/1000.0, lambda km: km*1000.0))
            ax_right.set_ylabel("Alongshore distance (km)")
            ax_right.tick_params(labelsize=9)

    else:
        im = ax.imshow(Z_map, origin="lower", cmap=CMAP, vmin=vmin, vmax=vmax,
                       interpolation="nearest", resample=False, rasterized=True, aspect="auto")
        ax.set_xlabel("Cross-shore sample index (0 = dune crest; landward to left)", labelpad=8)
        ax.set_ylabel(f"Alongshore transect index (Domain {DOMAIN_MAX_ID} at bottom → {DOMAIN_MIN_ID} at top)")
        ax.axvline(Z_map.shape[1] - 0.5, color="k", lw=0.8, alpha=0.35)
        draw_separators_here = True

        if SHOW_METRIC_AXES:
            if CROSS_SHORE_DX_M is not None:
                def idx_to_m_x(i):  return (Z_map.shape[1]-1 - i) * CROSS_SHORE_DX_M
                def m_to_idx_x(m):  return (Z_map.shape[1]-1) - (m / CROSS_SHORE_DX_M)
                ax_top = ax.secondary_xaxis("top", functions=(idx_to_m_x, m_to_idx_x))
                ax_top.set_xlabel("Cross-shore distance landward of dune crest (m)")
                ax_top.xaxis.set_major_locator(MaxNLocator(6))
                ax_top.tick_params(labelsize=9)

            if ALONGSHORE_DY_M is not None:
                ax_right = ax.secondary_yaxis(
                    "right",
                    functions=(lambda i: (i*ALONGSHORE_DY_M)/1000.0,
                               lambda km: (km*1000.0)/ALONGSHORE_DY_M),
                )
                ax_right.set_ylabel("Alongshore distance (km)")
                ax_right.tick_params(labelsize=9)

    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.01)
    cb.set_label(f"Elevation above MHW ({unit_label}); {MHW_NOTE}")
    ax.set_title(build_title(stitched, domain_stems), pad=6)

    # ----- domain separators & labels (INDEX view only) -----
    if (not metric_mode) and stitched and domain_widths and DRAW_SEPARATORS and draw_separators_here:
        cum = np.cumsum(domain_widths)
        for y in cum[:-1]:
            ax.hlines(y - 0.5, xmin=-0.5, xmax=Z_map.shape[1] - 0.5,
                      colors="k", lw=0.4, alpha=0.22)
        y0 = 0
        for idx, (stem, w) in enumerate(zip(domain_stems, domain_widths), start=1):
            dom_num = domain_num_from_name(stem)
            periodic = ((idx - 1) % max(LABEL_EVERY_N_DOMAINS, 1) == 0)
            if periodic and (dom_num not in EXCLUDE_DOMAIN_LABELS):
                y_mid_ax = (y0 + w/2.0) / max(Z_map.shape[0]-1, 1)
                ax.text(-0.08, y_mid_ax, nice_domain_label(stem, idx),
                        transform=ax.transAxes, va="center", ha="right",
                        fontsize=8, rotation=90, alpha=0.85, clip_on=False)
            y0 += w

    if metric_mode:
        ax.yaxis.set_major_locator(MaxNLocator(10))
    else:
        if AUTO_SIZE and Z_map.shape[0] > 1500:
            ax.yaxis.set_major_locator(MaxNLocator(10))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax.tick_params(axis="y", labelsize=8, pad=2)
        ax.xaxis.set_major_locator(MaxNLocator(8))

    ax.tick_params(direction="out", length=4)
    for spine in ("top","right"):
        ax.spines[spine].set_visible(False)

    # save (print absolute paths so you can click them in the console)
    SAVE_DIR.mkdir(parents=True, exist_ok=True)
    png_path = (SAVE_DIR / out_name).resolve()
    fig.savefig(png_path, dpi=DPI, pad_inches=0.02)
    print(f"[saved] {png_path}")
    if SAVE_PDF:
        pdf_path = png_path.with_suffix(".pdf")
        fig.savefig(pdf_path, pad_inches=0.02)
        print(f"[saved] {pdf_path}")
    plt.close(fig)

def main():
    apply_view_preset()

    print(f"[repo ] PROJECT_ROOT={PROJECT_ROOT}")
    print(f"[paths] TOPO_DIR={TOPO_DIR}")
    print(f"[paths] DUNE_DIR={DUNE_DIR}")
    print(f"[save ] SAVE_DIR={SAVE_DIR}")
    if DOMAIN_ID:
        print(f"[mode] single domain: {DOMAIN_ID}")
    else:
        print(f"[mode] stitched island, domains {DOMAIN_MIN_ID}–{DOMAIN_MAX_ID}")

    scale = 10.0 if SHOW_IN_METERS else 1.0

    if DOMAIN_ID:
        topo_dm, dune_dm, stems, widths = load_single(DOMAIN_ID)
        Z_map, D_map = reorient_map_like(topo_dm, dune_dm, scale)
        fname = f"{nice_domain_label(stems[0], 1).replace(' ', '_').lower()}_topography_{TITLE_YEAR}.png"
        draw_figure(Z_map, D_map, stems, widths, stitched=False, out_name=fname)
    else:
        topo_dm, dune_dm, stems, widths = load_all_stitched()
        Z_map, D_map = reorient_map_like(topo_dm, dune_dm, scale)
        fname = f"stitched_island_topography_{TITLE_YEAR}_domains_{DOMAIN_MIN_ID}-{DOMAIN_MAX_ID}.png"
        draw_figure(Z_map, D_map, stems, widths, stitched=True, out_name=fname)

if __name__ == "__main__":
    main()
