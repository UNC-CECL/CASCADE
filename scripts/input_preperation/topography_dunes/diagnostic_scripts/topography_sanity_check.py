# =============================================================================
# CASCADE Dune & Topography Visualizer (Benton-parity outputs, MAP-LIKE VIEW)
# =============================================================================

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# =============================================================
# === USER INPUTS: DATA PATH CONFIGURATION (CRITICAL FIX) ===
# =============================================================

DATA_ROOT_DIR = Path(r'/data/hatteras_init')

# STEP 2: Sub-directories (should not need changing)
# These paths are constructed relative to DATA_ROOT_DIR
TOPO_DIR = DATA_ROOT_DIR / "topography_dunes" / "2009"
DUNE_DIR = DATA_ROOT_DIR / "dunes" / "2009"
YEAR_SUFFIX = "2009"

# Single domain (set to None to stitch all)
DOMAIN_ID = 65  # e.g., 65 or None for whole island

# Display options
SHOW_IN_METERS = True
CMAP = "terrain"

# EXPECTED SHAPES
N_ROWS_EXPECT = 200
ALONG_EXPECT = 50
SENTINEL_DAM = -0.3

# Flip alongshore direction *inside each domain* to match your GIS map view
FLIP_DOMAIN_ALONGSHORE = True


# =============================================================
# ======= HELPERS (No change needed) =======
# =============================================================

def domain_num_from_name(stem: str) -> int:
    m = re.search(r"domain_(\d+)", stem)
    return int(m.group(1)) if m else 10 ** 9


def _apply_per_domain_flip(topo_dm: np.ndarray, dune_dm: np.ndarray):
    """
    Inputs: topo_dm shape (N_ROWS, alongshore), dune_dm shape (alongshore,)
    If FLIP_DOMAIN_ALONGSHORE, reverse alongshore axis inside the domain.
    """
    if not FLIP_DOMAIN_ALONGSHORE:
        return topo_dm, dune_dm
    return topo_dm[:, ::-1], dune_dm[::-1]


def mask_sentinel(a, sentinel=SENTINEL_DAM):
    out = a.astype(float).copy()
    out[out <= sentinel + 1e-6] = np.nan
    return out


def _shape_warn(stem, top, dun):
    if top.ndim != 2 or dun.ndim != 1:
        print(f"[warn] {stem}: unexpected dims top.ndim={top.ndim}, dune.ndim={dun.ndim}")
    if top.shape[1] != dun.shape[0]:
        print(f"[warn] {stem}: alongshore mismatch top={top.shape[1]} vs dune={dun.shape[0]}")
    if top.shape[0] != N_ROWS_EXPECT:
        print(f"[warn] {stem}: N_ROWS={top.shape[0]} (expected {N_ROWS_EXPECT}). "
              f"Plot will still work (auto), but verify extractor settings.")


# ======= LOADING FUNCTIONS (with checks) =======

def load_one(domain_id):
    domain_stem = f"domain_{domain_id}_resampled"
    topo_path = TOPO_DIR / f"{domain_stem}_topography_{YEAR_SUFFIX}.npy"
    dune_path = DUNE_DIR / f"{domain_stem}_dune_{YEAR_SUFFIX}.npy"

    if not topo_path.exists():
        print(f"CRITICAL ERROR: Topography file not found at {topo_path}")
        print("HINT: Check the DATA_ROOT_DIR path defined in the script.")
        sys.exit(1)
    if not dune_path.exists():
        print(f"CRITICAL ERROR: Dune file not found at {dune_path}")
        print("HINT: Check the DATA_ROOT_DIR path defined in the script.")
        sys.exit(1)

    topo = np.load(topo_path)
    dune = np.load(dune_path)
    _shape_warn(domain_stem, topo, dune)
    topo, dune = _apply_per_domain_flip(topo, dune)
    return topo, dune, [domain_stem], [topo.shape[1]]


def load_all_stitched():
    tops, dunes, stems, widths = [], [], [], []

    # Check if directories exist first
    if not TOPO_DIR.exists() or not DUNE_DIR.exists():
        print(f"CRITICAL ERROR: One or both data directories not found:")
        print(f"  Topography Dir: {TOPO_DIR}")
        print(f"  Dune Dir: {DUNE_DIR}")
        print("HINT: Check the DATA_ROOT_DIR path defined in the script.")
        sys.exit(1)

    for p in sorted(TOPO_DIR.glob(f"domain_*_resampled_topography_{YEAR_SUFFIX}.npy"),
                    key=lambda x: domain_num_from_name(x.stem)):
        stem = p.stem.replace(f"_resampled_topography_{YEAR_SUFFIX}", "")
        dune_p = DUNE_DIR / f"{stem}_resampled_dune_{YEAR_SUFFIX}.npy"

        if not dune_p.exists():
            print(f"[skip] Missing dune file for {stem} at {dune_p}")
            continue

        try:
            top = np.load(p)
            dun = np.load(dune_p)
        except Exception as e:
            print(f"[error] Failed to load files for {stem}: {e}")
            continue

        _shape_warn(stem, top, dun)
        top, dun = _apply_per_domain_flip(top, dun)

        tops.append(top)
        dunes.append(dun)
        stems.append(stem)
        widths.append(top.shape[1])

    if not tops:
        raise RuntimeError(f"No domain outputs found in {TOPO_DIR} to stitch. Check file naming and DATA_ROOT_DIR.")

    topo_all = np.hstack(tops)
    dune_all = np.concatenate(dunes)
    return topo_all, dune_all, stems, widths


# =============================================================
# ======= LOAD & PLOT EXECUTION =======
# =============================================================

if DOMAIN_ID:
    print(f"Loading single domain: domain_{DOMAIN_ID}")
    topo_dm, dune_dm, domain_stems, domain_widths = load_one(DOMAIN_ID)
    stitched = False
else:
    print("Loading and stitching ALL domains...")
    topo_dm, dune_dm, domain_stems, domain_widths = load_all_stitched()
    stitched = True

# ======= REORIENT TO MAP-LIKE =======
scale = 10.0 if SHOW_IN_METERS else 1.0
unit_label = "m" if SHOW_IN_METERS else "dam"

topo_dm_masked = mask_sentinel(topo_dm)
# Transpose to (alongshore, cross-shore), then flip columns so dune (row 0 inland cell) appears at RIGHT
Z_map = np.flip((topo_dm_masked * scale).T, axis=1)
D_map = (dune_dm * scale)

# ======= STATS =======
finite_topo = np.isfinite(Z_map)
print(f"\n--- PLOT STATS ---")
print(f"Topography valid fraction: {finite_topo.mean():.2%}  (shape: {Z_map.shape[0]} x {Z_map.shape[1]})")

if finite_topo.any():
    vmin = np.nanpercentile(Z_map, 1)
    vmax = np.nanpercentile(Z_map, 99)
else:
    vmin, vmax = 0, 1

# ======= PLOT (map-like) =======
fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(111)

im = ax.imshow(
    Z_map,
    origin="lower",
    cmap=CMAP,
    aspect="auto",
    vmin=vmin, vmax=vmax
)
cb = plt.colorbar(im, ax=ax, pad=0.01)
cb.set_label(f"Elevation ({unit_label}, MHW-relative)")

# Axes labels
ax.set_ylabel("Alongshore index (domain_1 at bottom → increasing upward)")
ax.set_xlabel("Cross-shore (landward distance from dune → LEFT)   |   Dune/shoreline on RIGHT")
ax.axvline(Z_map.shape[1] - 0.5, color="k", lw=0.8, alpha=0.6)

# Optional: x-ticks in meters (10 m per pixel)
if SHOW_IN_METERS:
    step_cols = 20  # 200 m tick spacing
    xticks = np.arange(0, Z_map.shape[1] + 1, step_cols)
    ax.set_xticks(xticks)
    xticklabels = (Z_map.shape[1] - xticks) * 10
    ax.set_xticklabels([f"{int(v)}" for v in xticklabels])

# Domain separators & labels (horizontal lines because alongshore is vertical)
if stitched and domain_widths:
    cum = np.cumsum(domain_widths)
    for y in cum[:-1]:
        ax.hlines(y - 0.5, xmin=-0.5, xmax=Z_map.shape[1] - 0.5, colors="k", lw=0.5, alpha=0.5)
    y0 = 0
    for stem, w in zip(domain_stems, domain_widths):
        ax.text(-5, y0 + w / 2.0, stem, va="center", ha="right", fontsize=8, rotation=90)
        y0 += w

plt.title((domain_stems[0] if not stitched else "Stitched island") + f" — map-like view ({YEAR_SUFFIX})")
plt.tight_layout()
plt.show()