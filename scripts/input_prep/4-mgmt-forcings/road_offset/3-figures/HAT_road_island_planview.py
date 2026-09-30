"""
Where NC-12 sits on the island CASCADE runs, one plan view per vintage: the model's road, not the measured one.

    python scripts/input_prep/4-mgmt-forcings/road_offset/3-figures/HAT_road_island_planview.py

Writes one figure per vintage (PNG, with a vector PDF beside it). Details: scripts/input_prep/4-mgmt-forcings/road_offset/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as pe
from matplotlib.colors import FuncNorm, LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch, Rectangle

REPO = next(_p for _p in Path(__file__).resolve().parents
            if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_topo_version as htv  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, INK, INK_MUTED, GRID_C, apply_style, caption,
    figsize, save)

# --- CONFIG ------------------------------------------------------------------
INIT_ROOT = REPO / "data" / "hatteras_init"
from site_layer import hat_topo_version as _tv  # noqa: E402
ROADS_ROOT = _tv.ROADS_ROOT
DUNESTART = _tv.ROAD_SETBACK_ROOT
ROAD_ELEV_CSV = _tv.ROAD_ELEVATION_FILE

# Config - mirrors HAT_dune_topo_extractor.py and roadway_manager.py

YEARS = [1984, 2004]

CELL_SIZE_M = 10.0                  # dx = dy
DAM_TO_M = 10.0
SENTINEL_M = -3.0
ISLAND_PAD_ROWS = 200
ISLAND_INCLUDE_DUNE = True
NUM_REAL_DOMAINS = 90
N_BUFFER_DOMAINS = 15
FIRST_GIS, LAST_GIS = 1, 90

ROAD_WIDTH_M = 20.0                 # roadway_manager default; 2 cells at dx=10
MHW_M = 0.36
BERM_ELEV_NAVD_M = 1.7
BERM_ELEV_MHW_M = BERM_ELEV_NAVD_M - MHW_M      # 1.34 m MHW

# Island elevation: the extractor's poster ramp, 0 m pinned at 0.35
ISLAND_ELEV_MIN_M = -1.0
ISLAND_ELEV_MAX_M = 4.0
ISLAND_OCEAN_COLOR = C["WATER"]     # the house colour for cells at or below MHW

# Which island ramp: terrain (default) or oleron; each pins sea level differently (README)
ISLAND_CMAP_NAME = os.environ.get("HAT_ISLAND_CMAP", "terrain").lower()
SEA_LEVEL_POS = {"terrain": 0.35, "oleron": 0.50}
# -----------------------------------------------------------------------------


# The island ramp and its norm, with 0 m pinned where that ramp needs it
def island_cmap_norm():
    if ISLAND_CMAP_NAME not in SEA_LEVEL_POS:
        raise SystemExit(f"\nHAT_ISLAND_CMAP={ISLAND_CMAP_NAME!r}; "
                         f"expected one of {sorted(SEA_LEVEL_POS)}\n")
    lo, hi = ISLAND_ELEV_MIN_M, ISLAND_ELEV_MAX_M
    pos = SEA_LEVEL_POS[ISLAND_CMAP_NAME]

    def fwd(x):
        out = np.where(x < 0.0, pos * (x - lo) / (0.0 - lo),
                       pos + (1.0 - pos) * x / hi)
        return np.where(np.isnan(x), np.nan, out)

    def inv(x):
        return np.where(x < pos, lo + (x / pos) * (0.0 - lo),
                        (x - pos) / (1.0 - pos) * hi)

    if ISLAND_CMAP_NAME == "oleron":
        from cmcrameri import cm as cmc          # optional dependency
        cmap = cmc.oleron.copy()
    else:
        cmap = plt.cm.terrain.copy()
    cmap.set_bad(color=ISLAND_OCEAN_COLOR)
    return cmap, FuncNorm((fwd, inv), vmin=lo, vmax=hi)


# Non-default ramps get their own filenames so a comparison keeps both
def out_suffix():
    return "" if ISLAND_CMAP_NAME == "terrain" else f"_{ISLAND_CMAP_NAME}"

# Road elevation in magenta-purple, the one region terrain does not use
ROAD_CMAP = LinearSegmentedColormap.from_list(
    "road_rdpu", plt.cm.RdPu(np.linspace(0.30, 0.97, 256)))
ROAD_VMIN, ROAD_VMAX = 0.3, 1.6     # covers RoadElevation.csv, both vintages
ROAD_EDGE = "#3f0d2e"

# Stroke width for the road, in points
ROAD_LW_PT = 1.6

apply_style()
plt.rcParams.update({
    # Embed fonts as TrueType in the PDF
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# Load

# Offset in metres per domain, {domain
def load_offsets(year):
    # The CURRENT build (2026-09-18), not the first sorted match
    root = _tv.BRIE_ROOT
    hits = [p for p in [_tv.offset_file(year, "input")] if p.is_file()]
    if not hits:
        raise SystemExit(f"\nno offset CSV for {year} under {root}\n")
    v = np.loadtxt(hits[0], skiprows=1, delimiter=",", ndmin=2).astype(float)[:, 0]
    if v.size == NUM_REAL_DOMAINS + 2 * N_BUFFER_DOMAINS:
        v = v[N_BUFFER_DOMAINS:N_BUFFER_DOMAINS + NUM_REAL_DOMAINS]
    print(f"  [offsets] {v.size} domains, {v.min():.0f}-{v.max():.0f} m "
          f"({hits[0].name})")
    return {i + 1: float(v[i]) for i in range(v.size)}


# Read a 2-row (GIS ids, values) CASCADE forcing file -> {domain
def read_two_row_csv(path):
    raw = np.loadtxt(path, delimiter=",")
    if raw.ndim != 2 or raw.shape[0] != 2:
        raise SystemExit(f"{path}: expected 2 rows, got shape {raw.shape}")
    return {int(k): float(v) for k, v in zip(raw[0], raw[1])}


# Pad landward to ISLAND_PAD_ROWS with the sentinel, or crop to it
def pad_or_crop(topo_m):
    n = topo_m.shape[0]
    if n < ISLAND_PAD_ROWS:
        pad = ISLAND_PAD_ROWS - n
        return np.vstack([topo_m,
                          np.full((pad, topo_m.shape[1]), SENTINEL_M)]), 0
    if n > ISLAND_PAD_ROWS:
        return topo_m[:ISLAND_PAD_ROWS], int((topo_m[ISLAND_PAD_ROWS:] > 0).sum())
    return topo_m, 0


# Build

# Everything one figure needs, for one vintage
def build(year):
    product = htv.product_for_year(year)
    topo_dir, dune_dir, version = htv.topo_dirs(product)
    print(f"\n--- {year} -> {product}/{version} " + "-" * 40)

    offsets = load_offsets(year)
    domains = [n for n in range(FIRST_GIS, LAST_GIS + 1) if n in offsets]
    off_cells = {n: int(round(offsets[n] / CELL_SIZE_M)) for n in domains}

    blocks, dune_rows, cropped = [], [], []
    for n in domains:
        topo = np.load(topo_dir / htv.array_name("topography", n)) * DAM_TO_M
        dune = np.load(dune_dir / htv.array_name("dune", n)) * DAM_TO_M
        topo_p, lost = pad_or_crop(topo)
        if lost:
            cropped.append((n, lost))
        blocks.append(topo_p)
        dune_rows.append(np.asarray(dune).reshape(-1, topo_p.shape[1])[0]
                         + BERM_ELEV_MHW_M)

    if cropped:
        print(f"  [canvas] ISLAND_PAD_ROWS={ISLAND_PAD_ROWS} crops real land "
              f"from {len(cropped)} domain(s): "
              + ", ".join(f"D{d}({v})" for d, v in cropped[:6])
              + (" ..." if len(cropped) > 6 else ""))

    n_along = blocks[0].shape[1]
    canvas_rows = max(off_cells.values()) + ISLAND_PAD_ROWS + 5
    total_cols = n_along * len(blocks)

    dem = np.full((canvas_rows, total_cols), np.nan)
    road = np.full((canvas_rows, total_cols), np.nan)

    for k, n in enumerate(domains):
        origin, g = off_cells[n], blocks[k]
        c0, c1 = k * n_along, (k + 1) * n_along
        end = min(origin + g.shape[0], canvas_rows)
        dem[origin:end, c0:c1] = g[:end - origin, :]
        if ISLAND_INCLUDE_DUNE and origin >= 1:
            dem[origin - 1, c0:c1] = dune_rows[k]

    # The road, exactly as roadway_manager builds it
    setback = read_two_row_csv(htv.road_setback_file(year))
    elev = read_two_row_csv(ROAD_ELEV_CSV)
    width_cells = int(ROAD_WIDTH_M / CELL_SIZE_M)

    drawn, missing_elev, segs, seg_z = 0, [], [], []
    for k, n in enumerate(domains):
        if n not in setback:
            continue
        z = elev.get(n, np.nan)
        if not np.isfinite(z):
            missing_elev.append(n)
            continue
        origin = off_cells[n]
        start = origin + int(setback[n] / CELL_SIZE_M)   # int(), as the model does
        stop = min(start + width_cells, canvas_rows)
        if start >= canvas_rows:
            continue
        c0, c1 = k * n_along, (k + 1) * n_along
        road[start:stop, c0:c1] = z
        segs.append((c0, c1, start, stop))
        seg_z.append(z)
        drawn += 1

    print(f"  [road]    {drawn} domains drawn as {width_cells}-cell bands "
          f"({int(ROAD_WIDTH_M)} m), one constant row and elevation each")
    if missing_elev:
        print(f"  [road]    no elevation for {missing_elev} - not drawn")
    ez = [elev[n] for n in setback if n in elev and np.isfinite(elev[n])]
    print(f"  [elev]    {min(ez):.2f} to {max(ez):.2f} m MHW "
          f"(RoadElevation.csv, shared by both vintages)")

    # A domain with no line is ambiguous on the figure
    no_road = [n for n in domains if n not in setback]
    if no_road:
        print(f"  [road]    no road in {len(no_road)} domain(s): "
              f"{min(no_road)}-{max(no_road)}")

    return dict(year=year, product=product, version=version, domains=domains,
                dem=dem, road=road, segs=segs, seg_z=seg_z, canvas_rows=canvas_rows,
                total_cols=total_cols, n_along=n_along, drawn=drawn,
                width_cells=width_cells, no_road=no_road,
                setback_path=htv.road_setback_file(year))


# Draw

# The extractor's poster styling, with the road carrying a second scale
def draw(D):
    year, domains, n_along = D["year"], D["domains"], D["n_along"]
    n_cs, n_al = D["canvas_rows"], D["total_cols"]

    # Drawn at the width it will be printed
    fig_w, fig_h = figsize("double", aspect=0.52)
    fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")
    # Narrower than the extractor's 0.88 to leave a gutter for the cross-shore distance axis
    ax_rect = [0.085, 0.195, 0.715, 0.665]
    ax = fig.add_axes(ax_rect)
    ax.set_facecolor(ISLAND_OCEAN_COLOR)

    cmap, norm = island_cmap_norm()
    im = ax.pcolormesh(np.ma.masked_invalid(D["dem"]), cmap=cmap, norm=norm,
                       shading="auto", rasterized=True)

    # The road is two cells on an 835-row canvas
    rnorm = Normalize(vmin=ROAD_VMIN, vmax=ROAD_VMAX)
    for (c0, c1, r0, r1), z in zip(D["segs"], D["seg_z"]):
        ax.plot([c0, c1], [(r0 + r1) / 2.0] * 2, color=ROAD_CMAP(rnorm(z)),
                lw=ROAD_LW_PT, solid_capstyle="butt", zorder=4,
                path_effects=[pe.Stroke(linewidth=ROAD_LW_PT + 0.8,
                                        foreground=ROAD_EDGE), pe.Normal()])
    ax.plot([], [], color=ROAD_CMAP(0.62), lw=3,
            label=f"NC-12 {year}, as the model builds it")

    # Mark the domains that carry NO road, as a hatched strip along the bottom of the frame
    handles = None
    if D["no_road"]:
        idx = {n: k for k, n in enumerate(domains)}
        bar_h = D["canvas_rows"] * 0.022
        runs, run = [], [D["no_road"][0]]
        for n in D["no_road"][1:]:
            (run.append(n) if n == run[-1] + 1 else (runs.append(run),
                                                     run := [n]))
        runs.append(run)
        for r in runs:
            ax.add_patch(Rectangle((idx[r[0]] * n_along, 0),
                                   len(r) * n_along, bar_h,
                                   facecolor=C["BASE_FILL"], edgecolor=INK_MUTED,
                                   lw=0.5, hatch="///", zorder=6))
        span = ", ".join(f"{r[0]}–{r[-1]}" if len(r) > 1 else f"{r[0]}"
                         for r in runs)
        handles = [*ax.get_legend_handles_labels()[0],
                   Patch(facecolor=C["BASE_FILL"], edgecolor=INK_MUTED,
                         hatch="///", label=f"no NC-12 in domain {span}")]

    ax.legend(handles=handles, loc="upper right", fontsize=7)

    ax.set_xlim(0, n_al)
    ax.set_ylim(0, n_cs)

    # Real distance, alongside the model's own indices

    # The primary axes carry domain number and canvas cell, which is what you need to trace a value back to a file
    km = lambda c: c * CELL_SIZE_M / 1000.0     # noqa: E731
    cell = lambda k: k * 1000.0 / CELL_SIZE_M   # noqa: E731
    sx = ax.secondary_xaxis("top", functions=(km, cell))
    sx.set_xlabel("alongshore distance (km)", labelpad=4)
    sy = ax.secondary_yaxis("right", functions=(km, cell))
    sy.set_ylabel("cross-shore distance (km)", labelpad=4)

    # Vertical exaggeration, computed from the axes actually drawn rather than assumed
    bb = ax.get_position()
    ve = ((bb.height * fig_h) / n_cs) / ((bb.width * fig_w) / n_al)
    # Reported in the footer rather than inside the axes

    # Two colorbars share the extractor's single-colorbar column.
    cax = fig.add_axes([0.895, 0.545, 0.013, 0.315])
    cbar = plt.colorbar(im, cax=cax)
    cbar.set_label("island elevation (m MHW)", color=INK, labelpad=8,
                   rotation=270)
    cbar.ax.yaxis.set_tick_params(color=INK, labelcolor=INK)
    cbar.outline.set_edgecolor(INK_MUTED)
    cbar.outline.set_linewidth(0.6)
    cbar.set_ticks([-1, 0, 1, 2, 3, 4])

    cax_r = fig.add_axes([0.895, 0.195, 0.013, 0.265])
    cbr = plt.colorbar(plt.cm.ScalarMappable(norm=rnorm, cmap=ROAD_CMAP), cax=cax_r)
    cbr.set_label("road elevation (m MHW)", color=INK, labelpad=8, rotation=270)
    cbr.ax.yaxis.set_tick_params(color=INK, labelcolor=INK)
    cbr.outline.set_edgecolor(INK_MUTED)
    cbr.outline.set_linewidth(0.6)
    cbr.ax.axhline(BERM_ELEV_MHW_M, color=INK, lw=1.0, ls=(0, (2.4, 1.6)))

    ticks, labels = [], []
    for k, n in enumerate(domains):
        if n % 5 == 0 or n == 1:
            ticks.append(k * n_along + n_along // 2)
            labels.append(str(n))
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    # The endpoints used to be named here as "Cape Hatteras to Rodanthe"
    ax.set_xlabel(DOMAIN_AXIS_LABEL, labelpad=5)
    ax.set_ylabel("cross-shore cell (1 cell = 10 m)")
    for k, n in enumerate(domains):
        if n % 10 == 0:
            ax.axvline(k * n_along - 0.5, color=GRID_C, lw=0.4, zorder=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("bottom", "left"):
        ax.spines[sp].set_color(INK)

    # NO in-figure title, and no footnote either

    extent = f"{ISLAND_PAD_ROWS} cells / {ISLAND_PAD_ROWS * CELL_SIZE_M:.0f} m"
    nr = D["no_road"]
    cap = (
        f"NC-12 on Hatteras Island as CASCADE initialises it, {year} vintage. "
        f"Domain 1 is at Cape Point in the south, domain 90 at Pea Island in "
        f"the north. "
        f"Barrier interior and dune from the {D['product']}/{D['version']} "
        f"extraction, assembled in the alongshore-offset frame with every "
        f"domain padded to {extent} cross-shore; colour is elevation relative "
        f"to mean high water. NC-12 is drawn as Barrier3D represents it "
        f"(roadway_manager.py:99): per domain a single flat band, "
        f"int(setback/{CELL_SIZE_M:.0f} m) cells landward of interior row 0, "
        f"{D['width_cells']} cells ({ROAD_WIDTH_M:.0f} m) wide, spanning all "
        f"{n_along} alongshore profiles at one constant elevation — so it steps "
        f"between domains rather than following the surveyed centreline. Road "
        f"colour is elevation from RoadElevation.csv, a single file used for "
        f"both vintages, so road colour is identical between the {YEARS[0]} and "
        f"{YEARS[1]} figures by construction and any difference between them is "
        f"the setback or the island, never the roadbed; the dashed mark on that "
        f"scale is the {BERM_ELEV_MHW_M:.2f} m berm. "
        f"{D['drawn']} of {len(domains)} domains carry road"
        + (f"; domains {min(nr)}–{max(nr)} carry none and are hatched along the "
           f"lower frame. " if nr else ". ")
        + f"The road band is stroked at {ROAD_LW_PT:.1f} pt against a true "
        f"width of {(bb.height * fig_h * 72.0) / n_cs * (ROAD_WIDTH_M / CELL_SIZE_M):.1f} pt "
        f"so that it can carry colour; cross-shore position and alongshore "
        f"extent are exact. Vertical exaggeration ×{ve:.1f}. "
        f"Sources: {D['product']}/{D['version']}; setback "
        f"{D['setback_path'].name}; elevation {ROAD_ELEV_CSV.name}; offsets "
        f"Island_Dune_Offsets_{year}_CASCADE_Input.csv; drawn by "
        f"HAT_road_island_planview.py on the {ISLAND_CMAP_NAME} ramp."
    )
    # The caption lands in CAPTIONS.md beside the PNG, keyed by file name
    caption(fig, cap)

    out = DUNESTART / f"HAT_road_island_planview_{year}{out_suffix()}.png"
    # 300 dpi PNG, with a vector PDF beside it
    written = save(fig, out, bbox_inches="tight")
    plt.close(fig)

    true_lw = (bb.height * fig_h * 72.0) / n_cs * (ROAD_WIDTH_M / CELL_SIZE_M)
    print(f"  [road]    stroked {ROAD_LW_PT:.1f} pt vs {true_lw:.1f} pt true "
          f"width -> x{ROAD_LW_PT / true_lw:.2f}")
    print(f"  [scale]   vertical exaggeration x{ve:.2f}")
    for w in written:
        print(f"  [out]     {w}")
    return out


# Run: one plan view per vintage
def main():
    print("=" * 78)
    print("NC-12 on the island CASCADE runs - the model's road, per vintage")
    print("=" * 78)
    for year in YEARS:
        draw(build(year))


if __name__ == "__main__":
    main()
