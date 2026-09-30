"""
The island plan view reduced to one question: where is the unsurveyed ground in what CASCADE runs on?

    python scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/HAT_plot_island_nodata.py

Reads a dune-topo version's topography arrays; writes the plan view beside the
elevation plan view it is compared with, plus zooms. Details: scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

REPO = next(
    _p for _p in Path(__file__).resolve().parents
    if (_p / "pyproject.toml").exists())   # 1-extraction/nodata_audit/ since 2026-09-09
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_topo_version as htv  # noqa: E402

# Config - mirrors HAT_dune_topo_extractor.py

# --- CONFIG ------------------------------------------------------------------
TOPO_PRODUCT = "1984-start"
VERSION_OVERRIDE = None
OFFSET_YEAR = 1984

CELL_SIZE_M = 10.0
DAM_TO_M = 10.0
SENTINEL_M = -3.0
ISLAND_PAD_ROWS = 200
ISLAND_INCLUDE_DUNE = True
NUM_REAL_DOMAINS = 90
N_BUFFER_DOMAINS = 15
FIRST_GIS, LAST_GIS = 1, 90

ZOOM_DOMAINS = (1, 8)              # second figure: inclusive GIS id range

# categories
UNCOVERED = -1
LAND, WATER, DUNE, GAP_OUT, GAP_IN, GAP_TRUNC, LAND_HIDDEN = 0, 1, 2, 3, 4, 5, 6
COLORS = ["#d7d0c0", "#e8eef3", "#8a8378", "#f6b0a8", "#d7191c", "#67000d",
          "#f0a202"]
LABELS = ["land the model sees ($z$ > 0 m MHW)",
          "water the model sees",
          "dune row",
          "unsurveyed, beyond the island (open sound)",
          "unsurveyed, INSIDE the island",
          "unsurveyed AND the cell that stops FindWidths",
          "measured land hidden behind that cell"]
CMAP = ListedColormap(COLORS)
NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5], CMAP.N)

plt.rcParams.update({
    "font.size": 10,
    "axes.linewidth": 0.7,
    "axes.titlesize": 11,
    "axes.labelsize": 10.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "legend.frameon": False,
})


# All outputs in one nodata-audit/ folder beside the extraction; audit_dir() names it
AUDIT_SUBDIR = "nodata-audit"
# -----------------------------------------------------------------------------


# <product>/dune-topo/<version>/nodata-audit/, created on demand
def audit_dir(topo_dir):
    d = topo_dir.parent / AUDIT_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d


# Canvas

# Offset in metres per domain, {domain
def load_offsets(year):
    # The CURRENT build (2026-09-18)
    from site_layer.hat_topo_version import offset_file
    hits = [offset_file(year, "input")]
    if not hits[0].is_file():
        raise SystemExit(f"\nno offset CSV for {year}: {hits[0]}\n")
    v = np.loadtxt(hits[0], skiprows=1, delimiter=",", ndmin=2).astype(float)
    v = v[:, 0]
    if v.size == NUM_REAL_DOMAINS + 2 * N_BUFFER_DOMAINS:
        v = v[N_BUFFER_DOMAINS:N_BUFFER_DOMAINS + NUM_REAL_DOMAINS]
    print(f"[offsets] {year}: {v.size} domains, {v.min():.0f}-{v.max():.0f} m "
          f"({hits[0].name})")
    return {i + 1: float(v[i]) for i in range(v.size)}


# Pad landward to ISLAND_PAD_ROWS with the sentinel, or crop to it
def pad_or_crop(topo_m, nodata):
    n = topo_m.shape[0]
    if n < ISLAND_PAD_ROWS:
        pad = ISLAND_PAD_ROWS - n
        topo_m = np.vstack([topo_m, np.full((pad, topo_m.shape[1]), SENTINEL_M)])
        nodata = np.vstack([nodata,
                            np.zeros((pad, nodata.shape[1]), dtype=bool)])
        return topo_m, nodata, 0
    if n > ISLAND_PAD_ROWS:
        lost = int((topo_m[ISLAND_PAD_ROWS:] > 0.0).sum())
        return topo_m[:ISLAND_PAD_ROWS], nodata[:ISLAND_PAD_ROWS], lost
    return topo_m, nodata, 0


# Category codes for one domain, plus per-profile truncation flags
def classify(topo_m, nodata):
    n_rows, n_cols = topo_m.shape
    water = topo_m <= 0.0
    land = ~water

    last_land = np.where(land.any(axis=0),
                         n_rows - 1 - land[::-1, :].argmax(axis=0), -1)
    idx = np.arange(n_rows)[:, None]
    inside = (idx <= last_land[None, :]) & (last_land[None, :] >= 0)

    cat = np.where(water, WATER, LAND).astype(np.int8)
    cat[nodata & ~inside] = GAP_OUT
    cat[nodata & inside] = GAP_IN

    first_water = np.where(water.any(axis=0), water.argmax(axis=0), n_rows)
    trunc = np.zeros(n_cols, dtype=bool)
    hidden = np.zeros(n_cols, dtype=int)
    for c in range(n_cols):
        r = first_water[c]
        if r < n_rows and nodata[r, c]:
            cat[r, c] = GAP_TRUNC
            trunc[c] = True
            # Measured land behind the truncating cell
            beyond = land[r + 1:, c]
            hidden[c] = int(beyond.sum())
            rows_beyond = np.nonzero(beyond)[0] + r + 1
            cat[rows_beyond, c] = LAND_HIDDEN
    return cat, trunc, first_water, hidden


# Run: classify every domain, draw the island plan view and the zooms
def main():
    topo_dir, dune_dir, version = htv.topo_dirs(TOPO_PRODUCT, VERSION_OVERRIDE)
    out_png = (audit_dir(topo_dir)
               / f"HAT_island_nodata_{version}_{OFFSET_YEAR}_padded.png")
    print(f"product {TOPO_PRODUCT}, version {version}")

    offsets = load_offsets(OFFSET_YEAR)
    domains = [n for n in range(FIRST_GIS, LAST_GIS + 1) if n in offsets]
    off_cells = {n: int(round(offsets[n] / CELL_SIZE_M)) for n in domains}

    blocks, dunes, per_domain, cropped = [], [], [], []
    fw_rows, hidden_all = [], []
    for n in domains:
        topo = np.load(topo_dir / htv.array_name("topography", n)) * DAM_TO_M
        nod = np.load(topo_dir / htv.array_name("nodata", n))
        dune = np.load(dune_dir / htv.array_name("dune", n)) * DAM_TO_M

        topo_p, nod_p, lost = pad_or_crop(topo, nod)
        if lost:
            cropped.append((n, lost))
        cat, trunc, first_water, hidden = classify(topo_p, nod_p)
        blocks.append(cat)
        dunes.append(dune)
        fw_rows.append(first_water)
        hidden_all.append(hidden)
        per_domain.append(dict(domain=n,
                               n_in=int((cat == GAP_IN).sum()
                                        + (cat == GAP_TRUNC).sum()),
                               n_out=int((cat == GAP_OUT).sum()),
                               n_trunc=int(trunc.sum()),
                               hidden=int(hidden.sum())))

    if cropped:
        print(f"[canvas] ISLAND_PAD_ROWS = {ISLAND_PAD_ROWS} crops real land "
              f"from {len(cropped)} domain(s): "
              f"{', '.join(f'D{d}({v})' for d, v in cropped[:8])}"
              + (" ..." if len(cropped) > 8 else ""))

    n_along = blocks[0].shape[1]
    canvas_rows = max(off_cells.values()) + ISLAND_PAD_ROWS + 5
    total_cols = n_along * len(blocks)
    canvas = np.full((canvas_rows, total_cols), UNCOVERED, dtype=np.int8)

    col = 0
    for k, n in enumerate(domains):
        g = blocks[k]
        origin = off_cells[n]
        end = min(origin + g.shape[0], canvas_rows)
        canvas[origin:end, col:col + n_along] = g[:end - origin, :]
        if ISLAND_INCLUDE_DUNE and origin >= 1:
            canvas[origin - 1, col:col + n_along] = DUNE
        col += n_along

    cm = np.ma.masked_where(canvas < 0, canvas)

    n_in = sum(d["n_in"] for d in per_domain)
    n_out = sum(d["n_out"] for d in per_domain)
    n_trunc = sum(d["n_trunc"] for d in per_domain)
    with_in = [d["domain"] for d in per_domain if d["n_in"]]

    # Draw

    fig = plt.figure(figsize=(20, 10.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.05, 1.0], hspace=0.16,
                          left=0.055, right=0.988, top=0.838, bottom=0.075)

    dom_centre = np.arange(len(domains)) * n_along + n_along / 2.0
    ticks = [i for i, n in enumerate(domains) if n % 5 == 0 or n == 1]

    axA = fig.add_subplot(gs[0])
    axA.imshow(cm, cmap=CMAP, norm=NORM, origin="lower", aspect="auto",
               interpolation="nearest",
               extent=[0, total_cols, 0, canvas_rows])
    axA.set_facecolor("white")
    axA.set_ylabel("Cross-shore cell (raw_offset frame)")
    axA.set_xlim(0, total_cols)
    axA.set_xticks(dom_centre[ticks])
    axA.set_xticklabels([str(domains[t]) for t in ticks])
    axA.tick_params(labelbottom=False)
    axA.set_title("(a)  Every cell CASCADE starts from. Elevation is discarded: "
                  "only the unsurveyed ground is coloured.",
                  loc="left", fontsize=11)

    axB = fig.add_subplot(gs[1], sharex=axA)
    a_in = np.array([d["n_in"] for d in per_domain], float)
    a_out = np.array([d["n_out"] for d in per_domain], float)
    t_all = np.array([d["n_trunc"] for d in per_domain], float)
    w = n_along * 0.88
    axB.bar(dom_centre, a_in, width=w, color=COLORS[GAP_IN], lw=0,
            label="inside the island")
    axB.bar(dom_centre, a_out, width=w, bottom=a_in, color=COLORS[GAP_OUT],
            lw=0, label="beyond it, in the open sound")
    axB.set_ylabel("unsurveyed cells")
    axB.set_xlabel("GIS domain (south → north)")
    axB.set_xticks(dom_centre[ticks])
    axB.set_xticklabels([str(domains[t]) for t in ticks])
    axB.grid(axis="y", lw=0.4, color="0.88")
    axB.set_axisbelow(True)
    axB.set_ylim(0, max((a_in + a_out).max() * 1.34, 1))

    # Truncated profiles on their own axis: a count of profiles, not cells
    axT = axB.twinx()
    axT.plot(dom_centre[t_all > 0], t_all[t_all > 0], ls="none", marker="o",
             ms=4.5, color=COLORS[GAP_TRUNC],
             label="profiles truncated by one (of 50)")
    axT.set_ylabel("profiles truncated\n(of 50)", color=COLORS[GAP_TRUNC])
    axT.tick_params(axis="y", colors=COLORS[GAP_TRUNC])
    axT.set_ylim(0, max(t_all.max() * 1.34, 1))
    axT.spines["right"].set_color(COLORS[GAP_TRUNC])

    h1, l1 = axB.get_legend_handles_labels()
    h2, l2 = axT.get_legend_handles_labels()
    axB.legend(h1 + h2, l1 + l2, loc="upper right", ncol=3, frameon=True,
               framealpha=0.95, edgecolor="0.8", fancybox=False, fontsize=9.5)
    axB.set_title("(b)  Counted per domain, because one cell is a third of a "
                  "pixel wide above and can be invisible there.",
                  loc="left", fontsize=11)

    for ax in (axA, axB):
        for s in ax.spines.values():
            s.set_color("0.35")

    fig.text(0.055, 0.982,
             f"Unsurveyed cells in the CASCADE input  |  {OFFSET_YEAR} offsets  "
             f"|  {TOPO_PRODUCT} dune-topo {version}, padded to "
             f"{ISLAND_PAD_ROWS} cells / "
             f"{ISLAND_PAD_ROWS * CELL_SIZE_M:.0f} m cross-shore",
             fontsize=13, fontweight="bold", va="top", ha="left")
    fig.text(0.055, 0.953,
             f"{n_in + n_out:,} unsurveyed cells are in the arrays CASCADE reads. "
             f"{n_in:,} lie INSIDE the island, in {len(with_in)} of "
             f"{len(domains)} domains; the other {n_out:,} are open sound "
             f"beyond the bay margin, where no lidar return is the expected "
             f"result and nothing about the barrier changes.\n"
             f"None is flagged - Barrier3D has no \"unknown\", so each is read "
             f"as an elevation of exactly $-$3.0 m MHW. The ones that "
             f"demonstrably change the run are the {n_trunc:,} of "
             f"{len(domains) * n_along:,} profiles where an unsurveyed cell is "
             f"the first water cell, because FindWidths stops there.",
             fontsize=10, color="0.28", va="top", ha="left", linespacing=1.5)

    handles = [Patch(facecolor=c, edgecolor="0.45", lw=0.5, label=l)
               for c, l in zip(COLORS, LABELS)]
    fig.legend(handles=handles, loc="upper left", ncol=7, frameon=False,
               fontsize=9, bbox_to_anchor=(0.055, 0.892))

    fig.savefig(out_png, dpi=190, facecolor="white")
    print(f"wrote {out_png}\n")

    print(f"  unsurveyed, inside the island  : {n_in:,}")
    print(f"  unsurveyed, open sound beyond  : {n_out:,}")
    print(f"  domains with any inside        : {len(with_in)} of {len(domains)}")
    print(f"  profiles truncated by one      : {n_trunc:,} of "
          f"{len(domains) * n_along:,}")
    print("  worst domains, by unsurveyed cells INSIDE the island:")
    for d in sorted(per_domain, key=lambda d: -d["n_in"])[:8]:
        print(f"    domain {d['domain']:>3}   {d['n_in']:>5,} inside   "
              f"{d['n_out']:>5,} beyond   {d['n_trunc']:>2} profiles truncated")
    clean = [d["domain"] for d in per_domain if not d["n_in"]]
    print(f"  domains with none inside       : {len(clean)}"
          + (f"  -> {clean}" if 0 < len(clean) <= 12 else ""))
    print()

    # Zoom

    z0, z1 = ZOOM_DOMAINS
    zi = [k for k, n in enumerate(domains) if z0 <= n <= z1]
    if zi:
        zoom_png = (audit_dir(topo_dir)
                    / f"HAT_island_nodata_{version}_{OFFSET_YEAR}"
                      f"_D{z0}-{z1}.png")
        draw_zoom(canvas, domains, zi, off_cells, fw_rows, hidden_all,
                  per_domain, n_along, version, zoom_png)


# The same canvas, cropped to a few domains, at one pixel per cell
def draw_zoom(canvas, domains, zi, off_cells, fw_rows, hidden_all, per_domain,
              n_along, version, out_png):
    c0 = zi[0] * n_along
    c1 = (zi[-1] + 1) * n_along
    sub = canvas[:, c0:c1]
    rows = np.nonzero((sub >= 0).any(axis=1))[0]
    r0, r1 = max(int(rows.min()) - 3, 0), min(int(rows.max()) + 4, sub.shape[0])
    sub = sub[r0:r1]
    cm = np.ma.masked_where(sub < 0, sub)

    # FindWidths boundary in cropped-canvas coordinates
    fw_x, fw_y = [], []
    for j, k in enumerate(zi):
        origin = off_cells[domains[k]]
        for p in range(n_along):
            fw_x += [j * n_along + p, j * n_along + p + 1]
            fw_y += [origin + fw_rows[k][p] - r0] * 2

    hid = np.concatenate([hidden_all[k] for k in zi]).astype(float) * CELL_SIZE_M

    fig = plt.figure(figsize=(17, 11))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.9, 1.0], hspace=0.16,
                          left=0.062, right=0.986, top=0.855, bottom=0.070)

    axA = fig.add_subplot(gs[0])
    axA.imshow(cm, cmap=CMAP, norm=NORM, origin="lower", aspect="auto",
               interpolation="nearest",
               extent=[0, sub.shape[1], r0, r0 + sub.shape[0]])
    axA.plot(fw_x, np.array(fw_y) + r0, color="#111111", lw=1.1, zorder=6,
             label="FindWidths boundary: the model's island ends here")
    axA.set_facecolor("white")
    axA.set_ylabel("Cross-shore cell (raw_offset frame)")
    axA.set_xlim(0, sub.shape[1])
    axA.legend(loc="lower left", fontsize=9.5, frameon=True, framealpha=0.95,
               edgecolor="0.8", fancybox=False)
    axA.set_title(f"(a)  Domains {domains[zi[0]]}-{domains[zi[-1]]} at one pixel "
                  f"per 10 m cell. Amber is measured land Barrier3D never sees, "
                  f"because an unsurveyed cell stopped the scan below it.",
                  loc="left", fontsize=11)

    axB = fig.add_subplot(gs[1], sharex=axA)
    axB.bar(np.arange(hid.size) + 0.5, hid, width=1.0, color=COLORS[LAND_HIDDEN],
            lw=0, label="land hidden behind an unsurveyed cell")
    axB.set_ylabel("hidden land (m)")
    axB.set_xlabel(f"Alongshore profile, 50 per domain "
                   f"(S $\\rightarrow$ N)")
    axB.grid(axis="y", lw=0.4, color="0.88")
    axB.set_axisbelow(True)
    axB.legend(loc="upper right", fontsize=9.5, frameon=True, framealpha=0.95,
               edgecolor="0.8", fancybox=False)
    axB.set_ylim(0, max(hid.max() * 1.25, 10))
    axB.set_title("(b)  Per profile, how much measured land the truncation "
                  "removes from the island width.",
                  loc="left", fontsize=11)

    for ax in (axA, axB):
        ax.set_xticks([(j + 0.5) * n_along for j in range(len(zi))])
        ax.set_xticklabels([str(domains[k]) for k in zi])
        for x in range(1, len(zi)):
            ax.axvline(x * n_along, color="0.55", lw=0.7, zorder=7)
        for s in ax.spines.values():
            s.set_color("0.35")
    axA.tick_params(labelbottom=False)
    axB.set_xlabel(f"Domain (S $\\rightarrow$ N), 50 profiles each")

    d_zoom = [per_domain[k] for k in zi]
    n_tr = sum(d["n_trunc"] for d in d_zoom)
    fig.text(0.062, 0.982,
             f"Domains {domains[zi[0]]}-{domains[zi[-1]]}: where the unsurveyed "
             f"cells actually cost the model island  |  {TOPO_PRODUCT} "
             f"dune-topo {version}",
             fontsize=13, fontweight="bold", va="top", ha="left")
    fig.text(0.062, 0.953,
             f"{n_tr} of {len(zi) * n_along} profiles here are truncated by an "
             f"unsurveyed cell, against {sum(d['n_trunc'] for d in per_domain)} "
             f"island-wide - so this reach carries "
             f"{100 * n_tr / max(sum(d['n_trunc'] for d in per_domain), 1):.0f}% "
             f"of the problem in {len(zi)} of {len(per_domain)} domains. "
             f"Hidden land: {hid[hid > 0].mean():.0f} m mean, "
             f"{hid.max():.0f} m worst.",
             fontsize=10, color="0.28", va="top", ha="left")

    handles = [Patch(facecolor=c, edgecolor="0.45", lw=0.5, label=l)
               for c, l in zip(COLORS, LABELS)]
    fig.legend(handles=handles, loc="upper left", ncol=4, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.062, 0.930))

    fig.savefig(out_png, dpi=190, facecolor="white")
    print(f"wrote {out_png}")
    print(f"  profiles truncated in D{domains[zi[0]]}-{domains[zi[-1]]} : "
          f"{n_tr} of {len(zi) * n_along}")
    for d in d_zoom:
        print(f"    domain {d['domain']:>2}   {d['n_trunc']:>2}/50 truncated   "
              f"{d['hidden'] * CELL_SIZE_M / max(d['n_trunc'], 1):>6.0f} m "
              f"hidden per truncated profile")
    print()


if __name__ == "__main__":
    main()
