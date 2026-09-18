"""
HAT_plot_island_nodata.py

The island plan view, stripped to one question: where is the unsurveyed ground
in what CASCADE actually runs on?

WHY A SECOND PLAN VIEW
----------------------
HAT_dune_topo_island_planview_<run>_<year>_padded.png shows elevation, and at
that colour scale an unsurveyed cell is indistinguishable from water: both sit
at the -3.0 m sentinel and both render as the same dark blue. That is not a
flaw in the figure - it is the honest consequence of Barrier3D having no
representation for "unknown" - but it means the elevation view cannot answer
"is any of this no-data affecting my model".

So this draws the same canvas, at the same offsets, with the same padding, and
throws the elevation away. Land is one flat grey, water another, and the only
thing with a colour is the no-data.

THE CANVAS IS THE SAME ONE, DELIBERATELY
-----------------------------------------
Every geometric rule here is copied from _build_island_canvas() in
HAT_dune_topo_extractor.py so the two figures overlay cell for cell:

    offsets     2-brie-offset/<year>/Island_Dune_Offsets_*.csv,
                metres, seaward positive, row 0 = domain 1 (Cape Point).
                A 120-row file is stripped of its 15 buffer domains per end.
    origin      round(offset_m / 10) - the canvas row interior row 0 lands on
    padding     every domain padded landward to ISLAND_PAD_ROWS = 200 cells,
                or cropped to it
    dune        written into canvas row origin - 1, one row, matching
                ISLAND_INCLUDE_DUNE
    columns     domains concatenated in ascending order, 50 profiles each,
                no per-domain flip - the arrays already run south to north

If those constants move in the extractor they must move here. The alternative -
importing the extractor - drags in its interactive picker and its own
TOPO_PRODUCT literal, which is how the figure scripts came to disagree with the
road scripts before.

TWO SHADES OF RED, AND THE DIFFERENCE MATTERS
----------------------------------------------
    unsurveyed              a cell CASCADE reads as -3.0 m water that was
                            never measured
    unsurveyed, truncating   the same, AND it is the first water cell on its
                            profile, so barrier3d.FindWidths stops there

The second is the one with a demonstrable effect. FindWidths measures the
island as the run of land from interior row 0 to the first water cell, and land
behind that cell is invisible to the model. A truncating unsurveyed cell
therefore deletes every real, measured cell behind it from the island width
Barrier3D uses. The rest of the red is inside the barrier and may or may not
matter, depending on what the run does with it.

Panel (b) counts both per domain, so nothing is missed at 45 km: a single
unsurveyed cell is a third of a pixel wide in panel (a) and can be invisible
there while still being a real bar below.

INPUT   <product>/dune-topo/<version>/topography/domain_<N>_topography.npy  dam
        <product>/dune-topo/<version>/topography/domain_<N>_nodata.npy     bool
        <product>/dune-topo/<version>/dunes/domain_<N>_dune.npy            dam
        2-brie-offset/<year>/Island_Dune_Offsets_*.csv            m

        Product and version resolve through scripts/site_layer/hat_topo_version.py.

OUTPUT  <product>/dune-topo/<version>/HAT_dune_topo_island_nodata_<version>_<year>_padded.png
        Written beside the elevation plan view it is meant to be compared with.
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

# =============================================================================
# CONFIG - mirrors HAT_dune_topo_extractor.py
# =============================================================================

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



# Every output of this folder lands under one directory beside the extraction it
# describes, rather than being scattered through the run folder it did not
# produce. audit_dir() is the only place that name is spelled.
AUDIT_SUBDIR = "nodata-audit"


def audit_dir(topo_dir):
    """<product>/dune-topo/<version>/nodata-audit/, created on demand."""
    d = topo_dir.parent / AUDIT_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d


# =============================================================================
# CANVAS
# =============================================================================

def load_offsets(year):
    """Offset in metres per domain, {domain: offset_m}. Mirrors load_offsets()."""
    # The CURRENT build (2026-09-18). This took the first sorted match under
    # 2-brie-offset/, which for 1984 and 2004 is superseded_20260915_flat/ --
    # a build that differs from the current one.
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


def pad_or_crop(topo_m, nodata):
    """Pad landward to ISLAND_PAD_ROWS with the sentinel, or crop to it.

    Returns the padded arrays and how many real land cells the crop discarded,
    which the extractor also reports - a domain wider than 2000 m loses its bay
    margin to this figure's frame, not to the model.
    """
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


def classify(topo_m, nodata):
    """Category codes for one domain, plus per-profile truncation flags.

    Two things are separated here, and the distinction is the whole point of
    the figure:

    INSIDE the island envelope - row 0 up to that profile's last cell above
    MHW - an unsurveyed cell is a hole in the barrier. Beyond it, the same cell
    is open sound the survey never flew over, which is the expected state of a
    lidar return over water and changes nothing about the barrier.

    first_water is barrier3d.FindWidths' stopping point: the first cell at or
    below sea level, scanning landward from interior row 0. Sea level is 0 in
    the Lagrangian frame, and these arrays are MHW-relative. When that cell is
    unsurveyed, the profile's island is truncated there and every measured cell
    behind it is invisible to the model.
    """
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
            # Measured land behind the truncating cell. Flagged ONLY on
            # profiles an unsurveyed cell truncated: land behind a genuine
            # water cell is also invisible to FindWidths, but that is a real
            # bay, not a data artefact, and colouring it here would blame the
            # survey for the island's actual shape.
            beyond = land[r + 1:, c]
            hidden[c] = int(beyond.sum())
            rows_beyond = np.nonzero(beyond)[0] + r + 1
            cat[rows_beyond, c] = LAND_HIDDEN
    return cat, trunc, first_water, hidden


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

    # =========================================================================
    # DRAW
    # =========================================================================
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

    # Truncated profiles are a count of profiles, not of cells, so they get
    # their own axis rather than being stacked onto a bar they do not belong on.
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

    # =========================================================================
    # ZOOM
    # =========================================================================
    z0, z1 = ZOOM_DOMAINS
    zi = [k for k, n in enumerate(domains) if z0 <= n <= z1]
    if zi:
        zoom_png = (audit_dir(topo_dir)
                    / f"HAT_island_nodata_{version}_{OFFSET_YEAR}"
                      f"_D{z0}-{z1}.png")
        draw_zoom(canvas, domains, zi, off_cells, fw_rows, hidden_all,
                  per_domain, n_along, version, zoom_png)


def draw_zoom(canvas, domains, zi, off_cells, fw_rows, hidden_all, per_domain,
              n_along, version, out_png):
    """The same canvas, cropped to a few domains, at one pixel per cell.

    The FindWidths boundary is drawn on top as a step line. Above it, on a
    truncated profile, is measured land the model cannot see - which is the
    whole reason this zoom exists. The step is drawn per profile rather than
    smoothed: it moves by whole cells, and interpolating it would suggest a
    precision the 10 m grid does not have.
    """
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
