"""
road_setback_figures.py
==============================================================================
How the NC-12 road setback the model reads is measured, on one domain, and
what the 1996 and 2010 runs are handed along the reach.

    python scripts/figure_making/pipeline/4-mgmt-forcings/road_setback_figures.py

Writes output/figures/3-model-inputs/4-management/:
    road_setback_measurement.png   one domain: raw grid + rasterised road,
                                   the straightened profiles, the per-profile
                                   setbacks and the domain value
    road_setback_inputs.png        the setbacks the 1996 and 2010 runs read,
                                   the measurements they come from, and the
                                   corrections applied on the way

Everything is read from the producers' saved products (nothing re-measured):
    raster/<line vintage>/masks/domain_N_road_<vintage>.npy
        HAT_rasterize_road_to_domains.py: the NC-12 centreline buffered to
        the road width and burned onto each domain's raw 10 m grid
    dunestart_offset/measured/<year>/RoadOffset_<year>_profiles.csv, _domains.csv
        HAT_road_offset_from_dune_start.py: the mask sheared with the same
        per-profile shear as the topography, then per profile the distance
        from interior row 0 (one cell landward of the picked dune crest) to
        the road's seaward edge; the domain value is the median, then the
        negative floor (ocean side) and the drowning-road move (bay side)
    dunestart_offset/derived/<1996|2010>/RoadSetback_*_dunestart.csv
        HAT_road_setback_derived_vintages.py: 1996 = the 1984 measurement +
        the 1989 Pea Island relocation; 2010 = the 2004 measurement.
Paths resolve through site_layer.hat_topo_version. Ocean on the right.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1984, C_1997, INK, INK_MUTED, CELL_M, DOMAIN_AXIS_LABEL, figsize, figure_dir,
    save, record_caption, _title, open_frame, town_bands,
)

OUT = figure_dir("inputs", "4-management")
EXAMPLE_GIS = 31
EXAMPLE_YEAR = 2004          # the measurement the 2010 run reads unchanged
PRODUCT = {1984: "1984-start", 2004: "2004-start"}


def two_row(path):
    a = np.loadtxt(path, delimiter=",")
    return pd.Series(a[1], index=a[0].astype(int))


def measured_dir(year):
    return tv.road_setback_file(year).parent


def fig_measurement():
    year, gis = EXAMPLE_YEAR, EXAMPLE_GIS
    vintage = tv.road_line_for_year(year)
    raw = np.load(tv.npy_dirs(PRODUCT[year])[0] / f"domain_{gis}.npy")
    mask = np.load(tv.road_mask_file(vintage, gis)).astype(bool)
    prof = pd.read_csv(measured_dir(year) / f"RoadOffset_{year}_profiles.csv")
    prof = prof[prof.domain == gis].sort_values("profile")
    dom = pd.read_csv(measured_dir(year) / f"RoadOffset_{year}_domains.csv").set_index("domain").loc[gis]

    fig = plt.figure(figsize=figsize("double", height=6.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1])

    # (a) the raw grid and the rasterised road
    ax = fig.add_subplot(gs[0, 0])
    z = np.ma.masked_less_equal(raw, -9.9)
    cols_land = np.nonzero((raw > -9.9).any(axis=0))[0]
    c0, c1 = max(cols_land.min(), mask.nonzero()[1].min() - 60), cols_land.max() + 3
    ext = (c0 * CELL_M - 0.5 * CELL_M, (c1 + 0.5) * CELL_M, 0, raw.shape[0] * CELL_M)
    ax.imshow(z[:, c0:c1 + 1], cmap="Greys", vmin=-1.5, vmax=6, origin="lower", aspect="auto",
              extent=ext, interpolation="nearest")
    ax.imshow(np.ma.masked_where(~mask[:, c0:c1 + 1], mask[:, c0:c1 + 1]), cmap=ListedColormap([C["ACCENT"]]),
              origin="lower", aspect="auto", extent=ext, interpolation="nearest", alpha=0.9)
    ax.set_facecolor(C["WATER"])
    ax.set_xlabel("raw grid column (m)  ·  ocean →")
    ax.set_ylabel("alongshore (m)")
    open_frame(ax)
    ax.legend(handles=[Line2D([], [], color=C["ACCENT"], lw=5, label=f"NC-12 ({vintage} line), rasterised")],
              loc="upper left", frameon=True, fontsize=7)
    _title(ax, 0, f"GIS {gis}: the raw 10 m grid")

    # (b) the straightened profiles: dune crest, row 0, road band
    ax = fig.add_subplot(gs[0, 1])
    y = prof.profile.to_numpy() * CELL_M
    ax.fill_betweenx(y, prof.road_seaward_cell * CELL_M, (prof.road_landward_cell + 1) * CELL_M,
                     color=C["ACCENT_FILL"], step="mid", label="road cells")
    ax.step(prof.dune_crest_cell * CELL_M, y, where="mid", color=C["ADDED"], lw=1.2, label="dune crest")
    ax.step(prof.interior_row0_cell * CELL_M, y, where="mid", color=INK, lw=1.2, label="interior row 0")
    ax.step(prof.road_seaward_cell * CELL_M, y, where="mid", color=C["ACCENT"], lw=1.2, label="road, seaward edge")
    k = len(prof) // 2
    r = prof.iloc[k]
    ax.annotate("", xy=(r.road_seaward_cell * CELL_M, y[k]), xytext=(r.interior_row0_cell * CELL_M, y[k]),
                arrowprops=dict(arrowstyle="<->", lw=0.8, color=INK))
    ax.text((r.road_seaward_cell + r.interior_row0_cell) / 2 * CELL_M, y[k] + 12, f"{r.setback_m:.0f} m",
            ha="center", va="bottom", fontsize=7.5)
    ax.invert_xaxis()
    ax.set_xlabel("along the straightened profile (m)  ·  ocean →")
    ax.set_ylabel("profile (m alongshore)")
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7, loc="lower center")
    _title(ax, 1, "Sheared like the topography")

    # (c) per-profile setbacks -> the domain value
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(prof.profile, prof.setback_m, "o", color=C["ACCENT"], ms=3)
    ax.axhline(dom.setback_dunestart_m, color=INK, lw=1.2, label=f"median {dom.setback_dunestart_m:.0f} m")
    ax.axhspan(dom.setback_p10_m, dom.setback_p90_m, color="0.92", lw=0, label="10th-90th percentile")
    ax.set_xlabel("profile (alongshore, 10 m each)")
    ax.set_ylabel("setback from interior row 0 (m)")
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7.5, loc="best")
    _title(ax, 2, "One value per profile, one per domain")

    # (d) the shear it removed: raw road columns vs straightened
    ax = fig.add_subplot(gs[1, 1])
    rows = np.arange(mask.shape[0])
    raw_c = [np.nonzero(mask[i])[0] for i in rows]
    # landward of the road's most seaward cell, in both frames (raw columns
    # grow toward the ocean; straightened cells grow away from it)
    lo = np.array([c.min() if len(c) else np.nan for c in raw_c])
    hi = np.array([c.max() if len(c) else np.nan for c in raw_c])
    ref = np.nanmax(hi)
    ax.fill_betweenx(rows * CELL_M, (ref - hi) * CELL_M, (ref - lo + 1) * CELL_M, color="0.8", step="mid",
                     label=f"raw grid: {(np.nanmax(hi) - np.nanmin(lo) + 1) * CELL_M:.0f} m of cross-shore")
    s0 = prof.road_seaward_cell.min()
    ax.fill_betweenx(y, (prof.road_seaward_cell - s0) * CELL_M, (prof.road_landward_cell + 1 - s0) * CELL_M,
                     color=C["ACCENT"], alpha=0.7, step="mid",
                     label=f"straightened: {dom.measured_road_width_m:.0f} m road")
    ax.invert_xaxis()
    ax.set_xlabel("landward of the road's seaward edge (m)")
    ax.set_ylabel("alongshore (m)")
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    _title(ax, 3, f"Obliquity {dom.obliquity_deg:.0f}°")

    out = save(fig, OUT / "road_setback_measurement.png")
    plt.close(fig)
    record_caption(out[0],
        f"How the NC-12 setback is measured on one domain: GIS {gis}, the {year} measurement (the "
        f"{vintage}-digitised centreline on the {PRODUCT[year]} extraction), which the 2010 run reads unchanged. "
        "(a) The domain's raw 10 m elevation grid (grey, ocean at the right, water blue) with the road "
        "centreline buffered to its width and burned onto it (purple; HAT_rasterize_road_to_domains.py). "
        f"NC-12 crosses this north-up box at {dom.obliquity_deg:.0f}°, so the road cells smear across many "
        "columns. (b) The same road after shearing each alongshore profile by the shear the topography "
        "extractor applied, drawn against the picked dune crest and interior row 0 (one cell landward of the "
        "crest), the row CASCADE counts the setback from (roadway_manager: road_start = setback / dy). The "
        "setback of a profile is row 0 to the road's seaward edge. (c) Every profile's setback (points) and "
        f"their median, {dom.setback_dunestart_m:.0f} m, which is the domain's value. (d) What the shear buys: the "
        "road occupies a wide cross-shore band in the raw grid and one road width once straightened. Two "
        "corrections follow in the producer (HAT_road_offset_from_dune_start.py) and are shown per domain in "
        "road_setback_inputs: negative setbacks floored to row 0, and roads that would drown at "
        "initialisation moved seaward.")
    return out


def fig_inputs():
    gis = np.arange(1, 91)
    s96 = two_row(tv.road_setback_file(1996)).reindex(gis)
    s10 = two_row(tv.road_setback_file(2010)).reindex(gis)
    d84 = pd.read_csv(measured_dir(1984) / "RoadOffset_1984_domains.csv").set_index("domain").reindex(gis)
    d04 = pd.read_csv(measured_dir(2004) / "RoadOffset_2004_domains.csv").set_index("domain").reindex(gis)

    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.6), sharex=True, constrained_layout=True)
    for ax, (label, model, meas, col, i) in zip(axes, [
            ("1996 run: the 1984 measurement + the 1989 relocation", s96, d84, C_1984, 0),
            ("2010 run: the 2004 measurement", s10, d04, C_1997, 1)]):
        ax.plot(gis, meas.setback_dunestart_m, "o", mfc="none", mec="0.55", ms=3.5, label="measured (median)")
        neg = meas.setback_dunestart_m < 0
        ax.plot(gis[neg], meas.setback_dunestart_m[neg], "v", color=C["ADDED"], ms=5,
                label="negative, floored to row 0")
        moved = meas.relocated_seaward_m.fillna(0) > 0
        ax.plot(gis[moved], meas.setback_dunestart_floored_m[moved], "^", color=C["ACCENT"], ms=5,
                label="drowned at start, moved seaward")
        ax.step(gis, model, where="mid", color=col, lw=1.4, label=f"model input, {(1996, 2010)[i]} run")
        ax.fill_between(gis, meas.setback_p10_m, meas.setback_p90_m, step="mid", color=col, alpha=0.15, lw=0,
                        label=f"10th-90th percentile of profiles, {(1984, 2004)[i]}")
        ax.axhline(0, color=INK_MUTED, lw=0.5)
        ax.set_ylabel("setback from row 0 (m)")
        town_bands(ax)
        open_frame(ax)
        _title(ax, i, label)
    diff = (s96 - two_row(tv.road_setback_file(1984)).reindex(gis))
    for g in gis[diff.fillna(0).abs() > 0]:
        axes[0].annotate("", xy=(g, s96[g]), xytext=(g, s96[g] - diff[g]),
                         arrowprops=dict(arrowstyle="->", lw=0.8, color=INK))
    axes[1].set_xlabel(DOMAIN_AXIS_LABEL)
    axes[1].set_xlim(0.5, 90.5)
    h, l_ = axes[1].get_legend_handles_labels()
    h0, l0 = axes[0].get_legend_handles_labels()
    keep = dict(zip(l0, h0)) | dict(zip(l_, h))
    fig.legend(list(keep.values()), list(keep.keys()), loc="outside lower center", ncol=3, frameon=False,
               fontsize=7.5)
    out = save(fig, OUT / "road_setback_inputs.png")
    plt.close(fig)
    moved_1996 = [int(g) for g in gis[diff.fillna(0).abs() > 0]]
    record_caption(out[0],
        "The NC-12 setbacks the two runs read, metres landward of interior row 0, and where they come from. "
        "(a) The 1996 run: no NC-12 line of 1996 vintage exists, so its file is the 1984 measurement (the "
        "1978-digitised line on the 1984-start extraction) with the one relocation between the two dates "
        f"applied, the 1989 Pea Island move (arrows, GIS {moved_1996[0]}-{moved_1996[-1]}; displacements read "
        "from HATTERAS_ROAD_EVENTS, the same values the model fires). (b) The 2010 run: the 2004 measurement "
        "(the 2008-digitised line on the 2004-start extraction), copied unchanged because no relocation falls "
        "between 2004 and 2010. Open circles are each domain's measured median over its alongshore profiles "
        "and the band their 10th-90th percentile; the line is the model-facing value after the producer's two "
        "flagged corrections, negative setbacks floored to 0 (down triangles) and roads drowning at start moved "
        "seaward (up triangles; none fire on the current topographies). Domains with no road carry no value, "
        "and GIS 8 is measured but excluded from the managed span (flag EXCLUDED_FROM_SPAN: a 750 m-wide, "
        "scattered road footprint). "
        "GIS 1 is Cape Point, GIS 90 Pea Island.")
    return out


def main():
    apply_style()
    print(fig_measurement()[0].relative_to(REPO))
    print(fig_inputs()[0].relative_to(REPO))


if __name__ == "__main__":
    main()
