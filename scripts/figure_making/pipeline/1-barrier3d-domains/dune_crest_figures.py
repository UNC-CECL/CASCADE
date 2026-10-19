"""
Dune crest height along the island at the two hindcast starts (1996, 2009), domain by domain.

    python scripts/figure_making/pipeline/1-barrier3d-domains/dune_crest_figures.py [--band min_max|p10_p90|both]

Reads each start's dune arrays through hat_topo_version (height above the berm,
dam) and draws, per GIS domain, the median crest of its 50 alongshore cells, a
band for their spread and the lowest cell, against the design height, the Dmaxel
ceiling and Isabel's run-up. One figure per band, each in its own subfolder of
dune_crest_heights/, written twice: beside the arrays in
data/hatteras_init/1-barrier3d-domains/figures/ and published in
output/figures/3-model-inputs/1-domains/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import yaml  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.transforms import blended_transform_factory  # noqa: E402


# --- CONFIG ------------------------------------------------------------------
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_topo_version as tv  # noqa: E402
from cascade_pipeline.coastsat_lowess import lowess_transect_values  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as A  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, C, C_1984, C_1997, C_1984_FILL, C_1997_FILL, INK, INK_MUTED,
    DOMAIN_AXIS_LABEL, figsize, figure_dir, save, record_caption, _title, town_bands,
)
OUT = figure_dir("inputs", "1-domains", "dune_crest_heights")
# The same figures beside the arrays; shared by both products, so not under either one
DATA_OUT = tv.DOMAIN_FIGURES_DIR / "dune_crest_heights"
PARAMS = REPO / "data" / "hatteras_init" / "Hatteras-CASCADE-parameters.yaml"

# Calibration start red, test start blue
STARTS = {1996: (C_1984, C_1984_FILL), 2009: (C_1997, C_1997_FILL)}
# The survey each start's dunes come from, for the panel titles; the product and version go in the caption
DUNE_SOURCE = {1996: "1996 ALACE survey", 2009: "2009–2014 lidar"}
# The spread band: percentiles of a domain's 50 cells, its name in the caption and its legend entry, by subfolder
BANDS = {"min_max": (0, 100, "min–max", "domain min–max of 50 cells"),
         "p10_p90": (10, 90, "10th–90th percentile", "middle 80% of cells (10th–90th percentile)")}

DESIGN_MHW = 3.94       # m MHW, NC-12 safe crest 4.3 m NAVD88 (Velasquez-Montoya et al. 2021)
ISABEL_RHIGH_MHW = 5.22  # m MHW, Isabel 2003 run-up (5.58 m NAVD88), the only storm over 5 m
N_GIS = 90
DOMAIN_M = 500.0        # alongshore length of one domain
LOWESS_DOMAINS = 7      # the project smoothing width
# Change bars, RdBu mid shades (Hannah, 2026-10-19): blue the crest grew, red it lost height;
# lighter than the dark year lines in (a) and (b)
GAIN_C, LOSS_C = "#4393c3", "#d6604d"
# -----------------------------------------------------------------------------


# Datum from the model's own parameter file, so the conversion follows the yaml
def datum():
    p = yaml.safe_load(PARAMS.read_text(encoding="utf-8"))
    mhw, berm, dmaxel = float(p["MHW"]), float(p["BermEl"]), float(p["Dmaxel"])
    return dict(mhw=mhw, berm=berm, berm_mhw=berm - mhw, ceil_mhw=dmaxel - mhw, ceil_navd=dmaxel)


# One start's crests, m MHW, shape (90 domains, 50 cells)
def crests(year, berm_mhw):
    product = tv.product_for_year(year)
    _, dune_dir, ver = tv.topo_dirs(product)
    a = np.array([np.load(dune_dir / tv.array_name("dune", g)).ravel()[:50] for g in range(1, N_GIS + 1)])
    return a * 10.0 + berm_mhw, product, ver


# Shoal zones hatched gold, as in the management ladder
def shoals(ax, label=True):
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for name, (lo, hi) in A.shoal_zones.items():
        ax.axvspan(lo - 0.5, hi + 0.5, facecolor=A.color_shoal, alpha=0.10, zorder=1, hatch="///",
                   edgecolor=A.color_shoal, linewidth=0)
        if label:
            ax.text((lo + hi) / 2.0, 0.03, name, transform=trans, ha="center", va="bottom",
                    fontsize=7, color="#7A5800", style="italic", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))


# A labelled horizontal reference, the label at the right edge
def ref_line(ax, y, text, ls, color, above=True):
    ax.axhline(y, color=color, ls=ls, lw=0.9, zorder=3)
    ax.annotate(text, (1.0, y), xycoords=("axes fraction", "data"), xytext=(-3, 2 if above else -2),
                textcoords="offset points", ha="right", va="bottom" if above else "top",
                fontsize=7, color=color, zorder=8,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))


def fig_crests(band, data, d):
    ceil_mhw = d["ceil_mhw"]
    lo_q, hi_q, band_name, band_legend = BANDS[band]
    gis = np.arange(1, N_GIS + 1)
    fig, axes = plt.subplots(3, 1, figsize=figsize("double", height=8.2), sharex=True,
                             gridspec_kw=dict(height_ratios=[1, 1, 0.75]), constrained_layout=True)

    for i, (ax, year) in enumerate(zip(axes[:2], STARTS)):
        c, product, ver = data[year]
        col, fill = STARTS[year]
        med = np.median(c, axis=1)
        ax.fill_between(gis, np.percentile(c, lo_q, axis=1), np.percentile(c, hi_q, axis=1),
                        step="mid", color=fill, alpha=0.55, lw=0, zorder=2)
        ax.step(gis, med, where="mid", color=col, lw=1.4, zorder=5)
        ax.step(gis, c.min(axis=1), where="mid", color=col, lw=0.8, ls=(0, (3, 1.5)), zorder=4)
        ax.set_ylim(0, 12)
        ax.set_ylabel("Dune crest (m MHW)")
        ref_line(ax, ceil_mhw, f"Dmaxel ceiling {ceil_mhw:.2f} m", "-", INK, above=False)
        ref_line(ax, DESIGN_MHW, f"design height {DESIGN_MHW:.2f} m", (0, (5, 2)), C["REF"], above=False)
        if year == 1996:
            ref_line(ax, ISABEL_RHIGH_MHW, f"Isabel 2003 run-up {ISABEL_RHIGH_MHW:.2f} m", ":", INK_MUTED)
        _title(ax, i, f"{year} start: dunes from the {DUNE_SOURCE[year]}")

    # Change in the domain median between the starts: a bar per domain, the LOWESS-7 trend over it
    ax = axes[2]
    (c96, *_), (c09, *_) = data[1996], data[2009]
    dmed = np.median(c09, axis=1) - np.median(c96, axis=1)
    trend, _ = lowess_transect_values(gis * DOMAIN_M, dmed, LOWESS_DOMAINS)
    ax.bar(gis, dmed, width=0.8, color=np.where(dmed >= 0, GAIN_C, LOSS_C), lw=0, zorder=4)
    ax.axhline(0, color=INK, lw=0.6, zorder=5)
    ax.plot(gis, trend, color=INK, lw=1.6, zorder=6)
    ax.legend(handles=[Patch(facecolor=GAIN_C, label="crest higher in 2009"),
                       Patch(facecolor=LOSS_C, label="crest lower in 2009"),
                       Line2D([], [], color=INK, lw=1.6, label=f"LOWESS {LOWESS_DOMAINS} domains")],
              loc="upper right", ncol=3, fontsize=7)
    ax.set_ylim(-3.5, 3.5)
    ax.set_ylabel("Median change (m)")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    _title(ax, 2, "Change in domain median crest, 2009 − 1996")
    ax.set_xlim(0.5, N_GIS + 0.5)
    ax.set_xticks([1, 10, 20, 30, 40, 50, 60, 70, 80, 90])

    for k, a_ in enumerate(axes):
        town_bands(a_, strip=0.06, label=(k == 0))
        shoals(a_, label=(k == 2))
        a_.grid(axis="y")
        a_.set_axisbelow(True)

    handles = [Line2D([], [], color=INK, lw=1.4, label="domain median"),
               Patch(facecolor="0.8", label=band_legend),
               Line2D([], [], color=INK, lw=0.8, ls=(0, (3, 1.5)), label="domain minimum")]
    fig.legend(handles=handles, loc="outside upper center", ncol=3, frameon=False)

    name = f"dune_crest_heights_1996_2009_{band}.png"
    paths = [OUT / band / name, DATA_OUT / band / name]
    for p in paths:
        save(fig, p)
    plt.close(fig)

    # Numbers the caption quotes
    def n_dom(year, f):
        return int(f(data[year][0]).sum())
    s = {y: dict(above=n_dom(y, lambda c: np.median(c, 1) > ceil_mhw),
                 below=n_dom(y, lambda c: np.median(c, 1) < DESIGN_MHW),
                 low=n_dom(y, lambda c: c.min(1) < DESIGN_MHW)) for y in STARTS}
    p96, v96 = data[1996][1:]
    p09, v09 = data[2009][1:]
    text = (
        f"Dune crest height along Hatteras Island at the two hindcast starts, by GIS domain (1 at Cape Point, "
        f"90 at Pea Island), in m MHW: the dune arrays CASCADE starts from (height above the berm, "
        f"converted with BermEl {d['berm']:.2f} m and MHW {d['mhw']:.2f} m NAVD88). "
        f"(a) 1996 start, {p96} {v96} (2009–2014 lidar with the 1996 ALACE survey grafted on the ocean side); "
        f"(b) 2009 start, {p09} {v09} (the gap-filled 2009–2014 lidar). Each domain is 50 alongshore cells (10 m each), and each cell's crest is the "
        f"highest point in its cross-shore search window (the extractor's find_dunes), so the spread is along the "
        f"dune ridge, not across it. Solid: the median of a domain's 50 crests; band: their {band_name}; "
        f"dashed: the lowest crest, where overwash gets through first. "
        f"Reference lines: the design height {DESIGN_MHW:.2f} m MHW (4.3 m NAVD88, Velasquez-Montoya et al. 2021), "
        f"the island-wide Dmaxel ceiling {ceil_mhw:.2f} m MHW ({d['ceil_navd']:.1f} m NAVD88; dunes do not grow above it), "
        f"and in (a) Isabel's 2003 run-up {ISABEL_RHIGH_MHW:.2f} m MHW. The domain median stands above the ceiling "
        f"in {s[1996]['above']} of 90 domains in 1996 and {s[2009]['above']} in 2009, and below the design height "
        f"in {s[1996]['below']} and {s[2009]['below']}; the lowest cell is below the design height in "
        f"{s[1996]['low']} and {s[2009]['low']} domains. (c) Change in the domain median, 2009 minus 1996: "
        f"bars per domain (blue higher in 2009, red lower; {int((dmed > 0).sum())} of 90 domains higher), "
        f"black the LOWESS fit over {LOWESS_DOMAINS} domains, all 90 smoothed (no raw southern end: "
        f"this is not a shoreline). GIS 84 in 2009 (Mirlo Beach S-curves) sits almost entirely at the extractor's floor "
        f"(0.1 m above the berm). Grey: villages; gold hatching: Avon and Wimble Shoals.")
    for p in paths:
        record_caption(p, text)
    return paths


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--band", choices=[*BANDS, "both"], default="both")
    args = ap.parse_args()
    apply_style()
    d = datum()
    data = {y: crests(y, d["berm_mhw"]) for y in STARTS}
    for band in (BANDS if args.band == "both" else [args.band]):
        print(fig_crests(band, data, d))


if __name__ == "__main__":
    main()
