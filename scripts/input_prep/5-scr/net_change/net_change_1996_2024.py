"""
net_change_1996_2024.py
==============================================================================
Net shoreline change (CoastSat) against net dune-line change, per GIS domain,
over 1996-2024 and its two halves, in METRES. Built 2026-09-18 (Hannah, by
interview: "the total long term shoreline position change and then compare
this to the total change in duneline position").

WHAT IS COMPARED
    Both sides are read from the stored endpoint products, never recomputed:
        shoreline   3-rates/coastsat/endpoint/<window>/   mean CoastSat position
                    within +/-6 months of each dune-line survey date, end minus
                    start
        dune line   3-rates/duneline/endpoint/<window>/   end line minus start
                    line
    Same windows, same survey dates (1997-10-12, 2009-05-30, 2023-07-01
    ASSUMED), seaward positive on both, so the gap is
        beach-width change = shoreline change - dune-line change
    positive where the beach widened (the waterline gained on the dune).
    Both products share the 2009 date, so for each of them the two halves add
    up to the whole exactly; the script checks it.

WHAT IS DRAWN
    Three stacked panels, one per window (1997-2023, 1997-2009, 2009-2023),
    one y axis in metres: the shoreline blue, the dune line red (the house
    pair for FEATURE in duneline_vs_coastsat), the gap between them shaded
    grey. The village spans are a strip along the top of each panel, not the
    usual full-height wash, because the grey gap already shades. Groin and
    piers as hairlines, the offshore shoals as faint hatched amber boxes, the
    model-input fills as bars above the top panel -- the same marks as the
    two halves figures.

OUTPUT   data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/net_change/chains/
    (was 4-comparisons/net_change_1996_2024/ until 2026-09-19)
    net_change_chain_1996_2010_2024.png   also published to
                                          output/figures/shoreline/
    supporting/
        domain_comparison.csv   one row per window x domain: shoreline, dune,
                                beach-width change, whether they agree in sign
        island_summary.csv      per window: means, r, slope, RMSE, sign
                                agreement
        net_change_shoreline_vs_dune.pdf, CAPTIONS.md
    The per-transect CoastSat change (with positions per end window) is in
    3-rates/coastsat/endpoint/<window>/transect_endpoint.csv.

USAGE
    python scripts/input_prep/5-scr/net_change/net_change_1996_2024.py
==============================================================================
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "CoastSat"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

import coastsat_lrr_windows as cw  # noqa: E402
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "duneline_vs_coastsat"))
from duneline_vs_coastsat import beach_width_handles, shade_beach_width  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C_1984, C_1997, DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style, caption,
    figsize, figure_dir, open_frame, save, structures, support_dir, town_bands,
)
from site_layer.hat_observed_rates import (  # noqa: E402
    NET_CHANGE_1996_2024, coastsat_endpoint_csv, dune_endpoint_csv,
)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402

WHOLE = (1996, 2024)
HALVES = [(1996, 2010), (2010, 2024)]
WINDOWS = [WHOLE] + HALVES
N_DOMAINS = 90
# The chain figure of shoreline_vs_duneline/net_change/ since 2026-09-19
# (was 4-comparisons/net_change_1996_2024/net_change_shoreline_vs_dune).
STEM = "net_change_chain_1996_2010_2024"

C_SHORE = C_1997          # "#2166ac", the CoastSat blue of duneline_vs_coastsat
C_DUNE = C_1984           # "#b2182b", its dune red
C_GAP = "0.86"            # the beach-width gap
TOWN_STRIP = 0.055        # village bands as a strip, the gap owns the grey
Y_STEP_M = 10.0
Y_TICK_M = 20.0


def load(window):
    """Per-domain change from both products, plus the survey metadata."""
    s, e = window
    cs = pd.read_csv(coastsat_endpoint_csv(s, e, "domain")).set_index("domain_number")
    du = pd.read_csv(dune_endpoint_csv(s, e, "domain")).set_index("domain_number")
    meta = pd.read_csv(dune_endpoint_csv(s, e, "transect")).iloc[0]
    df = pd.DataFrame({
        "shoreline_change_m": cs["mean_change_m"],
        "dune_change_m": du["mean_change_m"],
        "n_coastsat_transects": cs["n_transects"],
        "n_dune_transects": du["n_transects"],
    }).reindex(range(1, N_DOMAINS + 1))
    df["beach_width_change_m"] = df["shoreline_change_m"] - df["dune_change_m"]
    df["same_sign"] = np.sign(df["shoreline_change_m"]) == np.sign(df["dune_change_m"])
    return df, meta


def summarise(label, df, meta):
    x, y = df["shoreline_change_m"], df["dune_change_m"]
    ok = x.notna() & y.notna()
    slope, icpt = np.polyfit(x[ok], y[ok], 1)
    return {
        "window": label,
        "lines": f"{int(meta['start_vintage'])}-{int(meta['end_vintage'])}",
        "start_date": meta["start_date"], "end_date": meta["end_date"],
        "end_date_assumed": bool(meta["end_date_assumed"]),
        "n_domains": int(ok.sum()),
        "mean_shoreline_change_m": round(x[ok].mean(), 2),
        "mean_dune_change_m": round(y[ok].mean(), 2),
        "mean_beach_width_change_m": round((x - y)[ok].mean(), 2),
        "r_dune_vs_shoreline": round(float(np.corrcoef(x[ok], y[ok])[0, 1]), 3),
        "slope_dune_on_shoreline": round(float(slope), 3),
        "intercept_m": round(float(icpt), 2),
        "rmse_dune_minus_shoreline_m": round(float(np.sqrt(((y - x)[ok] ** 2).mean())), 2),
        "domains_same_sign": int(df.loc[ok, "same_sign"].sum()),
        "domains_shoreline_landward": int((x[ok] < 0).sum()),
        "domains_dune_landward": int((y[ok] < 0).sum()),
    }


def draw(ax, df, half, label):
    ax.set_xlim(0.5, N_DOMAINS + 0.5)
    ax.set_ylim(-half, half)
    town_bands(ax, label=label, strip=TOWN_STRIP)
    cw.draw_shoals(ax, label=False)
    ax.axhline(0, color=INK_MUTED, lw=0.6, zorder=2)
    x = df.index.to_numpy(float)
    ys, yd = df["shoreline_change_m"].to_numpy(float), df["dune_change_m"].to_numpy(float)
    shade_beach_width(ax, x, ys, yd)
    ax.plot(x, yd, color=C_DUNE, lw=1.1, zorder=5)
    ax.plot(x, ys, color=C_SHORE, lw=1.1, zorder=5)
    structures(ax, label)
    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_minor_locator(MultipleLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(Y_TICK_M))
    ax.yaxis.grid(True, zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)


def caption_text(summ, half, fills):
    by = {r["window"]: r for r in summ}
    w = by[f"{WHOLE[0]}_{WHOLE[1]}"]
    fill_txt = "; ".join(f"{y} at GIS {lo}–{hi}" for y, lo, hi in fills)
    shoal_txt = "; ".join(f"{n} GIS {lo}–{hi}" for n, (lo, hi)
                          in HATTERAS_ANNOTATIONS.shoal_zones.items())
    return (
        "Net change in shoreline and dune-line position by GIS domain (1 at Cape "
        "Point, 90 at Pea Island), seaward positive, in metres: (a) 1997–2023, "
        "(b) 1997–2009, (c) 2009–2023, the digitized dune lines standing in for "
        "the model years 1996, 2010 and 2024. Blue: the CoastSat shoreline, the "
        "mean satellite position within six months of each dune-line image date "
        "(1997-10-12, 2009-05-30, and 2023-07-01, assumed; the 2023 flight date "
        "is not known), end minus start, averaged over the ~10 CoastSat transects "
        "of each 500 m domain. Red: the dune line, end line minus start line "
        "along the 100 m transects, averaged over the ~5 of each domain. The gap "
        "between them is the change in beach width: solid grey where the beach "
        "widened (blue above red), hatched where it narrowed. Both are measured between the same "
        "dates, so (b) and (c) add up to (a) for each. Over 1997–2023 the "
        f"shoreline moved {w['mean_shoreline_change_m']:+.1f} m on average and "
        f"the dune line {w['mean_dune_change_m']:+.1f} m; r = "
        f"{w['r_dune_vs_shoreline']:.2f}, and the two agree in sign in "
        f"{w['domains_same_sign']} of {w['n_domains']} domains. Black bars above "
        f"(a) mark the beach fills inside the window at the footprint the "
        f"hindcast uses ({fill_txt}); hatched amber boxes mark the offshore "
        f"shoals ({shoal_txt}). Village spans are the grey strip along the top of "
        "each panel; the solid hairline is the Buxton groin and the dotted "
        "hairlines are the Avon and Rodanthe piers. All panels share a y axis of "
        f"±{half:g} m, the smallest multiple of {Y_STEP_M:g} m that holds every "
        "value.")


def main() -> int:
    frames, summ, rows = {}, [], []
    for w in WINDOWS:
        df, meta = load(w)
        label = f"{w[0]}_{w[1]}"
        frames[w] = (df, meta)
        summ.append(summarise(label, df, meta))
        rows.append(df.assign(window=label).rename_axis("domain_number").reset_index())

    # the halves add up to the whole on both sides (shared 2009 date)
    whole = frames[WHOLE][0]
    for col in ("shoreline_change_m", "dune_change_m"):
        resid = (frames[HALVES[0]][0][col] + frames[HALVES[1]][0][col] - whole[col]).abs().max()
        assert resid < 0.01, f"{col}: halves do not add up to the whole ({resid:.4f} m)"

    out = NET_CHANGE_1996_2024
    sup = support_dir(out)
    table = pd.concat(rows)[["window", "domain_number", "shoreline_change_m",
                             "dune_change_m", "beach_width_change_m", "same_sign",
                             "n_coastsat_transects", "n_dune_transects"]]
    table.round(3).to_csv(sup / "domain_comparison.csv", index=False)
    pd.DataFrame(summ).to_csv(sup / "island_summary.csv", index=False)

    apply_style()
    extreme = max(float(np.nanmax(np.abs(df[c]))) for df, _ in frames.values()
                  for c in ("shoreline_change_m", "dune_change_m"))
    half = float(math.ceil(extreme / Y_STEP_M) * Y_STEP_M)
    fills = cw.fills_in(int(frames[WHOLE][1]["start_vintage"]),
                        int(frames[WHOLE][1]["end_vintage"]))

    fig, axes = plt.subplots(3, 1, sharex=True, sharey=True, constrained_layout=True,
                             figsize=figsize("double", height=6.6))
    for i, (ax, w) in enumerate(zip(axes, WINDOWS)):
        df, meta = frames[w]
        draw(ax, df, half, label=(i == 0))
        _title(ax, i, f"{int(meta['start_vintage'])}–{int(meta['end_vintage'])}")
    cw.draw_fills(axes[0], fills, half)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel("Net change in position (m)", fontsize=9)
    fig.legend(handles=[Line2D([], [], color=C_SHORE, lw=1.2,
                               label="Shoreline change (CoastSat endpoint)"),
                        Line2D([], [], color=C_DUNE, lw=1.2,
                               label="Dune-line change (endpoint)")]
               + beach_width_handles(),
               loc="outside lower center", ncol=4, frameon=False)
    caption(fig, caption_text(summ, half, fills))
    written = save(fig, out / STEM)
    written += save(fig, figure_dir("shoreline") / STEM)
    plt.close(fig)

    print(pd.DataFrame(summ)[["window", "lines", "mean_shoreline_change_m",
                              "mean_dune_change_m", "mean_beach_width_change_m",
                              "r_dune_vs_shoreline", "slope_dune_on_shoreline",
                              "rmse_dune_minus_shoreline_m", "domains_same_sign"]]
          .to_string(index=False))
    print(f"y bounds  +/-{half:g} m")
    for p in written:
        print("wrote    ", p.relative_to(_REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
