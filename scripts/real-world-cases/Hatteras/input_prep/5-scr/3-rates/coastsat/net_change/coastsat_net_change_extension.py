"""
The observed 1996-2025 net shoreline change from GIS 1 to Oregon Inlet: the model domains plus the extension.

    python scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change_extension.py

GIS 1-90 come from coastsat_net_change.py's 1996_2025 output. GIS 91-120 (Pea Island,
transect_domain_lookup_ext.csv) have no window means on disk, so their 1996 and 2025
window means are computed here with coastsat_mean_shoreline.window_means, the same
averaging and 10-position minimum. Writes coastsat/net_change/1996_2025/ext/, kept
apart from the model-domain output as the extension is everywhere else.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-07
"""
from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in Path(__file__).resolve().parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as _REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)

import coastsat_mean_shoreline as ms  # noqa: E402
import coastsat_net_change as nc  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK_MUTED, apply_style, figsize, mark_offaxis, offaxis_clause,
    open_frame, record_caption, save, structures, town_bands)
from site_layer.hat_observed_rates import EXT_DIR, net_change_dir, transect_lookup_ext  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
PERIOD = nc.FULL                                      # 1996-2025
N_MODEL = nc.N_DOMAINS                                # GIS 1-90, the model domains
N_ALL = 120                                           # GIS 91-120, the extension to Oregon Inlet
MIN_OBS = ms.DEFAULT_MIN_OBS
# -----------------------------------------------------------------------------


# Window means for the extension transects, in the same columns as mean_shoreline's CSVs
def ext_window_means(lookup, geometry, window):
    df, _ = ms.window_means(lookup, geometry, ms.Window.from_dates(*window), MIN_OBS)
    df = df[df["included"].astype(bool)]
    return df.set_index("transect_id")


# Per transect: end mean minus start mean for the extension transects in both windows
def ext_transect_net_change():
    lookup = pd.read_csv(transect_lookup_ext())[["transect_id", "domain_number"]]
    geometry = ms.transect_geometry(set(lookup["transect_id"]))
    start_w, end_w = nc.WINDOWS[PERIOD]
    s, e = ext_window_means(lookup, geometry, start_w), ext_window_means(lookup, geometry, end_w)
    both = s.index.intersection(e.index)
    out = pd.DataFrame({
        "domain_number": s.loc[both, "domain_number"].astype(int),
        "start_mean_chainage_m": s.loc[both, "mean_chainage_m"],
        "end_mean_chainage_m": e.loc[both, "mean_chainage_m"],
        "start_n_obs": s.loc[both, "n_obs"].astype(int),
        "end_n_obs": e.loc[both, "n_obs"].astype(int),
    })
    out["net_change_m"] = out["end_mean_chainage_m"] - out["start_mean_chainage_m"]
    out["se_net_change_m"] = np.hypot(s.loc[both, "se_chainage_m"], e.loc[both, "se_chainage_m"])
    dropped = len(lookup) - len(out)
    return out.sort_values("domain_number").rename_axis("transect_id"), dropped, len(lookup)


# GIS 1-120: model-domain transects from disk plus the extension, smoothed together like the model side
def domain_net_change_all(ext_tr):
    base_tr = pd.read_csv(net_change_dir(*PERIOD) / "transect_net_change_1996_2025.csv",
                          index_col="transect_id")
    tr = pd.concat([base_tr, ext_tr])
    g = tr.groupby("domain_number")
    idx = pd.RangeIndex(1, N_ALL + 1, name="gis_domain")
    raw = pd.Series(g["net_change_m"].mean().reindex(idx).values, index=idx)
    se = (g["se_net_change_m"].apply(lambda v: np.sqrt((v ** 2).sum())) / g.size()).reindex(idx)
    dom = pd.DataFrame({"n_transects": g.size().reindex(idx).fillna(0).astype(int).values,
                        "net_change_m": raw.values, "se_net_change_m": se.values,
                        "net_change_lowess7_m": nc.smooth_like_model(raw).values}, index=idx)
    dom.insert(0, "reach", np.where(idx <= N_MODEL, "model", "extension"))
    return dom


def figure(dom, out):
    fig, ax = plt.subplots(figsize=figsize("double", height=3.4), constrained_layout=True)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1, N_ALL)
    ax.set_ylim(-nc.POS_HALF, nc.POS_HALF)
    ax.grid(axis="y")
    open_frame(ax)
    town_bands(ax, label=True)
    ax.axvspan(N_MODEL + 0.5, N_ALL, color="0.97", zorder=0)
    ax.axvline(N_MODEL + 0.5, color=INK_MUTED, lw=0.6, ls="--", zorder=1)
    ax.text(N_MODEL + 1.2, 0.03, "beyond the model domains", transform=ax.get_xaxis_transform(),
            ha="left", va="bottom", rotation=90, fontsize=7, color=INK_MUTED)
    x = dom.index
    ax.plot(x, dom["net_change_m"], color=nc.C_RAW, lw=nc.LW_RAW, zorder=5)
    ax.plot(x, dom["net_change_lowess7_m"], color=nc.C_SMOOTH, lw=nc.LW_SMOOTH, zorder=6)
    off = [("domain means", mark_offaxis(ax, x, dom["net_change_m"].values, nc.POS_HALF,
                                         color=nc.C_RAW)),
           ("LOWESS", mark_offaxis(ax, x, dom["net_change_lowess7_m"].values, nc.POS_HALF,
                                   color=nc.C_SMOOTH))]
    ax.set_ylabel("Net shoreline change (m)")
    for y, va, word in ((0.97, "top", "seaward ↑"), (0.03, "bottom", "landward ↓")):
        ax.text(0.995, y, word, transform=ax.transAxes, ha="right", va=va,
                fontsize=7, color=INK_MUTED, zorder=7)
    s, e = nc.WINDOWS[PERIOD]
    c0, c1 = nc.CENTRES[PERIOD]
    island = f"{dom['net_change_m'].mean():+.1f}".replace("-", "−")
    ax.legend(handles=[
        Line2D([], [], color=nc.C_RAW, lw=nc.LW_RAW, label=f"Domain mean (GIS 1–{N_ALL} mean {island} m)"),
        Line2D([], [], color=nc.C_SMOOTH, lw=nc.LW_SMOOTH, label=f"{nc.SMOOTH_DOMAINS}-domain LOWESS")],
        title=(f"Mean shoreline {c1[:4]} ({nc.half_width(e)}) − "
               f"mean shoreline {c0[:4]} ({nc.half_width(s)})"),
        loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6.5, title_fontsize=6.5,
        frameon=False)
    structures(ax, label=True)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    png = out / "coastsat_net_change_1996_2025_gis1-120.png"
    save(fig, png, dpi=300, close=True)
    model = dom.loc[:N_MODEL, "net_change_m"].mean()
    ext = dom.loc[N_MODEL + 1:, "net_change_m"].mean()
    record_caption(png, (
        f"{nc.WHAT} {nc.window_sentence(PERIOD)} GIS 1-{N_MODEL} are the model domains; GIS "
        f"{N_MODEL + 1}-{N_ALL} (shaded, right of the dashed line) extend the same 500 m domain "
        "frame north across Pea Island to Oregon Inlet and are not modelled. Context, not a "
        f"target. Domain means average {model:+.1f} m over GIS 1-{N_MODEL} and {ext:+.1f} m over "
        f"GIS {N_MODEL + 1}-{N_ALL}. Pale line: the mean over each domain's transects (transects "
        "that pass the 10-position minimum in both windows). Dark line: the same at "
        f"{nc.SMOOTH_DOMAINS}-domain LOWESS across all {N_ALL} domains, GIS 1-{nc.RAW_SOUTH} raw. "
        f"At GIS {N_ALL}, the last domain before Oregon Inlet ({int(dom.loc[N_ALL, 'n_transects'])} "
        f"transects), the LOWESS end ({dom.loc[N_ALL, 'net_change_lowess7_m']:+.1f} m) follows "
        f"GIS {N_ALL - 1} ({dom.loc[N_ALL - 1, 'net_change_m']:+.1f} m), not its own mean "
        f"({dom.loc[N_ALL, 'net_change_m']:+.1f} m). "
        "The legend title names the two means differenced. Seaward positive."
        + offaxis_clause(off, nc.POS_HALF, "m")))
    return png


def write_provenance(dom, ext_tr, dropped, n_lookup, out):
    s, e = nc.WINDOWS[PERIOD]
    c0, c1 = nc.CENTRES[PERIOD]
    built = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    ext = dom.loc[N_MODEL + 1:]
    text = f"""# coastsat/net_change/1996_2025/ext -- the full-period net change, GIS 1-{N_ALL}

Written by `scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change_extension.py` on {built}.

## What this is

Context, not a target: the 1996-2025 net change of `../` carried north past the model domains to
Oregon Inlet. GIS 1-{N_MODEL} are the parent folder's transects unchanged. GIS {N_MODEL + 1}-{N_ALL} are
the extension domains of `2-transect-frame/transect_domains/transect_domain_lookup_ext.csv`.

| | |
|---|---|
| start window | {s[0]} to {s[1]}, centred on {c0} |
| end window | {e[0]} to {e[1]}, centred on {c1} |
| extension transects in both windows | {len(ext_tr)} of {n_lookup} (dropped: {dropped}) |
| domain mean net change, GIS 1-{N_MODEL} | {dom.loc[:N_MODEL, 'net_change_m'].mean():+.1f} m |
| domain mean net change, GIS {N_MODEL + 1}-{N_ALL} | {ext['net_change_m'].mean():+.1f} m (range {ext['net_change_m'].min():+.1f} to {ext['net_change_m'].max():+.1f}) |
| domain mean net change, GIS 1-{N_ALL} | {dom['net_change_m'].mean():+.1f} m |

## How the extension was built

- No mean-shoreline window files exist for the extension transects. Their window means are
  computed in memory by `coastsat_mean_shoreline.window_means`, the function that writes
  `1-observations/mean_shoreline/`, with the same {MIN_OBS}-position minimum. Nothing is
  written to `mean_shoreline/`.
- `net_change_lowess7_m` is LOWESS at 7 domains over all {N_ALL} domain means (frac 7/{N_ALL}),
  GIS 1-{nc.RAW_SOUTH} raw. Near GIS 90 it therefore differs slightly from the parent folder's
  curve, which stops at GIS 90.
- Seaward positive. `se_net_change_m` as in the parent folder: sampling noise only.
- The CoastSat transects south of GIS 1 (Cape Point, `usa_NC_0032_0001`-`0020`) have no domain,
  and the Hatteras Inlet to Cape Point shore (`usa_NC_0030`, `usa_NC_0031`) was never
  downloaded, so neither is here.
"""
    (out / "PROVENANCE.md").write_text(text, encoding="utf-8")


def main():
    apply_style()
    out = net_change_dir(*PERIOD) / EXT_DIR
    out.mkdir(parents=True, exist_ok=True)
    ext_tr, dropped, n_lookup = ext_transect_net_change()
    dom = domain_net_change_all(ext_tr)
    ext_tr.round(4).to_csv(out / "transect_net_change_1996_2025_ext.csv")
    dom.round(4).to_csv(out / "domain_net_change_1996_2025_gis1-120.csv")
    png = figure(dom, out)
    write_provenance(dom, ext_tr, dropped, n_lookup, out)
    print(f"extension: {len(ext_tr)} of {n_lookup} transects, GIS 91-120 mean "
          f"{dom.loc[N_MODEL + 1:, 'net_change_m'].mean():+.1f} m; wrote {png.relative_to(_REPO)}")


if __name__ == "__main__":
    main()
