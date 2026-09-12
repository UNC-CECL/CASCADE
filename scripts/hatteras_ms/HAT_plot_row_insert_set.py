#!/usr/bin/env python3
r"""
HAT_plot_row_insert_set.py
==============================================================================
Compares the arms of the 1984 seaward-row insert set IN THE MODEL: island-wide
skill, every road domain's relocation record, and the GIS 80-90 reach in
detail. Reads the runs HAT_run_row_insert_set.py wrote under

    output/raw_runs/row-insert/<arm>/1984_2004/calibBE/<run_name>/

THE ARMS -- since 2026-09-07 only the two unmodified ones can exist
    original         v1                        the original pick set, no rows (reference)
    none             v2                        no rows - the re-pick base

    RETIRED 2026-09-07: measured-floor v4, median v5, platform v6,
    matched-crest v7, matched-nocrest v8. Hannah decided to keep only
    unmodified topography, so the layers, the set's runs (including the two
    controls) and this script's output folder were all deleted. The findings
    are recorded in 1-barrier3d-domains/LINEAGE.md (2026-09-04 entries). To
    use this script again, run HAT_run_row_insert_set.py first.

WHAT IS WRITTEN (output/experiments/row_insert_set/)
    set_rates_island.png        model OLS rate vs the CoastSat target, all
                                90 domains, six arms; skill in the caption
    set_relocations.png         first relocation year and count, every road
                                domain, six arms
    set_setback_GIS80_90.png    road setback through time, relocations marked
    set_dune_GIS80_90.png       max dune height above berm through time
    set_scrape_GIS80_90.png     cumulative sand bulldozed off NC-12
    set_domains.csv             per (domain, arm): N, setback, relocations,
                                dune, scrape, rates, target
    set_summary.txt             the island-wide table and the 80-90 table

TIME INDEXING
    RoadwayManager writes at time_index - 1 and Barrier3D's time_index is 1
    at t = 0, so index i is the END of model year i, calendar 1984 + i.

USAGE
    python HAT_plot_row_insert_set.py [--out DIR] [--arms a,b,c]
==============================================================================
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "scripts"))
from hat_topo_version import dune_topo_root  # noqa: E402
from hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402
from cascade_pipeline.run_layout import resolve  # noqa: E402
from cascade_pipeline.run_registry import preset_dir_for  # noqa: E402
from cascade_pipeline.coastsat_loess import (  # noqa: E402
    CoastSatDataset, LoessConfig, build_coastsat_series)
from cascade_pipeline.hindcast import build_target_table  # noqa: E402

SET = "row-insert"
RAW_RUNS = REPO / "output" / "raw_runs"
RUN_INDEX = RAW_RUNS / "run_index.csv"
# Moved out of the scripts tree 2026-09-12: the rate fits are DATA and
# the model reads them. Resolve through hat_observed_rates.py in new code.
COASTSAT_BASE_DIR = REPO / "data" / "hatteras_init" / "5-scr" / "coastsat_lrr"
OUT = REPO / "output" / "experiments" / "row_insert_set"
PERIOD, PRESET = "1984_2004", "calibBE"
START_YEAR = 1984
BUFFER = HATTERAS_DOMAINS.num_buffer_domains
REACH = tuple(range(80, 91))

# arm, version, label, colour, linestyle
ARMS = [
    ("original",        "v1", "original picks, no insert (reference)", "0.7", ":"),
    ("none",            "v2", "no insert (re-pick base)",      "0.35",    "-"),
    # measured-floor v4 / median v5 / platform v6 / matched-crest v7 /
    # matched-nocrest v8 retired 2026-09-07 with the layers (see docstring).
]


def coastsat_target() -> pd.Series:
    """The per-domain CoastSat LRR target, built as section 8 of the runner
    builds it. The run's rate CSV holds two MODEL estimators and no target."""
    datasets = [CoastSatDataset(
        label="CoastSat LRR (1984-2004)", period_start=1984,
        csv_path=str(COASTSAT_BASE_DIR / "1984_2004" / "transect_lrr_full.csv"))]
    loess = LoessConfig(window_domains=(7, 10), skip_southern_domains=10)
    cs = build_coastsat_series(datasets, active_period_start=START_YEAR,
                               loess_config=loess, domains=HATTERAS_DOMAINS)
    active = next(c for c in cs if c["active"])
    return build_target_table(active, loess, HATTERAS_DOMAINS, 10) \
        .set_index("gis_domain")["target_lrr_m_yr"]


def find_run(arm: str) -> Path:
    d = preset_dir_for(RAW_RUNS, PERIOD, PRESET, arm="{}/{}".format(SET, arm))
    # The arm folder also holds the relocation-ON partner (run name with the
    # `reloc` token) since the relocation comparison; this reads the set run,
    # prescribed relocations off.
    runs = sorted(p for p in d.glob("*") if p.is_dir() and list(p.glob("*.npz"))
                  and "_reloc_" not in p.name)
    if not runs:
        raise SystemExit(
            "\nno run with a model .npz under {}\nRun HAT_run_row_insert_set.py "
            "first (arm {!r}).\n".format(d, arm))
    if len(runs) > 1:
        raise SystemExit("\nmore than one run under {}: {}\n".format(
            d, [r.name for r in runs]))
    return runs[0]


def series(obj, name: str, nt: int) -> np.ndarray:
    v = np.asarray(getattr(obj, name, []), dtype=float).ravel()
    out = np.full(nt, np.nan)
    k = min(nt, v.size)
    out[:k] = v[:k]
    return out


def n_rows(version: str) -> dict:
    audit = dune_topo_root("1984-start") / version / "HAT_seaward_row_insert_audit.csv"
    if not audit.is_file():
        return {}
    with open(audit, newline="") as fh:
        return {int(r["domain"]): int(r["n_rows_inserted"])
                for r in csv.DictReader(fh)}


def index_row(arm: str):
    if not RUN_INDEX.is_file():
        return None
    df = pd.read_csv(RUN_INDEX)
    sel = df[df["arm"].astype(str) == "{}/{}".format(SET, arm)]
    return None if sel.empty else sel.sort_values("timestamp").iloc[-1]


def grid(n: int, title: str):
    ncol = 3
    nrow = int(np.ceil((n + 1) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(15, 3.4 * nrow), squeeze=False)
    fig.suptitle(title, fontsize=14, fontweight="bold", x=0.02, ha="left")
    return fig, axes.ravel()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--arms", default=",".join(a[0] for a in ARMS))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    wanted = [a.strip() for a in args.arms.split(",") if a.strip()]
    arms = [a for a in ARMS if a[0] in wanted]

    loaded = {}
    for arm, version, label, colour, ls in arms:
        run_dir = find_run(arm)
        npz = next(run_dir.glob("*.npz"))
        c = np.load(npz, allow_pickle=True)["cascade"][0]
        rates = pd.read_csv(
            resolve(run_dir, "rate_csv", run_dir.name, must_exist=True)) \
            .set_index("gis_domain")
        loaded[arm] = dict(c=c, run_dir=run_dir, n=n_rows(version),
                           rates=rates, idx=index_row(arm))
        print("  {:16s} {:26s} {}".format(arm, version, run_dir))

    target = coastsat_target()
    first = loaded[arms[0][0]]["c"]
    nt = len(np.asarray(first.barrier3d[BUFFER].x_s_TS))
    years = START_YEAR + np.arange(nt)
    all_gis = range(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1)

    # ------------------------------------------------- per-domain records
    rows = []
    for D in all_gis:
        i = HATTERAS_DOMAINS.gis_to_pad(D)
        for arm, version, label, colour, ls in arms:
            L = loaded[arm]
            c = L["c"]
            b3d = c.barrier3d[i]
            has_road = bool(np.asarray(c.roadway_management_module)[i])
            crest = np.array([np.max(b3d.DuneDomain[t]) * 10.0
                              if t < b3d.DuneDomain.shape[0] else np.nan
                              for t in range(nt)])
            rec = {"domain": D, "arm": arm, "version": version,
                   "rows_inserted": L["n"].get(D, 0), "roadway_module": has_road,
                   "crest_t0_m": round(float(crest[0]), 2),
                   "crest_t1_m": round(float(crest[1]), 2),
                   "crest_max_m": round(float(np.nanmax(crest)), 2),
                   "rate_model_lrr_m_yr": round(float(L["rates"].loc[D, "lrr_m_yr"]), 3),
                   "rate_model_endpoint_m_yr": round(float(L["rates"].loc[D, "change_rate_m_yr"]), 3),
                   "rate_coastsat_target_m_yr": round(float(target.loc[D]), 3),
                   "setback_t0_m": np.nan, "n_relocations": np.nan,
                   "first_relocation_year": np.nan,
                   "road_scrape_yr1_m3": np.nan, "road_scrape_total_m3": np.nan}
            if has_road:
                mgr = c.roadways[i]
                sb = series(mgr, "_road_setback_TS", nt)
                rel = np.nan_to_num(series(mgr, "_road_relocated_TS", nt))
                scrape = np.nan_to_num(series(mgr, "_road_overwash_volume", nt))
                jumps = np.flatnonzero(rel > 0)
                rec.update({
                    "setback_t0_m": float(sb[0]), "n_relocations": int(jumps.size),
                    "first_relocation_year": (START_YEAR + int(jumps[0])) if jumps.size else np.nan,
                    "road_scrape_yr1_m3": float(scrape[1]) if nt > 1 else 0.0,
                    "road_scrape_total_m3": float(scrape.sum())})
            rows.append(rec)
    df = pd.DataFrame(rows)
    df.to_csv(out / "set_domains.csv", index=False)

    # ------------------------------------------------- fig: island rates
    fig, ax = plt.subplots(figsize=(15, 5.5))
    ax.plot(list(all_gis), target.loc[list(all_gis)], "k", marker="s", ms=4,
            lw=0, label="CoastSat LRR target (LOESS window 10)")
    caption = []
    for arm, version, label, colour, ls in arms:
        r = loaded[arm]["rates"].loc[list(all_gis), "lrr_m_yr"]
        ax.plot(list(all_gis), r, ls, color=colour, lw=1.4, label="{} [{}]".format(label, version))
        ir = loaded[arm]["idx"]
        if ir is not None:
            caption.append("{:16s} {:26s} interior RMSE {:.3f}  bias {:+.3f} m/yr  "
                           "roads drowned {}".format(
                               arm, version, float(ir["rmse_interior_m_yr"]),
                               float(ir["mean_bias_interior_m_yr"]),
                               int(ir["roads_drowned"])))
    ax.axvspan(REACH[0] - .5, REACH[-1] + .5, color="0.92", zorder=0)
    ax.set_xlabel("GIS domain")
    ax.set_ylabel("shoreline change rate, OLS (m/yr)")
    ax.set_title("Shoreline change rate by fill, all domains "
                 "(GIS 80-90 shaded; island-wide interior skill below)", loc="left")
    ax.grid(alpha=.3)
    ax.legend(fontsize=8, ncol=2)
    fig.text(0.01, -0.02, "\n".join(caption), family="monospace", fontsize=8.5, va="top")
    fig.tight_layout()
    fig.savefig(out / "set_rates_island.png", dpi=130, bbox_inches="tight")

    # ------------------------------------------------- fig: relocations
    road = df[df.roadway_module].copy()
    doms = sorted(road.domain.unique())
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8), sharex=True)
    for k, (arm, version, label, colour, ls) in enumerate(arms):
        sub = road[road.arm == arm].set_index("domain")
        off = (k - (len(arms) - 1) / 2) * 0.12
        ax1.plot(np.array(doms) + off, sub.loc[doms, "first_relocation_year"],
                 "o", color=colour, ms=4, label=label)
        ax2.bar(np.array(doms) + off, sub.loc[doms, "n_relocations"], width=0.12,
                color=colour, label=label)
    ax1.set_ylabel("first relocation year")
    ax1.grid(alpha=.3)
    ax1.legend(fontsize=8, ncol=3)
    ax1.set_title("NC-12 relocations by fill, every road domain "
                  "(prescribed events off; blank = never relocated)", loc="left")
    ax2.set_ylabel("relocations in 20 yr")
    ax2.set_xlabel("GIS domain")
    ax2.grid(alpha=.3, axis="y")
    fig.tight_layout()
    fig.savefig(out / "set_relocations.png", dpi=130)

    # ------------------------------------------------- figs: the reach
    fig_s, ax_s = grid(len(REACH), "Road setback through time, GIS 80-90 "
                       "(relocations marked)")
    fig_c, ax_c = grid(len(REACH), "Dune height above berm (max of the dune "
                       "array), GIS 80-90")
    fig_v, ax_v = grid(len(REACH), "Cumulative sand bulldozed off NC-12, GIS 80-90")
    for k, D in enumerate(REACH):
        i = HATTERAS_DOMAINS.gis_to_pad(D)
        has_road = bool(np.asarray(first.roadway_management_module)[i])
        for arm, version, label, colour, ls in arms:
            c = loaded[arm]["c"]
            b3d = c.barrier3d[i]
            crest = np.array([np.max(b3d.DuneDomain[t]) * 10.0 for t in range(nt)])
            ax_c[k].plot(years, crest, ls, color=colour, lw=1.6, label=label)
            if has_road:
                mgr = c.roadways[i]
                sb = series(mgr, "_road_setback_TS", nt)
                rel = np.nan_to_num(series(mgr, "_road_relocated_TS", nt))
                cum = np.cumsum(np.nan_to_num(series(mgr, "_road_overwash_volume", nt)))
                ax_s[k].plot(years, sb, ls, color=colour, lw=1.6, label=label)
                for j in np.flatnonzero(rel > 0):
                    ax_s[k].plot(years[j], sb[j], "o", color=colour, ms=5,
                                 mec="k", mew=.6, zorder=5)
                ax_v[k].plot(years, cum / 1e3, ls, color=colour, lw=1.6, label=label)
        if not has_road:
            for ax in (ax_s[k], ax_v[k]):
                ax.text(0.5, 0.5, "beach-dune managed\n(no roadway module)",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=10, color="0.4")
        n_here = loaded[arms[-1][0]]["n"].get(D, 0)
        for ax, ylab in ((ax_s[k], "setback (m)"), (ax_c[k], "height (m)"),
                         (ax_v[k], "cumulative (1000 m3)")):
            ax.set_title("GIS {}   N = {} rows".format(D, n_here), fontsize=11, loc="left")
            ax.set_ylabel(ylab)
            ax.grid(alpha=.3)
        ax_s[k].axhline(0, color="k", lw=.6)
    for axes_, fig in ((ax_s, fig_s), (ax_c, fig_c), (ax_v, fig_v)):
        h, l = axes_[0].get_legend_handles_labels()
        axes_[len(REACH)].axis("off")
        axes_[len(REACH)].legend(h, l, loc="center", fontsize=9, frameon=False)
        for ax in axes_[len(REACH) + 1:]:
            ax.axis("off")
        fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig_s.savefig(out / "set_setback_GIS80_90.png", dpi=130)
    fig_c.savefig(out / "set_dune_GIS80_90.png", dpi=130)
    fig_v.savefig(out / "set_scrape_GIS80_90.png", dpi=130)

    # ------------------------------------------------- text summary
    lines = ["ROW-INSERT FILL SET  --  1984-2004 calibBE full management, groin on, "
             "prescribed relocations off", ""]
    lines.append("  {:16s} {:26s} {:32s} {}".format("arm", "version", "fill", "island-wide"))
    for arm, version, label, colour, ls in arms:
        ir = loaded[arm]["idx"]
        rd = road[road.arm == arm]
        skill = ("interior RMSE {:.3f}, bias {:+.3f} m/yr".format(
            float(ir["rmse_interior_m_yr"]), float(ir["mean_bias_interior_m_yr"]))
            if ir is not None else "no run_index row")
        lines.append("  {:16s} {:26s} {:32s} {}; relocations {} in {} domains; "
                     "road scrape {:.0f} m3".format(
                         arm, version, label, skill, int(rd.n_relocations.sum()),
                         int((rd.n_relocations > 0).sum()),
                         float(rd.road_scrape_total_m3.sum())))
    lines.append("")
    hdr = ("  dom  N | arm              setb0  1st reloc  n | dune t0   t1   max | "
           "scrape yr1     total  | rate model  target")
    lines.append(hdr)
    lines.append("  " + "-" * (len(hdr) - 2))
    for D in REACH:
        for _, r in df[df.domain == D].iterrows():
            if r.roadway_module:
                road_cols = "{:5.0f} {:>9} {:2d} | ".format(
                    r.setback_t0_m,
                    "" if pd.isna(r.first_relocation_year) else int(r.first_relocation_year),
                    int(r.n_relocations)) + "{:5.2f} {:5.2f} {:5.2f} | {:10.0f} {:9.0f}  | ".format(
                    r.crest_t0_m, r.crest_t1_m, r.crest_max_m,
                    r.road_scrape_yr1_m3, r.road_scrape_total_m3)
            else:
                road_cols = "{:>5} {:>9} {:>2} | {:5.2f} {:5.2f} {:5.2f} | {:>10} {:>9}  | ".format(
                    "", "", "", r.crest_t0_m, r.crest_t1_m, r.crest_max_m, "", "")
            lines.append("  {:3d} {:2d} | {:16s} ".format(D, int(r.rows_inserted), r.arm)
                         + road_cols + "{:+6.2f}  {:+6.2f}".format(
                             r.rate_model_lrr_m_yr, r.rate_coastsat_target_m_yr))
        lines.append("")
    (out / "set_summary.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print("\nwrote figures, CSV and summary to {}".format(out))


if __name__ == "__main__":
    main()
