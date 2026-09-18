#!/usr/bin/env python3
r"""
HAT_version_pair_report.py
==============================================================================
v2 against v3 in ONE report: the relocation comparison's report.txt, but with
the two dune-topo versions side by side in every section instead of one
version per file. The companion of HAT_version_pair_gif.py, which does the
same for the animations (Hannah, 2026-09-10: "a txt report directly comparing
the two versions, similar to v3/calibBE_groin/report.txt").

WHAT IS READ. The two per-version sets that HAT_relocation_comparison.py
wrote, v2/calibBE_groin/tables/ and v3/calibBE_groin/tables/, the run
metadata of the four runs behind them, and the v3 footprint audit
(dune-topo/v3/HAT_footprint_audit.csv: which domains got rows and how many).
Nothing is re-run and nothing is re-scored: every number here is one of
theirs, or a difference of two of theirs. A section whose numbers disagree
with the per-version reports means one of the three is stale.

WHAT IS WRITTEN, all under output/comparisons/relocation/versions/v2_vs_v3/
    report.txt          the console output of this run, with the provenance
                        of the four runs and the two table sets above it
    tables/*.csv        the side-by-side tables the report prints, in full
                        (the report prints the historical domains and the
                        domains that differ; the CSVs hold every road domain)

USAGE
    python HAT_version_pair_report.py
    python HAT_version_pair_report.py --preset calibBE --set calibBE_groin
==============================================================================
"""
from __future__ import annotations

import argparse
import contextlib
import datetime
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "hatteras_ms"))
# the relocation comparison moved into hatteras_ms/experiments/ on 2026-09-13
# (cfd0b475); this import broke silently until 2026-09-17
sys.path.insert(0, str(REPO / "scripts" / "hatteras_ms" / "experiments"))
import HAT_relocation_comparison as RC  # noqa: E402   _Tee, _arm_provenance, TOLERANCE_YEARS, windows

# raw_runs by purpose since 2026-09-16: the version pair is under versions/
RAW = REPO / "output" / "raw_runs" / "versions" / "version-pair"
COMP = REPO / "output" / "comparisons" / "relocation" / "1984_2004"
OUT = REPO / "output" / "comparisons" / "relocation" / "versions" / "v2_vs_v3"
from site_layer import hat_topo_version as _b3d  # noqa: E402
AUDIT = _b3d.dune_topo_root("1984-start") / "v3" / "HAT_footprint_audit.csv"
VERSIONS = ("v2", "v3")
LABEL = {"v2": "the extraction (today's setbacks)", "v3": "the 1984 reconstruction (1984 setbacks)"}
ARMS = {"free": "HAT_1984_2004_{preset}_road_bdm_groin",
        "prescribed": "HAT_1984_2004_{preset}_road_reloc_bdm_groin"}
TABLES = ("confusion", "first_relocation_year", "indexing_check", "near_miss_margin",
          "road_outcomes", "setback_by_year", "setback_summary")
CELL_M = 10.0
RULE = "-" * 74
BAR = "=" * 74

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_rows", 200)


# =============================================================================
# READING
# =============================================================================

def read_set(version: str, set_name: str) -> dict:
    """The seven CSVs of one per-version comparison set, plus the 'generated'
    stamp of the report beside them."""
    d = COMP / version / set_name
    if not d.is_dir():
        raise SystemExit(f"no comparison set at {d} -- run HAT_relocation_comparison.py on {version} first")
    tabs = {}
    for t in TABLES:
        p = d / "tables" / f"{t}.csv"
        if not p.is_file():
            raise SystemExit(f"{p} is missing; the {version} set is incomplete")
        tabs[t] = pd.read_csv(p)
    stamp = "?"
    rep = d / "report.txt"
    if rep.is_file():
        m = re.search(r"^generated\s+(\S+ \S+)", rep.read_text(encoding="utf-8"), re.M)
        stamp = m.group(1) if m else "?"
    return {"dir": d, "tables": tabs, "report_stamp": stamp}


def read_audit() -> pd.DataFrame:
    """v3's footprint audit: one row per GIS domain, n_cells signed (+ added,
    - removed, 0 unchanged)."""
    if not AUDIT.is_file():
        raise SystemExit(f"{AUDIT} is missing")
    a = pd.read_csv(AUDIT).rename(columns={"domain": "gis"})
    a["rows_v3_minus_v2"] = a["n_cells"].fillna(0).astype(int)
    return a[["gis", "rows_v3_minus_v2", "operation"]]


def _fmt_years(vals) -> str:
    return " ".join(str(int(v)) for v in vals)


def _dom_list(s) -> list[int]:
    if isinstance(s, str) and s.strip():
        return [int(x) for x in s.split()]
    return []


def _tab(df: pd.DataFrame, na_rep: str = "NaN") -> str:
    """A ruled table: pandas' own per-column formatting, columns separated by
    ' | ' and a rule under the header (Hannah, 2026-09-10: lines between the
    columns so the tables read more easily). Right-aligned like to_string."""
    cols = []
    for c in df.columns:
        vals = df[c].to_string(index=False, na_rep=na_rep).splitlines() if len(df) else []
        vals = [v.strip().replace("<NA>", na_rep) for v in vals]
        w = max([len(str(c))] + [len(v) for v in vals])
        cols.append((str(c).rjust(w), [v.rjust(w) for v in vals], w))
    head = " | ".join(c[0] for c in cols)
    rule = "-+-".join("-" * c[2] for c in cols)
    body = [" | ".join(c[1][i] for c in cols) for i in range(len(df))]
    return chr(10).join(["  " + head, "  " + rule] + ["  " + b for b in body])


# =============================================================================
# THE COMPARISON. Everything printed becomes report.txt.
# =============================================================================

def compare(sets: dict, audit: pd.DataFrame, preset: str, out_dir: Path) -> None:
    tdir = out_dir / "tables"
    tdir.mkdir(parents=True, exist_ok=True)
    T = {v: sets[v]["tables"] for v in VERSIONS}
    first = {v: T[v]["first_relocation_year"] for v in VERSIONS}
    hist = first["v2"]["gis"].astype(int).tolist()
    hist_year = dict(zip(first["v2"]["gis"].astype(int), first["v2"]["historical_year"].astype(int)))
    summ = {v: T[v]["setback_summary"] for v in VERSIONS}
    road = sorted(set(summ["v2"]["gis"]).union(summ["v3"]["gis"]))

    print(BAR)
    print(f"NC-12 RELOCATION: v2 against v3, each as emergent vs prescribed")
    print(BAR)
    print(f"  preset              {preset}")
    for v in VERSIONS:
        print(f"  {v}                  {LABEL[v]}")
    print(f"  managed road domains  v2: {len(summ['v2'])}   v3: {len(summ['v3'])}")
    print(f"  historical domains    {hist}")
    for yr in sorted(set(hist_year.values())):
        print(f"    {yr}  GIS {[g for g in hist if hist_year[g] == yr]}")
    print()
    print("  Sign convention throughout: a 'change' or 'v3 - v2' column is v3")
    print("  minus v2. Setbacks are metres landward of dune row 0, so a positive")
    print("  change means v3 starts the road FURTHER from the dune line.")

    # ---------------------------------------------------------------- A
    print()
    print(RULE)
    print("A. WHAT v3 CHANGED AT THE ROAD, before the model ran")
    print(RULE)
    print("  v3 = v2 + the symmetric 1984 footprint (rows added behind NC-12 as")
    print("  placed under the 1984 setback, rows removed in front of today's road)")
    print("  and the 1984 setback CSV in place of today's. The dune array is the")
    print("  same. So the two runs differ in (i) how much island sits behind the")
    print("  road and (ii) where the road starts.")
    print()
    n_add = int((audit["rows_v3_minus_v2"] > 0).sum())
    n_rem = int((audit["rows_v3_minus_v2"] < 0).sum())
    print(f"  island (90 domains)   rows added at {n_add} domains (+{int(audit.loc[audit['rows_v3_minus_v2'] > 0, 'rows_v3_minus_v2'].sum())}), "
          f"removed at {n_rem} ({int(audit.loc[audit['rows_v3_minus_v2'] < 0, 'rows_v3_minus_v2'].sum())}), "
          f"unchanged at {90 - n_add - n_rem}")
    start = pd.DataFrame({"gis": road})
    for v in VERSIONS:
        s = summ[v].set_index("gis")
        start[f"setback_1984_{v}_m"] = start["gis"].map(s["free_start_m"])
    start["change_m"] = start["setback_1984_v3_m"] - start["setback_1984_v2_m"]
    start = start.merge(audit, on="gis", how="left")
    start["historical"] = start["gis"].isin(hist)
    start["historical_year"] = start["gis"].map(hist_year).astype("Int64")
    ra = start[start["rows_v3_minus_v2"] > 0]
    rr = start[start["rows_v3_minus_v2"] < 0]
    print(f"  road domains ({len(road)})     rows added at {len(ra)} (GIS {ra['gis'].tolist()}),")
    print(f"                        removed at {len(rr)} (GIS {rr['gis'].tolist()})")
    moved = start[start["change_m"].abs() >= 5]
    print(f"  setback at 1984       differs by >= 5 m at {len(moved)} of {len(road)} road domains; "
          f"median change {start['change_m'].median():+.1f} m, range {start['change_m'].min():+.1f} to {start['change_m'].max():+.1f} m")
    print(f"                        v2 is quantised to whole cells; v3 carries the measured value unrounded, no floor "
          f"(min {start['setback_1984_v3_m'].min():.1f} m at GIS {int(start.loc[start['setback_1984_v3_m'].idxmin(), 'gis'])})")
    print()
    print("  HISTORICAL domains:")
    cols = ["gis", "historical_year", "rows_v3_minus_v2", "operation", "setback_1984_v2_m", "setback_1984_v3_m", "change_m"]
    print(_tab(start.loc[start["historical"], cols]))
    start_out = start[["gis", "historical", "historical_year"] + cols[2:]]
    start_out.to_csv(tdir / "start_conditions.csv", index=False)
    print()
    print("  saved all road domains -> start_conditions.csv")

    # ---------------------------------------------------------------- 0
    print()
    print(RULE)
    print("0. VALIDITY CHECKS, per version")
    print(RULE)
    for v in VERSIONS:
        ic = T[v]["indexing_check"]
        print(f"  {v}  prescribed jumps land on the event year: {int(ic['matches'].sum())}/{len(ic)} domains")
    j = T["v2"]["indexing_check"].merge(T["v3"]["indexing_check"], on="gis", suffixes=("_v2", "_v3"))
    dj = (j["largest_jump_m_v3"] - j["largest_jump_m_v2"]).abs().max()
    print(f"  prescribed displacements identical between versions: {'YES' if dj < 0.01 else 'NO'} "
          f"(largest difference {dj:.3f} m) -- the same measured 1989/1999 moves are applied to both")

    # ---------------------------------------------------------------- 1
    print()
    print(RULE)
    print("1. FIRST MODELLED RELOCATION, free-running arm")
    print(RULE)
    f = first["v2"][["gis", "historical_year"]].copy()
    for v in VERSIONS:
        s = first[v].set_index("gis")
        f[f"first_{v}"] = f["gis"].map(s["modelled_first_year"])
        f[f"err_{v}"] = f["gis"].map(s["error_years"])
        f[f"n_reloc_{v}"] = f["gis"].map(s["n_relocations"]).astype(int)
        f[f"outcome_{v}"] = f["gis"].map(s["outcome"])
    f["abs_err_change"] = f["err_v3"].abs() - f["err_v2"].abs()

    def verdict(r):
        a, b = np.isfinite(r["err_v2"]), np.isfinite(r["err_v3"])
        if a and b:
            d = abs(r["err_v3"]) - abs(r["err_v2"])
            return "same" if d == 0 else ("v3 closer" if d < 0 else "v2 closer")
        if a and not b:
            return "lost in v3"
        if b and not a:
            return "gained in v3"
        return "neither"
    f["verdict"] = f.apply(verdict, axis=1)
    print(_tab(f))
    print()
    for v in VERSIONS:
        d = first[v]
        rel = d["modelled_first_year"].notna()
        print(f"  {v}  relocated at all {int(rel.sum())}/{len(d)}   signed error median {d['error_years'].median():+.1f} yr, "
              f"mean {d['error_years'].mean():+.1f} yr   mean absolute error {d['error_years'].abs().mean():.1f} yr   "
              f"relocations in these domains {int(d['n_relocations'].sum())}")
    vc = f["verdict"].value_counts()
    print("  per domain: " + ", ".join(f"{k} {int(n)}" for k, n in vc.items()))
    print("  a domain that never relocates has no error and is not in the median;")
    print("  compare the 'relocated at all' counts alongside it.")
    f.to_csv(tdir / "first_relocation_year.csv", index=False)
    print("  saved -> first_relocation_year.csv")

    # ---------------------------------------------------------------- 1b
    print()
    print(RULE)
    print("1b. NEAR MISSES: how much more dune migration would have fired it")
    print(RULE)
    print("  'needed' is the extra landward migration that would have tripped")
    print("  `setback < 0`: closest approach + one 10 m cell. 'relocated' where")
    print("  the domain did fire in that version.")
    nm = {v: T[v]["near_miss_margin"].set_index("gis") for v in VERSIONS}
    print()
    print("  HISTORICAL domains that never relocated in at least one version:")
    rows = []
    for g in hist:
        r = {"gis": g, "hist": hist_year[g]}
        for v in VERSIONS:
            if g in nm[v].index:
                r[f"closest_{v}_m"] = nm[v].loc[g, "min_setback_m"]
                r[f"needed_{v}_m"] = nm[v].loc[g, "migration_needed_m"]
                r[f"at_year_{v}"] = int(nm[v].loc[g, "closest_year"])
            else:
                r[f"closest_{v}_m"] = np.nan
                r[f"needed_{v}_m"] = np.nan
                r[f"at_year_{v}"] = np.nan
        rows.append(r)
    nh = pd.DataFrame(rows)
    nh = nh[nh[["needed_v2_m", "needed_v3_m"]].notna().any(axis=1)]
    nh["needed_change_m"] = nh["needed_v3_m"] - nh["needed_v2_m"]
    for v in VERSIONS:
        nh[f"at_year_{v}"] = nh[f"at_year_{v}"].astype("Int64")
    print(_tab(nh, na_rep="relocated"))
    print()
    print("  CONTROL domains (history never relocated these):")
    ctrl = {}
    for v in VERSIONS:
        c = nm[v][nm[v]["kind"] == "control"]
        ctrl[v] = c
        one = c[c["cells_remaining"] <= 1].index.astype(int).tolist()
        print(f"  {v}  closest approach median {c['min_setback_m'].median():.0f} m, min {c['min_setback_m'].min():.1f} m "
              f"(the per-version reports' 'control margin'); within ONE cell of firing: {len(one)}  GIS {one}")
    cc = ctrl["v2"][["migration_needed_m", "closest_year"]].join(
        ctrl["v3"][["migration_needed_m", "closest_year"]], lsuffix="_v2", rsuffix="_v3", how="outer")
    cc["needed_change_m"] = cc["migration_needed_m_v3"] - cc["migration_needed_m_v2"]
    cc = cc.reset_index().sort_values("migration_needed_m_v3")
    print()
    print("  closest 10 controls in v3, with the v2 margin beside:")
    print(_tab(cc.head(10).rename(columns={"migration_needed_m_v2": "needed_v2_m", "migration_needed_m_v3": "needed_v3_m"})))
    tight = cc[cc["needed_change_m"] <= -20].sort_values("gis")
    print(f"  controls whose margin SHRANK by 20 m or more in v3: {len(tight)}  GIS {tight['gis'].astype(int).tolist()}")
    wide = cc[cc["needed_change_m"] >= 20].sort_values("gis")
    print(f"  controls whose margin GREW by 20 m or more in v3:   {len(wide)}  GIS {wide['gis'].astype(int).tolist()}")
    print("  -- v3's 1984 setbacks are the measured ones, no floor, so a control")
    print("     that today sits a few cells back can start almost on the dune line.")
    pd.concat([nh.assign(kind="historical"), cc.assign(kind="control")]).to_csv(tdir / "near_miss_margin.csv", index=False)
    print("  saved -> near_miss_margin.csv")

    # ---------------------------------------------------------------- 2
    print()
    print(RULE)
    print("2. HIT / MISS, with false positives")
    print(RULE)
    conf_rows = []
    for tol in RC.TOLERANCE_YEARS:
        c = {v: T[v]["confusion"].set_index("tolerance_years").loc[tol] for v in VERSIONS}
        hits = {v: set(_dom_list(c[v]["hit_domains"])) for v in VERSIONS}
        fps = {v: set(_dom_list(c[v]["false_positive_domains"])) for v in VERSIONS}
        print()
        print(f"  +/-{tol} yr")
        for v in VERSIONS:
            print(f"    {v}  hits {int(c[v]['hits'])}/{int(c[v]['historical_domains'])}  (recall {c[v]['recall']:.2f})   "
                  f"hit GIS {sorted(hits[v])}   false pos {int(c[v]['false_positives'])}/{int(c[v]['control_domains'])} "
                  f"(rate {c[v]['false_positive_rate']:.2f})")
        print(f"    hit in both {sorted(hits['v2'] & hits['v3'])}   only v2 {sorted(hits['v2'] - hits['v3'])}   "
              f"only v3 {sorted(hits['v3'] - hits['v2'])}")
        conf_rows.append({"tolerance_years": tol,
                          **{f"{k}_{v}": c[v][k] for v in VERSIONS for k in ("hits", "recall", "false_positives", "false_positive_rate")},
                          "hit_both": _fmt_years(sorted(hits["v2"] & hits["v3"])),
                          "hit_only_v2": _fmt_years(sorted(hits["v2"] - hits["v3"])),
                          "hit_only_v3": _fmt_years(sorted(hits["v3"] - hits["v2"])),
                          "false_positive_domains_v2": _fmt_years(sorted(fps["v2"])),
                          "false_positive_domains_v3": _fmt_years(sorted(fps["v3"]))})
    print()
    cr = pd.DataFrame(conf_rows)
    if (cr["recall_v2"] == cr["recall_v3"]).all():
        print("  the recall is the same in both versions at every tolerance; what")
        print("  changes is WHICH domains are hit:")
    else:
        print("  recall differs between the versions; the domains behind it:")
    swapped = sorted(set().union(*[set(_dom_list(r["hit_only_v2"])) | set(_dom_list(r["hit_only_v3"])) for r in conf_rows]))
    if swapped:
        sw = f[f["gis"].isin(swapped)][["gis", "historical_year", "first_v2", "err_v2", "first_v3", "err_v3"]].copy()
        st = start.set_index("gis")
        sw["setback_1984_v2_m"] = sw["gis"].map(st["setback_1984_v2_m"])
        sw["setback_1984_v3_m"] = sw["gis"].map(st["setback_1984_v3_m"])
        print(_tab(sw))
    cr.to_csv(tdir / "confusion.csv", index=False)
    print("  saved -> confusion.csv")

    # ---------------------------------------------------------------- 3
    print()
    print(RULE)
    print("3. SETBACK TRAJECTORIES")
    print(RULE)
    s = summ["v2"][["gis", "historical", "event_year", "measured_2004_m"]].copy()
    for v in VERSIONS:
        m = summ[v].set_index("gis")
        for k in ("free_start_m", "free_last_m", "prescribed_last_m", "rmse_m"):
            s[f"{k[:-2]}_{v}_m"] = s["gis"].map(m[k])
    s["free_last_change_m"] = s["free_last_v3_m"] - s["free_last_v2_m"]
    s["prescribed_last_change_m"] = s["prescribed_last_v3_m"] - s["prescribed_last_v2_m"]
    for v in VERSIONS:
        s[f"free_net_{v}_m"] = s[f"free_last_{v}_m"] - s[f"free_start_{v}_m"]
    s["free_net_change_m"] = s["free_net_v3_m"] - s["free_net_v2_m"]
    print("  HISTORICAL domains, 2004 setback under each arm and version:")
    cols = ["gis", "event_year", "free_start_v2_m", "free_start_v3_m", "free_last_v2_m", "free_last_v3_m",
            "prescribed_last_v2_m", "prescribed_last_v3_m", "rmse_v2_m", "rmse_v3_m", "measured_2004_m"]
    print(_tab(s.loc[s["historical"], cols].round(1)))
    print()
    print("  rmse is the free-vs-prescribed gap within a version (the per-version")
    print("  reports' column), not a gap between versions.")
    for v in VERSIONS:
        h = s[s["historical"]]
        print(f"  {v}  |prescribed 2004 - measured 2004| over the 10 historical domains: "
              f"median {(h[f'prescribed_last_{v}_m'] - h['measured_2004_m']).abs().median():.1f} m, "
              f"max {(h[f'prescribed_last_{v}_m'] - h['measured_2004_m']).abs().max():.1f} m")
    print()
    print("  ALL road domains, free arm: net 20-yr change of the setback (2004 - 1984)")
    print("  per version, and where the two versions differ by a cell or more:")
    dif = s[s["free_net_change_m"].abs() >= CELL_M]
    cols2 = ["gis", "historical", "free_start_v2_m", "free_start_v3_m", "free_net_v2_m", "free_net_v3_m", "free_net_change_m"]
    if len(dif):
        print(_tab(dif[cols2].round(1)))
    else:
        print("    none")
    print(f"  {len(dif)} of {len(s)} road domains differ by >= {CELL_M:.0f} m in net change; "
          f"{int((s['free_net_change_m'].abs() < 0.05).sum())} are identical to 0.05 m")
    print("  -- the module's relocation resets the setback to the run's STARTING")
    print("     value (see HAT_relocation_comparison.py, 'what is deliberately not")
    print("     done'), so a domain's net change is the migration since its last")
    print("     relocation, and a larger v3 start does not by itself move 2004.")

    # per-year gap between versions, free arm
    by = {}
    for v in VERSIONS:
        by[v] = T[v]["setback_by_year"][["gis", "year", "setback_free_m", "setback_prescribed_m"]].rename(
            columns={"setback_free_m": f"free_{v}_m", "setback_prescribed_m": f"prescribed_{v}_m"})
    yy = by["v2"].merge(by["v3"], on=["gis", "year"])
    yy["free_v3_minus_v2_m"] = yy["free_v3_m"] - yy["free_v2_m"]
    yy["prescribed_v3_minus_v2_m"] = yy["prescribed_v3_m"] - yy["prescribed_v2_m"]
    yy = yy[["gis", "year", "free_v2_m", "free_v3_m", "free_v3_minus_v2_m",
             "prescribed_v2_m", "prescribed_v3_m", "prescribed_v3_minus_v2_m"]]
    g = yy.groupby("year")["free_v3_minus_v2_m"]
    print()
    print("  free arm, v3 - v2 setback across the road domains, by year:")
    yr_tab = pd.DataFrame({"year": g.median().index, "median_m": g.median().values,
                           "p10_m": g.quantile(0.1).values, "p90_m": g.quantile(0.9).values})
    print(_tab(yr_tab.round(1)))
    s.to_csv(tdir / "setback_summary.csv", index=False)
    yy.to_csv(tdir / "setback_by_year.csv", index=False)
    print("  saved -> setback_summary.csv, setback_by_year.csv (per year, both arms)")

    # ---------------------------------------------------------------- 4
    print()
    print(RULE)
    print("4. ROAD OUTCOMES, per arm and version")
    print(RULE)
    ro = {v: T[v]["road_outcomes"] for v in VERSIONS}
    for arm in ("free", "prescribed"):
        for v in VERSIONS:
            r = ro[v][ro[v]["arm"] == arm]
            print(f"  {arm:<11} {v}  drowned {int(r['drowned'].sum()):>3}   blocked {int(r['relocation_blocked'].sum()):>3}   "
                  f"relocations {int(r['relocations'].sum()):>3}   dunes rebuilt {int(r['dunes_rebuilt'].sum()):>4}   "
                  f"overwash removed {r['overwash_removed_m3'].sum() / 1e3:>8.1f} x10^3 m3   of {len(r)} managed domains")
    merged = []
    for arm in ("free", "prescribed"):
        a = ro["v2"][ro["v2"]["arm"] == arm].set_index("gis")
        b = ro["v3"][ro["v3"]["arm"] == arm].set_index("gis")
        m = a[["drowned", "relocation_blocked", "relocations", "dunes_rebuilt", "overwash_removed_m3"]].join(
            b[["drowned", "relocation_blocked", "relocations", "dunes_rebuilt", "overwash_removed_m3"]],
            lsuffix="_v2", rsuffix="_v3", how="outer")
        m["relocations_change"] = m["relocations_v3"] - m["relocations_v2"]
        m["overwash_removed_change_m3"] = m["overwash_removed_m3_v3"] - m["overwash_removed_m3_v2"]
        m = m.reset_index()
        m.insert(0, "arm", arm)
        merged.append(m)
    mm = pd.concat(merged)
    print()
    for arm in ("free", "prescribed"):
        d = mm[(mm["arm"] == arm) & (mm["relocations_change"] != 0)]
        print(f"  {arm}: domains whose relocation count changed: {len(d)}")
        if len(d):
            print(_tab(d[["gis", "relocations_v2", "relocations_v3", "relocations_change"]]))
    fate = mm[(mm["drowned_v2"] != mm["drowned_v3"]) | (mm["relocation_blocked_v2"] != mm["relocation_blocked_v3"])]
    print(f"  domains whose fate (drowned / blocked) differs between versions: {len(fate)}"
          + (f"  GIS {fate['gis'].tolist()}" if len(fate) else ""))
    mm.to_csv(tdir / "road_outcomes.csv", index=False)
    print("  saved -> road_outcomes.csv")

    # ---------------------------------------------------------------- 5
    print()
    print(RULE)
    print("5. ANIMATIONS (v2 left, v3 right; written by HAT_version_pair_gif.py)")
    print(RULE)
    gifs = sorted(out_dir.glob("*/*/*.gif"))
    if not gifs:
        print("  none found under this folder -- run HAT_version_pair_gif.py")
    for p in gifs:
        print(f"  {p.relative_to(out_dir)}")
    print()
    print("  the island-wide geometry (interior width, elevation, overwash) of the")
    print("  same four runs is in 2-domain-reconstruction-1984/6-result/HAT_compare_versions.txt")


# =============================================================================
# HEADER AND MAIN
# =============================================================================

def header(preset: str, runs: dict, sets: dict) -> str:
    L = [BAR,
         f"generated   {datetime.datetime.now():%Y-%m-%d %H:%M:%S} by {Path(__file__).name}",
         f"preset      {preset}"]
    for v in VERSIONS:
        L.append(f"{v}          {LABEL[v]}")
        for arm in ("free", "prescribed"):
            L.append(f"  {arm:<11}{RC._arm_provenance(runs[v][arm])}")
    L += ["",
          "Read from the two per-version comparison sets written by",
          "HAT_relocation_comparison.py (the report.txt and tables/ in each):"]
    for v in VERSIONS:
        L.append(f"  {v}  {sets[v]['dir'].relative_to(REPO)}   (its report generated {sets[v]['report_stamp']})")
    L += [f"and the v3 footprint audit, {AUDIT.relative_to(REPO)}.",
          "Nothing is re-run and nothing is re-scored: every number below is one",
          "of theirs, or a difference of two of theirs.",
          "",
          "This file is the console output of the run that wrote the CSVs in",
          "tables/ beside it. If the run identities above do not match the runs on",
          "disk, or a per-version report is newer than this one, this report is",
          "stale -- re-run HAT_relocation_comparison.py on both versions, then this.",
          BAR, ""]
    return "\n".join(L) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description="v2 against v3: the relocation comparison side by side")
    ap.add_argument("--preset", default="calibBE")
    ap.add_argument("--set", default=None, help="per-version set folder name (default <preset>_groin)")
    a = ap.parse_args()
    set_name = a.set or f"{a.preset}_groin"
    runs = {v: {arm: RAW / v / "1984_2004" / a.preset / ARMS[arm].format(preset=a.preset) for arm in ARMS}
            for v in VERSIONS}
    for v in VERSIONS:
        for arm, d in runs[v].items():
            if not d.is_dir():
                raise SystemExit(f"no {arm} run for {v} at {d}")
    sets = {v: read_set(v, set_name) for v in VERSIONS}
    # the per-version sets must have been built from THESE runs
    for v in VERSIONS:
        for arm in ARMS:
            if RC._topo_version(runs[v][arm]) != v:
                raise SystemExit(f"{runs[v][arm]} is not a {v} run")
    audit = read_audit()

    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "report.txt"
    if report_path.exists():                       # a stale report is worse than none
        report_path.unlink()
    tee = RC._Tee()
    failure = None
    try:
        with contextlib.redirect_stdout(tee):
            compare(sets, audit, a.preset, OUT)
    except BaseException as exc:                   # noqa: BLE001
        failure = exc
        raise
    finally:
        body = tee.getvalue()
        if failure is not None:
            body += ("\n" + "!" * 74 + "\n"
                     f"RUN FAILED before completing: {type(failure).__name__}: {failure}\n"
                     "This report is PARTIAL.\n" + "!" * 74 + "\n")
        report_path.write_text(header(a.preset, runs, sets) + body
                               + "\n" + BAR + f"\nartifacts -> {OUT}\n" + BAR + "\n", encoding="utf-8")
        print(f"report -> {report_path}")


if __name__ == "__main__":
    main()
