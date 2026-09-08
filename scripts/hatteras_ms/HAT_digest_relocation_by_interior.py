#!/usr/bin/env python3
r"""
HAT_digest_relocation_by_interior.py
==============================================================================
One table across the interiors: what HAT_relocation_comparison.py found for
each arm of the row-insert set. Reads the per-arm reports under

    output/comparisons/relocation_1984_2004/row-insert/<arm>/

and writes output/experiments/row_insert_set/relocation/
    relocation_by_interior.csv    per (arm, historical domain): first year,
                                  error, outcome, plus each arm's recall and
                                  false positives
    relocation_by_interior.txt    the same, readable

The seven reports stay authoritative; this only puts their headline rows side
by side so the interiors can be compared without opening seven folders.

USAGE
    python HAT_digest_relocation_by_interior.py
==============================================================================
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE))
from HAT_run_row_insert_set import ARMS  # noqa: E402  the arm -> version map

CMP = REPO / "output" / "comparisons" / "relocation_1984_2004" / "row-insert"
OUT = REPO / "output" / "experiments" / "row_insert_set" / "relocation"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    head, rows = [], []
    for arm, version in ARMS.items():
        d = CMP / arm
        if not (d / "confusion.csv").is_file():
            print("  (no report for {!r} under {})".format(arm, d))
            continue
        conf = pd.read_csv(d / "confusion.csv").set_index("tolerance_years")
        first = pd.read_csv(d / "first_relocation_year.csv")
        outc = pd.read_csv(d / "road_outcomes.csv")
        free = outc[outc.arm == "free"]
        head.append({
            "arm": arm, "version": version,
            "recall_2yr": float(conf.loc[2, "recall"]),
            "recall_5yr": float(conf.loc[5, "recall"]),
            "hits_5yr": conf.loc[5, "hit_domains"] if isinstance(conf.loc[5, "hit_domains"], str) else "",
            "false_pos": int(conf.loc[5, "false_positives"]),
            "false_pos_domains": conf.loc[5, "false_positive_domains"] if isinstance(conf.loc[5, "false_positive_domains"], str) else "",
            "relocations_free": int(free.relocations.sum()),
            "domains_relocating_free": int((free.relocations > 0).sum()),
            "roads_drowned_free": int(free.drowned.sum()),
        })
        for _, r in first.iterrows():
            rows.append({"arm": arm, "version": version, "gis": int(r.gis),
                         "historical_year": int(r.historical_year),
                         "modelled_first_year": r.modelled_first_year,
                         "error_years": r.error_years,
                         "n_relocations": int(r.n_relocations),
                         "outcome": r.outcome,
                         "min_setback_m": r.min_setback_m,
                         "closest_year": r.closest_year})
    H = pd.DataFrame(head)
    R = pd.DataFrame(rows)
    H.to_csv(OUT / "relocation_by_interior_summary.csv", index=False)
    R.to_csv(OUT / "relocation_by_interior.csv", index=False)

    lines = ["NC-12 RELOCATION, EMERGENT, BY INTERIOR   "
             "(1984-2004 calibBE full management, groin on; arm A of each pair)", ""]
    lines.append("  {:16s} {:4s} {:>7s} {:>7s} {:>6s} {:>7s} {:>6s}  {}".format(
        "arm", "ver", "rec+-2", "rec+-5", "FP/45", "n_reloc", "doms", "hits at +-5 yr"))
    for _, h in H.iterrows():
        lines.append("  {:16s} {:4s} {:7.2f} {:7.2f} {:6d} {:7d} {:6d}  {}".format(
            h.arm, h.version, h.recall_2yr, h.recall_5yr, h.false_pos,
            h.relocations_free, h.domains_relocating_free, h.hits_5yr))
    lines.append("")
    lines.append("  first emergent relocation year per historical domain (hist. year in "
                 "brackets); '-' = never relocated by 2004")
    doms = sorted(R.gis.unique())
    lines.append("  {:16s} ".format("arm") + " ".join(
        "{:>9s}".format("{}[{}]".format(g, int(R[R.gis == g].historical_year.iloc[0]) % 100))
        for g in doms))
    for arm in H.arm:
        sub = R[R.arm == arm].set_index("gis")
        cells = []
        for g in doms:
            y = sub.loc[g, "modelled_first_year"]
            cells.append("{:>9s}".format("-" if pd.isna(y) else "{:d}".format(int(y))))
        lines.append("  {:16s} ".format(arm) + " ".join(cells))
    lines.append("")
    lines.append("  signed error, years (modelled - historical); blank = never relocated")
    for arm in H.arm:
        sub = R[R.arm == arm].set_index("gis")
        cells = []
        for g in doms:
            e = sub.loc[g, "error_years"]
            cells.append("{:>9s}".format("" if pd.isna(e) else "{:+d}".format(int(e))))
        lines.append("  {:16s} ".format(arm) + " ".join(cells))
    text = "\n".join(lines)
    (OUT / "relocation_by_interior.txt").write_text(text, encoding="utf-8")
    print(text)
    print("\nwrote {}".format(OUT))


if __name__ == "__main__":
    main()
