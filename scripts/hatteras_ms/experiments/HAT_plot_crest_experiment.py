#!/usr/bin/env python3
r"""
HAT_plot_crest_experiment.py
==============================================================================
Compares the three arms written by HAT_run_crest_experiment.py.

    pea1989base    v1 as shipped -- GIS 85 setback floored to 0
    pea1989keep    +N rows, the 1996 dune crest left standing in the interior
    pea1989lower   +N rows, that crest shaved to the backdune platform

    FROZEN 2026-09-07. The keep/lower run outputs were deleted with every run
    on modified topography (their topography had gone on 2026-09-03), so this
    script can no longer be re-run; output/experiments/pea1989_crest/ is the
    record. Only pea1989base (v1) still exists under output/raw_runs/.

WHAT THE FIGURES ARE FOR
    Figure 1 is the test. It draws the road setback through time at GIS 84, 85
    and 86 with the prescribed 1989 relocation marked. The claim being checked
    is that the baseline's setback hits zero YEARS BEFORE 1989 and triggers an
    emergent relocation, so the prescribed displacement is then added to a
    synthetic base rather than to the evolved 1984 position. If the inserts
    work, their traces reach 1989 still positive.

    Figure 2 is the control. The insert touches three domains out of ninety, so
    the island-wide shoreline change rate should be unmoved everywhere else. A
    difference out at GIS 40 would mean the insert leaked through the alongshore
    coupling and the experiment is not isolating what it claims to.

READING THE SETBACK TRACE
    roadway_manager keeps `_road_setback_TS` in metres and rewrites it every
    year as `setback += dune_migrated`. A relocation shows as an upward jump: an
    EMERGENT one to the fixed `relocation_setback_m`, a PRESCRIBED one by that
    domain's own measured displacement. Which is which is the point, so both are
    marked rather than left to the eye.

USAGE
    python HAT_plot_crest_experiment.py
==============================================================================
"""

from __future__ import annotations

import argparse
import csv
import glob
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
# Anchored by SEARCHING UPWARD for the project root rather than by
# counting parent directories (2026-09-13). A counted depth is correct
# only while the file stays where it was written, and these moved into
# subfolders of hatteras_ms. Six files here already did it this way.
REPO = next(_p for _p in HERE.parents if (_p / 'pyproject.toml').exists())
sys.path.insert(0, str(REPO / "scripts"))

ARMS = [("pea1989base", "baseline (v1, setback floored to 0)", "0.35"),
        ("pea1989keep", "insert, 1996 crest kept", "#2166ac"),
        ("pea1989lower", "insert, crest shaved", "#b2182b")]

DOMAINS = (84, 85, 86)
CONTROL = (40, 50, 60)
BUFFER = 15          # padded domains each side; npz index = GIS + BUFFER - 1
START_YEAR = 1984
EVENT_YEAR = {84: 1989, 85: 1989, 86: 1989}
OUT = REPO / "output" / "experiments" / "pea1989_crest"


def load_arm(arm: str):
    # experiments/2026-09-02-pea1989/<member>/ since 2026-09-16; the two
    # older layouts (arms/<arm>/ and the loose <arm>/) are tried after it.
    member = arm.replace("pea1989", "", 1)
    roots = [REPO / "output" / "raw_runs" / "experiments" / "2026-09-02-pea1989" / member,
             REPO / "output" / "raw_runs" / "arms" / arm,
             REPO / "output" / "raw_runs" / arm]
    hits = []
    for root in roots:
        hits = glob.glob(str(root / "1984_2004" / "calibBE" / "*" / "*.npz"))
        if hits:
            break
    if not hits:
        raise SystemExit(
            "\nno run found for arm {!r}. pea1989keep/pea1989lower and every "
            "other insert arm were deleted 2026-09-07 (only unmodified "
            "topography is kept); see 1-barrier3d-domains/"
            "archive_purge_20260907.csv. Only pea1989base can be re-run, with "
            "HAT_run_crest_experiment.py.\n".format(arm))
    return np.load(hits[0], allow_pickle=True)["cascade"][0], Path(hits[0]).parent.name


def series(mgr, name, nt):
    v = np.asarray(getattr(mgr, name, []), dtype=float)
    if v.size == 0:
        return np.full(nt, np.nan)
    out = np.full(nt, np.nan)
    out[:min(nt, v.size)] = v[:min(nt, v.size)]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--arm-suffix", default="",
                    help="\"noreloc\" reads the arms run with the "
                         "prescribed HATTERAS_ROAD_EVENTS switched off")
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    loaded = {}
    for arm, label, colour in ARMS:
        c, run_name = load_arm(arm + args.arm_suffix)
        loaded[arm] = c
        print("  {:14s} {}".format(arm, run_name))

    nt = len(np.asarray(loaded[ARMS[0][0]].barrier3d[BUFFER].x_s_TS))
    years = START_YEAR + np.arange(nt)

    # ---------------------------------------------------------------- figure 1
    fig, axes = plt.subplots(len(DOMAINS), 2, figsize=(15, 4.2 * len(DOMAINS)),
                             squeeze=False)
    rows = []
    for r, D in enumerate(DOMAINS):
        i = D + BUFFER - 1
        ax_s, ax_d = axes[r]
        for arm, label, colour in ARMS:
            c = loaded[arm]
            mgr = c.roadways[i]
            b3d = c.barrier3d[i]
            if mgr is None:
                continue
            sb = series(mgr, "_road_setback_TS", nt)
            rel = series(mgr, "_road_relocated_TS", nt)
            ax_s.plot(years[:sb.size], sb, "-", color=colour, lw=2, label=label)
            jumps = np.flatnonzero(np.nan_to_num(rel) > 0)
            for j in jumps:
                ax_s.plot(years[j], sb[j], "o", color=colour, ms=7,
                          mec="k", mew=.8, zorder=5)
            crest = np.array([np.max(b3d.DuneDomain[t, :, :]) * 10.0
                              if t < b3d.DuneDomain.shape[0] else np.nan
                              for t in range(nt)])
            ax_d.plot(years, crest, "-", color=colour, lw=2, label=label)

            first = int(jumps[0]) if jumps.size else None
            rows.append({
                "domain": D, "arm": arm,
                "setback_t0_m": round(float(sb[0]), 1),
                "n_relocations": int(jumps.size),
                "first_relocation_year": (START_YEAR + first) if first is not None else "",
                "prescribed_year": EVENT_YEAR[D],
                "pre_empts_event": (first is not None
                                    and START_YEAR + first < EVENT_YEAR[D]),
                "mean_dune_crest_m": round(float(np.nanmean(crest)), 3),
                "final_dune_crest_m": round(float(crest[-1]), 3),
            })

        ax_s.axvline(EVENT_YEAR[D], color="k", ls="--", lw=1.4)
        ax_s.text(EVENT_YEAR[D] + 0.3, ax_s.get_ylim()[1] * 0.92,
                  ("prescribed\n1989 relocation" if not args.arm_suffix
                   else "1989: when NC-12\nwas ACTUALLY moved"),
                  fontsize=8, va="top")
        ax_s.axhline(0, color="#3b6ea5", lw=1, ls=":")
        ax_s.set_ylabel("road setback (m)")
        ax_s.set_xlabel("year")
        ax_s.set_title("GIS {} | setback through time\n"
                       "markers = a relocation fired that year".format(D),
                       fontsize=11)
        ax_s.legend(fontsize=8)
        ax_s.grid(alpha=.3)

        ax_d.axvline(EVENT_YEAR[D], color="k", ls="--", lw=1.4)
        ax_d.set_ylabel("max dune crest (m above berm)")
        ax_d.set_xlabel("year")
        ax_d.set_title("GIS {} | dune crest".format(D), fontsize=11)
        ax_d.legend(fontsize=8)
        ax_d.grid(alpha=.3)

    fig.suptitle(("Does the row insert let NC-12 reach its 1989 relocation?\n"
                  "baseline vs insert, GIS 84-86, relocations ON in all three arms")
                 if not args.arm_suffix else
                 ("Left to itself, when does the module relocate NC-12?\n"
                  "prescribed events OFF -- the relocation year is a PREDICTION, "
                  "to be scored against 1989"),
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    tag = ("_" + args.arm_suffix) if args.arm_suffix else ""
    f1 = out_dir / "HAT_crest_experiment_setback{}.png".format(tag)
    fig.savefig(f1, dpi=130)
    print("  wrote {}".format(f1))

    # ---------------------------------------------------------------- figure 2
    fig2, ax2 = plt.subplots(2, 1, figsize=(14, 9))
    gis = np.arange(1, 91)
    base_rate = None
    for arm, label, colour in ARMS:
        c = loaded[arm]
        rate = []
        for D in gis:
            xs = np.asarray(c.barrier3d[D + BUFFER - 1].x_s_TS, dtype=float) * 10.0
            t = np.arange(xs.size)
            rate.append(np.polyfit(t, xs, 1)[0] if xs.size > 2 else np.nan)
        rate = np.array(rate)
        ax2[0].plot(gis, rate, "-", color=colour, lw=1.8, label=label)
        if base_rate is None:
            base_rate = rate
        else:
            ax2[1].plot(gis, rate - base_rate, "-", color=colour, lw=1.8,
                        label="{} - baseline".format(label))
    for a in ax2:
        for D in DOMAINS:
            a.axvspan(D - .5, D + .5, color="#ffd27f", alpha=.45, zorder=0)
        a.grid(alpha=.3)
        a.set_xlabel("GIS domain")
        a.legend(fontsize=8.5)
    ax2[0].set_ylabel("shoreline change rate (m/yr, LRR)")
    ax2[0].set_title("Island-wide rate, all three arms  "
                     "(shaded = the three inserted domains)", fontsize=11)
    ax2[1].set_ylabel("difference from baseline (m/yr)")
    ax2[1].axhline(0, color="k", lw=1)
    ax2[1].set_title("CONTROL: the insert touches 3 domains of 90, so everything "
                     "outside the shaded bands should sit on zero", fontsize=11)
    fig2.tight_layout()
    f2 = out_dir / "HAT_crest_experiment_island{}.png".format(tag)
    fig2.savefig(f2, dpi=130)
    print("  wrote {}".format(f2))

    csv_p = out_dir / "HAT_crest_experiment_summary{}.csv".format(tag)
    with open(csv_p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("  wrote {}".format(csv_p))

    print()
    hdr = "{:>4} {:>14} {:>10} {:>6} {:>10} {:>10}".format(
        "GIS", "arm", "setback t0", "n_rel", "first rel", "pre-empts")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print("{:>4} {:>14} {:>9.0f}m {:>6} {:>10} {:>10}".format(
            r["domain"], r["arm"].replace("pea1989", ""), r["setback_t0_m"],
            r["n_relocations"], r["first_relocation_year"] or "-",
            "YES" if r["pre_empts_event"] else "no"))


if __name__ == "__main__":
    main()
