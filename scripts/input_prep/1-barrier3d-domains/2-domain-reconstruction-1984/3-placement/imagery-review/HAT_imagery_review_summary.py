#!/usr/bin/env python3
r"""
HAT_imagery_review_summary.py
==============================================================================
What the imagery review says, once the sheet is filled: the verdicts along the
island against the rows the footprint changes, and the reviewer's toe picks
against the shift the footprint measured from the digitized lines.

Reads 2-domain-reconstruction-1984/3-placement/imagery-review/imagery_review_1984.csv (written by
HAT_imagery_review_1984.py, filled by hand or through HAT_imagery_review_gui.py)
and footprint_1984_by_domain.csv. Writes nothing model-facing. Runs on a
half-filled or empty sheet and says so: unjudged domains are drawn hollow and
counted as such.

FIGURE  figures/3-placement/imagery-review/island/HAT_imagery_review_summary.png
    (a) every reviewed domain along the island: the signed rows N as a bar,
        coloured by the verdict `extra_width_was` (where the 1984 island was
        wider / narrower than in 1997): behind the road (v3's placement),
        between the crest and the road, seaward of the crest, or no real
        change; hollow where not yet judged. Controls (N = 0) as markers,
        filled where their `dune_field_change` says "same". A cross marks a
        domain whose `placement_ok` is "no".
    (b) THE PLACEMENT QUESTION, from the picks: per measured domain, where the
        1984 width sat, as stacked bars of N x 10 m - the part lost or gained
        in the dune band (toe to back of dune), in the strip from the back of
        the dune to the road, and the remainder behind the road - with the
        island-wide medians of the three shares by sign printed in the panel.
        The rule reads off the medians.
    (c) which feature moved: the shift of the toe, the back of the dune and
        the road edge (1997 - 1984, + = the 1984 feature seaward) against
        N x 10 m, per measured domain.

REPORT  2-domain-reconstruction-1984/3-placement/imagery-review/HAT_imagery_review_summary.txt
    counts by verdict for the changed domains (adds and removals apart),
    placement_ok, confidence, the controls, the domains where v3's placement
    is contradicted, and the agreement of the picks with N.

USAGE
    python HAT_imagery_review_summary.py
    python HAT_imagery_review_summary.py --sheet <other.csv> --out-dir <dir>   # a copy, for testing
==============================================================================
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "2-extent"))   # HAT_footprint_1984, for the community bands
import HAT_imagery_review_1984 as R  # noqa: E402
from site_layer.hat_topo_version import insert_figures_dir  # noqa: E402

off = R.off
CATEGORY = {                       # extra_width_was -> colour, label
    "behind_road": ("#4d4d4d", "behind the road (v3's placement)"),
    "crest_to_road": ("#e08214", "between the crest and the road"),
    "seaward_of_crest": (R.C_SEA_ALT, "seaward of the crest (the alternative)"),
    "none": ("#bdbdbd", "no real change (artefact of the lines)"),
}
COMMUNITY_C = "#dfe9f3"


def _s(v) -> str:
    return v.strip() if isinstance(v, str) else ""


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return np.nan


def load(sheet: Path) -> pd.DataFrame:
    df = pd.read_csv(sheet, dtype=str).fillna("")
    df["domain"] = df["domain"].astype(int)
    df = df.set_index("domain").sort_index()
    for c in R.VERDICT_COLS + R.MEASURED_COLS:
        if c not in df:
            df[c] = ""
    df["n_cells"] = df["n_cells"].astype(float).astype(int)
    for c in ("shift_m_median", "toe_shift_m", "back_shift_m", "road_shift_m", "toe_shift_cells",
              "lost_dune_band_m", "lost_back_to_road_m", "lost_behind_road_m"):
        df[c + "_f"] = df[c].map(_f)
    return df


def fig_summary(df: pd.DataFrame, out: Path) -> Path:
    off.apply_style()
    fig = plt.figure(figsize=(15.0, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1.0])
    ax = fig.add_subplot(gs[0, :])
    ab = fig.add_subplot(gs[1, 0])
    ac = fig.add_subplot(gs[1, 1])

    # ---- (a) verdict along the island
    try:
        import HAT_footprint_1984 as F
        F._community_bands(ax)
    except Exception:
        pass
    ch = df[df["n_cells"] != 0]
    ct = df[df["n_cells"] == 0]
    for d, r in ch.iterrows():
        v = _s(r["extra_width_was"])
        col = CATEGORY.get(v, (None,))[0]
        ax.bar(d, r["n_cells"], width=0.8, facecolor=col if col else "white",
               edgecolor=col if col else "0.5", lw=0.8 if col else 1.0, zorder=3)
        if _s(r["placement_ok"]) == "no":
            ax.plot(d, r["n_cells"] + (0.5 if r["n_cells"] > 0 else -0.5), marker="x", color="#b2182b",
                    ms=7, mew=1.8, zorder=5, ls="none")
    for d, r in ct.iterrows():
        same = _s(r["dune_field_change"]) == "same"
        judged = bool(_s(r["dune_field_change"]))
        ax.plot(d, 0, marker="o", ms=6, mfc="0.3" if same else ("white" if not judged else "#e08214"),
                mec="0.3", mew=1.0, zorder=4, ls="none")
    ax.axhline(0, color=off.INK, lw=0.8, zorder=2)
    ax.set_xlim(0, 91)
    ax.set_xlabel("GIS domain (south at left)")
    ax.set_ylabel("rows in the footprint, N (signed)")
    ax.set_xticks(range(5, 91, 5))
    ax.grid(True, axis="y", alpha=0.4)
    n_j = int(ch["extra_width_was"].map(_s).astype(bool).sum())
    off._title(ax, 0, f"where the 1984 width was, by the reviewer: {n_j} of {len(ch)} changed domains judged, "
                      f"{int(ct['dune_field_change'].map(_s).astype(bool).sum())} of {len(ct)} controls")
    handles = [Patch(facecolor=c, edgecolor=c, label=lab) for c, lab in CATEGORY.values()]
    handles += [Patch(facecolor="white", edgecolor="0.5", label="not yet judged"),
                Line2D([0], [0], marker="o", mfc="0.3", mec="0.3", ls="none", label="control, judged 'same'"),
                Line2D([0], [0], marker="o", mfc="#e08214", mec="0.3", ls="none", label="control, judged changed"),
                Line2D([0], [0], marker="o", mfc="white", mec="0.3", ls="none", label="control, not judged"),
                Line2D([0], [0], marker="x", color="#b2182b", mew=1.8, ls="none", label="placement_ok = no")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=5, fontsize=8, frameon=False)

    # ---- (b) where the 1984 width sat, from the picks, as shares of N x 10 m
    BAND_C = {"lost_dune_band_m": (R.C_SEA_ALT, "dune band (toe to back of dune)"),
              "lost_back_to_road_m": ("#e08214", "back of dune to road"),
              "lost_behind_road_m": ("#4d4d4d", "behind the road (remainder)")}
    m = df[np.isfinite(df["lost_behind_road_m_f"]) & (df["n_cells"] != 0)]
    if len(m):
        m = m.sort_values("n_cells")
        xs = np.arange(len(m))
        for j, (d, r) in enumerate(m.iterrows()):
            pos_base, neg_base = 0.0, 0.0
            for c, (col, _) in BAND_C.items():
                v = r[c + "_f"]
                if v >= 0:
                    ab.bar(j, v, bottom=pos_base, width=0.8, facecolor=col, edgecolor="white", lw=0.5, zorder=3)
                    pos_base += v
                else:
                    ab.bar(j, v, bottom=neg_base, width=0.8, facecolor=col, edgecolor="white", lw=0.5, zorder=3)
                    neg_base += v
            ab.plot([j - 0.4, j + 0.4], [r["n_cells"] * R.CELL_M] * 2, color=off.INK, lw=1.4, zorder=5)
        ab.set_xticks(xs)
        ab.set_xticklabels([str(d) for d in m.index], fontsize=7, rotation=90)
        ab.set_xlabel("GIS domain, sorted by N")
        ab.set_ylabel("1984 width relative to 1997, m (+ = 1984 wider)")
        ab.axhline(0, color=off.INK, lw=0.6)
        txt = []
        for name, sub in (("adds", m[m["n_cells"] > 0]), ("removals", m[m["n_cells"] < 0])):
            if len(sub):
                tot = (sub["n_cells"] * R.CELL_M).abs()
                shares = {lab.split(" (")[0]: float(np.median(sub[c + "_f"].abs() / tot))
                          for c, (_, lab) in BAND_C.items()}
                txt.append(f"{name} ({len(sub)}): " + ", ".join(f"{k} {v:.0%}" for k, v in shares.items()))
        ab.text(0.02, 0.97, "median share of |N| x 10 m\n" + "\n".join(txt), transform=ab.transAxes,
                ha="left", va="top", fontsize=8, color=off.INK,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.2"))
        handles = [Patch(facecolor=col, label=lab) for col, lab in BAND_C.values()]
        handles.append(Line2D([0], [0], color=off.INK, lw=1.4, label="N x 10 m (the footprint)"))
        ab.legend(handles=handles, loc="lower right", fontsize=7.5)
        off._title(ab, 1, f"where the 1984 width sat, from your picks: {len(m)} domains")
    else:
        ab.text(0.5, 0.5, "no domain has all six picks yet", transform=ab.transAxes, ha="center", va="center",
                color=off.INK_MUTED)
        off._title(ab, 1, "where the 1984 width sat, from your picks")
    ab.grid(True, axis="y", alpha=0.4)

    # ---- (c) which feature moved, against N x 10
    KIND_C = {"toe_shift_m": (R.C_SEA_ALT, "o", "toe"), "back_shift_m": ("#e08214", "s", "back of dune"),
              "road_shift_m": ("#4d4d4d", "D", "road edge")}
    w = df[(df["n_cells"] != 0) & (np.isfinite(df["toe_shift_m_f"]) | np.isfinite(df["back_shift_m_f"])
                                   | np.isfinite(df["road_shift_m_f"]))]
    if len(w):
        lim = max(60.0, float(np.nanmax(np.abs(np.concatenate(
            [w[c + "_f"].to_numpy() for c in KIND_C] + [(w["n_cells"] * R.CELL_M).to_numpy()])))) + 10)
        ac.fill_between([-lim, lim], [-lim - R.CELL_M, lim - R.CELL_M], [-lim + R.CELL_M, lim + R.CELL_M],
                        color="0.9", lw=0, zorder=1)
        ac.plot([-lim, lim], [-lim, lim], color="0.4", lw=0.8, zorder=2, label="feature moved by all of N x 10 m")
        rng = np.random.default_rng(0)
        for c, (col, mk, lab) in KIND_C.items():
            ok = np.isfinite(w[c + "_f"])
            jitter = rng.uniform(-2.0, 2.0, size=int(ok.sum()))
            ac.plot(w.loc[ok, "n_cells"] * R.CELL_M + jitter, w.loc[ok, c + "_f"], marker=mk, ms=6, color=col,
                    ls="none", zorder=3, label=lab, alpha=0.85)
        ac.set_xlim(-lim, lim)
        ac.set_ylim(-lim, lim)
        ac.set_aspect("equal")
        ac.legend(loc="upper left", fontsize=7.5)
        off._title(ac, 2, f"which feature moved: {len(w)} domains")
    else:
        ac.text(0.5, 0.5, "no picks yet", transform=ac.transAxes, ha="center", va="center", color=off.INK_MUTED)
        off._title(ac, 2, "which feature moved")
    ac.axhline(0, color=off.INK, lw=0.6)
    ac.axvline(0, color=off.INK, lw=0.6)
    ac.set_xlabel("N x 10 m (+ = rows added)")
    ac.set_ylabel("feature shift 1997 - 1984, m (+ = 1984 seaward)")
    ac.grid(True, alpha=0.4)

    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def write_report(df: pd.DataFrame, fig: Path, out: Path, sheet: Path) -> Path:
    ch = df[df["n_cells"] != 0]
    ct = df[df["n_cells"] == 0]
    L = [f"HAT_imagery_review_summary.txt - what the imagery review says ({datetime.now():%Y-%m-%d %H:%M})",
         f"sheet: {sheet}", ""]

    def counts(sub: pd.DataFrame, col: str) -> str:
        v = sub[col].map(_s)
        vals = [x.strip() for x in R.VERDICT_VOCAB[col].split("|")]
        parts = [f"{k} {int((v == k).sum())}" for k in vals if k]
        return ", ".join(parts) + f", blank {int((v == '').sum())}"

    L.append(f"CHANGED DOMAINS ({len(ch)}: {int((ch['n_cells'] > 0).sum())} add, {int((ch['n_cells'] < 0).sum())} remove)")
    for name, sub in (("rows added", ch[ch["n_cells"] > 0]), ("rows removed", ch[ch["n_cells"] < 0])):
        L.append(f"  {name} ({len(sub)})")
        for c in ("extra_width_was", "dune_field_change", "edge_moved", "placement_ok", "confidence"):
            L.append(f"    {c:18s} {counts(sub, c)}")
    L.append("")
    L.append(f"CONTROLS ({len(ct)})")
    for c in ("dune_field_change", "edge_moved", "confidence"):
        L.append(f"    {c:18s} {counts(ct, c)}")
    L.append("")
    no = ch[ch["placement_ok"].map(_s) == "no"]
    L.append(f"V3'S PLACEMENT CONTRADICTED (placement_ok = no): {len(no)} domain(s)")
    for d, r in no.iterrows():
        L.append(f"  GIS {d:>2}  N {r['n_cells']:+d}  width was: {_s(r['extra_width_was']) or '-':18s} "
                 f"confidence {_s(r['confidence']) or '-':7s} {_s(r['notes'])}")
    L.append("")
    m = df[np.isfinite(df["lost_behind_road_m_f"]) & (df["n_cells"] != 0)]
    L.append(f"PICKS: {len(m)} changed domain(s) with all six picks (toe, back of dune, road; 1984 and 1997)")
    if len(m):
        for name, sub in (("rows added", m[m["n_cells"] > 0]), ("rows removed", m[m["n_cells"] < 0])):
            if not len(sub):
                continue
            tot = (sub["n_cells"] * R.CELL_M).abs()
            L.append(f"  {name} ({len(sub)}): median share of |N| x 10 m in the dune band "
                     f"{np.median(sub['lost_dune_band_m_f'].abs() / tot):.0%}, back of dune to road "
                     f"{np.median(sub['lost_back_to_road_m_f'].abs() / tot):.0%}, behind the road "
                     f"{np.median(sub['lost_behind_road_m_f'].abs() / tot):.0%}; bands suggest: "
                     + ", ".join(f"{k} {v}" for k, v in sub["band_suggests"].map(_s).value_counts().items()))
        t = df[np.isfinite(df["toe_shift_m_f"])]
        diff = t["toe_shift_m_f"] - t["shift_m_median_f"]
        L.append(f"  toe pick shift minus digitized-line shift ({len(t)}): median {np.median(diff):+.1f} m, "
                 f"p10 {np.percentile(diff, 10):+.1f}, p90 {np.percentile(diff, 90):+.1f}; "
                 f"{int((np.abs(diff) <= R.CELL_M).sum())} within one cell")
        L.append("")
        L.append("  domain   N   toe shift  back shift  road shift | lost: dune band  dune->road  behind road "
                 "| suggests / verdict")
        for d, r in m.iterrows():
            L.append(f"  {d:>6}  {r['n_cells']:+3d}   {r['toe_shift_m_f']:+8.1f}   {r['back_shift_m_f']:+8.1f}   "
                     f"{r['road_shift_m_f']:+8.1f} |   {r['lost_dune_band_m_f']:+8.1f}   "
                     f"{r['lost_back_to_road_m_f']:+8.1f}   {r['lost_behind_road_m_f']:+8.1f} | "
                     f"{_s(r['band_suggests']):13s} {_s(r['extra_width_was'])}")
    L.append("")
    L.append(f"figure: {fig}")
    out.write_text("\n".join(L) + "\n", encoding="utf-8")
    return out


def run(sheet: Path | None = None, out_dir: Path | None = None) -> tuple[Path, Path]:
    sheet = sheet or R.SHEET
    if not sheet.is_file():
        raise SystemExit(f"{sheet} not found - run HAT_imagery_review_1984.py first")
    df = load(sheet)
    if out_dir is None:
        fig = fig_summary(df, insert_figures_dir(R.PRODUCT, "3-placement", "imagery-review/island")
                          / "HAT_imagery_review_summary.png")
        rep = write_report(df, fig, R.STEP_DIR / "HAT_imagery_review_summary.txt", sheet)
        R.upsert_caption(
            "## `HAT_imagery_review_summary.png` (3-placement/imagery-review/island)",
            "What the imagery review says so far. (a) Every reviewed domain along the island, south at "
            "left: the signed rows of the footprint as a bar, coloured by the reviewer's verdict on where "
            "the 1984 island was wider (rows added) or narrower (rows removed) than in 1997: behind the "
            "road (v3's placement, dark), between the crest and the road (orange), seaward of the crest "
            "(purple, the alternative placement), or no real change (light grey, the two digitized lines "
            "differ but the photographs do not); hollow bars are not yet judged, a red cross marks a "
            "domain whose verdict says v3's placement is wrong. The ten controls, unchanged neighbours "
            "of the changed runs, are circles at zero: filled dark where the reviewer saw no change, "
            "orange where they saw change, hollow where not judged. Communities banded. (b) The placement question, from the reviewer's picks: for every domain with all six picks (toe, "
            "back of the dune band and road edge, in 1984 and 1997), where the 1984 width sat, as stacked bars "
            "of N × 10 m: the part lost or gained in the dune band (purple), in the strip from the back of the "
            "dune to the road (orange), and the remainder behind the road (dark), with N × 10 m as a tick; "
            "domains sorted by N; the island-wide median share of each band by sign printed in the panel. "
            "The rule reads off the medians. (c) Which feature moved: the shift of the toe (circles), the back "
            "of the dune (squares) and the road edge (diamonds), 1997 minus 1984, positive where the 1984 "
            "feature lay seaward, against N × 10 m; the 1:1 line and its one-cell band mark a feature that "
            "moved by all of N. "
            "Counts and the per-domain table in "
            "`HAT_imagery_review_summary.txt`; the picks and verdicts in `imagery_review_1984.csv`.")
    else:
        out_dir.mkdir(parents=True, exist_ok=True)
        fig = fig_summary(df, out_dir / "HAT_imagery_review_summary.png")
        rep = write_report(df, fig, out_dir / "HAT_imagery_review_summary.txt", sheet)
    print(f"wrote {fig}\nwrote {rep}")
    return fig, rep


def main() -> None:
    ap = argparse.ArgumentParser(description="summarize the imagery review sheet")
    ap.add_argument("--sheet", default="", help="a sheet other than the live one (testing)")
    ap.add_argument("--out-dir", default="", help="where to write when --sheet is given")
    a = ap.parse_args()
    run(Path(a.sheet) if a.sheet else None, Path(a.out_dir) if a.out_dir else None)


if __name__ == "__main__":
    main()
