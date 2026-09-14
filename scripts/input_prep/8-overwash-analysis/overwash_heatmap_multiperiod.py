"""
overwash_heatmap_multiperiod.py
==============================================================================
Observed overwash on Hatteras Island, per image and per CASCADE domain, with
the storm record beside it. One figure per period; all three are written in
one run.

    python overwash_heatmap_multiperiod.py            # period1, period2, combined
    python overwash_heatmap_multiperiod.py period2    # one of them
    python overwash_heatmap_multiperiod.py combined --uniform-years

PANELS
    (a) the observation matrix: one row per image (labelled year and date),
        one column per domain; red = overwash present, white = assessed and
        absent, dotted = the image does not reach that domain, hatched rows =
        no image that year. Runs of years without an image are collapsed to
        one thin row (--uniform-years keeps a row per year).
    (b) how many domains each image shows overwashed, against how many it
        assessed (grey).
    (c) the named storms from the reference sheet, at their year. Bold with a
        filled dot: the image on that row was taken after the storm, so what
        the row shows includes it. Light italic: the image predates the storm
        or there is none that year; the thin line then leads down to the
        first image taken after it, which is where its effects can appear.
    (d) how many images show each domain overwashed.

STYLE
    hat_figure_style; no in-image title or footnote. The words are in
    data/hatteras_init/8-overwash-analysis/CAPTIONS.md, written by this script.

OUTPUT   data/hatteras_init/8-overwash-analysis/
    figures/heatmaps/overwash_heatmap_<period>.png
    tables/overwash_observations.csv   one row per image and domain
    tables/storms_by_image.csv         which image first shows each storm
    CAPTIONS.md

The 2026-05 version of this script had a comparison mode for a modelled
overwash matrix that was never produced. It is gone; a model comparison
should align to tables/overwash_observations.csv.
==============================================================================
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = next(
    _p for _p in HERE.parents
    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))

from hat_figure_style import (C, DOMAIN_AXIS_LABEL, FIG_H_MAX, INK,   # noqa: E402
                              INK_MUTED, _title, apply_style, figsize,
                              open_frame, save, spines_for_image,
                              town_bands)
from overwash_data import (OUT_DIR, PERIODS, SECTIONS, assign_capture,   # noqa: E402
                           load_observations, load_storms,
                           observations_long, storms_table, upsert_caption)

FIG_DIR = OUT_DIR / "figures" / "heatmaps"
TAB_DIR = OUT_DIR / "tables"

# ---------------------------------------------------------------- colours
CLR_OW = C["ACCENT"]          # overwash present
CLR_NONE = "#ffffff"          # assessed, absent
CLR_GAP_FACE = "#f2f2f2"      # no image that year
CLR_GAP_EDGE = "#b8b8b8"
CLR_PART_FACE = "#f7f7f7"     # image does not cover the domain
CLR_PART_EDGE = "#9a9a9a"
CLR_ASSESSED = C["BASE_FILL"]
CLR_INK = INK
CLR_MUTED = INK_MUTED
CLR_PERIOD = ("#dcdcdc", "#bdbdbd")

# Storm intensity: hurricanes on a red ramp, nor'easters on a blue ramp,
# tropical/extratropical storms grey. (size, colour)
CAT_STYLE = {
    "H5": (8.0, "#67000d"), "H4": (7.2, "#a50f15"), "H3": (6.4, "#cb181d"),
    "H2": (5.6, "#ef3b2c"), "H1": (4.8, "#fb6a4a"),
    "NE 5": (8.0, "#08306b"), "NE 4": (7.2, "#2166ac"), "NE 3": (6.4, "#4292c6"),
    "TS": (4.5, "#737373"), "ET": (5.6, "#737373"),
}

ROW_IN = 0.30          # height of one image row, inches, at most
STORM_LINE_IN = 0.115  # one line of storm type, inches
TICK_EVERY = 10        # domain tick spacing


# ============================================================ row building
def build_rows(obs, period, collapse):
    """
    Display rows, top to bottom: dicts with kind ('obs' | 'gap'), label,
    years (list), obs (index into `obs` or None), height (row units, set
    later once the storms are placed).
    """
    lo, hi = period
    have = obs[obs["Year"].between(lo, hi)]
    rows = []
    for yr in range(lo, hi + 1):
        sub = have[have["Year"] == yr]
        if len(sub):
            for i, r in sub.iterrows():
                rows.append(dict(kind="obs", years=[yr], obs=int(i),
                                 label=f"{yr}  {r['Imagery_Date']:%d %b}"
                                       + ("*" if r["poor"] else "")))
        elif collapse and rows and rows[-1]["kind"] == "gap":
            rows[-1]["years"].append(yr)
        else:
            rows.append(dict(kind="gap", years=[yr], obs=None, label=""))
    for r in rows:
        if r["kind"] == "gap":
            y = r["years"]
            r["label"] = str(y[0]) if len(y) == 1 else f"{y[0]}–{y[-1]}"
    return rows


def place_storms(rows, storms, obs, period):
    """
    Attach each storm in the period to a display row and decide whether that
    row's image shows it. Returns {row_index: [storm, ...]} and sets
    storm['matched'], storm['capture_row'].
    """
    lo, hi = period
    obs_to_row = {r["obs"]: k for k, r in enumerate(rows) if r["kind"] == "obs"}
    year_rows = {}
    for k, r in enumerate(rows):
        for y in r["years"]:
            year_rows.setdefault(y, []).append(k)

    placed = {}
    for s in storms:
        if not (lo <= s["year"] <= hi):
            continue
        cap = s["capture"]
        cap_row = obs_to_row.get(cap) if cap is not None else None
        cands = year_rows[s["year"]]
        own = cands[0]
        if rows[own]["kind"] == "obs":
            # A year with several images: sit on the first image after the
            # storm, else on the last image of the year.
            after = [k for k in cands
                     if obs.loc[rows[k]["obs"], "Imagery_Date"] >= s["end"]]
            own = after[0] if after else cands[-1]
        s["row"] = own
        s["matched"] = (cap_row is not None and cap_row == own)
        s["capture_row"] = cap_row
        # the first image after the storm exists but is outside this period
        s["capture_date"] = (obs.loc[cap, "Imagery_Date"]
                             if cap is not None and cap_row is None else None)
        placed.setdefault(own, []).append(s)
    return placed


def set_heights(rows, placed, slot=0.34):
    """Row heights in row units. A gap row holds `slot` units per storm named
    on it; `slot` is set by make_figure so that one storm gets one line of
    type however far the rows have had to be squeezed."""
    for k, r in enumerate(rows):
        if r["kind"] == "obs":
            r["height"] = 1.0
        else:
            n = len(placed.get(k, []))
            r["height"] = max(0.55, slot * n + 0.22)
    edges = np.concatenate([[0.0], np.cumsum([r["height"] for r in rows])])
    for k, r in enumerate(rows):
        r["y0"], r["y1"] = edges[k], edges[k + 1]
        r["yc"] = 0.5 * (edges[k] + edges[k + 1])
    return edges


# ================================================================= drawing
def draw_matrix(ax, rows, domains, matrix):
    n_d = len(domains)
    x0, x1 = domains[0] - 0.5, domains[-1] + 0.5
    ax.set_facecolor(CLR_NONE)
    for r in rows:
        h = r["y1"] - r["y0"]
        if r["kind"] == "gap":
            ax.add_patch(mpatches.Rectangle(
                (x0, r["y0"]), n_d, h, facecolor=CLR_GAP_FACE,
                edgecolor=CLR_GAP_EDGE, hatch="////", linewidth=0, zorder=1))
            continue
        v = matrix[r["obs"]]
        for j, d in enumerate(domains):
            if np.isnan(v[j]):
                ax.add_patch(mpatches.Rectangle(
                    (d - 0.5, r["y0"]), 1, h, facecolor=CLR_PART_FACE,
                    edgecolor=CLR_PART_EDGE, hatch="....", linewidth=0,
                    zorder=2))
            elif v[j] >= 0.5:
                ax.add_patch(mpatches.Rectangle(
                    (d - 0.5, r["y0"]), 1, h, facecolor=CLR_OW,
                    edgecolor="none", linewidth=0, zorder=2))
    for r in rows:
        ax.axhline(r["y1"], color="white", lw=0.5, zorder=3)
    for _, lo, _, _ in SECTIONS[1:]:
        ax.axvline(lo - 0.5, color="#bfbfbf", lw=0.8, zorder=4)
    # the villages as light bands, named once under panel (d)
    town_bands(ax, label=False)

    ax.set_xlim(x0, x1)
    ax.set_ylim(rows[-1]["y1"], 0)
    ticks = [d for d in domains if d == 1 or d % TICK_EVERY == 0]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.tick_params(axis="x", length=2.5, pad=2)
    ax.set_yticks([r["yc"] for r in rows])
    ax.set_yticklabels([r["label"] for r in rows], fontsize=7.5)
    for lab, r in zip(ax.get_yticklabels(), rows):
        if r["kind"] == "gap":
            lab.set_color(CLR_MUTED)
            lab.set_fontstyle("italic")
            lab.set_fontsize(7.0)
    ax.tick_params(axis="y", length=0, pad=4)
    spines_for_image(ax)
    _title(ax, 0, "overwash observed, per image and domain")


def draw_period_bars(ax, rows):
    """Two thin columns, Period 1 over 1984–2004 and Period 2 over 2004–2024.
    They overlap at the 2004 row, which is the last image of one and the
    first of the other."""
    ax.set_xlim(0, 2)
    ax.set_ylim(rows[-1]["y1"], 0)
    ax.axis("off")
    spans = [(1984, 2004, "Period 1  1984–2004", 0),
             (2004, 2024, "Period 2  2004–2024", 1)]
    for lo, hi, text, col in spans:
        ys = [r for r in rows if any(lo <= y <= hi for y in r["years"])]
        if not ys:
            continue
        y0, y1 = ys[0]["y0"], ys[-1]["y1"]
        ax.add_patch(mpatches.Rectangle((col + 0.15, y0), 0.7, y1 - y0,
                                        facecolor=CLR_PERIOD[col],
                                        edgecolor="none"))
        ax.text(col + 0.5, 0.5 * (y0 + y1), text, rotation=90,
                ha="center", va="center", fontsize=7, color=CLR_INK)


def draw_counts_per_image(ax, rows, matrix, n_d):
    ax.set_ylim(rows[-1]["y1"], 0)
    ax.set_xlim(0, n_d)
    for r in rows:
        if r["kind"] != "obs":
            continue
        v = matrix[r["obs"]]
        n_ass = int(np.sum(~np.isnan(v)))
        n_ow = int(np.nansum(v))
        h = 0.62 * (r["y1"] - r["y0"])
        ax.barh(r["yc"], n_ass, height=h, color=CLR_ASSESSED, zorder=1)
        ax.barh(r["yc"], n_ow, height=h, color=CLR_OW, zorder=2)
        if n_ow:
            ax.text(n_ow + 2.5, r["yc"], str(n_ow), va="center", ha="left",
                    fontsize=7, color=CLR_INK, zorder=3)
    ax.set_xticks([0, n_d // 2, n_d])
    ax.tick_params(axis="x", labelsize=7, length=2.5, pad=2)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("domains" + chr(10) + "overwashed", fontsize=7.5, labelpad=2)
    open_frame(ax)
    _title(ax, 1, "")          # the panel is one bar wide: letter only


def draw_storms(ax, rows, placed, obs, slot=0.34):
    ax.set_xlim(0, 1)
    ax.set_ylim(rows[-1]["y1"], 0)
    ax.axis("off")
    AXIS_X, LABEL_X = 0.17, 0.215
    for r in rows:
        ax.axhline(r["y1"], color="#ececec", lw=0.5, zorder=0)
    ax.axvline(AXIS_X, color="#cfcfcf", lw=1.0, zorder=1)

    # A storm's label position, and a small horizontal offset per connector
    # so the lines leading from several storms to the same image stay apart.
    conn_k = {}
    labels = []
    for k, storms in placed.items():
        r = rows[k]
        n = len(storms)
        if r["kind"] == "gap":
            ys = [r["y0"] + (r["y1"] - r["y0"]) * (i + 0.5) / n for i in range(n)]
        else:
            ys = [r["yc"] + (i - (n - 1) / 2) * slot for i in range(n)]
        for s, y in zip(storms, ys):
            labels.append((s, r, y))

    for s, r, y in labels:
        size, col = CAT_STYLE.get(s["cat"], (5.0, "#555555"))
        matched = s["matched"]
        name = s["name"].upper() if s["cat"].startswith("H") else s["name"]
        when = s["month"] if r["kind"] == "obs" else f"{s['month']} {s['year']}"
        # the category is the marker (see the key beneath), so the label
        # carries name and month only: the column is one text width wide.
        text = f"{name}  ({when})"
        ax.plot(AXIS_X, y, "o", ms=size, color=col, zorder=5,
                mec="white", mew=0.6, alpha=1.0 if matched else 0.45)
        ax.text(LABEL_X, y, text, va="center", ha="left", fontsize=7,
                color=col if matched else "#7a7a7a",
                fontweight="bold" if matched else "normal",
                fontstyle="normal" if matched else "italic", zorder=6)
        cap = s["capture_row"]
        if not matched and cap is not None:
            k = conn_k.get(cap, 0)
            conn_k[cap] = k + 1
            x = AXIS_X - 0.035 - 0.013 * k
            yc = rows[cap]["yc"]
            ax.plot([AXIS_X, x], [y, y], color=col, lw=0.6, alpha=0.5, zorder=2)
            ax.plot([x, x], [y, yc], color=col, lw=0.6, alpha=0.5, zorder=2)
            ax.plot([x, AXIS_X - 0.012], [yc, yc], color=col, lw=0.6,
                    alpha=0.5, zorder=2)
            ax.plot(AXIS_X - 0.012, yc, marker=">", ms=3.2, color=col,
                    alpha=0.6, zorder=3, mec="none")
        elif not matched:
            note = ("no image since" if s["capture_date"] is None else
                    f"next image {s['capture_date']:%b %Y}")
            ax.text(0.985, y, note, va="center", ha="right",
                    fontsize=7, color=CLR_MUTED, fontstyle="italic")

    targets = {s["capture_row"] for s, _, _ in labels
               if not s["matched"] and s["capture_row"] is not None}
    for k, r in enumerate(rows):
        if k in placed:
            continue
        if r["kind"] == "obs":
            ax.plot(AXIS_X, r["yc"], "o", ms=3.4, color="#c8c8c8", zorder=4,
                    mec="white", mew=0.5)
            note = ("first image after them" if k in targets
                    else "no named storm")
            ax.text(LABEL_X, r["yc"], note, va="center", ha="left",
                    fontsize=7, color="#a8a8a8", fontstyle="italic")
        else:
            ax.plot(AXIS_X, r["yc"], "o", ms=3.0, color="none", zorder=4,
                    mec="#c8c8c8", mew=0.6)
    _title(ax, 2, "named storms")


def draw_storm_legend(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    items = [("H5", "hurricane 5"), ("H3", "hurricane 3"), ("H1", "hurricane 1"),
             ("NE 5", "nor'easter 5"), ("NE 3", "nor'easter 3"),
             ("TS", "TS / extratropical")]
    y = 0.94
    for cat, lab in items:
        size, col = CAT_STYLE[cat]
        ax.plot(0.06, y, "o", ms=size, color=col, mec="white", mew=0.6)
        ax.text(0.17, y, lab, va="center", ha="left", fontsize=7, color=CLR_INK)
        y -= 0.105
    y -= 0.05
    for line, colour in (("bold, filled dot: the image on", CLR_INK),
                         ("that row postdates the storm", CLR_INK),
                         ("faint italic: the image predates", CLR_MUTED),
                         ("it, and the line leads to the", CLR_MUTED),
                         ("first image that postdates it", CLR_MUTED)):
        ax.text(0.0, y, line, fontsize=7, color=colour, va="center")
        y -= 0.085


def draw_counts_per_domain(ax, rows, domains, matrix):
    idx = [r["obs"] for r in rows if r["kind"] == "obs"]
    sub = matrix[idx]
    n_ow = np.nansum(sub, axis=0)
    ax.bar(domains, n_ow, width=1.0, color=CLR_OW, zorder=2)
    for _, lo, _, _ in SECTIONS[1:]:
        ax.axvline(lo - 0.5, color="#bfbfbf", lw=0.8, zorder=1)
    ax.set_xlim(domains[0] - 0.5, domains[-1] + 0.5)
    top = int(np.nanmax(n_ow)) if n_ow.size else 1
    ax.set_ylim(0, max(top, 1) + 1.4)   # headroom for the village names
    ax.set_yticks(range(0, max(top, 1) + 1, 1 if top <= 6 else 2))
    ticks = [d for d in domains if d == 1 or d % TICK_EVERY == 0]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.tick_params(axis="both", length=2.5, pad=2)
    town_bands(ax)
    open_frame(ax)
    ax.set_xlabel(DOMAIN_AXIS_LABEL, labelpad=3)
    ax.set_ylabel("images with\noverwash")
    _title(ax, 3, "images with overwash, per domain")


def matrix_legend(fig, x, y):
    handles = [
        mpatches.Patch(fc=CLR_OW, ec="none", label="overwash present"),
        mpatches.Patch(fc=CLR_NONE, ec="#999999", lw=0.6,
                       label="assessed, no overwash"),
        mpatches.Patch(fc=CLR_PART_FACE, ec=CLR_PART_EDGE, hatch="....",
                       lw=0.0, label="image does not reach it"),
        mpatches.Patch(fc=CLR_GAP_FACE, ec=CLR_GAP_EDGE, hatch="////",
                       lw=0.0, label="no image that year"),
        mpatches.Patch(fc=CLR_ASSESSED, ec="none", label="domains assessed (b)"),
    ]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(x, y),
               ncol=2, frameon=False, handlelength=1.5, handleheight=1.0,
               columnspacing=1.4, borderaxespad=0, labelspacing=0.45)


# =================================================================== figure
def make_figure(tag, period, obs, domains, matrix, storms, collapse, out):
    rows = build_rows(obs, period, collapse)
    placed = place_storms(rows, storms, obs, period)
    combined = tag == "combined"
    n_d = len(domains)

    # ---- layout in inches at the printed width, then fractions ----
    fig_w = figsize("double")[0]
    PB = 0.34 if combined else 0.0          # period bars
    LM, G1, CW, G2, TW, RM = 1.00 + PB, 0.10, 0.55, 0.12, 1.95, 0.08
    HW = fig_w - LM - G1 - CW - G2 - TW - RM
    TM = 0.42
    XT, SB, G3, BH, BX, LG, BM = 0.22, 0.10, 0.34, 0.85, 0.34, 0.50, 0.10
    # One image row is ROW_IN tall unless the period has too many rows for a
    # page, when every row shrinks together and the figure still fits. The
    # storms named on a gap row each need a line of type whatever the rows
    # come out at, so the two are solved for together rather than fixed.
    avail = FIG_H_MAX - (TM + XT + SB + G3 + BH + BX + LG + BM)
    row_in, slot = ROW_IN, 0.34
    for _ in range(12):
        units = set_heights(rows, placed, slot)[-1]
        row_in = min(ROW_IN, avail / units)
        slot = max(0.34, STORM_LINE_IN / row_in)
    units = set_heights(rows, placed, slot)[-1]
    row_in = min(ROW_IN, avail / units)
    HH = units * row_in
    fig_h = TM + HH + XT + SB + G3 + BH + BX + LG + BM

    def ax_at(x, y_from_bottom, w, h):
        return fig.add_axes([x / fig_w, y_from_bottom / fig_h, w / fig_w, h / fig_h])

    fig = plt.figure(figsize=figsize("double", height=fig_h))
    y_hm = fig_h - TM - HH
    ax_hm = ax_at(LM, y_hm, HW, HH)
    ax_cnt = ax_at(LM + HW + G1, y_hm, CW, HH)
    ax_st = ax_at(LM + HW + G1 + CW + G2, y_hm, TW, HH)
    y_bh = y_hm - XT - SB - G3 - BH
    ax_bh = ax_at(LM, y_bh, HW, BH)
    ax_sl = ax_at(LM + HW + G1 + CW + G2, y_bh - BX, TW, BH + BX + G3 + SB)
    if combined:
        ax_pb = ax_at(0.08, y_hm, PB - 0.06, HH)
        draw_period_bars(ax_pb, rows)

    draw_matrix(ax_hm, rows, domains, matrix)
    draw_counts_per_image(ax_cnt, rows, matrix, n_d)
    draw_storms(ax_st, rows, placed, obs, slot)
    draw_counts_per_domain(ax_bh, rows, domains, matrix)
    draw_storm_legend(ax_sl)
    matrix_legend(fig, LM / fig_w, (y_bh - BX + 0.02) / fig_h)

    save(fig, out, close=True)
    print(f"  wrote {out.relative_to(REPO)} (+ .pdf)")
    return rows, placed


# ================================================================= captions
def caption_text(tag, period, rows, placed, obs, storms):
    lo, hi = period
    n_img = sum(1 for r in rows if r["kind"] == "obs")
    gaps = [r["label"] for r in rows if r["kind"] == "gap"]
    src = ("Period 1 rows are the Hapke and Henderson (2007) delineations, "
           "Period 2 rows are read from Google Earth imagery"
           if tag == "combined" else
           "Rows are the Hapke and Henderson (2007) delineations"
           if tag == "period1" else
           "Rows are read from Google Earth imagery")
    unmatched = [s for s in storms if lo <= s["year"] <= hi and not s["matched"]]
    ex = "; ".join(
        f"{s['name']} ({s['month']} {s['year']}) → "
        + (obs.loc[s['capture'], 'Imagery_Date'].strftime('%d %b %Y')
           if s['capture'] is not None else "no image yet")
        for s in unmatched[:6])
    return (
        f"Observed overwash on Hatteras Island, {lo}–{hi}. "
        f"(a) One row per image assessed ({n_img} images, labelled by year and "
        f"date; * poor or partial image), one column per CASCADE domain, 1 at "
        f"Cape Point to 90 at Pea Island; the light bands are the villages, "
        f"named under (d). Purple "
        f"where overwash was present, white where the image was assessed and "
        f"showed none, dotted where the image does not reach the domain, "
        f"hatched rows where no image exists; runs of such years are drawn as "
        f"one thin row ({', '.join(gaps)}). {src}. "
        f"(b) Domains overwashed in each image, over the number assessed (grey). "
        f"(c) The named storms of the reference sheet at their year, the marker "
        f"giving the peak category, bold with "
        f"a filled dot where the image on that row was taken after the storm "
        f"(so the row can show its effects), faint italic where the image "
        f"predates it or there is none that year, with a line to the first "
        f"image taken after the storm; e.g. {ex}. The rule is the first image "
        f"on or after the storm's last listed day, with a 7-day grace; the "
        f"May 2022 nor'easter is routed to the Oct 2023 image by hand. "
        f"(d) Number of images showing each domain overwashed. 2004 belongs to "
        f"both periods.")


def write_captions(entries):
    for name, text in entries:
        p = upsert_caption(name, "heatmaps", text)
    print(f"  wrote {p.relative_to(REPO)}")


# ===================================================================== main
def main(argv):
    apply_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    uniform = "--uniform-years" in argv
    tags = [a for a in argv if a in PERIODS] or list(PERIODS)

    obs, domains, matrix = load_observations()
    storms = assign_capture(load_storms(), obs)
    TAB_DIR.mkdir(parents=True, exist_ok=True)
    observations_long(obs, domains, matrix).to_csv(
        TAB_DIR / "overwash_observations.csv", index=False)
    storms_table(storms, obs).to_csv(TAB_DIR / "storms_by_image.csv", index=False)
    print(f"  wrote {TAB_DIR.relative_to(REPO)}/overwash_observations.csv, storms_by_image.csv")

    entries = []
    for tag in tags:
        period = PERIODS[tag]
        suffix = "_uniform" if uniform else ""
        out = FIG_DIR / f"overwash_heatmap_{tag}{suffix}.png"
        rows, placed = make_figure(tag, period, obs, domains, matrix, storms,
                                   collapse=not uniform, out=out)
        entries.append((out.name, caption_text(tag, period, rows, placed, obs, storms)))
    if not uniform:
        write_captions(entries)


if __name__ == "__main__":
    main(sys.argv[1:])
