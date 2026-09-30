#!/usr/bin/env python3
"""
The management rules the hindcast applies, as a table: a manuscript version and a slide version.

    python scripts/figure_making/management/management_rules_table_figure.py

Rows are built from hatteras_site_config (nourishment projects, road-rule
exceptions); writes rules_table.png and rules_table_slide.png to
output/figures/3-model-inputs/4-management/. Details: scripts/figure_making/management/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# House style (site_layer/hat_figure_style.py), applied at import
import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(next(_q for _q in _P(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED,
                              FIG_W_DOUBLE, save, record_caption)
apply_style()

from site_layer.hatteras_site_config import (HATTERAS_NOURISHMENT_PROJECTS,
                                  HATTERAS_COMMUNITY_ZONES,
                                  HATTERAS_ROAD_EVENTS,
                                  HATTERAS_FIRST_ROAD_DOMAIN,
                                  HATTERAS_DOMAINS,
                                  HATTERAS_ANNOTATIONS)

# Paths resolved from the repo root (ORGANIZATION.md rule 5)
_PATH_REPO = next(_p for _p in _P(__file__).resolve().parents
                  if (_p / "pyproject.toml").exists())
from site_layer import hat_figure_style as _hs  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
OUT_DIR = _hs.figure_dir("inputs", "4-management")

Y0, Y1 = 1984, 2024          # the modelled span, as the timeline figure draws it
EN = "\u2013"                # en dash: ranges only
SPACING = HATTERAS_DOMAINS.domain_spacing_m
# -----------------------------------------------------------------------------


# A domain range as text
def _span(first, last):
    return f"D{first}" if first == last else f"D{first}{EN}{last}"


# A year range as text
def _years(first, last):
    return str(first) if first == last else f"{first}{EN}{last}"


# The only editorial text: each project's description, keyed by its config name
_FILL_PROSE = {
    "Rodanthe emergency fill": (
        "Rodanthe",
        "Emergency fill at the Mirlo Beach S-curves following storm erosion, "
        "placed north of Rodanthe village along the NC-12 corridor"),
    "Buxton shore protection": (
        "Buxton",
        "USACE shore-protection project, extending north from the lighthouse "
        "groin field past the village into the road corridor"),
    "Avon shore protection": (
        "Avon",
        "USACE shore-protection project, Due East Road south to Askins Creek "
        "North Drive, within the Avon community zone"),
}


# One row per configured nourishment project, in order of placement
def nourishment_rows():
    rows = []
    for p in sorted(HATTERAS_NOURISHMENT_PROJECTS,
                    key=lambda q: (q.year, q.gis_domains[0])):
        if not p.enabled:
            continue
        loc, prose = _FILL_PROSE.get(p.name, (p.name, p.note))
        per_m = p.volume_m3_per_m(SPACING)
        rows.append((
            loc,
            _span(p.gis_domains[0], p.gis_domains[-1]),
            _years(p.year, p.year),
            f"{prose}. Reported total {p.volume_cubic_yards:,.0f} yd³, spread "
            f"evenly over the footprint at {per_m:.0f} m\u00b3/m."))
    return rows


# The domains where the relocate-or-abandon rule does not run, and why
def road_rows():
    rows = [(
        "Cape Point",
        _span(1, HATTERAS_FIRST_ROAD_DOMAIN - 1),
        _years(Y0, Y1),
        "No NC-12 in the modelled span; the cape terminus carries no "
        "through-road infrastructure to manage.")]

    # Village footprints: the road is a maintained street network, so the rule is off
    _label = {tuple(v): k for k, v in HATTERAS_ANNOTATIONS.town_spans.items()}
    _long = {"Tri-Village": "Tri-Village\n(Salvo / Waves / Rodanthe)"}
    for lo, hi in HATTERAS_COMMUNITY_ZONES:
        name = _label.get((lo, hi), f"GIS {lo}{EN}{hi}")
        rows.append((
            _long.get(name, name),
            _span(lo, hi),
            _years(Y0, Y1),
            "Permanent community zone: the road is a maintained street "
            "network, so landward relocation does not describe it."))

    for e in HATTERAS_ROAD_EVENTS:
        if type(e).__name__ != "BridgeEvent":
            continue
        d = e.gis_domains
        rows.append((
            "N. Rodanthe /\nPea Island NWR",
            _span(d[0], d[-1]),
            _years(e.year, Y1),
            f"Jug Handle Bridge ({e.year}): NC-12 is carried off the barrier "
            f"surface, so these domains are unmanaged from {e.year} onward."))
    return rows


SECTIONS = [
    ("Beach nourishment", C["ADDED"],
     "sediment added to the named domains in the named year",
     ["Location", "Domains", "Year", "Placement and volume"],
     nourishment_rows()),
    ("Road relocation disabled", C["ROAD"],
     "NC-12 held at its starting position",
     ["Location", "Domains", "Years", "Why the rule is off"],
     road_rows()),
]


# width_in(text, size, bold, italic) -> rendered width in inches, from a scratch canvas
def _measurer(dpi=100):
    fig = plt.figure(figsize=(1, 1), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    renderer = fig.canvas.get_renderer()
    cache = {}

    def width_in(text, fs, bold=False, italic=False):
        key = (text, fs, bold, italic)
        if key not in cache:
            t = ax.text(0, 0, text, fontsize=fs,
                        fontweight="bold" if bold else "normal",
                        fontstyle="italic" if italic else "normal")
            cache[key] = t.get_window_extent(renderer=renderer).width / dpi
            t.remove()
        return cache[key]

    return width_in, fig


# Greedy wrap to `avail` inches, keeping explicit newlines
def wrap(width_in, text, avail, fs, **kw):
    lines = []
    for para in str(text).split("\n"):
        words = para.split()
        if not words:
            lines.append("")
            continue
        line = words[0]
        for w in words[1:]:
            if width_in(f"{line} {w}", fs, **kw) <= avail:
                line = f"{line} {w}"
            else:
                lines.append(line)
                line = w
        lines.append(line)
    return lines


# Draw one variant ('paper' or 'slide'); returns the saved paths
def build(variant):
    slide = variant == "slide"

    FW = 10.0 if slide else FIG_W_DOUBLE
    FS_BODY = 10.0 if slide else 7.2
    FS_HEAD = 10.0 if slide else 7.2
    FS_SEC = 11.0 if slide else 7.8
    FS_TITLE = 13.5
    FS_NOTE = 8.5 if slide else 7.0

    ML = MR = 0.30 if slide else 0.04
    TL, TR = ML, FW - MR
    TW = TR - TL

    # Location, Domains and Years sized to their widest entry; Description gets the rest
    width_in, scratch = _measurer()
    GUT = 0.26 if slide else 0.18
    fixed = []
    for i in range(3):
        w = max(width_in(h[i], FS_HEAD, bold=True) for _, _, _, h, _ in SECTIONS)
        for _, _, _, _, rows in SECTIONS:
            for r in rows:
                for ln in str(r[i]).split("\n"):
                    w = max(w, width_in(ln, FS_BODY))
        fixed.append(w + GUT)
    CW = fixed + [TW - sum(fixed)]
    CX = [TL + sum(CW[:i]) for i in range(4)]

    LINE = FS_BODY * 1.32 / 72.0             # baseline-to-baseline, inches
    PAD_V = 0.075 if slide else 0.050        # above and below a cell's text
    PAD_SEC = 0.090 if slide else 0.058

    # Lay the rows out: wrap every cell, row height from the tallest
    blocks = []
    for name, accent, gloss, heads, rows in SECTIONS:
        head_cells = [wrap(width_in, h, CW[i] - GUT, FS_HEAD, bold=True)
                      for i, h in enumerate(heads)]
        body = []
        for r in rows:
            cells = [wrap(width_in, r[i], CW[i] - GUT, FS_BODY) for i in range(4)]
            body.append((cells, max(len(c) for c in cells) * LINE + 2 * PAD_V))
        blocks.append(dict(
            name=name, accent=accent, gloss=gloss,
            head=head_cells,
            head_h=max(len(c) for c in head_cells) * LINE + 2 * PAD_V,
            sec_h=FS_SEC * 1.32 / 72.0 + 2 * PAD_SEC,
            body=body))

    table_h = sum(b["sec_h"] + b["head_h"] + sum(h for _, h in b["body"])
                  for b in blocks)
    GAP = 0.14 if slide else 0.10            # air between the two sections
    table_h += GAP * (len(blocks) - 1)

    title_h = (FS_TITLE * 1.32 / 72.0 + 0.26) if slide else 0.0
    note_h = 0.34 if slide else 0.0
    FH = 0.10 + title_h + table_h + note_h + 0.10

    fig = plt.figure(figsize=(FW, FH), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, FW)
    ax.set_ylim(0, FH)
    ax.axis("off")

    def rule(y, lw, color=INK):
        ax.plot([TL, TR], [y, y], color=color, lw=lw, zorder=5,
                solid_capstyle="butt")

    def put(lines, x, y_top, w, fs, ha="left", **kw):
        for k, ln in enumerate(lines):
            tx = (x + GUT / 2 if ha == "left"
                  else x + w - GUT / 2 if ha == "right" else x + w / 2)
            ax.text(tx, y_top - PAD_V - (k + 1) * LINE + 0.26 * LINE, ln,
                    fontsize=fs, ha=ha, va="baseline", zorder=6,
                    clip_on=False, **kw)

    ALIGN = ["left", "center", "center", "left"]

    y = FH - 0.10
    if slide:
        # A slide has no caption beside it, so the title lives on the canvas.
        ax.text(TL, y - FS_TITLE * 1.32 / 72.0 * 0.80,
                "Management rules applied in the Hatteras Island hindcast, "
                f"{Y0}{EN}{Y1}",
                fontsize=FS_TITLE, fontweight="bold", color=INK,
                ha="left", va="baseline", zorder=6)
        y -= title_h

    rule(y, 1.1)                                    # booktabs \toprule

    for bi, b in enumerate(blocks):
        if bi:
            y -= GAP
            rule(y, 0.6)                            # \midrule between sections

        # Section heading: an accent tick in the rule's colour, the name, a one-line gloss
        y_s = y - b["sec_h"]
        if slide:
            ax.add_patch(Rectangle((TL, y_s), TW, b["sec_h"], fc="0.965",
                                   ec="none", zorder=1))
        tick = 0.055
        ax.add_patch(Rectangle((TL, y_s + 0.014), tick, b["sec_h"] - 0.028,
                               fc=b["accent"], ec="none", zorder=4))
        base = y_s + b["sec_h"] / 2 - FS_SEC / 72.0 * 0.34
        ax.text(TL + tick + 0.10, base, b["name"], fontsize=FS_SEC,
                fontweight="bold", color=INK, ha="left", va="baseline", zorder=6)
        ax.text(TL + tick + 0.10 + width_in(b["name"], FS_SEC, bold=True) + 0.13,
                base, f"{EN} {b['gloss']}", fontsize=FS_SEC * 0.92,
                color=INK_MUTED, ha="left", va="baseline", zorder=6)
        y = y_s

        # COLUMN HEADS
        for i, cell in enumerate(b["head"]):
            put(cell, CX[i], y, CW[i], FS_HEAD, ha=ALIGN[i],
                fontweight="bold", color=INK)
        y -= b["head_h"]
        rule(y, 0.5)                                # under the column heads

        # BODY
        for ri, (cells, h) in enumerate(b["body"]):
            for i, cell in enumerate(cells):
                put(cell, CX[i], y, CW[i], FS_BODY, ha=ALIGN[i], color=INK)
            y -= h
            if ri < len(b["body"]) - 1:
                rule(y, 0.3, color="0.86")          # hairline between entries

    rule(y, 1.1)                                    # \bottomrule

    if slide:
        ax.text(TL, y - 0.15,
                f"Inter-village domains (D{HATTERAS_FIRST_ROAD_DOMAIN}{EN}20, "
                f"D32{EN}67) keep road relocation throughout. "
                f"D1 = Cape Point (south), D90 = Oregon Inlet (north).",
                fontsize=FS_NOTE, color=INK_MUTED, fontstyle="italic",
                ha="left", va="top", zorder=6)

    plt.close(scratch)
    stem = "rules_table_slide" if slide else "rules_table"
    return save(fig, OUT_DIR / stem, close=True, bbox_inches="tight",
                pad_inches=0.06 if slide else 0.02, facecolor="white")


CAPTION_PAPER = (
    "The management rules the hindcast applies. Beach nourishment: the three "
    "projects in the record, each placed in the named year across the named "
    "domains, with the reported project total and the volume per unit "
    "alongshore length the model receives after spreading it evenly over the "
    "footprint. Road relocation: the domains where the relocate-or-abandon "
    "rule does not run, so NC-12 is held at its starting position. "
    f"Inter-village domains (D{HATTERAS_FIRST_ROAD_DOMAIN}\u201320, "
    "D32\u201367) retain relocation throughout. D1 is Cape Point at the "
    "south end, D90 Oregon Inlet at the north. The two lists overlap at "
    "D7\u20138, which lies both outside the modelled extent of NC-12 and "
    "inside the Buxton community zone; either condition alone switches the "
    "rule off there. Every domain span, year and volume is read from "
    "hatteras_site_config.py rather than transcribed."
)

CAPTION_SLIDE = (
    "The presentation rendering of rules_table.png: the same rows, set larger, "
    "carrying the title and the note on the canvas because a slide has no "
    "caption beside it. Use rules_table.png in a manuscript."
)

if __name__ == "__main__":
    for variant, caption in (("paper", CAPTION_PAPER), ("slide", CAPTION_SLIDE)):
        paths = build(variant)
        record_caption(paths[0], caption)
        print("Saved: " + str(paths[0]))
