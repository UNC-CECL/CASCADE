"""
Two 1984-start versions drawn against each other in a common frame (default: the two extractions, v1 and v2).

    python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/6-result/HAT_plot_seaward_insert_compare.py
    python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/6-result/HAT_plot_seaward_insert_compare.py --domains 85,86

Profiles and road position per domain; --versions compares any other pair. Details: scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Walk up until a directory holds data/hatteras_init
def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer.hat_topo_version import array_name, dune_topo_root  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
PRODUCT = "1984-start"
YEAR = 1984
LAND_DAM = 0.0
ROAD_WIDTH_M = 20.0

# Drawn thick to thin, since the versions coincide over most of the profile
VARIANTS = [
    ("v1", "v1 (original picks, as extracted)", "0.55", "-", 6.0),
    ("v2", "v2 (re-pick, as extracted; CURRENT)", "#b2182b", "-", 2.0),
]
# -----------------------------------------------------------------------------


# Island cells per alongshore column, BARRIER3D'S definition
def b3d_width(topo: np.ndarray) -> np.ndarray:
    out = np.empty(topo.shape[1])
    for c in range(topo.shape[1]):
        col = topo[:, c]
        first = next((i for i, v in enumerate(col) if v <= LAND_DAM), col.size)
        out[c] = max(first - 1, 0)
    return out


# One version's topography for a domain, with its insert audit
def load(version: str, domain: int):
    root = dune_topo_root(PRODUCT) / version
    topo = np.load(root / "topography" / array_name("topography", domain))
    audit = {}
    ap = root / "HAT_seaward_row_insert_audit.csv"
    if ap.is_file():
        for r in csv.DictReader(open(ap)):
            audit[int(r["domain"])] = r
    return topo, audit.get(domain)


# v1's own raw setbacks, from the road-offset measurement
def baseline_setbacks():
    from site_layer.hat_topo_version import road_setback_dir
    p = road_setback_dir(YEAR) / "RoadOffset_{}_domains.csv".format(YEAR)
    out = {}
    for r in csv.DictReader(open(p)):
        if r["setback_dunestart_m"] not in ("", "nan"):
            out[int(r["domain"])] = float(r["setback_dunestart_m"])
    return out


# Run: the profile comparison for each domain
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domains", default="85,86")
    ap.add_argument("--out", default=None)
    ap.add_argument("--title", default=None,
                    help="suptitle; the default describes the built-in VARIANTS "
                         "and would mislabel a --versions run")
    ap.add_argument("--versions", default=None,
                    help="semicolon list of version:label pairs, "
                         "overriding the built-in VARIANTS")
    args = ap.parse_args()
    domains = [int(x) for x in args.domains.split(",")]
    base = baseline_setbacks()
    global VARIANTS
    if args.versions:
        cols = ["0.55", "#2166ac", "#b2182b", "#1b7837"]
        lws = [6.0, 3.0, 1.6, 1.6]
        VARIANTS = []
        for k, item in enumerate(args.versions.split(";")):
            v, _, lab = item.partition(":")
            VARIANTS.append((v, lab or v, cols[k % 4], "-", lws[k % 4]))

    fig, axes = plt.subplots(len(domains), 3, figsize=(18, 5.0 * len(domains)),
                             squeeze=False,
                             gridspec_kw={"width_ratios": [2.1, 1.25, 1.15]})

    for r, D in enumerate(domains):
        ax_p, ax_w, ax_t = axes[r]
        rows = []

        road_x = np.nan
        n_max = 0
        v1_rows = load("v1", D)[0].shape[0]
        for version, label, colour, ls, lw in VARIANTS:
            try:
                topo, aud = load(version, D)
            except FileNotFoundError:
                print("  [skip] {} has no domain {}".format(version, D))
                continue
            n = int(aud["n_rows_inserted"]) if aud else 0
            # The origin shift is read off the ARRAY, not guessed from the version name
            shift = n if (n and topo.shape[0] == v1_rows + n) else 0
            n_max = max(n_max, shift)
            if aud is not None and road_x != road_x:
                # The road does not move; what changes is where row 0 sits relative to it
                road_x = float(aud["setback_raw_before_m"])

            med = np.median(topo, axis=1) * 10.0            # dam -> m MHW
            x = (np.arange(topo.shape[0]) - shift) * 10.0   # m landward of v1 row 0
            keep = x <= 500
            ax_p.plot(x[keep], med[keep], ls, color=colour, lw=lw, label=label,
                      solid_capstyle="butt")
            if shift:
                ax_p.axvline(-shift * 10.0, color=colour, lw=1.0, ls="--", alpha=.7)

            land = topo > LAND_DAM
            width = b3d_width(topo) * 10.0
            ax_w.plot(np.arange(topo.shape[1]), width, ls, color=colour, lw=lw,
                      label=label, solid_capstyle="butt")

            setb = (base.get(D, np.nan) if version == "v1"
                    else float(aud["setback_model_after_m"]) if aud else np.nan)
            rows.append((label, n, setb,
                         float(np.median(width)),
                         float(np.mean(topo[land]) * 10.0)))

        if road_x == road_x:
            ax_p.add_patch(plt.Rectangle((road_x, -0.85), ROAD_WIDTH_M, 0.5,
                                         color="k", zorder=6))
            ax_p.annotate("NC-12 (fixed)", xy=(road_x + 10, -0.35),
                          xytext=(70, -0.72), fontsize=8.5,
                          arrowprops=dict(arrowstyle="->", lw=1.0))

        ax_p.axvline(0, color="k", lw=1.2)
        ax_p.text(4, 5.15, "v1 interior row 0", fontsize=8, rotation=90, va="top")
        pad_x = max(n_max * 10.0, 30.0)
        ax_p.axvspan(-pad_x, 0, color="0.9", zorder=0)
        ax_p.text(-pad_x + 4, 5.3, "fabricated", fontsize=8, color="0.35")
        ax_p.set_xlim(-pad_x - 15, 500)
        ax_p.axhline(0, color="#3b6ea5", lw=0.8, ls=":")
        ax_p.set_ylim(-1.0, 5.6)
        ax_p.set_xlabel("metres landward of v1 interior row 0")
        ax_p.set_ylabel("median elevation (m MHW)")
        ax_p.set_title("GIS {} | cross-shore section, common frame\n"
                       "black bar = NC-12, fixed; dashed = each variant's row 0"
                       .format(D), fontsize=11)
        ax_p.legend(fontsize=8.5, loc="upper right")
        ax_p.grid(alpha=.3)

        ax_w.set_xlabel("alongshore column")
        ax_w.set_ylabel("land width (m, Barrier3D FindWidths)")
        ax_w.set_title("island width, alongshore", fontsize=11)
        ax_w.grid(alpha=.3)
        ax_w.legend(fontsize=7.5)

        ax_t.axis("off")
        txt = ["GIS {}".format(D), ""]
        txt.append("{:<24}{:>4}{:>8}{:>9}{:>8}".format(
            "variant", "N", "setb", "width", "elev"))
        txt.append("-" * 53)
        for label, n, setb, w, h in rows:
            txt.append("{:<24}{:>4}{:>7.0f}m{:>8.0f}m{:>7.2f}m".format(
                label[:24], n, setb, w, h))
        txt += ["", "N       rows prepended",
                "setb    model-facing setback",
                "width   median land width",
                "elev    mean elevation of land cells"]
        ax_t.text(0.0, 0.98, "\n".join(txt), family="monospace", fontsize=9,
                  va="top", transform=ax_t.transAxes)

    fig.suptitle(args.title or
                 ("Fixing the negative 1984 setback: three variants, same "
                  "setback, different islands\nfabricated rows filled at the "
                  "median of interior rows 1-3 (backdune platform)"),
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out = Path(args.out) if args.out else (
        dune_topo_root(PRODUCT) / "HAT_seaward_insert_compare.png")
    fig.savefig(out, dpi=130)
    print("wrote {}".format(out))


if __name__ == "__main__":
    main()
