#!/usr/bin/env python3
"""Does the dune line ever reach NC-12? Original vs observed vs modelled.

WHY THIS FIGURE EXISTS

    The 1984-2004 relocation comparison reports that almost no domain relocates
    unaided once the road setbacks are measured against `1984-start` row 0. The
    obvious objection is that the roads are still only 3-6 cells behind the
    dune, so something must be wrong. This figure is the check: it puts the
    three cross-shore positions on one axis, per domain, so the claim can be
    read off rather than argued.

        original   the 1984 dune line -- the run's own year-0 position, and the
                   datum every other quantity here is measured from
        observed   the surveyed 2004 dune line, from
                   2-brie-offset/raw_offsets/2004_duneline_offset_raw.csv minus
                   the 1984 file. This is the SAME target the run scores its
                   misfit against (cascade_pipeline.hindcast.
                   build_shoreline_target), not a second opinion.
        modelled   where the run put the dune line in 2004

    and draws NC-12 as the 20 m band it occupies, at the setback the model was
    initialised with.

WHAT "CORRECT BEHAVIOUR" LOOKS LIKE HERE

    `road_relocation_checks` relocates when the setback goes STRICTLY negative,
    and the setback is driven by

        dune_migration = barrier3d.ShorelineChangeTS[t-1] * 10       # m

    -- whole 10 m cells, because ShorelineChangeTS counts cells. So the road is
    overrun only when the dune line travels PAST the near edge of the road band.
    A domain whose observed and modelled 2004 dune lines both stop short of that
    edge SHOULD not relocate, and a model that agrees with the survey about
    where the dune line got to is behaving correctly even though it misses the
    historical relocation. That is the distinction this figure is drawn to make
    visible: a miss caused by the TRIGGER (geometric overrun) is not the same as
    a miss caused by the PHYSICS (dune line in the wrong place).

SIGN CONVENTION
    +x is LANDWARD throughout, matching x_s_TS and the raw offset files. The
    1984 dune line is 0 by construction on every domain.

USAGE
    python scripts/hatteras_ms/HAT_relocation_dune_position_check.py
    python scripts/hatteras_ms/HAT_relocation_dune_position_check.py --preset calibBE

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _p in (SCRIPTS_DIR, _HERE.parent):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from cascade_pipeline.hindcast import load_absolute_dune_distance  # noqa: E402
from hatteras_site_config import (                                 # noqa: E402
    HATTERAS_DOMAINS,
    HATTERAS_ROAD_EVENTS,
    resolve_be_preset,
)
from cascade_pipeline.roadway import RelocationEvent               # noqa: E402

START_YEAR, END_YEAR = 1984, 2004
RAW_OFFSET_DIR = PROJECT_BASE_DIR / "data" / "hatteras_init" / "2-brie-offset" / "raw_offsets"
RUN_ROOT = PROJECT_BASE_DIR / "output" / "raw_runs" / "1984_2004"
OUT_ROOT = PROJECT_BASE_DIR / "output" / "comparisons" / "relocation_1984_2004"

ROAD_WIDTH_M = 20.0          # roadway_manager default, and what the runs use

# Palette. Deliberately not a rainbow: the three positions are one family
# (where is the dune line) and the road is the thing they are compared against,
# so the road is the only warm colour on the figure.
C_ORIGINAL = "#1f4e79"       # 1984 dune line, the datum
C_OBSERVED = "#2E7D32"       # surveyed 2004
C_MODELLED = "#B36AE2"       # modelled 2004
C_ROAD = "#C1440E"
C_RELOC = "#C1440E"


def arm_a_dir(preset):
    """The free-running arm for a preset. Same name rule as the comparison."""
    return RUN_ROOT / preset / f"HAT_{START_YEAR}_{END_YEAR}_{preset}_road_bdm_nogroin"


def load_run(run_dir):
    """(cascade, shoreline matrix) for a finished run."""
    npz = sorted(glob.glob(os.path.join(str(run_dir), "*.npz")))
    if not npz:
        raise SystemExit(
            f"\nno .npz in {run_dir}\n"
            f"This figure needs the saved model state -- re-run with "
            f"HAT_SAVE_MODEL_STATE=1.\n")
    with np.load(npz[0], allow_pickle=True) as h:
        cascade = h["cascade"][0]
    mat = sorted(glob.glob(os.path.join(str(run_dir), "*_shoreline_matrix.npy")))
    return cascade, (np.load(mat[0]) if mat else None)


def historical_targets():
    """{gis: event year} for the relocations inside this period."""
    return {g: e.year for e in HATTERAS_ROAD_EVENTS
            if isinstance(e, RelocationEvent) and e.enabled
            and START_YEAR <= e.year <= END_YEAR
            for g in e.displacement_m}


def collect(preset):
    """Per-domain positions, all metres landward of the 1984 dune line."""
    cascade, shoreline = load_run(arm_a_dir(preset))

    # Observed: difference two ABSOLUTE surveyed distances. The padded offset
    # files each subtract their own year's minimum, so differencing THOSE is
    # not a shoreline change -- this is the same call build_shoreline_target
    # makes for the run's own misfit line.
    d0 = load_absolute_dune_distance(START_YEAR, HATTERAS_DOMAINS, RAW_OFFSET_DIR)
    d1 = load_absolute_dune_distance(END_YEAR, HATTERAS_DOMAINS, RAW_OFFSET_DIR)
    observed = d1 - d0                                    # + = landward

    real = slice(HATTERAS_DOMAINS.start_real_index, HATTERAS_DOMAINS.end_real_index)
    modelled = ((shoreline[-1] - shoreline[0])[real]
                if shoreline is not None else np.full(90, np.nan))

    roadways = cascade.roadways
    mgmt = cascade.roadway_management_module
    rows = []
    for gis in range(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1):
        pad = HATTERAS_DOMAINS.gis_to_pad(gis)
        managed = (roadways[pad] is not None
                   and (mgmt is None or bool(mgmt[pad])))
        if not managed:
            continue
        m = roadways[pad]
        sb = np.asarray(m._road_setback_TS, dtype=float)
        # `margin_m` mirrors HAT_relocation_comparison.relocation_margin():
        # the trigger is `setback < 0` STRICT and the setback moves in whole
        # 10 m cells, so the extra migration needed is the closest approach
        # plus one cell. Meaningless where the road did relocate -- the reset
        # in _apply_relocation puts an artificial minimum in the series.
        closest = float(sb.min())
        relocations = int(np.asarray(m._road_relocated_TS).sum())
        rows.append(dict(
            gis=gis,
            setback0=float(sb[0]),
            setback_end=float(sb[-1]),
            setback_ts=sb,
            relocations=relocations,
            min_setback_m=None if relocations else closest,
            margin_m=None if relocations else closest + 10.0,
            observed=float(observed[gis - 1]),
            modelled=float(modelled[gis - 1]),
        ))
    return rows


def draw(rows, preset, out_path):
    targets = historical_targets()
    gis = np.array([r["gis"] for r in rows])
    sb0 = np.array([r["setback0"] for r in rows])
    obs = np.array([r["observed"] for r in rows])
    mod = np.array([r["modelled"] for r in rows])
    reloc = np.array([r["relocations"] > 0 for r in rows])

    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=(15.5, 10.2), height_ratios=[1.45, 1],
        constrained_layout=True)

    # ---- panel A: every managed domain -------------------------------------
    # The road band is drawn as a bar from its near (seaward) edge landward,
    # because that near edge is the thing the dune line has to cross.
    ax.bar(gis, ROAD_WIDTH_M, bottom=sb0, width=0.85, color=C_ROAD,
           alpha=0.85, zorder=2, label=f"NC-12 (20 m wide, at its 1984 setback)")
    ax.axhline(0, color=C_ORIGINAL, lw=2.0, zorder=3)
    ax.plot(gis, obs, "o", ms=4.5, color=C_OBSERVED, zorder=4,
            label=f"observed {END_YEAR} dune line (survey)")
    ax.plot(gis, mod, "^", ms=4.5, color=C_MODELLED, zorder=5,
            label=f"modelled {END_YEAR} dune line ({preset})")

    for g in sorted(targets):
        ax.axvline(g, color="0.75", lw=4.0, alpha=0.35, zorder=0)
    for g in gis[reloc]:
        ax.plot(g, sb0[list(gis).index(g)] + ROAD_WIDTH_M + 6, "v",
                ms=7, color=C_RELOC, zorder=6)

    ax.set_ylabel(f"metres landward of the {START_YEAR} dune line", fontsize=10)
    ax.set_xlabel("GIS domain   (S | Cape Hatteras  →  Pea Island | N)", fontsize=10)
    ax.set_title(
        f"Does the dune line reach NC-12?   {START_YEAR}–{END_YEAR}, "
        f"{preset}, full management, groin off",
        fontsize=13, fontweight="bold", loc="left")
    ax.set_xlim(gis.min() - 1, gis.max() + 1)
    ax.grid(axis="y", alpha=0.25)

    handles = [
        Line2D([], [], color=C_ORIGINAL, lw=2.0,
               label=f"{START_YEAR} dune line (datum, = 0)"),
        Line2D([], [], marker="o", ls="none", color=C_OBSERVED,
               label=f"observed {END_YEAR} dune line (survey)"),
        Line2D([], [], marker="^", ls="none", color=C_MODELLED,
               label=f"modelled {END_YEAR} dune line"),
        Patch(facecolor=C_ROAD, alpha=0.85, label="NC-12, 20 m wide"),
        Line2D([], [], marker="v", ls="none", color=C_RELOC,
               label="relocated at least once"),
        Patch(facecolor="0.75", alpha=0.35, label="historical relocation domain"),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8.5, ncol=3,
              framealpha=0.95)

    sel_gis = set(targets)
    o10 = np.array([r["observed"] for r in rows if r["gis"] in sel_gis])
    m10 = np.array([r["modelled"] for r in rows if r["gis"] in sel_gis])
    ceiling = int((o10 > np.array([r["setback0"] for r in rows
                                   if r["gis"] in sel_gis])).sum())
    note = (f"The road is overrun only where a marker sits ABOVE the bottom of its red bar.\n"
            f"The SURVEYED dune line reaches the road in only {ceiling} of the "
            f"{len(o10)} historical domains (grey) -- that is the ceiling on recall "
            f"from geometric overrun alone.\n"
            f"Here the model's dune line at those domains sits "
            f"{(m10 - o10).mean():+.0f} m against the survey.")
    ax.text(0.995, 0.03, note, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=8.5, color="0.25",
            bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="0.8", alpha=0.92))

    # ---- panel B: the ten scored domains, zoomed ---------------------------
    sel = [r for r in rows if r["gis"] in targets]
    x = np.arange(len(sel))
    s0 = np.array([r["setback0"] for r in sel])
    ob = np.array([r["observed"] for r in sel])
    md = np.array([r["modelled"] for r in sel])
    lab = [f"{r['gis']}\n{targets[r['gis']]}" for r in sel]

    bx.bar(x, ROAD_WIDTH_M, bottom=s0, width=0.55, color=C_ROAD, alpha=0.85,
           zorder=2)
    bx.axhline(0, color=C_ORIGINAL, lw=2.0, zorder=3)
    bx.plot(x, ob, "o", ms=9, color=C_OBSERVED, zorder=4)
    bx.plot(x, md, "^", ms=9, color=C_MODELLED, zorder=5)
    # Both labels go ABOVE the road bar. An earlier version put the relocation
    # count below the lowest marker, where it collided with the tick labels on
    # exactly the two domains that relocate.
    #
    # For a domain that never fired, the useful number is not "it did not
    # relocate" but HOW CLOSE it came: the extra landward dune migration that
    # would have fired the trigger. The setback moves in whole 10 m cells and
    # the test is `< 0` strict, so a road sitting at setback 0 still needs one
    # more full cell -- hence min_setback + 10, not min_setback.
    for i, r in enumerate(sel):
        tag = f"{r['setback0']:.0f} m"
        if r["relocations"]:
            tag += f"   relocated ×{r['relocations']}"
        else:
            tag += f"\nneeded +{r['margin_m']:.0f} m"
        bx.annotate(tag, (i, r["setback0"] + ROAD_WIDTH_M),
                    textcoords="offset points", xytext=(0, 6), ha="center",
                    fontsize=8, color=C_ROAD, fontweight="bold")
    bx.set_xticks(x, lab, fontsize=9)
    bx.set_xlabel("GIS domain / historical relocation year", fontsize=10)
    bx.set_ylabel(f"metres landward of the {START_YEAR} dune line", fontsize=10)
    bx.set_title("The ten domains history relocated", fontsize=11,
                 fontweight="bold", loc="left")
    bx.grid(axis="y", alpha=0.25)

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return sel


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", default="zeroBE")
    args = ap.parse_args()
    preset, _ = resolve_be_preset(args.preset)

    out_dir = OUT_ROOT / "dune_position_check"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"dune_position_vs_road_{preset}.png"

    rows = collect(preset)
    sel = draw(rows, preset, out_path)

    print("=" * 86)
    print(f"DUNE LINE vs NC-12  --  {preset}, {START_YEAR}-{END_YEAR}, "
          f"metres landward of the {START_YEAR} dune line")
    print("=" * 86)
    print(f"{'GIS':>4} {'road near edge':>15} {'observed':>10} {'modelled':>10} "
          f"{'obs-mod':>8} {'overrun?':>20} {'relocs':>7} {'needed':>9}")
    for r in sel:
        edge = r["setback0"]
        obs_over = r["observed"] > edge
        mod_over = r["modelled"] > edge
        verdict = ("both" if obs_over and mod_over else
                   "observed only" if obs_over else
                   "modelled only" if mod_over else "neither")
        need = "-" if r["margin_m"] is None else f"+{r['margin_m']:.0f} m"
        print(f"{r['gis']:>4} {edge:>15.0f} {r['observed']:>10.1f} "
              f"{r['modelled']:>10.1f} {r['observed'] - r['modelled']:>8.1f} "
              f"{verdict:>20} {r['relocations']:>7} {need:>9}")

    allr = np.array([(r["observed"], r["modelled"]) for r in rows])
    ok = ~np.isnan(allr).any(axis=1)
    d = allr[ok, 0] - allr[ok, 1]
    print(f"\n  all {ok.sum()} managed domains: observed mean "
          f"{allr[ok, 0].mean():+.1f} m, modelled mean {allr[ok, 1].mean():+.1f} m")
    print(f"  observed - modelled: mean {d.mean():+.1f} m, "
          f"RMSE {np.sqrt((d ** 2).mean()):.1f} m")
    n_obs = sum(1 for r in rows if r["observed"] > r["setback0"])
    n_mod = sum(1 for r in rows if r["modelled"] > r["setback0"])
    print(f"  domains where the dune line passes the road's near edge by {END_YEAR}:"
          f"  observed {n_obs}/{len(rows)},  modelled {n_mod}/{len(rows)}")

    # ---- the decomposition this figure was drawn to produce ----------------
    # Splitting the misses by CAUSE is the whole point. An island-wide misfit
    # near zero hides it: the model tracks the survey on average and still
    # under-predicts badly on exactly the domains under test.
    obs_over = [r for r in sel if r["observed"] > r["setback0"]]
    obs_short = [r for r in sel if r["observed"] <= r["setback0"]]
    got = [r for r in obs_over if r["relocations"]]
    missed = [r for r in obs_over if not r["relocations"]]

    print("\n" + "=" * 86)
    print("WHY EACH HISTORICAL DOMAIN IS OR IS NOT RECOVERED")
    print("=" * 86)
    print(f"  the SURVEYED dune line overruns the road in "
          f"{len(obs_over)}/{len(sel)} domains: GIS {[r['gis'] for r in obs_over]}")
    print(f"    -> that is the CEILING on recall from geometric overrun alone.")
    print(f"       Even a model that reproduced the dune line exactly could not")
    print(f"       relocate the other {len(obs_short)}: GIS {[r['gis'] for r in obs_short]}")
    print(f"       -- history moved NC-12 there for reasons this trigger does not")
    print(f"       represent (storm damage, burial, maintenance cost).")
    print(f"  of those {len(obs_over)}, the model relocates "
          f"{len(got)}: GIS {[r['gis'] for r in got]}")
    if missed:
        print(f"  and MISSES {len(missed)}, by under-predicting dune retreat:")
        for r in missed:
            short = r["observed"] - r["modelled"]
            need = r["setback0"] - r["modelled"]
            print(f"    GIS {r['gis']:>3}  observed {r['observed']:+.1f} m, "
                  f"modelled {r['modelled']:+.1f} m  ({short:.0f} m short of the "
                  f"survey; {need:+.0f} m from firing)")

    o = np.array([r["observed"] for r in sel])
    m = np.array([r["modelled"] for r in sel])
    local = float((m - o).mean())
    island = float(d.mean() * -1)          # d is observed-modelled; flip to modelled-observed
    print(f"\n  at these 10 domains: observed {o.mean():+.1f} m, "
          f"modelled {m.mean():+.1f} m, misfit {local:+.1f} m")
    print(f"  across all {ok.sum()} managed domains the misfit is {island:+.1f} m.")
    if abs(local) > abs(island) + 5:
        print(f"  -> the error is CONCENTRATED on the domains under test: the "
              f"island-wide\n     figure understates it by "
              f"{abs(local) - abs(island):.0f} m, so recall here is limited by the "
              f"DUNE LINE,\n     not only by the trigger.")
    else:
        print(f"  -> the domains under test are no worse than the island as a "
              f"whole, so what\n     limits recall here is the TRIGGER's "
              f"definition, not the modelled dune line.")
    print(f"\n  RECALL AGAINST THE ACHIEVABLE CEILING: "
          f"{len(got)}/{len(obs_over)}   "
          f"(against all 10 historical domains: {len(got)}/{len(sel)})")
    print(f"\nfigure -> {out_path}")


if __name__ == "__main__":
    main()
