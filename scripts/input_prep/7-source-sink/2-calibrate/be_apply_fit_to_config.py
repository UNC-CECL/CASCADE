"""
Write be_zone_residual_fit.py's fitted rates into hatteras_site_config.py.

    python scripts/input_prep/7-source-sink/2-calibrate/be_apply_fit_to_config.py
    python scripts/input_prep/7-source-sink/2-calibrate/be_apply_fit_to_config.py --check
    python scripts/input_prep/7-source-sink/2-calibrate/be_apply_fit_to_config.py --add

Keeps the locked end domains and the other periods, backs the config up
first; --add adds an iteration pass to the field in place. Details: scripts/input_prep/7-source-sink/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

import argparse
import os
import re
import shutil
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve()
_REPO_ROOT = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())

CONFIG = _REPO_ROOT / "scripts" / "site_layer" / "hatteras_site_config.py"
NL = chr(10)

# Resolved the way the generator resolves its OUTPUT_DIR, not typed again

# THE PAIR BEING APPLIED -- read from the SAME variable the generator reads
sys.path.insert(0, str(_REPO_ROOT / "scripts"))
from site_layer import hat_source_sink as _be  # noqa: E402
from site_layer.hatteras_site_config import (HATTERAS_BE_EDGE_D90,     # noqa: E402
                                  HATTERAS_BE_RATES_EDGE,
                                  HATTERAS_PERIODS)

# THE ONE WAY TO SEED A NEW PERIOD'S GIS 90
_GIS90_ENV = os.environ.get("HAT_BE_GIS90", "").strip()
_GIS90_SEED = {}
for _item in filter(None, (x.strip() for x in _GIS90_ENV.split(","))):
    _k, _, _v = _item.partition("=")
    _GIS90_SEED[int(_k)] = float(_v)

DEFAULT_PERIOD_STARTS = (1984, 2004)
_PERIOD_ENV = os.environ.get("HAT_BE_PERIODS", "").strip()
PERIOD_STARTS = (tuple(int(x) for x in _PERIOD_ENV.split(","))
                 if _PERIOD_ENV else DEFAULT_PERIOD_STARTS)
if len(PERIOD_STARTS) != 2:
    raise SystemExit(f"{NL}HAT_BE_PERIODS must name exactly two period "
                     f"starts, got {PERIOD_STARTS!r}.{NL}")
for _st in PERIOD_STARTS:
    if _st not in HATTERAS_PERIODS:
        raise SystemExit(f"{NL}period start {_st} is not in "
                         f"HATTERAS_PERIODS ({sorted(HATTERAS_PERIODS)}).{NL}")
_P1, _P2 = PERIOD_STARTS
_ENDS = {p: HATTERAS_PERIODS[p]["end_year"] for p in PERIOD_STARTS}

# Every pair keeps its products in its own folder, so applying one never reads another's
_PAIR_TAG = _be.pair_tag(_P1, _ENDS[_P1], _P2, _ENDS[_P2])
# --- CONFIG ------------------------------------------------------------------
RATES_TXT = (Path(os.environ.get("HAT_BE_OUTPUT_DIR", "").strip()
                  or _be.calibrate_dir(_PAIR_TAG))
             / "DOMAIN_BE_RATES.txt")
LOCKED_GIS = (1, 90)

# The generator writes this file through the Windows console encoding, not UTF-8
RATES_ENCODING = "cp1252"

_ROW = re.compile(r"^\s*(\d+):\s*([+-]?\d+\.?\d*),\s*(?:#\s*(.*))?$", re.M)
# -----------------------------------------------------------------------------


# Reads the P1 and P2 blocks out of the calibration's output
def parse_generated(path=RATES_TXT):
    if not path.exists():
        raise SystemExit(
            f"No calibration output at {path}{NL}"
            f"Run 2-calibrate/be_zone_residual_fit.py first. If that "
            f"pass set HAT_BE_OUTPUT_DIR, set it here too -- this reads "
            f"whichever directory the generator wrote.")
    text = path.read_text(encoding=RATES_ENCODING)
    out = {}
    for block in re.split(r"^# ", text, flags=re.M):
        head = block.split(NL, 1)[0]
        if head.startswith("P1"):
            period = _P1
        elif head.startswith("P2"):
            period = _P2
        else:
            continue
        body = block[block.index("{"):block.index("}") + 1]
        out[period] = {int(m.group(1)): (float(m.group(2)),
                                         (m.group(3) or "").strip())
                       for m in _ROW.finditer(body)}
    missing = set(PERIOD_STARTS) - set(out)
    if missing:
        raise ValueError(f"{path} has no block for period(s) {sorted(missing)}")
    return out


# Locates the calibrated table in the config source
def existing_block(source):
    start = source.index("HATTERAS_BE_RATES_CALIBRATED = {")
    end = source.index(NL + "}" + NL, start) + len(NL + "}" + NL)
    return start, end, source[start:end]


# Renders one period's dict, keeping locked lines and zone labels
def render(period, generated, locked_lines, old_labels, old_rates, add,
           write=True):
    lines = ["    %d: {" % period]
    for gis in range(1, 91):
        if gis in LOCKED_GIS:
            lines.append(locked_lines[(period, gis)])
            continue
        if not write:
            # Period not selected: re-emit exactly what is there already.
            rate = old_rates[period][gis]
            label = old_labels[period].get(gis, "")
            line = "        %3d: %+.1f," % (gis, rate)
            if label and not label.startswith("LOCKED"):
                line += "  # " + label
            lines.append(line)
            continue
        rate, label = generated[period][gis]
        if add:
            rate += old_rates[period][gis]
        label = label or old_labels[period].get(gis, "")
        line = "        %3d: %+.1f," % (gis, rate)
        if label and not label.startswith("LOCKED"):
            line += "  # " + label
        lines.append(line)
    lines.append("    },")
    return NL.join(lines)


# Run: parse the fit, diff against the config, back up and write
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="report the diff and write nothing")
    parser.add_argument("--periods",
                        default=",".join(str(p) for p in PERIOD_STARTS),
                        help="comma-separated start years to write (default "
                             "both). Use this when one period has converged "
                             "and the other has not: re-writing a converged "
                             "period spends runs re-fitting an irreducible "
                             "residual, and under --add it would keep adding "
                             "to a field that is already done.")
    parser.add_argument("--add", action="store_true",
                        help="ADD the generated values to the ones already in "
                             "the config, instead of replacing them. Use this "
                             "for iteration passes, where the calibration was "
                             "run against a calibBE base "
                             "(HAT_BE_BASE_PRESET=calibBE) and so measured "
                             "what is LEFT rather than the whole correction. "
                             "Replacing there would discard the field the base "
                             "run was carrying and undo every earlier pass.")
    args = parser.parse_args()

    periods = tuple(int(x) for x in args.periods.split(",") if x.strip())
    unknown = [p for p in periods if p not in PERIOD_STARTS]
    if unknown:
        parser.error(f"unknown period(s) {unknown}; this pass covers "
                     f"{PERIOD_STARTS} (set HAT_BE_PERIODS to change it)")

    print(f"reading {RATES_TXT}")
    generated = parse_generated()
    source = CONFIG.read_text(encoding="utf-8")
    start, end, block = existing_block(source)

    # EVERY PERIOD THE CONFIG ALREADY HOLDS
    config_periods = [int(m) for m in
                      re.findall(r"^    (\d{4}): \{", block, flags=re.M)]
    new_periods = [p for p in periods if p not in config_periods]

    locked_lines, old_labels, old_rates = {}, {}, {}
    for period in new_periods:
        # Hard blocker, checked first: a period without a D90 entry breaks the config's import
        if period not in HATTERAS_BE_EDGE_D90:
            raise SystemExit(
                f"{NL}{period} is not in HATTERAS_BE_EDGE_D90 "
                f"({sorted(HATTERAS_BE_EDGE_D90)}).{NL}{NL}"
                f"HATTERAS_BE_RATES_EDGE is built by slicing "
                f"HATTERAS_BE_RATES_CALIBRATED and injecting "
                f"HATTERAS_BE_EDGE_D90[period], so adding {period} to the "
                f"calibrated table without a D90 entry makes "
                f"hatteras_site_config raise KeyError on IMPORT -- breaking "
                f"every script in the project, not just this one.{NL}{NL}"
                f"Solve edgeBE's GIS 90 for {period} and add it to "
                f"HATTERAS_BE_EDGE_D90 first.{NL}")

        # A period with no block yet needs its two END domains, which are solved by buffer-cell reproduction
        edge = HATTERAS_BE_RATES_EDGE.get(period, {})
        # NOT "does edge have a GIS 90" -- it does, for every period
        seeded = _GIS90_SEED.get(period)
        if seeded is None:
            raise SystemExit(
                f"{NL}{period} has no calibrated block yet, and GIS 90 cannot "
                f"be seeded.{NL}{NL}"
                f"GIS 1 is SHARED between the edge and calibrated presets -- "
                f"the config slices it out of the calibrated fit, and the two "
                f"agree at every solved period -- so it can come from "
                f"HATTERAS_BE_RATES_EDGE. GIS 90 is NOT shared: it reads 13.0 "
                f"in edge against 32.8 in calibrated for 1984-2004, because an "
                f"isolated forced cell is diffused by its neighbours and needs "
                f"about ten times the misfit it closes.{NL}{NL}"
                f"Solve GIS 90 for {period} under calibBE, then pass it "
                f"deliberately:{NL}"
                f'    HAT_BE_GIS90="{period}=<value>" python '
                f"2-calibrate/be_apply_fit_to_config.py{NL}")
        old_rates[period] = {g: 0.0 for g in range(1, 91)}
        old_labels[period] = {}
        locked_lines[(period, 1)] = (
            "        %3d: %+.1f,  # LOCKED - shared with the edge preset"
            % (1, edge.get(1, 0.0)))
        locked_lines[(period, 90)] = (
            "        %3d: %+.1f,  # LOCKED - solved separately for calibBE"
            % (90, seeded))
        print(f"  {period}: new block, GIS 1 seeded from the edge preset "
              f"({edge.get(1, 0.0):+.1f})")

    for period in config_periods:
        b0 = block.index("    %d: {" % period)
        b1 = block.index("    },", b0)
        sub = block[b0:b1]
        rows = {int(m.group(1)): (float(m.group(2)), (m.group(3) or "").strip())
                for m in _ROW.finditer(sub)}
        old_rates[period] = {g: v[0] for g, v in rows.items()}
        old_labels[period] = {g: v[1] for g, v in rows.items()}
        for gis in LOCKED_GIS:
            hit = re.search(r"^(\s*%d: [+-]?\d+\.\d+,\s*# LOCKED.*)$" % gis,
                            sub, flags=re.M)
            if hit is None:
                raise ValueError(
                    f"{period}: GIS {gis} is not marked LOCKED in the config. "
                    f"This script refuses to overwrite an end domain it cannot "
                    f"identify as separately solved.")
            locked_lines[(period, gis)] = hit.group(1).rstrip()

    for period in periods:
        final = {g: (old_rates[period][g] + generated[period][g][0] if args.add
                     else generated[period][g][0]) for g in range(2, 90)}
        moved = [(g, old_rates[period][g], final[g])
                 for g in range(2, 90)
                 if abs(old_rates[period][g] - final[g]) > 1e-9]
        deltas = [abs(a - b) for _, a, b in moved]
        print(f"  {period}: {len(moved)} of 88 interior domains move  "
              f"(mean {sum(deltas)/max(len(deltas),1):.3f}, "
              f"max {max(deltas, default=0):.2f} m/yr)")
        for gis in LOCKED_GIS:
            print(f"    GIS {gis:2d} locked, kept: "
                  f"{locked_lines[(period, gis)].strip().split(',')[0]}")

    if args.check:
        print("\n  --check: nothing written")
        return 0

    # The snapshot is this step's OUTPUT -- the BE field as it stood going in
    prebe_dir = _be.PREBE_DIR
    prebe_dir.mkdir(parents=True, exist_ok=True)
    backup = prebe_dir / (
        f"hatteras_site_config_prebe_{time.strftime('%Y%m%d_%H%M%S')}.py")
    shutil.copy2(CONFIG, backup)

    # The union, in period order: periods not selected are re-emitted as they were
    all_periods = sorted(set(config_periods) | set(periods))
    new_block = ("HATTERAS_BE_RATES_CALIBRATED = {" + NL
                 + NL.join(render(p, generated, locked_lines, old_labels,
                                  old_rates, args.add, p in periods)
                           for p in all_periods) + NL
                 + "}" + NL)
    print(f"  block holds {all_periods}; "
          f"{sorted(periods)} updated, "
          f"{sorted(set(config_periods) - set(periods))} untouched")
    CONFIG.write_text(source[:start] + new_block + source[end:],
                      encoding="utf-8")
    print(f"\n  wrote   {CONFIG}")
    print(f"  backup  {backup}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
