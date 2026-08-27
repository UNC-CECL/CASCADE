"""Writes HAT_be_zone_LOESS_analysis.py's fitted rates into the site config.

The calibration prints a DOMAIN_BE_RATES block "ready to paste". Pasting it by
hand loses two things every time, so this does it instead:

  - THE LOCKED END DOMAINS. The generator writes GIS 1 and 90 as 0.0 with a
    "use your solved value" comment, because they are solved separately by
    buffer-cell reproduction rather than fit from a residual. Pasting the block
    verbatim silently zeroes both, and HATTERAS_BE_RATES_EDGE is SLICED from
    this table -- so it would zero the edge preset too.
  - THE ZONE LABEL ON A ZERO-VALUED DOMAIN. The generator omits the trailing
    comment where the correction solved to 0.0, so a domain's physical zone
    disappears from the config purely because its rate came out zero.

Both are preserved here: locked lines are copied verbatim from the config, and
a missing zone label falls back to the one already there.

Usage:
    python scripts/input_prep/7-source-sink/apply_be_fit.py [--check] [--add]
"""

import argparse
import re
import shutil
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve()
_REPO_ROOT = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())

CONFIG = _REPO_ROOT / "scripts" / "hatteras_site_config.py"
RATES_TXT = (_HERE.parent / "loess_smooth" / "output" / "DOMAIN_BE_RATES.txt")
LOCKED_GIS = (1, 90)
NL = chr(10)

# The generator writes this file through the Windows console encoding, not
# UTF-8, so its em-dashes arrive as cp1252 bytes.
RATES_ENCODING = "cp1252"

_ROW = re.compile(r"^\s*(\d+):\s*([+-]?\d+\.?\d*),\s*(?:#\s*(.*))?$", re.M)


def parse_generated(path=RATES_TXT):
    """Reads the P1 and P2 blocks out of the calibration's output.

    Args:
        path: DOMAIN_BE_RATES.txt.

    Returns:
        {period: {gis: (rate, zone_label)}} for 1984 and 2004. The file also
        holds three forecast blocks; those are ignored -- this writes the
        hindcast presets only.
    """
    text = path.read_text(encoding=RATES_ENCODING)
    out = {}
    for block in re.split(r"^# ", text, flags=re.M):
        head = block.split(NL, 1)[0]
        if head.startswith("P1"):
            period = 1984
        elif head.startswith("P2"):
            period = 2004
        else:
            continue
        body = block[block.index("{"):block.index("}") + 1]
        out[period] = {int(m.group(1)): (float(m.group(2)),
                                         (m.group(3) or "").strip())
                       for m in _ROW.finditer(body)}
    missing = {1984, 2004} - set(out)
    if missing:
        raise ValueError(f"{path} has no block for period(s) {sorted(missing)}")
    return out


def existing_block(source):
    """Locates the calibrated table in the config source.

    Returns:
        (start, end, block_text) character offsets into `source`.
    """
    start = source.index("HATTERAS_BE_RATES_CALIBRATED = {")
    end = source.index(NL + "}" + NL, start) + len(NL + "}" + NL)
    return start, end, source[start:end]


def render(period, generated, locked_lines, old_labels, old_rates, add,
           write=True):
    """Renders one period's dict, keeping locked lines and zone labels.

    Args:
        add: If True, the generated value is a DELTA measured against a base
            run that already carries `old_rates`, so it is added rather than
            replacing. See --add in main().
    """
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true",
                        help="report the diff and write nothing")
    parser.add_argument("--periods", default="1984,2004",
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
    unknown = [p for p in periods if p not in (1984, 2004)]
    if unknown:
        parser.error(f"unknown period(s) {unknown}; expected 1984 and/or 2004")

    generated = parse_generated()
    source = CONFIG.read_text(encoding="utf-8")
    start, end, block = existing_block(source)

    locked_lines, old_labels, old_rates = {}, {}, {}
    for period in (1984, 2004):
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

    for period in (1984, 2004):
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

    backup = CONFIG.with_name(
        f"hatteras_site_config_prebe_{time.strftime('%Y%m%d_%H%M%S')}.py")
    shutil.copy2(CONFIG, backup)

    new_block = ("HATTERAS_BE_RATES_CALIBRATED = {" + NL
                 + render(1984, generated, locked_lines, old_labels, old_rates,
                          args.add, 1984 in periods) + NL
                 + render(2004, generated, locked_lines, old_labels, old_rates,
                          args.add, 2004 in periods) + NL
                 + "}" + NL)
    CONFIG.write_text(source[:start] + new_block + source[end:],
                      encoding="utf-8")
    print(f"\n  wrote   {CONFIG}")
    print(f"  backup  {backup.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
