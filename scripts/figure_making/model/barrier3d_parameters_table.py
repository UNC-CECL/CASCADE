"""
A table of the Barrier3D parameter yaml: the CASCADE default beside the value each hindcast ran with.

    python scripts/figure_making/model/barrier3d_parameters_table.py

Reads the run's own copy of Hatteras-CASCADE-parameters.yaml for each period
(CASCADE rewrites it per domain from the Cascade() arguments, so the copy, not
the live file, holds what ran), then lays the current settings over it: Dmaxel
from the live yaml, rmin/rmax and beta from the runner. Defaults come from the CASCADE template
data/default_cascade_variables/barrier3d-default-parameters.yaml. Elevations are
shown in m NAVD88 and m MHW. Writes
output/figures/4-model-mechanics/barrier3d/barrier3d_parameters_table.png.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import yaml  # noqa: E402

REPO = next(_p for _p in Path(__file__).resolve().parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1997, INK, INK_MUTED, apply_style, caption, figsize, figure_dir, save)

# --- CONFIG ------------------------------------------------------------------
DEFAULT_YAML = REPO / "data" / "default_cascade_variables" / "barrier3d-default-parameters.yaml"
RUN_ROOT = REPO / "output" / "raw_runs" / "experiments" / "code-checks" / "2026-10-20-groin-cleanup-check" / "new"
RUN_YAMLS = {   # period label -> that run's copy of the parameter yaml
    "1996–2009": RUN_ROOT / "1996_2009" / "edgeBE" / "HAT_1996_2009_edgeBE_offsetmetres_road_bdm_groin"
                 / "HAT_1996_2009_edgeBE_offsetmetres_road_bdm_groin-parameters.yaml",
    "2009–2025": RUN_ROOT / "2009_2025" / "edgeBE" / "HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_groin"
                 / "HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_groin-parameters.yaml",
}
LIVE_YAML = REPO / "data" / "hatteras_init" / "Hatteras-CASCADE-parameters.yaml"
RUNNER = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
OUT = figure_dir("mechanics", "barrier3d") / "barrier3d_parameters_table"
# Grouped rows: (key, meaning, kind). kind picks the formatting:
# elev = m NAVD88 converted to m MHW too; bay = depth below MHW; else a unit string
GROUPS = (
    ("Elevations", (
        ("MHW", "Mean high water", "mhw"),
        ("BermEl", "Berm elevation", "elev"),
        ("Dmaxel", "Maximum dune elevation", "elev"),
        ("BayDepth", "Bay floor (depth below MHW)", "bay"),
        ("Sandbag_elevation", "Sandbag height (sandbags off)", "m"),
    )),
    ("Dunes and beach", (
        ("rmin", "Minimum dune growth rate", "/yr"),
        ("rmax", "Maximum dune growth rate", "/yr"),
        ("DuneParamStart", "Dune heights from file", ""),
        ("GrowthParamStart", "Growth rates from file", ""),
        ("beta", "Beach slope for run-up", "–"),
    )),
    ("Shoreface (from BRIE)", (
        ("DShoreface", "Shoreface height", "m"),
        ("LShoreface", "Shoreface length", "m"),
        ("ShorefaceToe", "Shoreface toe position", "m"),
        ("s_sf_eq", "Equilibrium shoreface slope", "–"),
        ("k_sf", "Shoreface flux constant", "m³/m/yr"),
    )),
    ("Run and forcing", (
        ("TMAX", "Model years + 1", "yr"),
        ("BarrierLength", "Domain length alongshore", "m"),
        ("RSLR_const", "Relative sea-level rise", "m/yr"),
        ("RSLR_Constant", "Constant sea-level rise", ""),
        ("Rat", "Background alongshore retreat", "m/yr"),
        ("Shrub_ON", "Shrubs", ""),
        ("enable_sandbags", "Sandbags", ""),
        ("storm_file", "Storm series", "file"),
        ("elevation_file", "Interior topography", "file"),
        ("dune_file", "Dune heights", "file"),
    )),
)
SET_IN_YAML = {"Dmaxel"}      # the only key CASCADE does not overwrite per domain
COL_X = (0.0, 0.205, 0.475, 0.63, 0.815)   # left edge of each column, axes fraction
# -----------------------------------------------------------------------------


# One yaml as a dict
def load(path):
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


# The current settings the run copies predate: Dmaxel from the live yaml, rmin/rmax and beta from the runner
def current_settings():
    src = RUNNER.read_text(encoding="utf-8")
    now = {"Dmaxel": load(LIVE_YAML)["Dmaxel"]}
    for key in ("rmin", "rmax"):
        m = re.search(rf"^{key.upper()} = \[([0-9.]+)\]", src, re.M)
        if m is None:
            raise SystemExit(f"{key.upper()} not found in {RUNNER.name}")
        now[key] = float(m.group(1))
    m = re.search(r"^BEACH_SLOPE = ([0-9.]+)", src, re.M)
    if m is None:
        raise SystemExit(f"BEACH_SLOPE not found in {RUNNER.name}")
    now["beta"] = float(m.group(1))
    return now


# A number in as few digits as it needs
def num(v):
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, int):
        return f"{v:d}"
    a = abs(v)
    if a >= 1000:
        return f"{v:,.0f}"
    if a >= 10:
        return f"{v:.1f}"
    if a >= 0.1 or a == 0:
        return f"{v:.2f}"
    return f"{v:.4f}".rstrip("0")


# One cell: the value in the units the row asks for; elevations in both datums
def cell(prm, key, kind):
    v = prm.get(key)
    if v is None:
        return "—"
    if kind == "mhw":
        return f"{v:.2f} m NAVD88\n0.00 m MHW"
    if kind == "elev":
        return f"{v:.2f} m NAVD88\n{v - prm['MHW']:.2f} m MHW"
    if kind == "bay":
        return f"{-v + prm['MHW']:.2f} m NAVD88\n{-v:.2f} m MHW"
    if kind == "file":
        name = Path(str(v).replace("\\", "/")).name
        if len(name) <= 18:
            return name
        cut = name.rfind("_", 0, len(name) // 2 + 4)
        return name[:cut + 1] + "\n" + name[cut + 1:]
    if key == "Shrub_ON":
        return "on" if v else "off"
    unit = f" {kind}" if kind and kind not in ("–",) else ""
    return num(v) + unit


# Whether the run's value differs from the default, in the units the row shows
def changed(default, run, key, kind):
    d, r = default.get(key), run.get(key)
    if kind == "file":
        return Path(str(d)).name != Path(str(r).replace("\\", "/")).name
    if kind in ("elev", "bay"):
        return cell(default, key, kind) != cell(run, key, kind)
    return d != r


# Draw the table on one bare axes
def draw(ax, default, runs):
    n_rows = sum(len(rows) for _, rows in GROUPS)
    n_lines = n_rows + len(GROUPS) * 0.9 + 1.2
    dy = 1.0 / n_lines
    y = 1.0
    heads = ("Parameter", "Meaning", "CASCADE default", *runs)
    colours = (INK, INK, INK, C_1984, C_1997)
    for x, h, col in zip(COL_X, heads, colours):
        ax.text(x, y - 0.5 * dy, h, weight="bold", color=col, va="center", fontsize=7.5)
    y -= 1.2 * dy
    ax.plot([0, 1], [y, y], color=INK, lw=0.8)
    for group, rows in GROUPS:
        y -= 0.9 * dy
        ax.text(0, y + 0.35 * dy, group, weight="bold", color=C["ACCENT"], va="center", fontsize=7)
        for key, meaning, kind in rows:
            mid = y - 0.5 * dy
            tall = kind in ("elev", "bay", "mhw", "file")
            fs = 6.3 if tall else 6.8
            ax.text(COL_X[0], mid, key + ("*" if key in SET_IN_YAML else ""), va="center",
                    fontsize=6.8, family="monospace", color=INK)
            ax.text(COL_X[1], mid, meaning, va="center", fontsize=6.8, color=INK)
            ax.text(COL_X[2], mid, cell(default, key, kind), va="center", fontsize=fs,
                    color=INK_MUTED, linespacing=1.1)
            for x, prm in zip(COL_X[3:], runs.values()):
                diff = changed(default, prm, key, kind)
                ax.text(x, mid, cell(prm, key, kind), va="center", fontsize=fs, linespacing=1.1,
                        color=INK, weight="bold" if diff else "normal")
            y -= dy
            ax.plot([0, 1], [y, y], color=C["BASE_FILL"], lw=0.4)
    ax.plot([0, 1], [y, y], color=INK, lw=0.8)


# Run: read, draw, caption, save
def main():
    apply_style()
    default = load(DEFAULT_YAML)
    now = current_settings()
    runs = {label: {**load(path), **now} for label, path in RUN_YAMLS.items()}

    fig = plt.figure(figsize=figsize("double", height=8.6))
    ax = fig.add_axes((0.02, 0.02, 0.96, 0.96))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    draw(ax, default, runs)

    run_names = "; ".join(f"{p}: {path.relative_to(REPO).as_posix()}" for p, path in RUN_YAMLS.items())
    caption(fig, (
        "Every key in the Barrier3D parameter yaml, with the CASCADE default (its template "
        "barrier3d-default-parameters.yaml, which differs from Barrier3D's built-in defaults in places, e.g. Rat) beside the value each hindcast period ran with. Bold marks a "
        "value that differs from the default. Elevations are given in m NAVD88 and m MHW; the default "
        f"converts with its own MHW ({default['MHW']:.2f} m NAVD88), the Hatteras runs with "
        f"{runs[next(iter(runs))]['MHW']:.2f} m NAVD88 (NOAA Duck, 8651370). Barrier3D works internally "
        "in decametres above MHW and makes this conversion itself. CASCADE rewrites every key but Dmaxel "
        "(*) from its constructor arguments for each domain before the run starts, so the values are read "
        "from each run's own copy of the yaml, not the live data/hatteras_init/Hatteras-CASCADE-parameters.yaml; "
        "the per-domain keys (ShorefaceToe, LShoreface, file names) show the last domain written, a buffer "
        "domain. Sandbags are off, and Sandbag_elevation is shown as written, without a datum conversion. "
        "Keys not in the yaml take Barrier3D's built-in defaults (dune erosion, overwash, shrub and "
        f"vegetation constants; MODEL_VARIABLES.md). Runs read: {run_names}. Those runs predate the "
        f"2026-10-20 dune and beach-slope settings, so Dmaxel ({now['Dmaxel']} m NAVD88) is taken from the live yaml and "
        f"rmin/rmax ({now['rmin']}/{now['rmax']}) and beta ({now['beta']}, BEACH_SLOPE) from the runner, "
        "HAT_hindcast_1984_2024.py."))
    save(fig, OUT, vector=True, close=True)
    print(OUT.with_suffix(".png"))


if __name__ == "__main__":
    main()
