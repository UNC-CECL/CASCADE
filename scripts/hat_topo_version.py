# ==============================================================================
# hat_topo_version.py
#
# Which dune/topography run do the road scripts measure against?
#
# WHY THIS EXISTS
#   Four scripts here used to hardcode ".../2009-dune-topo/2009_v3/topography":
#
#       1-produce/HAT_road_placement_on_domains.py
#       2-audit/HAT_road_setback_audit.py
#       3-figures/HAT_road_domain_views.py
#       4-compare/HAT_road_method_diagnostic.py
#
#   HAT_road_offset_from_dune_start.py does NOT -- it takes the path off the
#   extractor, and says why: "Paths come off the extractor rather than being
#   hardcoded, so bumping VERSION does not silently point this at stale arrays."
#   The other four did not follow that, so when the dune windows were re-picked
#   into 2009_v4 (2026-08-19) they kept reading v3 interiors while consuming v4
#   setbacks. Nothing errored. 18 domains had different interiors, 10 of them a
#   different SHAPE -- including D79 and D80, two of the three roadways the
#   relocation logic acts on, and D11/D85/D86 among the floored negatives. Every
#   drown verdict and placement number for those was computed on the wrong grid.
#
#   So the version is resolved ONCE, here, from the extractor that produced the
#   arrays. Bumping VERSION in the extractor now moves the whole road tree with
#   it, and a version that does not exist on disk is an immediate, loud error
#   rather than a silently stale read.
#
# WHY IT PARSES RATHER THAN IMPORTS
#   Importing HAT_dune_topo_extractor.py pulls in matplotlib and a windowing
#   backend, which is a heavy and fragile dependency for an audit that otherwise
#   needs neither. Only two literals are needed, so they are read out of the
#   source. HAT_road_offset_from_dune_start.py still imports the module properly,
#   because it needs the actual functions.
#
# USAGE
#     from hat_topo_version import topo_dirs
#     TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs()
#
#   To pin a road script to a specific run instead (e.g. to reproduce an old
#   figure), pass it explicitly and the extractor is not consulted:
#     TOPO_DIR, DUNE_DIR, RUN_NAME = topo_dirs(override="2009_v3")
# ==============================================================================

from __future__ import annotations

import re
from pathlib import Path

_HERE = Path(__file__).resolve()
# scripts/hat_topo_version.py -> repo root. MOVED here 2026-08-20 from
# input_prep/4-mgmt-forcings/road_offset/, because the model runner now
# resolves its topography through this too and a runner importing out of an
# input-prep subfolder is backwards. It sits next to hatteras_site_config.py:
# both answer "what does this site use", for every consumer.
PROJECT_ROOT = _HERE.parents[1]
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
DUNE_TOPO_ROOT = INIT_ROOT / "1-barrier3d-domains" / "2009-dune-topo"

# The one that has ALONGSHORE_FLIP = True. Three other copies of this file exist
# in the repo and all are unflipped -- see the note in
# HAT_road_offset_from_dune_start.py.
EXTRACTOR = (PROJECT_ROOT / "scripts" / "input_prep" / "1-barrier3d-domains"
             / "topography_dunes" / "HAT_dune_topo_extractor.py")


def _literal(name: str, text: str) -> str:
    """Read `NAME = "value"` out of the extractor source."""
    m = re.search(rf'^{name}\s*=\s*["\']([^"\']+)["\']', text, re.MULTILINE)
    if not m:
        raise SystemExit(
            f"\nCannot find {name} in:\n    {EXTRACTOR}\n"
            f"hat_topo_version.py reads it to locate the topography. If the "
            f"extractor was restructured, pass override= explicitly instead.\n")
    return m.group(1)


def run_name() -> str:
    """RUN_NAME of the dune/topo run currently produced by the extractor."""
    if not EXTRACTOR.is_file():
        raise SystemExit(f"\nExtractor not found:\n    {EXTRACTOR}\n")
    text = EXTRACTOR.read_text(encoding="utf-8")
    # RUN_NAME = f"{DEM_YEAR}_{VERSION}" in the extractor; rebuilt rather than
    # parsed, because it is an f-string there and not a plain literal.
    return f"{_literal('DEM_YEAR', text)}_{_literal('VERSION', text)}"


def topo_dirs(override: str | None = None) -> tuple[Path, Path, str]:
    """(topography dir, dunes dir, run name), checked to exist.

    Raises with the list of runs actually on disk rather than returning a path
    that will later read as "no arrays for this domain".
    """
    name = override or run_name()
    run = DUNE_TOPO_ROOT / name
    topo, dune = run / "topography", run / "dunes"

    if not topo.is_dir():
        avail = (sorted(p.name for p in DUNE_TOPO_ROOT.iterdir() if p.is_dir())
                 if DUNE_TOPO_ROOT.is_dir() else [])
        raise SystemExit(
            f"\nTopography directory does not exist:\n    {topo}\n"
            f"runs present: {avail}\n"
            f"This came from {'override=' + name if override else 'VERSION in ' + EXTRACTOR.name}. "
            f"Re-run the extractor for that version, or pass a different "
            f"override to topo_dirs().\n")
    return topo, dune, name
