"""
Where the Hatteras code and the Hatteras workspace are: the two roots every script reads.

    from hatteras_paths import HATTERAS_DIR, WORKSPACE, BARRIER3D_DIR, workspace_twin

HATTERAS_DIR is this folder, scripts/real-world-cases/Hatteras/ in the
obx-cascade clone: the scripts, and nothing else (the lab rule for the clone).
WORKSPACE is everything the scripts read and write -- data/, output/,
hard-structures/ results -- laid out at the paths the old repository had them,
D:\\hatteras_ch1\\model_runs, unless the environment variable HATTERAS_WORKSPACE
names another folder (a second drive under another letter, say).

BARRIER3D_DIR is the Barrier3D checkout beside the obx-cascade clone (../Barrier3D,
on hatteras/adopted); its worktrees (Barrier3D-overwashfix, ...) sit beside it, so
BARRIER3D_DIR.with_name("Barrier3D-overwashfix") names one. In the old repository
they sat beside the workspace, which scripts reached as WORKSPACE.parent.

workspace_twin(folder) is where a folder of scripts keeps its non-code files:
the same place in the workspace. In the old repository a study wrote its
figures and tables beside its script; here the script stays in Hatteras/ and
they go to the twin. Hatteras/hard-structures/<x> pairs with
WORKSPACE/hard-structures/<x>, and every other Hatteras/<x> with
WORKSPACE/scripts/<x>, where the old repository's scripts/ folder kept them.

A script finds this file by walking up from its own location to the folder
that holds hatteras_paths.py, puts that folder on sys.path, and imports it.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""
from __future__ import annotations

import os
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
DEFAULT_WORKSPACE = Path(r"D:\hatteras_ch1\model_runs")
# -----------------------------------------------------------------------------

HATTERAS_DIR = Path(__file__).resolve().parent
WORKSPACE = Path(os.environ.get("HATTERAS_WORKSPACE", DEFAULT_WORKSPACE))
BARRIER3D_DIR = next(p for p in HATTERAS_DIR.parents if (p / ".git").exists()).parent / "Barrier3D"



def workspace_twin(folder) -> Path:
    """The workspace folder that holds what used to sit beside a Hatteras/ script folder.

    A path already in the workspace is returned as it is, so a name built
    from a twin can be passed again.
    """
    folder = Path(folder)
    if folder.is_relative_to(WORKSPACE):
        return folder
    rel = folder.resolve().relative_to(HATTERAS_DIR)
    if rel.parts and rel.parts[0] == "hard-structures":
        return WORKSPACE / rel
    return WORKSPACE / "scripts" / rel


if not (WORKSPACE / "data" / "hatteras_init").is_dir():
    raise RuntimeError(
        f"Hatteras workspace not found: {WORKSPACE} has no data/hatteras_init. "
        f"Plug in the D: drive, or set HATTERAS_WORKSPACE to the folder that "
        f"holds data/ and output/.")
