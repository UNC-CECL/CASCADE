"""
Where the Hatteras GIS drive is: the one setting every script reads it from.

    from site_layer.hat_gis_root import GIS_ROOT

D:\\Hatteras_GIS unless the environment variable HATTERAS_GIS_ROOT names another
folder (a second drive under another letter, say). What is in the folder:
D:\\Hatteras_GIS\\README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import os
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
DEFAULT_GIS_ROOT = Path(r"D:\Hatteras_GIS")
# -----------------------------------------------------------------------------

GIS_ROOT = Path(os.environ.get("HATTERAS_GIS_ROOT", DEFAULT_GIS_ROOT))
