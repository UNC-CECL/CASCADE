"""
Where the Hatteras GIS drive is: the one setting every script reads it from.

    from site_layer.hat_gis_root import GIS_ROOT, GIS_PROJECT

GIS_ROOT is the raw GIS sources, D:\\hatteras_ch1\\gis_sources, unless the
environment variable HATTERAS_GIS_ROOT names another folder (a second drive
under another letter, say). GIS_PROJECT is the ArcGIS Pro project beside it,
D:\\hatteras_ch1\\gis_project. What is in them: D:\\hatteras_ch1\\README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import os
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
DEFAULT_GIS_ROOT = Path(r"D:\hatteras_ch1\gis_sources")
# -----------------------------------------------------------------------------

GIS_ROOT = Path(os.environ.get("HATTERAS_GIS_ROOT", DEFAULT_GIS_ROOT))
GIS_PROJECT = GIS_ROOT.parent / "gis_project"
