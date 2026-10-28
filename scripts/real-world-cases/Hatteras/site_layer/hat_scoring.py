"""
How a Hatteras run is compared with CoastSat: the target, the scored interior and the model's smoothing.

    from site_layer.hat_scoring import coastsat_target, interior, smooth_like_target

Moved unchanged on 2026-10-23 out of hatteras_ms/experiments/offset/HAT_metres_1_offset_units.py, where
the matrix analysis, the hindcast figures and the groin studies imported them from as `common`.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""

import numpy as np

# The group's smoothing window, in domains (2026-09-28)
SMOOTH_DOMAINS = 7


# The CoastSat LRR target for a window, by GIS domain, as the runner builds it
def coastsat_target(start=1996, end=None, lowess_domains=None):
    lowess_domains = lowess_domains or SMOOTH_DOMAINS
    from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_PERIODS
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.coastsat_lowess import (CoastSatDataset, LowessConfig,
                                                 build_coastsat_series)
    # `end` scores a window inside a period; its CoastSat table must exist
    window = f"{start}_{end or HATTERAS_PERIODS[start]['end_year']}"
    ds = CoastSatDataset(label=f"CoastSat LRR ({window.replace('_', '-')})",
                         period_start=start,
                         csv_path=str(COASTSAT_LRR_ROOT / window / "transect_lrr_full.csv"))
    cfg = LowessConfig(window_domains=(lowess_domains,), skip_southern_domains=10)
    cs = build_coastsat_series([ds], active_period_start=start, lowess_config=cfg,
                               domains=HATTERAS_DOMAINS)[0]
    return build_target_table(cs, cfg, HATTERAS_DOMAINS, lowess_domains).set_index(
        "gis_domain")["target_lrr_m_yr"]


# The scored interior, GIS 2-89: the two locked end domains left out
def interior(series):
    from site_layer.hatteras_site_config import SCORE_INTERIOR_GIS
    lo, hi = SCORE_INTERIOR_GIS
    return series.loc[lo:hi]


# A model rate by GIS domain smoothed as the CoastSat target is: 7-domain LOWESS, GIS 1-10 left raw
def smooth_like_target(series):
    import pandas as pd
    from statsmodels.nonparametric.smoothers_lowess import lowess
    from cascade_pipeline.coastsat_lowess import DEFAULT_LOWESS
    x = series.index.to_numpy(dtype=float)
    y = series.to_numpy(dtype=float)
    ok = np.isfinite(y)
    sm = lowess(y[ok], x[ok], frac=SMOOTH_DOMAINS / len(x), return_sorted=False)
    out = pd.Series(np.nan, index=series.index)
    out[ok] = sm
    raw = series.index <= DEFAULT_LOWESS.skip_southern_domains
    out[raw] = series[raw]
    return out
