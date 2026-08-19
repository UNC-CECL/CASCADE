"""Shoreline time-series extraction from a completed CASCADE run.

Pure data extraction -- nothing here touches matplotlib. Feeds both the
shoreline GIF and the rate-comparison figures in cascade_pipeline.plotting.
"""

import numpy as np


def get_x_s_TS(b3d):
    """Extract the shoreline position time series from a Barrier3D object.

    Args:
        b3d: A single-domain Barrier3D instance (element of cascade.barrier3d).

    Returns:
        1-D float array, one value per model year, in the raw x_s_TS
        convention (dam, increasing landward with erosion).

    Raises:
        AttributeError: Neither x_s_TS nor _x_s_TS is present.
    """
    if hasattr(b3d, "x_s_TS"):
        return np.asarray(b3d.x_s_TS, dtype=float)
    if hasattr(b3d, "_x_s_TS"):
        return np.asarray(b3d._x_s_TS, dtype=float)
    raise AttributeError(
        "No shoreline time series found (x_s_TS / _x_s_TS) on Barrier3D object."
    )


def build_shoreline_matrix(cascade, to_meters=True):
    """Build a [time x domain] shoreline matrix from a completed CASCADE run.

    Args:
        cascade: A Cascade instance after cascade.update() has finished.
        to_meters: Convert dam -> m (CASCADE's native unit is dam).

    Returns:
        2-D float array, shape (n_years, n_domains), raw x_s_TS convention
        (increases landward / with erosion). No sign flip is applied here;
        see compute_change_rate for the plotting convention used elsewhere
        in cascade_pipeline.
    """
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt = len(get_x_s_TS(b3d_list[0]))
    shoreline = np.zeros((nt, ndom), dtype=float)
    for j in range(ndom):
        shoreline[:, j] = get_x_s_TS(b3d_list[j])
    if to_meters:
        shoreline *= 10.0  # dam -> m
    return shoreline


def compute_change_rate(shoreline_m, span_years=None, flip_sign=True):
    """Modeled shoreline change rate (m/yr) from a [time x domain] matrix.

    change_rate = (shoreline_m[-1, :] - shoreline_m[0, :]) / span_years

    Args:
        shoreline_m: Array from build_shoreline_matrix(cascade, to_meters=True).
        span_years: Denominator for the rate, e.g. END_YEAR - START_YEAR.
            Defaults to (n_years - 1) if not given.
        flip_sign: CASCADE's x_s_TS increases landward (erosion). True
            (default) flips the sign so a positive rate means
            seaward/accreting -- the convention used throughout
            cascade_pipeline's plotting functions. Pass False to keep the raw
            x_s_TS sign convention.

    Returns:
        1-D float array, length n_domains, m/yr.
    """
    shoreline_m = np.asarray(shoreline_m, dtype=float)
    nt = shoreline_m.shape[0]
    denom = span_years if span_years is not None else max(nt - 1, 1)
    rate = (shoreline_m[-1, :] - shoreline_m[0, :]) / float(denom)
    if flip_sign:
        rate *= -1.0
    return rate
