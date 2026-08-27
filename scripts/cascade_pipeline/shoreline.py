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


def compute_lrr(shoreline_m, span_years=None, flip_sign=True):
    """Modeled shoreline linear regression rate (m/yr), with fit quality.

    The ordinary-least-squares slope through EVERY annual state, which is the
    estimator the observational target is defined by: CoastSat's
    transect_lrr_full.csv holds a per-transect OLS slope fit through the full
    set of satellite shoreline positions in the period (`lrr_m_yr`, with
    `r_squared` and `unc_m_yr` beside it). compute_change_rate's endpoint
    difference is a different quantity -- a net displacement divided by a
    span -- so plotting one against the other compares two estimators rather
    than a model against an observation.

    The two also differ in what they are sensitive to. An endpoint difference
    reads only years 0 and N, so any single-year excursion still present in
    the final state lands in the result at full amplitude. That is not
    hypothetical here: a nourishment fill enters as an instantaneous step in
    x_s, and BRIE's Crank-Nicolson alongshore solve answers a step with a
    grid-scale (2*dy) mode whose amplification factor is (1-4r)/(1+4r) for
    r = D*dt/(2*dy^2). At Hatteras r ~= 1.05, so that factor is about -0.61:
    the mode alternates sign along the coast and decays only ~39% per year.
    A fill late in a run therefore leaves a sawtooth in the final state, and
    an endpoint rate reports all of it. The OLS slope spreads the same
    excursion over every year and reports about a quarter of it.

    r_squared is returned rather than left implicit because a slope is only a
    summary of a trend where a trend is what the domain has. A nourished
    domain sits flat for most of the period and then steps, which fits a line
    badly (r_squared near 0) however the slope is computed -- that is a fact
    about the trajectory worth carrying next to the number, not a defect in
    the estimator. Note also that r_squared is a variance ratio, so it goes
    small wherever the total signal is small, quite apart from linearity.

    Args:
        shoreline_m: Array from build_shoreline_matrix(cascade, to_meters=True),
            shape (n_states, n_domains).
        span_years: Elapsed years the states span, e.g. END_YEAR - START_YEAR.
            Used only to space the regressor, so the slope comes out per
            calendar year. Defaults to (n_states - 1), which is the same
            thing whenever the states are annual.
        flip_sign: CASCADE's x_s_TS increases landward (erosion). True
            (default) flips the sign so a positive rate means
            seaward/accreting -- matching compute_change_rate and the
            convention used throughout cascade_pipeline's plotting.

    Returns:
        An (lrr, r_squared) tuple of 1-D float arrays, length n_domains. lrr
        is in m/yr; r_squared is NaN for a domain that never moved, where the
        ratio is 0/0 rather than a perfect fit.
    """
    shoreline_m = np.asarray(shoreline_m, dtype=float)
    n_states = shoreline_m.shape[0]
    span = span_years if span_years is not None else max(n_states - 1, 1)

    # Evenly spaced in calendar years, so the slope is per year rather than
    # per state. The two agree for annual states; being explicit means a
    # sub-annual save spacing would not silently rescale the rate.
    years = np.linspace(0.0, float(span), n_states)

    slope, intercept = np.polyfit(years, shoreline_m, 1)
    fitted = years[:, None] * slope[None, :] + intercept[None, :]

    total = ((shoreline_m - shoreline_m.mean(axis=0)) ** 2).sum(axis=0)
    residual = ((shoreline_m - fitted) ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r_squared = np.where(total > 0, 1.0 - residual / total, np.nan)

    if flip_sign:
        slope = -slope
    return slope, r_squared
