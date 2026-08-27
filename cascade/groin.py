
"""Simplified groin representation for the CASCADE barrier-island model.

This module provides a minimal parameterization of a groin for use with CASCADE. It is designed
to be attached to a :class:`Cascade` instance via ``cascade._groin_callback`` and
invoked once per model year, immediately before BRIE's alongshore-transport solve.

Overview
--------
A groin blocks alongshore sediment transport, accreting the updrift beach and
starving the downdrift beach. This module represents that effect with the smallest
possible signal: a source/sink pair straddling the groin boundary. Each year it
adds ``-M`` to the updrift domain (seaward advance) and ``+M`` to the immediately
downdrift domain (landward retreat), where the two domains are adjacent and share
the blocked boundary.

Because BRIE's alongshore transport is a diffusion solver, this localized pair does
not remain localized: BRIE spreads it on its own, growing the updrift fillet and
downdrift notch through time. The fillet's extent, taper, and profile are therefore
*emergent* -- they are not prescribed. The only free parameter is the amplitude M.

Design intent
-------------
The parameterization is deliberately minimal: no decay, no bypass/leakage, and no
multi-era schedule. A single tunable knob (``trapping_rate_m_yr``) keeps the model
falsifiable rather than over-fit, in keeping with the project's preference for
physically motivated forms with amplitude-only tuning.

Why one knob is a genuine test (not circular)
---------------------------------------------
For BRIE's implicit diffusion with a ``+M/-M`` dipole, the steady-state updrift
offset scales with amplitude,

    A  ~=  M / (4 * r_ipl),        r_ipl = D * dt / (2 * dy**2),

while the alongshore extent of the fillet grows with time alone,

    L  ~=  dy * sqrt(2 * r_ipl * t),    (M does not appear).

So M sets amplitude only; extent is fixed by diffusivity. One tunes M to match the
observed updrift amplitude -- guaranteed achievable, since it is a single linear
knob -- and then checks whether the emergent extent matches observations *for
free*. That independent extent check is the actual scientific test.

Conventions (verified against ``brie_coupler``)
-----------------------------------------------
Units
    ``x_s_dt`` is passed to BRIE in metres (``brie_coupler`` converts back to
    decametres via ``/10`` on return). M is therefore in metres per step directly;
    no decametre conversion happens in this module.
Sign
    BRIE / Barrier3D ``x_s`` increases landward. A negative ``x_s_dt`` is a seaward
    advance (accretion); a positive ``x_s_dt`` is a landward retreat (erosion). The
    :data:`ACCRETION` and :data:`EROSION` constants encode this.
Indexing
    The runner resolves 1-based GIS domain IDs (D1..D90) to padded array indices
    via its own ``gis_to_pad`` translator and passes those indices in, so the
    padding width ``N_PAD`` is defined in exactly one place (the runner), never
    here.

Author
------
Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

from typing import Dict, List

import numpy as np

__all__ = ["GroinCallback", "predict_fillet", "ACCRETION", "EROSION"]


# ---------------------------------------------------------------------------
# Sign convention (BRIE / Barrier3D: x_s increases landward)
# ---------------------------------------------------------------------------
ACCRETION: float = -1.0   # seaward advance  -> x_s decreases
EROSION:   float = +1.0   # landward retreat -> x_s increases


# ===========================================================================
# Groin callback
# ===========================================================================
class GroinCallback:
    """Per-year source/sink injection representing a single groin.

    An instance is attached to a CASCADE run via ``cascade._groin_callback`` and
    called once per model year from inside ``Cascade.update()``, immediately before
    the alongshore-transport solve, with the signature::

        x_s_dt = callback(cascade, x_s_dt)

    Each active year the callback adds ``-M`` to the updrift domain and ``+M`` to
    the downdrift domain of the ``x_s_dt`` array, then returns it. Before the
    install year it returns ``x_s_dt`` unchanged. Per-year diagnostics are recorded
    for the scientific record.

    Parameters
    ----------
    updrift_pad, downdrift_pad : int
        Padded array indices of the two domains flanking the groin. They must be
        adjacent (differ by exactly one) because they share the blocked boundary.
        The runner resolves these from GIS IDs, so ``N_PAD`` lives only there.
    trapping_rate_m_yr : float
        M, the shoreline change applied per year at each flank (metres). This is
        the single tunable amplitude.
    start_year : int
        Model year of the run's first update step (e.g. 1967).
    install_year : int
        The groin is inert until ``year >= install_year``.
    n_domains : int
        Padded domain count; used for bounds checks and diagnostic sizing.

    Attributes
    ----------
    year_TS : list of int
        Model year recorded at each call.
    active_TS : list of bool
        Whether the groin was active (installed) at each call.
    applied_dx_updrift_TS, applied_dx_downdrift_TS : list of float
        Signed shoreline change (metres) applied updrift / downdrift each year.

    Raises
    ------
    ValueError
        If the two domains are not adjacent, or if either index is out of range
        for ``n_domains``.
    """

    def __init__(
        self,
        updrift_pad: int,
        downdrift_pad: int,
        trapping_rate_m_yr: float,
        start_year: int,
        install_year: int,
        n_domains: int,
        deterioration_delay_years: float = None,
        deterioration_mode: str = "instant",
        deterioration_fraction: float = 1.0,
        deterioration_ramp_years: float = 0.0,
        sink_fraction: float = 1.0,
    ) -> None:
        if abs(updrift_pad - downdrift_pad) != 1:
            raise ValueError(
                "groin domains must be adjacent (share the blocked boundary); "
                f"got updrift_pad={updrift_pad}, downdrift_pad={downdrift_pad}."
            )
        if not (0 <= updrift_pad < n_domains and 0 <= downdrift_pad < n_domains):
            raise ValueError(
                f"groin pads out of range for n_domains={n_domains}: "
                f"updrift_pad={updrift_pad}, downdrift_pad={downdrift_pad}."
            )

        # Configuration
        self.updrift_pad: int = int(updrift_pad)
        self.downdrift_pad: int = int(downdrift_pad)
        self.M: float = float(trapping_rate_m_yr)
        self.start_year: int = int(start_year)
        self.install_year: int = int(install_year)
        self.n_domains: int = int(n_domains)

        # Deterioration (optional): onset is specified as a delay relative to
        # install_year (structure age), so the same config generalizes across
        # runs with different install years. Two explicit, mutually exclusive
        # modes -- no silent parameter-collapsing between them:
        #   "instant"     -> M steps to M*deterioration_fraction the year the
        #                    delay elapses (e.g. a specific storm/abandonment).
        #   "linear_ramp" -> M declines linearly to M*deterioration_fraction
        #                    over deterioration_ramp_years (gradual bypass
        #                    increase / slow structural failure).
        if deterioration_mode not in ("instant", "linear_ramp"):
            raise ValueError(
                f"deterioration_mode must be 'instant' or 'linear_ramp', "
                f"got {deterioration_mode!r}."
            )
        if deterioration_mode == "instant" and deterioration_ramp_years != 0.0:
            raise ValueError(
                "deterioration_ramp_years must be 0 when deterioration_mode "
                "is 'instant' -- pass deterioration_mode='linear_ramp' for a "
                "gradual decline instead."
            )
        if deterioration_mode == "linear_ramp" and deterioration_ramp_years <= 0:
            raise ValueError(
                "deterioration_ramp_years must be > 0 when deterioration_mode "
                "is 'linear_ramp'."
            )
        if deterioration_delay_years is not None and deterioration_delay_years < 0:
            raise ValueError("deterioration_delay_years must be >= 0.")
        if not (0.0 <= deterioration_fraction <= 1.0):
            raise ValueError("deterioration_fraction must be in [0, 1].")

        self.deterioration_delay_years = (
            None if deterioration_delay_years is None else float(deterioration_delay_years)
        )
        self.deterioration_year = (
            None if self.deterioration_delay_years is None
            else self.install_year + self.deterioration_delay_years
        )
        self.deterioration_mode: str = deterioration_mode
        self.deterioration_fraction: float = float(deterioration_fraction)
        self.deterioration_ramp_years: float = float(deterioration_ramp_years)

        # ASYMMETRY. The original pair is strictly volume neutral: every metre
        # of advance imposed updrift is a metre of retreat imposed downdrift.
        # At Buxton that is falsified. The emergent extent measured on the
        # 2026-08-24 sweep runs 2,000 m updrift (against 2,250 m observed --
        # a good match) but 2,500 m DOWNDRIFT, where the observed extent is
        # ZERO: the real structure accretes updrift without a measurable
        # deficit on the other side.
        #
        # `sink_fraction` scales the downdrift sink only. 1.0 keeps the
        # volume-neutral pair, so every run predating this is reproducible.
        # 0.0 makes the groin a pure source, which is the physical statement
        # that the trapped sand arrives from OUTSIDE the pair -- at this site,
        # from Cape Point -- rather than from the downdrift beach.
        #
        # It is a MODEL FORM, not a fitted knob: the observed downdrift extent
        # of zero picks 0.0 directly. Sweeping it would be fitting a parameter
        # the data already determines.
        if not (0.0 <= sink_fraction <= 1.0):
            raise ValueError("sink_fraction must be in [0, 1].")
        self.sink_fraction: float = float(sink_fraction)

        # Per-year diagnostics (appended on each call): the scientific record of
        # how much shoreline change the groin imposed, where, and when.
        self.year_TS: List[int] = []
        self.active_TS: List[bool] = []
        self.applied_dx_updrift_TS: List[float] = []
        self.applied_dx_downdrift_TS: List[float] = []
        self.trapping_rate_applied_TS: List[float] = []
        self._call_count: int = 0

    # -- main hook ----------------------------------------------------------
    def __call__(self, cascade, x_s_dt):
        """Inject the source/sink into ``x_s_dt`` for the current model year.

        Parameters
        ----------
        cascade : Cascade
            The calling CASCADE instance (unused here; present for API symmetry
            and to allow future state-aware behaviour).
        x_s_dt : sequence of float
            Per-domain shoreline change increment (metres) for this step, prior to
            the alongshore-transport solve. Modified in place and returned.

        Returns
        -------
        x_s_dt : sequence of float
            The same object, with the groin source/sink applied when active.
        """
        year = self.start_year + self._call_count
        self._call_count += 1

        active = year >= self.install_year
        M_eff = self._effective_trapping_rate(year) if active else 0.0
        dx_updrift = ACCRETION * M_eff
        dx_downdrift = EROSION * M_eff * self.sink_fraction

        if active:
            x_s_dt[self.updrift_pad] += dx_updrift      # source: updrift accretes
            x_s_dt[self.downdrift_pad] += dx_downdrift  # sink:   downdrift erodes

        self.year_TS.append(year)
        self.active_TS.append(bool(active))
        self.applied_dx_updrift_TS.append(dx_updrift)
        self.applied_dx_downdrift_TS.append(dx_downdrift)
        self.trapping_rate_applied_TS.append(M_eff)

        return x_s_dt

    def _effective_trapping_rate(self, year: int) -> float:
        """Return M for ``year``, applying post-deterioration decline.

        Full M before ``deterioration_year`` (or always, if no deterioration
        was configured). From ``deterioration_year`` onward:
          "instant"     -- steps immediately to M * deterioration_fraction.
          "linear_ramp" -- ramps linearly to M * deterioration_fraction over
                           deterioration_ramp_years, then holds.
        """
        if self.deterioration_year is None or year < self.deterioration_year:
            return self.M

        floor = self.M * self.deterioration_fraction
        if self.deterioration_mode == "instant":
            return floor

        years_since = year - self.deterioration_year
        taper = min(1.0, years_since / self.deterioration_ramp_years)
        return self.M - taper * (self.M - floor)

    # -- diagnostics --------------------------------------------------------
    def summary(self) -> Dict[str, object]:
        """Return a compact run-metadata dictionary for logging.

        Returns
        -------
        dict
            Configuration plus aggregate totals (years active and cumulative
            shoreline change applied updrift and downdrift).
        """
        return dict(
            updrift_pad=self.updrift_pad,
            downdrift_pad=self.downdrift_pad,
            trapping_rate_m_yr=self.M,
            install_year=self.install_year,
            deterioration_delay_years=self.deterioration_delay_years,
            deterioration_year=self.deterioration_year,
            deterioration_mode=self.deterioration_mode,
            deterioration_fraction=self.deterioration_fraction,
            deterioration_ramp_years=self.deterioration_ramp_years,
            start_year=self.start_year,
            years_active=int(np.sum(self.active_TS)),
            cumulative_updrift_m=float(np.sum(self.applied_dx_updrift_TS)),
            cumulative_downdrift_m=float(np.sum(self.applied_dx_downdrift_TS)),
        )

    def diagnostics_frame(self) -> List[Dict[str, object]]:
        """Return the per-year record as a list of dictionaries.

        Each row carries the model year, active flag, the signed change applied
        updrift and downdrift that year, and the running cumulative totals. The
        runner writes these rows to CSV alongside each run.

        Returns
        -------
        list of dict
            One dictionary per recorded model year.
        """
        cum_up = np.cumsum(self.applied_dx_updrift_TS) if self.applied_dx_updrift_TS else []
        cum_down = np.cumsum(self.applied_dx_downdrift_TS) if self.applied_dx_downdrift_TS else []

        rows: List[Dict[str, object]] = []
        for i, year in enumerate(self.year_TS):
            rows.append(dict(
                model_year=year,
                groin_active=self.active_TS[i],
                trapping_rate_applied_m_yr=self.trapping_rate_applied_TS[i],
                applied_dx_updrift_m=self.applied_dx_updrift_TS[i],
                applied_dx_downdrift_m=self.applied_dx_downdrift_TS[i],
                cumulative_updrift_m=float(cum_up[i]),
                cumulative_downdrift_m=float(cum_down[i]),
            ))
        return rows


# ===========================================================================
# On-paper prediction (run before any model run)
# ===========================================================================
def predict_fillet(
    trapping_rate_m_yr: float,
    r_ipl: float,
    run_years: float,
    dy_m: float = 500.0,
):
    """Predict the fillet amplitude and extent analytically, before running CASCADE.

    Implements the scaling relations in the module docstring: amplitude depends on
    the trapping rate M, while extent depends only on the diffusion number and
    elapsed time. Comparing the predicted extent against the emergent model extent
    is the module's independent test.

    Parameters
    ----------
    trapping_rate_m_yr : float
        M, the amplitude knob (metres per step).
    r_ipl : float
        BRIE's dimensionless diffusion number at the groin face,
        ``r_ipl = D * dt / (2 * dy**2)``. Read it from a base run or estimate D.
    run_years : float
        Elapsed model years over which the fillet develops.
    dy_m : float, optional
        Alongshore domain width in metres (default 500.0), used to convert the
        extent from domains to metres.

    Returns
    -------
    amplitude_m : float
        Predicted steady-state updrift offset (metres).
    extent_domains : float
        Predicted alongshore extent of the fillet (number of domains).
    extent_m : float
        Predicted extent in metres (``extent_domains * dy_m``).

    Raises
    ------
    ValueError
        If ``r_ipl`` is not positive.
    """
    r_ipl = float(r_ipl)
    if r_ipl <= 0:
        raise ValueError("r_ipl must be > 0")

    amplitude_m = trapping_rate_m_yr / (4.0 * r_ipl)
    extent_domains = np.sqrt(2.0 * r_ipl * run_years)
    return amplitude_m, extent_domains, extent_domains * dy_m


# ===========================================================================
# Illustrative self-test (run directly: `python -m cascade.groin`)
# ===========================================================================
def _demo() -> None:
    """Print two short worked examples. Not used when the module is imported.

    1. A minimal install-only example (no deterioration), as before.
    2. The real Buxton Groin Field timeline: installed 1969, last repaired
       1996, damaged by a major storm in 2003 -- modeled as a linear ramp
       bridging the two documented dates, printed across 1990-2010 so the
       decline is visible.
    """
    class _FakeCascade:
        """Placeholder standing in for a real Cascade instance in the demo."""

    fake = _FakeCascade()

    # -- Example 1: install only, no deterioration ---------------------------
    callback = GroinCallback(
        updrift_pad=20,
        downdrift_pad=19,
        trapping_rate_m_yr=40.0,
        start_year=1967,
        install_year=1970,
        n_domains=45,   # 15 real + 15 + 15 padding
    )
    print("Example 1: install only (no deterioration)")
    print("year  active  x_s_dt[D6=20]  x_s_dt[D5=19]")
    for _ in range(6):  # 1967..1972, crossing the 1970 install date
        x_s_dt = [0.0] * callback.n_domains
        x_s_dt = callback(fake, x_s_dt)
        year = callback.year_TS[-1]
        print(f"{year}  {callback.active_TS[-1]!s:<6}  "
              f"{x_s_dt[20]:+12.1f}  {x_s_dt[19]:+12.1f}")
    print("\nsummary:", callback.summary())

    amplitude, extent_domains, extent_m = predict_fillet(
        trapping_rate_m_yr=40.0, r_ipl=0.2, run_years=30,
    )
    print(f"\npredicted updrift amplitude = {amplitude:.1f} m")
    print(f"predicted fillet extent     = {extent_domains:.1f} domains "
          f"({extent_m:.0f} m)")

    # -- Example 2: real Buxton Groin Field deterioration timeline -----------
    install_year = 1969
    last_repair_year = 1996
    storm_year = 2003
    groin_buxton = GroinCallback(
        updrift_pad=20,
        downdrift_pad=19,
        trapping_rate_m_yr=60.0,
        start_year=1984,
        install_year=install_year,
        n_domains=120,
        deterioration_delay_years=last_repair_year - install_year,   # 27
        deterioration_mode="linear_ramp",
        deterioration_ramp_years=storm_year - last_repair_year,      # 7
        deterioration_fraction=0.2,
    )
    print("\n\nExample 2: Buxton Groin Field (install 1969, "
          f"deterioration {last_repair_year}->{storm_year})")
    print("year  M_applied  x_s_dt[updrift]  x_s_dt[downdrift]")
    x_s_dt = [0.0] * groin_buxton.n_domains
    for _ in range(1990 - 1984 + 21):  # print through 2010
        x_s_dt = [0.0] * groin_buxton.n_domains
        x_s_dt = groin_buxton(fake, x_s_dt)
        year = groin_buxton.year_TS[-1]
        if 1990 <= year <= 2010:
            m_eff = groin_buxton.trapping_rate_applied_TS[-1]
            print(f"{year}  {m_eff:7.2f}   {x_s_dt[20]:+8.2f}         {x_s_dt[19]:+8.2f}")
    print("\nsummary:", groin_buxton.summary())


if __name__ == "__main__":
    _demo()
