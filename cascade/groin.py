
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

__all__ = ["GroinCallback", "BlockingGroinCallback", "predict_fillet",
           "ACCRETION", "EROSION"]


# ---------------------------------------------------------------------------
# Sign convention (BRIE / Barrier3D: x_s increases landward)
# ---------------------------------------------------------------------------
ACCRETION: float = -1.0   # seaward advance  -> x_s decreases
EROSION:   float = +1.0   # landward retreat -> x_s increases


# ===========================================================================
# Deterioration schedule, shared by both groin forms
# ===========================================================================
def _validate_deterioration(deterioration_mode, deterioration_ramp_years,
                            deterioration_delay_years, deterioration_fraction):
    """Raise ValueError on an inconsistent deterioration configuration.

    Two explicit, mutually exclusive modes -- no silent parameter-collapsing
    between them:
      "instant"     -> strength steps to strength * deterioration_fraction the
                       year the delay elapses (e.g. a specific storm).
      "linear_ramp" -> strength declines linearly to that floor over
                       deterioration_ramp_years (gradual structural failure).
    """
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


def _scheduled_strength(value, year, deterioration_year, deterioration_mode,
                        deterioration_fraction, deterioration_ramp_years):
    """Return ``value`` for ``year`` after the deterioration schedule.

    Full value before ``deterioration_year`` (or always, if None). From then
    on: "instant" steps to value * fraction; "linear_ramp" ramps linearly to
    that floor over the ramp years, then holds. The arithmetic is the one
    GroinCallback has always used, so dipole runs are bit-for-bit unchanged.
    """
    if deterioration_year is None or year < deterioration_year:
        return value

    floor = value * deterioration_fraction
    if deterioration_mode == "instant":
        return floor

    years_since = year - deterioration_year
    taper = min(1.0, years_since / deterioration_ramp_years)
    return value - taper * (value - floor)


def _check_pads(updrift_pad, downdrift_pad, n_domains):
    """Raise ValueError unless the two pads are adjacent and in range."""
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
    kind: str = "dipole"

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
        _check_pads(updrift_pad, downdrift_pad, n_domains)

        # Configuration
        self.updrift_pad: int = int(updrift_pad)
        self.downdrift_pad: int = int(downdrift_pad)
        self.M: float = float(trapping_rate_m_yr)
        self.start_year: int = int(start_year)
        self.install_year: int = int(install_year)
        self.n_domains: int = int(n_domains)

        # Deterioration (optional): onset is specified as a delay relative to
        # install_year (structure age), so the same config generalizes across
        # runs with different install years. See _validate_deterioration for
        # the two modes.
        _validate_deterioration(deterioration_mode, deterioration_ramp_years,
                                deterioration_delay_years, deterioration_fraction)

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
        return _scheduled_strength(
            self.M, year, self.deterioration_year, self.deterioration_mode,
            self.deterioration_fraction, self.deterioration_ramp_years)

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
# Blocking groin callback
# ===========================================================================
class BlockingGroinCallback:
    """A groin that blocks a fraction ``b`` of alongshore transport at its face.

    Where :class:`GroinCallback` imposes a fixed ``-M`` / ``+M`` every year,
    this one intercepts transport that actually arrives. Each year it reads
    BRIE's current shoreline, estimates the shoreline change BRIE's
    alongshore solve is about to move across the face between the two
    domains, and cancels a fraction ``b`` of it through ``x_s_dt``:

        x_s_dt[lo] -= b * 2 * r_lo * (x_s[hi] - x_s[lo])
        x_s_dt[hi] -= b * 2 * r_hi * (x_s[lo] - x_s[hi])

    where ``lo < hi`` are the two pads and ``r_i`` is BRIE's diffusion number
    for cell ``i``, computed exactly as ``brie.py`` does (forward-difference
    shoreline angle, ``_coast_diff`` lookup, clipped at zero). BRIE's step is
    Crank-Nicolson and row-scaled, so each cell's coupling to the face is
    ``r_i`` in the explicit half and ``r_i`` in the implicit half; the
    explicit term doubled stands in for both. That is an approximation --
    the implicit half uses next year's shoreline -- and was measured against
    scaling the face coupling inside the solve itself at under ~1 m over 14
    years at Buxton (hard-structures/groin/2-module-tests/3-real-planform/blocking_groin_emulator.py).

    Why this form. ``b`` is a trapping fraction, 0 (no structure) to 1 (a
    wall), comparable with published groin trapping efficiencies; the
    trapped volume is emergent and bounded by the transport arriving, so the
    structure cannot overshoot, run away, or exceed the sediment budget the
    way a fixed M can. It needs no change to BRIE: it uses the same pre-AST
    hook as the dipole.

    The equivalent trapping rate each year (the |shoreline change| applied at
    the updrift flank) is recorded as ``trapping_rate_applied_TS`` so budget
    and comparison with the dipole's M are read after the run.

    Parameters
    ----------
    updrift_pad, downdrift_pad : int
        Padded indices of the two domains sharing the blocked face. Adjacent,
        and not the periodic wrap (0 and n_domains - 1).
    blocking_fraction : float
        b, the fraction of the face's transport intercepted, in [0, 1].
    start_year, install_year, n_domains
        As :class:`GroinCallback`.
    deterioration_delay_years, deterioration_mode, deterioration_fraction,
    deterioration_ramp_years
        As :class:`GroinCallback`, applied to b instead of M.
    """
    kind: str = "blocking"

    def __init__(
        self,
        updrift_pad: int,
        downdrift_pad: int,
        blocking_fraction: float,
        start_year: int,
        install_year: int,
        n_domains: int,
        deterioration_delay_years: float = None,
        deterioration_mode: str = "instant",
        deterioration_fraction: float = 1.0,
        deterioration_ramp_years: float = 0.0,
    ) -> None:
        _check_pads(updrift_pad, downdrift_pad, n_domains)
        if not (0.0 <= blocking_fraction <= 1.0):
            raise ValueError("blocking_fraction must be in [0, 1].")
        _validate_deterioration(deterioration_mode, deterioration_ramp_years,
                                deterioration_delay_years, deterioration_fraction)

        self.updrift_pad: int = int(updrift_pad)
        self.downdrift_pad: int = int(downdrift_pad)
        self.blocking_fraction: float = float(blocking_fraction)
        self.start_year: int = int(start_year)
        self.install_year: int = int(install_year)
        self.n_domains: int = int(n_domains)
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

        self._lo: int = min(self.updrift_pad, self.downdrift_pad)
        self._hi: int = max(self.updrift_pad, self.downdrift_pad)

        # Per-year diagnostics, named as GroinCallback's where they mean the
        # same thing so the runner's CSV writer handles both.
        self.year_TS: List[int] = []
        self.active_TS: List[bool] = []
        self.blocking_applied_TS: List[float] = []
        self.applied_dx_updrift_TS: List[float] = []
        self.applied_dx_downdrift_TS: List[float] = []
        self.trapping_rate_applied_TS: List[float] = []
        self.face_offset_m_TS: List[float] = []
        self.r_ipl_updrift_TS: List[float] = []
        self.r_ipl_downdrift_TS: List[float] = []
        self._call_count: int = 0

    def _effective_trapping_rate(self, year: int) -> float:
        """Return b for ``year`` after the deterioration schedule.

        Named as GroinCallback's so schedule reporting reads either form; for
        this class the value is the blocking fraction, not a rate.
        """
        return _scheduled_strength(
            self.blocking_fraction, year, self.deterioration_year,
            self.deterioration_mode, self.deterioration_fraction,
            self.deterioration_ramp_years)

    @staticmethod
    def _r_ipl(brie, x_s, i):
        """BRIE's diffusion number for cell ``i``, as brie.py computes it."""
        ny = len(x_s)
        theta = 180.0 * np.arctan2(x_s[(i + 1) % ny] - x_s[i], brie._dy) / np.pi
        index = int(np.maximum(1, np.minimum(brie._wave_climl,
                                             np.round(90 - theta).astype(int))))
        return max(0.0, float(brie._coast_diff[index] * brie._dt / 2 / brie._dy ** 2))

    def __call__(self, cascade, x_s_dt):
        """Cancel a fraction b of this year's transport across the face.

        Parameters
        ----------
        cascade : Cascade
            The calling CASCADE instance; BRIE's shoreline and diffusivity are
            read from ``cascade._brie_coupler._brie``.
        x_s_dt : sequence of float
            Per-domain shoreline change (metres) before the alongshore solve.
            Modified in place and returned.

        Returns
        -------
        x_s_dt : sequence of float
        """
        year = self.start_year + self._call_count
        self._call_count += 1

        active = year >= self.install_year
        b = self._effective_trapping_rate(year) if active else 0.0

        brie = cascade._brie_coupler._brie
        x_s = np.asarray(brie.x_s, dtype=float)
        lo, hi = self._lo, self._hi
        r_lo = self._r_ipl(brie, x_s, lo)
        r_hi = self._r_ipl(brie, x_s, hi)
        offset = float(x_s[hi] - x_s[lo])

        dx_lo = -b * 2.0 * r_lo * offset
        dx_hi = b * 2.0 * r_hi * offset
        if active and b > 0.0:
            x_s_dt[lo] += dx_lo
            x_s_dt[hi] += dx_hi
        else:
            dx_lo = dx_hi = 0.0

        dx_up, dx_down = ((dx_lo, dx_hi) if self.updrift_pad == lo
                          else (dx_hi, dx_lo))
        self.year_TS.append(year)
        self.active_TS.append(bool(active))
        self.blocking_applied_TS.append(b)
        self.applied_dx_updrift_TS.append(dx_up)
        self.applied_dx_downdrift_TS.append(dx_down)
        self.trapping_rate_applied_TS.append(abs(dx_up))
        self.face_offset_m_TS.append(offset)
        self.r_ipl_updrift_TS.append(r_hi if self.updrift_pad == hi else r_lo)
        self.r_ipl_downdrift_TS.append(r_lo if self.updrift_pad == hi else r_hi)
        return x_s_dt

    @property
    def mean_trapping_rate_m_yr(self) -> float:
        """Mean equivalent trapping rate over the active years, m/yr."""
        active = [m for m, a in zip(self.trapping_rate_applied_TS, self.active_TS) if a]
        return float(np.mean(active)) if active else float("nan")

    def summary(self) -> Dict[str, object]:
        """Return a compact run-metadata dictionary for logging."""
        return dict(
            kind=self.kind,
            updrift_pad=self.updrift_pad,
            downdrift_pad=self.downdrift_pad,
            blocking_fraction=self.blocking_fraction,
            install_year=self.install_year,
            deterioration_delay_years=self.deterioration_delay_years,
            deterioration_year=self.deterioration_year,
            deterioration_mode=self.deterioration_mode,
            deterioration_fraction=self.deterioration_fraction,
            deterioration_ramp_years=self.deterioration_ramp_years,
            start_year=self.start_year,
            years_active=int(np.sum(self.active_TS)),
            mean_trapping_rate_m_yr=self.mean_trapping_rate_m_yr,
            cumulative_updrift_m=float(np.sum(self.applied_dx_updrift_TS)),
            cumulative_downdrift_m=float(np.sum(self.applied_dx_downdrift_TS)),
        )

    def diagnostics_frame(self) -> List[Dict[str, object]]:
        """Return the per-year record, GroinCallback's columns plus b and r."""
        cum_up = np.cumsum(self.applied_dx_updrift_TS) if self.applied_dx_updrift_TS else []
        cum_down = np.cumsum(self.applied_dx_downdrift_TS) if self.applied_dx_downdrift_TS else []
        rows: List[Dict[str, object]] = []
        for i, year in enumerate(self.year_TS):
            rows.append(dict(
                model_year=year,
                groin_active=self.active_TS[i],
                blocking_fraction_applied=self.blocking_applied_TS[i],
                trapping_rate_applied_m_yr=self.trapping_rate_applied_TS[i],
                applied_dx_updrift_m=self.applied_dx_updrift_TS[i],
                applied_dx_downdrift_m=self.applied_dx_downdrift_TS[i],
                cumulative_updrift_m=float(cum_up[i]),
                cumulative_downdrift_m=float(cum_down[i]),
                face_offset_m=self.face_offset_m_TS[i],
                r_ipl_updrift=self.r_ipl_updrift_TS[i],
                r_ipl_downdrift=self.r_ipl_downdrift_TS[i],
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
