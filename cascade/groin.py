"""Hold back a share of the alongshore drift at a groin

This module provides functions for adding a groin to the alongshore transport step in
CASCADE. Each model year, after every Barrier3D domain has run and before BRIE's
alongshore solve, the GroinManager adds to the per-domain shoreline change passed to
BRIE (``x_s_dt``), including:

1. the drift groin: a share (the strength) of BRIE's net drift arriving at the groin
   face, moved from the downdrift to the updrift domain,
2. the leak correction: BRIE's diffusive term at the groin's own bends, cancelled
   within ``half_width`` cells of the face so that BRIE does not leak the trapped sand
   back across the groin,
3. optionally, the spread: the updrift domain's yearly loss across the groin face
   shared with the two domains beyond it, for a fillet wider than one domain.

The groin stops the run with GroinOutsideTestedRange when it leaves the conditions it
was tested in. How it works, its tested range and its fit at Buxton, NC:
scripts/real-world-cases/Hatteras/hard-structures/groin/README.md.

References
----------

.. [1] Nienhuis, J. H., & Lorenzo Trueba, J. (2019). Simulating barrier island
       response to sea level rise with the barrier island and inlet environment
       (BRIE) model v1.0. Geoscientific Model Development, 12(9), 4013-4030.

Notes
-----
Units are metres and ``x_s`` is landward positive. Domain indices are padded indices,
i.e., indices into the full list of CASCADE domains, including any buffer domains.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
"""

from collections import namedtuple

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

__all__ = ["GroinManager", "GroinOutsideTestedRange", "strength_in_year"]


class GroinOutsideTestedRange(ValueError):
    """The run has left the conditions the groin module was tested in; its results would not hold."""


# Last year's state, kept to work out this year's footprint (each an array over the padded domains):
# shoreline at the start of the year, the year's shoreline change passed to BRIE (groin included), the
# groin's own part of that change, and the footprint used that year
_LastYear = namedtuple("_LastYear", ["x_s_m", "dx_m", "groin_dx_m", "footprint_m"])


def strength_in_year(strength_history, year):
    """
    The groin strength in a model year: the value of the last ``(from_year, strength)``
    pair at or before ``year``; 0 before the first pair.

    Parameters
    ----------
    strength_history: list of (int, float)
        ``(from_year, strength)`` pairs, sorted by year
    year: int
        Model year

    Returns
    -------
    float
        Strength in ``year`` [0-1]
    """
    value = 0.0
    for from_year, strength in strength_history:
        if year >= from_year:
            value = strength
    return value


def face_angles(brie, x_s_m):
    """
    Shoreline angle of each cell's forward face [deg]

    Parameters
    ----------
    brie: Brie
        BRIE instance (read only)
    x_s_m: array of float
        Shoreline position of each domain [m]
    """
    return 180.0 * np.arctan2(np.roll(x_s_m, -1) - x_s_m, brie._dy) / np.pi


def table_index(brie, angle_deg):
    """Index into BRIE's per-angle tables for a shoreline angle, clipped to the table"""
    return np.clip(np.round(90 - angle_deg).astype(int), 1, brie._wave_climl - 1)


class GroinManager:
    """Hold back a share of the drift arriving at a groin, and cancel BRIE's leak of it.

    Examples
    --------
    # >>> from cascade.groin import GroinManager
    # >>> groin = GroinManager(updrift_pad=21, downdrift_pad=20,
    # ...     strength_history=[(1969, 0.25), (2017, 0.0)], start_year=1996, n_domains=40)
    # >>> x_s_dt = groin.update(brie, x_s_dt)
    """

    kind = "wave-driven"
    # The tested range: outside it the run stops (groin README "The module's bounds")
    SAND_TOLERANCE_M = 10.0  # summed shoreline the groin may make or lose, m ...
    SAND_TOLERANCE_FRACTION = 0.05  # ... or this share of its footprint's gross, whichever is larger
    MAX_FACE_ANGLE_DEG = 45.0  # steepest tested groin face, deg
    NO_SPREAD = (1.0, 0.0, 0.0)  # the whole loss across the face on the updrift pad

    def __init__(
        self,
        updrift_pad,
        downdrift_pad,
        strength_history,
        start_year,
        n_domains,
        half_width=10,
        spread_weights=NO_SPREAD,
    ):
        """The GroinManager module

        Parameters
        ----------
        updrift_pad: int, required
            Padded index of the domain updrift of the groin face
        downdrift_pad: int, required
            Padded index of the domain downdrift of the groin face; must be adjacent
            to updrift_pad
        strength_history: list of (int, float), required
            ``(from_year, strength)`` pairs; each strength [0-1], the share of the
            arriving drift the groin holds back, holds until the next pair. The first
            year is the install year; the strength is 0 before it.
        start_year: int, required
            Model year of the first call to update
        n_domains: int, required
            Number of padded domains
        half_width: int, optional
            Cells either side of the groin face the leak correction acts on
        spread_weights: sequence of 3 floats, optional
            Shares of the updrift pad's yearly loss across the face taken from the
            updrift pad and the first and second pads beyond it (away from the face);
            non-negative, summing to 1. The default (1, 0, 0) leaves the loss on the
            updrift pad: no spread. Buxton's fit is (0, 2/3, 1/3).
        """

        history = sorted(
            (int(year), float(strength)) for year, strength in strength_history
        )
        if not history or any(not (0.0 <= strength <= 1.0) for _, strength in history):
            raise ValueError(
                "strength_history needs (year, strength) pairs with strength in [0, 1]."
            )
        self._updrift_pad = int(updrift_pad)
        self._downdrift_pad = int(downdrift_pad)
        self._strength_history = history
        self._install_year = history[0][0]
        self._start_year = int(start_year)
        self._n_domains = int(n_domains)
        self._half_width = int(half_width)
        self._lo_pad = min(self._updrift_pad, self._downdrift_pad)
        self._hi_pad = max(self._updrift_pad, self._downdrift_pad)
        self._check_zone(self._lo_pad, self._hi_pad, self._half_width, self._n_domains)
        self._spread_weights = tuple(float(weight) for weight in spread_weights)
        if (
            len(self._spread_weights) != 3
            or any(weight < 0.0 for weight in self._spread_weights)
            or abs(sum(self._spread_weights) - 1.0) > 1e-9
        ):
            raise ValueError("spread_weights needs 3 non-negative shares summing to 1.")
        # the updrift pad and the two beyond it, away from the face
        away = 1 if self._updrift_pad == self._hi_pad else -1
        self._spread_pads = tuple(self._updrift_pad + away * step for step in range(3))
        if min(self._spread_pads) < 0 or max(self._spread_pads) > self._n_domains - 1:
            raise ValueError(
                f"the spread pads {self._spread_pads} run past the end of the "
                f"{self._n_domains} domains."
            )

        self._record = []  # one row per model year
        self._call_count = 0
        self._coast_diff_at_start = None  # BRIE's diffusivity table in the first year
        self._drift_table_m3_yr = None  # BRIE's net drift per shoreline-angle bin, read once
        self._footprint_m = None  # how far the groin has moved each domain's shoreline
        self._last_year = None  # _LastYear, for the footprint

    def strength(self, year):
        """The strength in ``year``: the history's value, 0 before its first year."""
        return strength_in_year(self._strength_history, year)

    @staticmethod
    def _drift_table(brie):
        """BRIE's net alongshore drift for its wave climate, per shoreline-angle bin."""
        bin_width_deg = 180.0 / (brie._wave_climl - 1)
        angle_pdf = brie._angles.pdf(brie._angle_array) * np.deg2rad(bin_width_deg)
        convolved = np.convolve(angle_pdf, brie._coast_qs, mode="full")
        n_pad = len(brie._coast_qs) - 1
        first_bin = n_pad - n_pad // 2
        return -convolved[first_bin : first_bin + len(angle_pdf)]

    @staticmethod
    def _diffusion_number(brie, x_s_m):
        """BRIE's diffusion number r = D dt / (2 dy^2) at each cell's forward face, clipped at zero."""
        table_rows = table_index(brie, face_angles(brie, x_s_m))
        return np.maximum(0.0, brie._coast_diff[table_rows] * brie._dt / 2 / brie._dy**2)

    @classmethod
    def _leak(cls, brie, x_s_m):
        """The leak term 2 (r_i - r_{i-1}) (x_i - x_{i-1}) for shoreline x, from BRIE's diffusivity."""
        step_to_next_m = np.roll(x_s_m, -1) - x_s_m
        diffusion_number = cls._diffusion_number(brie, x_s_m)
        return (
            2 * (diffusion_number - np.roll(diffusion_number, 1)) * np.roll(step_to_next_m, 1)
        )

    @classmethod
    def _brie_step(cls, brie, x_s_m, dx_m):
        """BRIE's implicit alongshore step (brie.py:1293-1326) for shoreline x_s_m and change dx_m."""
        diffusion_number = cls._diffusion_number(brie, x_s_m)
        diagonals = np.r_[
            -diffusion_number[-1],
            -diffusion_number[1:],
            1 + 2 * diffusion_number,
            -diffusion_number[0:-1],
            -diffusion_number[0],
        ]
        right_hand_side = (
            x_s_m
            + diffusion_number * (np.roll(x_s_m, -1) - 2 * x_s_m + np.roll(x_s_m, 1))
            + dx_m
        )
        return spsolve(csr_matrix((diagonals, (brie._di, brie._dj))), right_hand_side)

    @staticmethod
    def _check_zone(lo_pad, hi_pad, half_width, n_domains):
        """Raise unless the pads are adjacent and the correction zone fits inside the array (it does not wrap)."""
        if hi_pad - lo_pad != 1:
            raise ValueError(
                f"groin domains must be adjacent (share the blocked boundary); got pads "
                f"{lo_pad} and {hi_pad}."
            )
        if lo_pad + 1 - half_width < 0 or hi_pad + half_width - 1 > n_domains - 1:
            raise GroinOutsideTestedRange(
                f"the correction zone ({half_width} cells either side of the face "
                f"{lo_pad}|{hi_pad}) runs past the end of the {n_domains} domains; it does "
                "not wrap, so the groin would leak sand there."
            )

    def _check_tested_range(self, brie, x_s_m):
        """Raise GroinOutsideTestedRange if the wave climate changed, sand is unbalanced, or the face is too steep."""
        # wave climate: the drift table is read once, so it must not change
        if self._coast_diff_at_start is None:
            self._coast_diff_at_start = np.array(brie._coast_diff, dtype=float, copy=True)
        elif not np.array_equal(brie._coast_diff, self._coast_diff_at_start):
            raise GroinOutsideTestedRange(
                "the wave climate changed during the run; the groin module is tested only "
                "with one wave climate for the whole run "
                "(hard-structures/groin/3-module-tests/6-stress/)."
            )
        # sand balance: the footprint should sum to about zero
        net_m = float(np.sum(self._footprint_m))
        gross_m = float(np.sum(np.abs(self._footprint_m)))
        if abs(net_m) > max(
            self.SAND_TOLERANCE_M, self.SAND_TOLERANCE_FRACTION * gross_m
        ):
            raise GroinOutsideTestedRange(
                f"the groin has {'made' if net_m < 0 else 'lost'} {abs(net_m):.1f} m of "
                f"summed shoreline (footprint gross {gross_m:.1f} m), more than the "
                "module's tolerance; BRIE's diffusion at the groin is outside the range "
                "the module was tested in (hard-structures/groin/3-module-tests/6-stress/)."
            )
        # face angle: past this BRIE's diffusivity is clipped to zero
        face_angle_deg = abs(face_angles(brie, x_s_m)[self._lo_pad])
        if face_angle_deg > self.MAX_FACE_ANGLE_DEG:
            raise GroinOutsideTestedRange(
                f"the groin face is at {face_angle_deg:.1f} deg, steeper than "
                f"{self.MAX_FACE_ANGLE_DEG:g} deg; past this BRIE's diffusion no longer "
                "carries sand back across the groin and the fillet grows without limit "
                "(hard-structures/groin/3-module-tests/6-stress/)."
            )

    def update(
        self,
        brie,
        x_s_dt,
    ):
        """Move a share of the arriving drift across the face, and cancel BRIE's leak of it.

        Parameters
        ----------
        brie: Brie
            BRIE instance, before its alongshore solve (read only)
        x_s_dt: list or array of float
            Per-domain shoreline change [m, landward +] before the alongshore solve.
            Modified in place and returned.

        Returns
        -------
        x_s_dt: list or array of float
            Per-domain shoreline change with the groin's contribution [m]
        """
        year = self._start_year + self._call_count
        self._call_count += 1
        active = year >= self._install_year
        strength = self.strength(year)

        x_s_m = np.asarray(brie.x_s, dtype=float)
        n_cells, lo_pad, hi_pad = len(x_s_m), self._lo_pad, self._hi_pad
        if self._drift_table_m3_yr is None:
            self._drift_table_m3_yr = self._drift_table(brie)
        if self._footprint_m is None:
            self._footprint_m = np.zeros(n_cells)
        if self._last_year is not None:
            # footprint: last year's BRIE step with and without the groin's own contribution
            last = self._last_year
            self._footprint_m = self._brie_step(
                brie, last.x_s_m, last.dx_m
            ) - self._brie_step(
                brie, last.x_s_m - last.footprint_m, last.dx_m - last.groin_dx_m
            )

        # drift groin: a share of the drift arriving at the face, from whichever side, as
        # shoreline change
        face_angles_deg = face_angles(brie, x_s_m)
        net_drift_m3_yr = self._drift_table_m3_yr[table_index(brie, face_angles_deg)]
        arrival_face = hi_pad if net_drift_m3_yr[hi_pad] < 0 else lo_pad - 1
        arriving_drift_m3_yr = float(net_drift_m3_yr[arrival_face])
        trapped_m = (
            strength
            * arriving_drift_m3_yr
            * brie._dt
            / (brie._dy * (brie._h_b_crit + brie._d_sf))
        )  # m3/yr to m of shoreline
        trap_lo_m, trap_hi_m = (-trapped_m, trapped_m) if strength > 0.0 else (0.0, 0.0)
        x_s_dt[lo_pad] += trap_lo_m
        x_s_dt[hi_pad] += trap_hi_m

        # leak correction near the face: this shoreline's leak term less the groin-free
        # shoreline's
        in_correction_zone = (
            np.abs(np.arange(n_cells) - (lo_pad + 0.5)) <= self._half_width
        )
        leak_correction_m = np.where(
            in_correction_zone,
            self._leak(brie, x_s_m) - self._leak(brie, x_s_m - self._footprint_m),
            0.0,
        )
        x_s_dt += leak_correction_m
        if active:
            self._check_tested_range(brie, x_s_m)

        # spread: the updrift pad's coming loss across the face (BRIE's term), shared with
        # the two pads beyond it, only while the groin works and the updrift pad stands
        # seaward of the downdrift one
        spread_loss_m = 0.0
        if strength > 0.0 and self._spread_weights != self.NO_SPREAD:
            spread_loss_m = max(
                0.0,
                2
                * self._diffusion_number(brie, x_s_m)[self._updrift_pad]
                * (x_s_m[self._downdrift_pad] - x_s_m[self._updrift_pad]),
            )
        self._apply_spread(x_s_dt, spread_loss_m)

        # keep this year's state and the groin's own contribution for next year's footprint
        groin_dx_m = leak_correction_m.copy()
        groin_dx_m[lo_pad] += trap_lo_m
        groin_dx_m[hi_pad] += trap_hi_m
        self._apply_spread(groin_dx_m, spread_loss_m)
        self._last_year = _LastYear(
            x_s_m.copy(),
            np.asarray(x_s_dt, dtype=float).copy(),
            groin_dx_m,
            self._footprint_m,
        )

        # record the year
        trap_updrift_m, trap_downdrift_m = (
            (trap_lo_m, trap_hi_m)
            if self._updrift_pad == lo_pad
            else (trap_hi_m, trap_lo_m)
        )
        self._record.append(
            dict(
                model_year=year,
                groin_active=active,
                strength_applied=strength,
                applied_dx_updrift_m=trap_updrift_m,
                applied_dx_downdrift_m=trap_downdrift_m,
                face_offset_m=float(x_s_m[hi_pad] - x_s_m[lo_pad]),
                net_drift_m3_yr=arriving_drift_m3_yr,
                approach_angle_deg=float(face_angles_deg[arrival_face]),
                leak_correction_m=float(leak_correction_m.sum()),
                footprint_max_m=float(np.abs(self._footprint_m).max()),
                spread_loss_m=float(spread_loss_m),
            )
        )
        return x_s_dt

    def _apply_spread(self, dx_m, spread_loss_m):
        """Give the updrift pad back the shares of the loss the two pads beyond it take instead (landward +)."""
        if spread_loss_m == 0.0:
            return
        updrift, first_beyond, second_beyond = self._spread_pads
        dx_m[updrift] -= (1 - self._spread_weights[0]) * spread_loss_m
        dx_m[first_beyond] += self._spread_weights[1] * spread_loss_m
        dx_m[second_beyond] += self._spread_weights[2] * spread_loss_m

    def summary(self):
        """Return a compact run-metadata dictionary for logging."""
        return dict(
            kind=self.kind,
            updrift_pad=self._updrift_pad,
            downdrift_pad=self._downdrift_pad,
            strength_history=[list(pair) for pair in self._strength_history],
            install_year=self._install_year,
            start_year=self._start_year,
            half_width=self._half_width,
            spread_weights=list(self._spread_weights),
            years_active=sum(row["groin_active"] for row in self._record),
            mean_trapping_rate_m_yr=self.mean_trapping_rate_m_yr,
            cumulative_updrift_m=float(
                np.sum([row["applied_dx_updrift_m"] for row in self._record])
            ),
            cumulative_downdrift_m=float(
                np.sum([row["applied_dx_downdrift_m"] for row in self._record])
            ),
        )

    def diagnostics_frame(self):
        """Return the per-year record, with the shoreline change moved so far."""
        cumulative_updrift_m = np.cumsum(
            [row["applied_dx_updrift_m"] for row in self._record]
        )
        cumulative_downdrift_m = np.cumsum(
            [row["applied_dx_downdrift_m"] for row in self._record]
        )
        return [
            dict(
                model_year=row["model_year"],
                groin_active=row["groin_active"],
                strength_applied=row["strength_applied"],
                trapping_rate_applied_m_yr=abs(row["applied_dx_updrift_m"]),
                applied_dx_updrift_m=row["applied_dx_updrift_m"],
                applied_dx_downdrift_m=row["applied_dx_downdrift_m"],
                cumulative_updrift_m=float(cumulative_updrift_m[year_index]),
                cumulative_downdrift_m=float(cumulative_downdrift_m[year_index]),
                face_offset_m=row["face_offset_m"],
                net_drift_m3_yr=row["net_drift_m3_yr"],
                approach_angle_deg=row["approach_angle_deg"],
                leak_correction_m=row["leak_correction_m"],
                footprint_max_m=row["footprint_max_m"],
                spread_loss_m=row["spread_loss_m"],
            )
            for year_index, row in enumerate(self._record)
        ]

    @property
    def mean_trapping_rate_m_yr(self):
        """Mean shoreline change moved across the face over the active years [m/yr]"""
        moved_m = [
            abs(row["applied_dx_updrift_m"]) for row in self._record if row["groin_active"]
        ]
        return float(np.mean(moved_m)) if moved_m else float("nan")

    @property
    def updrift_pad(self):
        return self._updrift_pad

    @property
    def downdrift_pad(self):
        return self._downdrift_pad

    @property
    def strength_history(self):
        return self._strength_history

    @property
    def install_year(self):
        return self._install_year

    @property
    def start_year(self):
        return self._start_year

    @property
    def n_domains(self):
        return self._n_domains

    @property
    def half_width(self):
        return self._half_width

    @property
    def spread_weights(self):
        return self._spread_weights

    @property
    def record(self):
        return self._record

    @property
    def footprint(self):
        return self._footprint_m
