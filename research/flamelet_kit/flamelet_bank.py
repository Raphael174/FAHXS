"""
FlameletBank: minimal manager pattern for a set of representative flamelets
keyed by operating condition, each advanced through its own SteadyCache, with
extraction of physical fields back onto a caller-supplied mixture fraction.

This is intentionally generic and NOT tied to any rocket/HX geometry -- it
shows the intended usage pattern: one Flamelet per "representative operating
condition" (e.g. one per heat-exchanger zone, or one per (T_in, chi_st) bin),
each updated only when its own condition has moved outside cache tolerance.

The "map back to physical space" step (`field_at`) is the same Bilger-Z
pattern the rocket solver uses: the flamelet solution lives on a 1-D Z grid;
the host CFD/1-D flow field carries its own local mixture fraction Z(x) (from
transporting Z as a conserved scalar, or computing it from local species via
the same Bilger formula); T(x) and Y(x) are then just T_at_Z(Z(x)) /
Y_at_Z(Z(x)) -- a lookup, not a re-solve.

SINGLE-FLAMELET NOTE: the rocket source kept an axial 4-group bank because
a 5 m rocket grain has real axial variation in local mixing/strain
conditions. A single non-premixed mixing zone (one fuel jet meeting one
oxidizer stream, e.g. a diesel/O2 injector) needs exactly ONE representative
flamelet -- there is nothing to bank. For that case, skip this module
entirely and use `Flamelet` + `SteadyCache` directly (see README.md
quickstart / example_run.py); reach for `FlameletBank` only if your design
genuinely has multiple distinct mixing zones/conditions to track at once.
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np

try:
    from .flamelet import Flamelet
    from .steady_cache import SteadyCache
except ImportError:  # executed as a top-level script, not as part of the package
    from flamelet import Flamelet
    from steady_cache import SteadyCache


class FlameletBank:
    """
    Holds Flamelet instances keyed by an arbitrary hashable `key` (e.g. a
    zone name or operating-condition id), each with its own SteadyCache.

    Parameters
    ----------
    mechanism : str
        Cantera mechanism passed to every Flamelet member.
    n_z : int
        Z-grid resolution for every member.
    cache_kwargs : dict, optional
        Keyword arguments forwarded to each member's SteadyCache
        (tol_dT, tol_p, tol_chi, tol_T_fuel, p_ema_tau, chi_ema_tau).
    flamelet_kwargs : dict, optional
        Extra keyword arguments forwarded to each Flamelet's constructor
        (e.g. diff_mask, cvode_rtol/atol).
    """

    def __init__(
        self,
        mechanism: str,
        n_z: int = 65,
        cache_kwargs: Optional[dict] = None,
        flamelet_kwargs: Optional[dict] = None,
    ):
        self._mechanism = mechanism
        self._n_z = n_z
        self._cache_kwargs = dict(cache_kwargs or {})
        self._flamelet_kwargs = dict(flamelet_kwargs or {})
        self._members: dict = {}   # key -> Flamelet
        self._caches: dict = {}    # key -> SteadyCache

    def members(self):
        return dict(self._members)

    def get(self, key) -> Optional[Flamelet]:
        return self._members.get(key)

    def cache_for(self, key) -> Optional[SteadyCache]:
        return self._caches.get(key)

    def ensure_member(
        self,
        key,
        T_ox: float,
        Y_ox: np.ndarray,
        T_fuel: float,
        Y_fuel: np.ndarray,
        p: float,
    ) -> Flamelet:
        """Create-and-initialize the member for `key` if it does not exist yet."""
        if key not in self._members:
            fl = Flamelet(self._mechanism, n_z=self._n_z, **self._flamelet_kwargs)
            fl.init_mixing(T_ox, Y_ox, T_fuel, Y_fuel, p)
            self._members[key] = fl
            self._caches[key] = SteadyCache(**self._cache_kwargs)
        return self._members[key]

    def advance(
        self,
        key,
        dt: float,
        p: float,
        T_ox: float,
        Y_ox: np.ndarray,
        T_fuel: float,
        Y_fuel: np.ndarray,
        chi_st: float,
        dt_elapsed_for_ema: float = 0.0,
    ) -> bool:
        """
        Advance the member `key` by dt, through its SteadyCache.

        Returns True if a real flamelet advance was performed (cache MISS),
        False if the cache HIT and the advance was skipped.

        Raises KeyError if `ensure_member` was never called for this key.
        """
        fl = self._members[key]
        cache = self._caches[key]

        p_key = cache.key_p(p, dt_elapsed_for_ema)
        chi_key = cache.key_chi(chi_st, dt_elapsed_for_ema)

        if cache.is_hit(p_key, chi_key, T_fuel):
            return False

        T_before = fl.T_max
        fl.step(dt, p, T_ox, Y_ox, T_fuel, Y_fuel, chi_st)
        dT = abs(fl.T_max - T_before)
        cache.record_advance(p_key, chi_key, T_fuel, dT)
        return True

    def field_at(
        self,
        key,
        Z_query,
        species_names: Optional[list] = None,
    ) -> dict:
        """
        Interpolate the member's converged flamelet solution onto a
        caller-supplied mixture fraction (array or scalar).

        Returns {"T": ..., "Y": {species: ...}} with the same shape as
        Z_query. If `species_names` is None, all species in the mechanism
        are returned.
        """
        fl = self._members[key]
        names = species_names if species_names is not None else fl.species_names
        return {
            "T": fl.T_at_Z(Z_query),
            "Y": {sp: fl.Y_at_Z(Z_query, sp) for sp in names},
        }
