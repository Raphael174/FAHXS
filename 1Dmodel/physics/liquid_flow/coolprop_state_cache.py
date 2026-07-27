"""Cached low-level CoolProp state objects for the liquid coolant march.

Every property lookup in ``correlations.py``/``dispatch.py`` used to go
through the high-level ``CP.PropsSI(...)`` call, which re-parses the fluid
string and rebuilds a throwaway state object on every single call. Under
repeated march-node evaluation (and especially under the counter-flow
shooting method, which re-runs the whole march many times) this dominates
runtime. Reusing one persistent ``CP.AbstractState`` per (backend, fluid)
pair - measured this session - cuts the per-call cost by roughly an order of
magnitude for the exact "HEOS" backend alone, with zero change in the
underlying equation of state (same physics, same answer to floating-point
precision).

``TTSE`` and ``BICUBIC`` request CoolProp's tabulated interpolation backends
instead of the exact HEOS evaluation - measured this session at roughly
300-1000x faster than a single HEOS PropsSI call once the (one-time, cached
to disk by CoolProp) table is built, at the cost of interpolation error.
These are opt-in (default stays "HEOS") because interpolation error is
largest exactly where correctness matters most: near the saturation dome,
where quality and CHF margin are evaluated. See
``validation/liquid_ttse_backend_validation.py`` for measured error bounds
before using TTSE/BICUBIC for anything but exploratory speed.

Important: CoolProp's high-level ``PropsSI("BICUBIC::Water", ...)`` string
form does NOT work for TTSE/BICUBIC in this CoolProp build (confirmed this
session - it raises "cannot be used in the high-level interface"). Only the
low-level ``AbstractState`` API supports these backends, which is exactly why
this module exists instead of just prefixing the fluid string.
"""

from __future__ import annotations

import CoolProp.CoolProp as CP

_SUPPORTED_BACKENDS = ("HEOS", "TTSE", "BICUBIC")


def parse_backend(fluid: str) -> tuple[str, str]:
    """Split a possibly backend-tagged fluid string, e.g. "BICUBIC::Water".

    Returns (backend, plain_fluid_name). Untagged strings default to "HEOS".
    """
    if "::" in fluid:
        backend, name = fluid.split("::", 1)
        if backend not in _SUPPORTED_BACKENDS:
            raise ValueError(
                f"unsupported CoolProp property backend: {backend!r} "
                f"(expected one of {_SUPPORTED_BACKENDS})"
            )
        return backend, name
    return "HEOS", fluid


def coolprop_fluid_string(fluid: str, backend: str = "HEOS") -> str:
    """Tag a plain fluid name with a backend selector for the liquid march.

    ``backend="HEOS"`` returns ``fluid`` unchanged (the exact, validated
    default). ``"TTSE"``/``"BICUBIC"`` return e.g. ``"BICUBIC::Water"``,
    which every ``correlations.py``/``dispatch.py`` property call accepts
    transparently (they all route through ``get_cached_state`` below).
    """
    if backend == "HEOS":
        return fluid
    if backend not in _SUPPORTED_BACKENDS:
        raise ValueError(
            f"unsupported CoolProp property backend: {backend!r} "
            f"(expected one of {_SUPPORTED_BACKENDS})"
        )
    return f"{backend}::{fluid}"


class _CachedFluidState:
    """Three persistent AbstractState objects for one (backend, fluid) pair.

    Three, not one: a single AbstractState only holds ONE state at a time, and
    saturated-liquid (Q=0) and saturated-vapor (Q=1) properties are needed
    simultaneously (e.g. for two-phase mixture properties) - a shared object
    would require re-``update()``-ing between reads, which is exactly the
    per-call overhead this cache exists to avoid.
    """

    __slots__ = (
        "_general", "_sat_l", "_sat_v", "_probe", "_wall", "_aux",
        "molar_mass_kg_mol", "p_crit_Pa", "T_crit_K",
    )

    def __init__(self, backend: str, fluid: str) -> None:
        self._general = CP.AbstractState(backend, fluid)
        self._sat_l = CP.AbstractState(backend, fluid)
        self._sat_v = CP.AbstractState(backend, fluid)
        # (T, p) states, separate from the (p, h) hot-path _general so wall
        # property lookups and the pseudo-critical-temperature search don't
        # clobber the bulk state mid-evaluation. _wall: wall props at (T_w, p)
        # for property-ratio corrections; _aux: cp scan for T_pc (rare, cached).
        self._wall = CP.AbstractState(backend, fluid)
        self._aux = CP.AbstractState(backend, fluid)
        # Plain HEOS probe, used for (a) molar_mass/p_crit (backend-independent
        # constants) and (b) surface tension - confirmed this session as a
        # CoolProp quirk: for TTSE/BICUBIC specifically, only the FIRST
        # AbstractState instance of a given tabulated backend+fluid pair can
        # compute surface_tension(); every other coexisting instance of that
        # same tabulated backend raises "only defined within the two-phase
        # region" even at an exact Q=0 point. Routing surface tension through
        # a dedicated HEOS instance (which has no such issue, verified with
        # multiple coexisting instances) sidesteps it entirely, and costs
        # nothing extra since it isn't a hot-path property.
        self._probe = CP.AbstractState("HEOS", fluid)
        self.molar_mass_kg_mol = float(self._probe.molar_mass())
        self.p_crit_Pa = float(self._probe.p_critical())
        self.T_crit_K = float(self._probe.T_critical())

    def flash_ph(self, p_Pa: float, h_J_kg: float):
        """Update the general-purpose state to (p, h) and return it."""
        self._general.update(CP.HmassP_INPUTS, h_J_kg, p_Pa)
        return self._general

    def wall_state_tp(self, T_K: float, p_Pa: float):
        """Update and return the wall (T, p) state (for property-ratio terms)."""
        self._wall.update(CP.PT_INPUTS, p_Pa, T_K)
        return self._wall

    def cp_at_tp(self, T_K: float, p_Pa: float) -> float:
        """Isobaric specific heat at (T, p), via the dedicated aux state.

        Used only by the pseudo-critical-temperature search (which is itself
        cached per fluid/pressure), so it never runs on the property hot path.
        """
        self._aux.update(CP.PT_INPUTS, p_Pa, T_K)
        return float(self._aux.cpmass())

    def saturated_liquid(self, p_Pa: float):
        self._sat_l.update(CP.PQ_INPUTS, p_Pa, 0.0)
        return self._sat_l

    def saturated_vapor(self, p_Pa: float):
        self._sat_v.update(CP.PQ_INPUTS, p_Pa, 1.0)
        return self._sat_v

    def surface_tension_liquid(self, p_Pa: float) -> float:
        self._probe.update(CP.PQ_INPUTS, p_Pa, 0.0)
        return float(self._probe.surface_tension())


_CACHE: dict[tuple[str, str], _CachedFluidState] = {}


def get_cached_state(fluid: str) -> _CachedFluidState:
    """Return the cached low-level state wrapper for a (possibly
    backend-tagged) fluid string, creating and caching it on first use.

    The cache is process-global and grows by at most a handful of entries in
    practice (one or two fluids, one backend each, for the lifetime of a
    Python process) - not a memory concern.
    """
    backend, name = parse_backend(fluid)
    key = (backend, name)
    cached = _CACHE.get(key)
    if cached is None:
        cached = _CachedFluidState(backend, name)
        _CACHE[key] = cached
    return cached
