"""Correlation registry: fluid/validity/geometry-tagged closures + selection.

The registry is the mechanism behind the project's fluid-agnostic goal: choose
the heat-transfer closure "most adapted to the fluid used and thermodynamic
conditions", and ALWAYS report when the chosen closure is being used outside its
validated envelope (mirroring the honest CHF/ONB/HTD warning channel elsewhere
in this package).

Two kinds of entry live here:

- **Executable closures** (supercritical family, Phase 1): a real ``callable``
  taking a :class:`ClosureContext` and returning an HTC. ``select_supercritical``
  filters by regime/geometry, ranks by fluid-specificity and in-validity, and
  returns the pick plus an :class:`ExtrapolationReport`.
- **Metadata-only records** (subcritical family): ``callable=None``. The
  subcritical dispatch path in ``dispatch.py`` still executes its existing inline
  hot path (kept byte-for-byte for bit-identical results), but consults the
  registry via ``validity_report_for`` to surface the SAME honest extrapolation
  reporting and to make previously-silent fluid restrictions explicit (e.g. the
  Groeneveld CHF LUT and Bergles-Rohsenow ONB are water-only fits).

Closure modules self-register at import (see supercritical.py); this module has
no dependency on them, avoiding an import cycle.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

# Tiers, ordered worst-to-best confidence (higher = prefer).
TIER_CONSERVATIVE_BOUND = "conservative_bound"
TIER_STRUCTURAL_EXTRAPOLATION = "structural_extrapolation"
TIER_VALIDATED_IN_RANGE = "validated_in_range"
_TIER_RANK = {
    TIER_CONSERVATIVE_BOUND: 0,
    TIER_STRUCTURAL_EXTRAPOLATION: 1,
    TIER_VALIDATED_IN_RANGE: 2,
}

FLUID_ANY = "any"


@dataclass(frozen=True)
class ClosureContext:
    """Everything a heat-transfer closure might need at one march node.

    Closures pull what they need and ignore the rest. ``wall_temp_K`` is the
    lagged wall temperature (one node / one sweep behind, same pattern as the
    boiling Bo term); None when not yet available (first node/sweep), in which
    case property-ratio-correction closures fall back to a ratio of 1.0.
    """

    fluid: str
    p_Pa: float
    h_J_kg: float
    T_bulk_K: float
    rho_b: float
    mu_b: float
    k_b: float
    cp_b: float
    Pr_b: float
    mass_flux_kg_m2_s: float
    diameter_m: float
    heat_flux_W_m2: float
    wall_temp_K: float | None = None
    # Flow-length / diameter ratio (Taylor's entrance-effect term). None ->
    # closures that use it fall back to their long-channel (fully developed)
    # limit.
    x_over_D: float | None = None
    # Calibration knobs (a CorrelationCoefficients instance) for closures that
    # need them -- e.g. the gas-side forced-convection closures registered in
    # gas_closures.py. None for closures that don't use calibration.
    corrCoeffs: object | None = None
    # Escape hatch for closure-specific scalars with no natural home in the
    # common bulk-property fields above (e.g. raw axial position "x_m" for a
    # developing-length correction -- distinct from x_over_D above, which
    # supercritical closures consume differently; "roughness_m";
    # "corrugation_thickness_m"/"corrugation_pitch_m"). Closures pull what
    # they need by key and raise a clear KeyError if it's missing, same as
    # any other required argument.
    extra: dict = field(default_factory=dict)


@dataclass(frozen=True)
class ExtrapolationReport:
    """Which validity limits the selected closure violates at this point.

    ``in_range`` is True when nothing is violated. ``violations`` maps a
    parameter name to (value, low, high, relative_overshoot) so callers can
    print a specific, quantified warning rather than a bare "out of range".
    """

    closure_name: str
    in_range: bool
    violations: dict[str, tuple[float, float | None, float | None, float]] = field(
        default_factory=dict
    )

    def message(self) -> str | None:
        if self.in_range:
            return None
        parts = []
        for param, (val, lo, hi, rel) in self.violations.items():
            bound = f"[{lo}, {hi}]"
            parts.append(f"{param}={val:.4g} outside {bound} (by {rel*100:.0f}%)")
        return f"closure {self.closure_name!r} extrapolated: " + "; ".join(parts)


@dataclass(frozen=True)
class ClosureRecord:
    """``callable`` returns an HTC [W/m2K] for every regime EXCEPT
    ``"gas_forced_convection_friction"`` (see gas_closures.py), where it
    returns a dimensionless Darcy friction factor instead -- a deliberate,
    documented broadening of this field's meaning, not silent reuse. Records
    are never ranked against each other across that regime boundary (the
    regime tag is a hard filter in both ``select_supercritical`` and
    ``get_record``), so an h-returning and f-returning record can never be
    compared as if they were interchangeable.
    """

    name: str
    regime_tags: frozenset
    geometry_tags: frozenset       # e.g. straight_tube, helical_coil, shell_crossflow
    orientation_tags: frozenset    # e.g. vertical, horizontal, any
    fluid_scope: frozenset         # explicit fluid names, or {FLUID_ANY}
    validity: dict                 # param -> (lo|None, hi|None)
    provenance: str                # docs/reference stem or citation
    tier: str
    callable: Callable[[ClosureContext], float] | None = None
    priority: int = 0              # final deterministic tiebreak (higher wins)

    def covers_fluid(self, fluid: str) -> bool:
        return FLUID_ANY in self.fluid_scope or fluid in self.fluid_scope


REGISTRY: dict[str, ClosureRecord] = {}


def register(record: ClosureRecord) -> None:
    if record.name in REGISTRY:
        raise ValueError(f"closure {record.name!r} already registered")
    REGISTRY[record.name] = record


def get_record(name: str) -> ClosureRecord:
    return REGISTRY[name]


def check_validity(record: ClosureRecord, operating_point: dict) -> ExtrapolationReport:
    """Compare an operating point against a record's validity ranges.

    ``operating_point`` maps parameter names to values; only parameters that
    appear in BOTH the operating point and the record's validity dict are
    checked (a range with no corresponding operating-point value is simply not
    evaluated -- the caller supplies whatever it can measure).
    """
    violations: dict[str, tuple[float, float | None, float | None, float]] = {}
    for param, (lo, hi) in record.validity.items():
        if param not in operating_point:
            continue
        val = float(operating_point[param])
        if lo is not None and val < lo:
            rel = (lo - val) / abs(lo) if lo != 0 else float("inf")
            violations[param] = (val, lo, hi, rel)
        elif hi is not None and val > hi:
            rel = (val - hi) / abs(hi) if hi != 0 else float("inf")
            violations[param] = (val, lo, hi, rel)
    return ExtrapolationReport(
        closure_name=record.name, in_range=not violations, violations=violations
    )


def validity_report_for(name: str, operating_point: dict) -> ExtrapolationReport:
    """Extrapolation report for an already-chosen closure (subcritical inline
    path uses this to get honest reporting without registry-driven selection)."""
    return check_validity(REGISTRY[name], operating_point)


def select_supercritical(
    *,
    regime: str,
    geometry: str,
    orientation: str,
    fluid: str,
    operating_point: dict,
) -> tuple[ClosureRecord, ExtrapolationReport]:
    """Pick the best supercritical closure for the given context.

    Ranking (best first), **in-validity-range is the PRIMARY key**:
      1. in range AND fluid-specific
      2. in range AND generic (any-fluid)
      3. out of range AND fluid-specific
      4. out of range AND generic
    Then geometry match, orientation match, tier (validated > structural >
    conservative), and finally the record's ``priority`` as a deterministic
    tiebreak between equally-ranked closures.

    Rationale for in-range-first (changed 2026-07-17 after the supercritical-N2
    literature review): a fluid-specific correlation used far outside its
    validated Reynolds/pressure envelope can be *less* trustworthy than an
    in-range generic property-ratio form -- Locke & Landrum (2008) show the
    bulk-reference correlations overpredict (spike) around the pseudo-critical
    line when extrapolated. So an in-range generic closure now beats an
    out-of-range fluid-specific one, while an in-range fluid-specific closure
    still wins outright when it exists (e.g. Cheng2020 inside its own
    7000-27000 Re window). Geometry/orientation are soft-matched; an
    ``any``-tagged closure is always eligible so there is a fallback. Raises
    LookupError if nothing matches the regime (should not happen given the
    registered conservative-bound fallback).
    """
    candidates = [
        r for r in REGISTRY.values()
        if r.callable is not None and regime in r.regime_tags
    ]
    if not candidates:
        raise LookupError(f"no executable closure registered for regime {regime!r}")

    def score(r: ClosureRecord):
        report = check_validity(r, operating_point)
        fluid_specific = FLUID_ANY not in r.fluid_scope and fluid in r.fluid_scope
        geom_match = geometry in r.geometry_tags
        orient_match = orientation in r.orientation_tags or "any" in r.orientation_tags
        return (
            int(report.in_range),        # PRIMARY: in validity range
            # SECOND: tier (validated_in_range > structural_extrapolation >
            # conservative_bound) -- this must outrank geometry/orientation
            # tag matching, or a broadly-tagged always-in-range fallback (e.g.
            # gnielinski_bulk_bound, tagged for shell_crossflow specifically
            # because nothing better is validated there) would out-score a
            # genuinely validated-in-range correlation just because the
            # validated one isn't tagged for this exact geometry (found via
            # the shell-and-tube LN2 case, 2026-07-19: this ordering bug
            # picked Gnielinski over an in-range McCarthy-Wolf).
            _TIER_RANK.get(r.tier, 0),
            int(fluid_specific),         # then fluid-specificity
            int(geom_match),
            int(orient_match),
            r.priority,                  # final deterministic tiebreak
        )

    best = max(candidates, key=score)
    return best, check_validity(best, operating_point)
