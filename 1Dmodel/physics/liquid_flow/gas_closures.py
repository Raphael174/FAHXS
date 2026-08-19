"""Gas-side forced-convection closures, registered into the same registry as
the liquid/supercritical family (Stage D, Slice 1 of
``docs/solver_design/FV_CORE_REWORK_PLAN.md``).

Scope, confirmed by reading ``transient_core/adapters_shelltube.py::
shelltube_tube_gas_film`` (the only function that computes shell-and-tube's
hot-gas film today): shell-and-tube's tube-side hot gas uses exactly four
correlation functions, selected by ``shellTubeProp.inside_tube_choice ∈
{"smooth", "grooved"}``. This module registers all four AS DELEGATING
WRAPPERS around the existing, validated functions in
``physics/heat_transfer_correlations.py`` / ``physics/friction_correlations.py``
-- the physics is not reimplemented here, only given honest
registry-selection and (future) extrapolation reporting. Bit-identical
equivalence to the direct legacy calls is asserted in
``tests/test_core_closures.py``.

Deliberately NOT registered here: ``dispatch_nu_coil``/``dispatch_nu_shell``/
``dispatch_friction_coil`` (helical-only -- shell-and-tube never calls them)
and Bell-Delaware (returns ``(h, dp)`` from one call against a whole geometry
dict, not a single scalar from a ``ClosureContext``; needs its own two-output
closure protocol, not a shoehorn here).

Two new regime tags, since none of ``regime.py``'s six labels fit "ideal-gas
forced convection": ``"gas_forced_convection"`` (Nu-type closures, return an
HTC) and ``"gas_forced_convection_friction"`` (friction-type closures, return
a dimensionless Darcy factor -- see ``ClosureRecord``'s docstring for why
that's a deliberate, documented broadening of the "callable returns an HTC"
contract, not silent reuse).
"""

from __future__ import annotations

from hps_combustor.physics.friction_correlations import (
    dispatch_friction_tube_straight,
    friction_corrugated_tube_vicente,
)
from hps_combustor.physics.heat_transfer_correlations import (
    dispatch_nu_tube_straight,
    nu_corrugated_tube_vicente,
)
from hps_combustor.physics.liquid_flow.registry import (
    FLUID_ANY,
    TIER_VALIDATED_IN_RANGE,
    ClosureContext,
    ClosureRecord,
    register,
)

GAS_FORCED_CONVECTION = "gas_forced_convection"
GAS_FORCED_CONVECTION_FRICTION = "gas_forced_convection_friction"


def _reynolds(ctx: ClosureContext) -> float:
    return ctx.mass_flux_kg_m2_s * ctx.diameter_m / ctx.mu_b


def _x_m(ctx: ClosureContext) -> float:
    """Raw axial position [m] for the developing-length correction --
    distinct from ``x_over_D`` (used differently by supercritical closures).
    ``10e10`` (the legacy fully-developed sentinel) if not supplied."""
    return float(ctx.extra.get("x_m", 10e10))


def _tube_straight_friction(ctx: ClosureContext) -> float:
    roughness_m = float(ctx.extra.get("roughness_m", 0.0))
    return dispatch_friction_tube_straight(
        _reynolds(ctx),
        roughness_m,
        ctx.diameter_m,
        x=_x_m(ctx),
        Re_lo=getattr(ctx.corrCoeffs, "Re_transition_lo", 2300.0),
        Re_hi=getattr(ctx.corrCoeffs, "Re_transition_hi", 4000.0),
    )


def _tube_straight_htc(ctx: ClosureContext) -> float:
    # Same dependency order as the legacy per-cell loop
    # (adapters_shelltube.py:465-502): friction is computed first and fed
    # into Nu as f_fd (Gnielinski needs it).
    f_fd = _tube_straight_friction(ctx)
    Nu = dispatch_nu_tube_straight(
        "gnielinski_blended",
        Re=_reynolds(ctx),
        Pr=ctx.Pr_b,
        d=ctx.diameter_m,
        x=_x_m(ctx),
        f_fd=f_fd,
        T_bulk=ctx.T_bulk_K,
        T_wall=ctx.wall_temp_K,
        error_factor=1.0,
        corrCoeffs=ctx.corrCoeffs,
    )
    return Nu * ctx.k_b / ctx.diameter_m


def _grooved_phi(ctx: ClosureContext) -> float:
    thickness_m = float(ctx.extra.get("corrugation_thickness_m", 0.0))
    pitch_m = float(ctx.extra.get("corrugation_pitch_m", 1.0))
    e = max(thickness_m, 0.0)
    p = max(pitch_m, 1e-12)
    return (e**2) / (p * ctx.diameter_m)


def _tube_grooved_friction(ctx: ClosureContext) -> float:
    # Re_lo/Re_hi from corrCoeffs (default 2300/4000), NOT this function's own
    # defaults (2000/4000) -- matches adapters_shelltube.py:460-461/467-472
    # exactly; the two default sets differ and using the wrong one silently
    # shifts the laminar/turbulent blend band.
    return friction_corrugated_tube_vicente(
        _reynolds(ctx),
        _grooved_phi(ctx),
        Re_lo=getattr(ctx.corrCoeffs, "Re_transition_lo", 2300.0),
        Re_hi=getattr(ctx.corrCoeffs, "Re_transition_hi", 4000.0),
    )


def _tube_grooved_htc(ctx: ClosureContext) -> float:
    Nu = nu_corrugated_tube_vicente(
        _reynolds(ctx),
        ctx.Pr_b,
        _grooved_phi(ctx),
        D_i=ctx.diameter_m,
        x=_x_m(ctx),
        Re_lo=getattr(ctx.corrCoeffs, "Re_transition_lo", 2300.0),
        Re_hi=getattr(ctx.corrCoeffs, "Re_transition_hi", 4000.0),
    )
    return Nu * ctx.k_b / ctx.diameter_m


register(
    ClosureRecord(
        name="tube_straight_gnielinski_blended",
        regime_tags=frozenset({GAS_FORCED_CONVECTION}),
        geometry_tags=frozenset({"straight_tube"}),
        orientation_tags=frozenset({"any"}),
        fluid_scope=frozenset({FLUID_ANY}),
        validity={},
        provenance=(
            "Gnielinski (1976) turbulent + laminar-entrance composite, "
            "linear-blended across the transition band -- see "
            "heat_transfer_correlations.py::dispatch_nu_tube_straight"
        ),
        tier=TIER_VALIDATED_IN_RANGE,
        callable=_tube_straight_htc,
    )
)

register(
    ClosureRecord(
        name="tube_straight_friction_blended",
        regime_tags=frozenset({GAS_FORCED_CONVECTION_FRICTION}),
        geometry_tags=frozenset({"straight_tube"}),
        orientation_tags=frozenset({"any"}),
        fluid_scope=frozenset({FLUID_ANY}),
        validity={},
        provenance=(
            "Hagen-Poiseuille + Colebrook (1939), linear-blended -- see "
            "friction_correlations.py::dispatch_friction_tube_straight"
        ),
        tier=TIER_VALIDATED_IN_RANGE,
        callable=_tube_straight_friction,
    )
)

register(
    ClosureRecord(
        name="tube_grooved_vicente",
        regime_tags=frozenset({GAS_FORCED_CONVECTION}),
        geometry_tags=frozenset({"straight_tube"}),
        orientation_tags=frozenset({"any"}),
        fluid_scope=frozenset({FLUID_ANY}),
        validity={},
        provenance=(
            "Vicente et al., summarized by Cruz et al. (2021) -- see "
            "heat_transfer_correlations.py::nu_corrugated_tube_vicente"
        ),
        tier=TIER_VALIDATED_IN_RANGE,
        callable=_tube_grooved_htc,
    )
)

register(
    ClosureRecord(
        name="tube_grooved_friction_vicente",
        regime_tags=frozenset({GAS_FORCED_CONVECTION_FRICTION}),
        geometry_tags=frozenset({"straight_tube"}),
        orientation_tags=frozenset({"any"}),
        fluid_scope=frozenset({FLUID_ANY}),
        validity={},
        provenance=(
            "Vicente et al., summarized by Cruz et al. (2021) -- see "
            "friction_correlations.py::friction_corrugated_tube_vicente"
        ),
        tier=TIER_VALIDATED_IN_RANGE,
        callable=_tube_grooved_friction,
    )
)
