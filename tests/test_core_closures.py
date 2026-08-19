"""Stage D, Slice 1 of docs/solver_design/FV_CORE_REWORK_PLAN.md.

Proves the four new registered gas-side closures (physics/liquid_flow/
gas_closures.py, wired through core/closures.py) are bit-identical to the
direct legacy calls `transient_core/adapters_shelltube.py::
shelltube_tube_gas_film` makes today -- the whole point of "delegate, don't
reimplement". Covers laminar, transitional, and turbulent Reynolds regimes
for both `inside_tube_choice` values, plus the friction-Nu dependency order
(friction computed first, fed into Nu as f_fd) and the corrCoeffs-sourced
Re_lo/Re_hi (2300/4000, NOT the grooved functions' own different defaults of
2000/4000 -- a real bug caught while writing this closure, see
gas_closures.py's comments).
"""
from __future__ import annotations

import pytest

from hps_combustor.core.closures import tube_friction_closure, tube_htc_closure
from hps_combustor.input_data import CorrelationCoefficients
from hps_combustor.physics.friction_correlations import (
    dispatch_friction_tube_straight,
    friction_corrugated_tube_vicente,
)
from hps_combustor.physics.heat_transfer_correlations import (
    dispatch_nu_tube_straight,
    nu_corrugated_tube_vicente,
)
from hps_combustor.physics.liquid_flow.registry import ClosureContext

CORR = CorrelationCoefficients()

# (Re, Pr, T_bulk, T_wall) points spanning laminar / transitional / turbulent
# -- Re_lo=2300, Re_hi=4000 per CorrelationCoefficients defaults.
POINTS = [
    ("laminar", 1500.0, 0.72, 1200.0, 900.0),
    ("transitional", 3000.0, 0.70, 1100.0, 850.0),
    ("turbulent", 20000.0, 0.68, 900.0, 700.0),
]

D = 0.0135  # representative shell-and-tube tube inner diameter [m]
K = 0.08    # representative hot-gas conductivity [W/m/K]
MU = 4.5e-5


def _ctx(Re, Pr, T_bulk, T_wall, *, extra=None):
    mass_flux = Re * MU / D
    return ClosureContext(
        fluid="combustion_products",
        p_Pa=5e5,
        h_J_kg=0.0,
        T_bulk_K=T_bulk,
        rho_b=1.0,
        mu_b=MU,
        k_b=K,
        cp_b=1500.0,
        Pr_b=Pr,
        mass_flux_kg_m2_s=mass_flux,
        diameter_m=D,
        heat_flux_W_m2=0.0,
        wall_temp_K=T_wall,
        corrCoeffs=CORR,
        extra=extra or {},
    )


@pytest.mark.parametrize("label,Re,Pr,T_bulk,T_wall", POINTS)
def test_tube_straight_friction_matches_legacy(label, Re, Pr, T_bulk, T_wall):
    ctx = _ctx(Re, Pr, T_bulk, T_wall, extra={"x_m": 0.5, "roughness_m": 1.5e-6})
    record = tube_friction_closure("smooth")
    got = record.callable(ctx)
    expected = dispatch_friction_tube_straight(
        Re, 1.5e-6, D, x=0.5,
        Re_lo=CORR.Re_transition_lo, Re_hi=CORR.Re_transition_hi,
    )
    assert got == expected


@pytest.mark.parametrize("label,Re,Pr,T_bulk,T_wall", POINTS)
def test_tube_straight_htc_matches_legacy(label, Re, Pr, T_bulk, T_wall):
    ctx = _ctx(Re, Pr, T_bulk, T_wall, extra={"x_m": 0.5, "roughness_m": 1.5e-6})
    record = tube_htc_closure("smooth")
    got = record.callable(ctx)

    f_fd = dispatch_friction_tube_straight(
        Re, 1.5e-6, D, x=0.5,
        Re_lo=CORR.Re_transition_lo, Re_hi=CORR.Re_transition_hi,
    )
    Nu = dispatch_nu_tube_straight(
        "gnielinski_blended", Re=Re, Pr=Pr, d=D, x=0.5, f_fd=f_fd,
        T_bulk=T_bulk, T_wall=T_wall, error_factor=1.0, corrCoeffs=CORR,
    )
    expected = Nu * K / D
    assert got == expected


@pytest.mark.parametrize("label,Re,Pr,T_bulk,T_wall", POINTS)
def test_tube_grooved_friction_matches_legacy(label, Re, Pr, T_bulk, T_wall):
    ctx = _ctx(
        Re, Pr, T_bulk, T_wall,
        extra={"x_m": 0.5, "corrugation_thickness_m": 3e-4, "corrugation_pitch_m": 4e-3},
    )
    record = tube_friction_closure("grooved")
    got = record.callable(ctx)

    phi = (3e-4**2) / (4e-3 * D)
    expected = friction_corrugated_tube_vicente(
        Re, phi, Re_lo=CORR.Re_transition_lo, Re_hi=CORR.Re_transition_hi,
    )
    assert got == expected


@pytest.mark.parametrize("label,Re,Pr,T_bulk,T_wall", POINTS)
def test_tube_grooved_htc_matches_legacy(label, Re, Pr, T_bulk, T_wall):
    ctx = _ctx(
        Re, Pr, T_bulk, T_wall,
        extra={"x_m": 0.5, "corrugation_thickness_m": 3e-4, "corrugation_pitch_m": 4e-3},
    )
    record = tube_htc_closure("grooved")
    got = record.callable(ctx)

    phi = (3e-4**2) / (4e-3 * D)
    Nu = nu_corrugated_tube_vicente(
        Re, Pr, phi, D_i=D, x=0.5,
        Re_lo=CORR.Re_transition_lo, Re_hi=CORR.Re_transition_hi,
    )
    expected = Nu * K / D
    assert got == expected


def test_grooved_re_thresholds_come_from_corrcoeffs_not_function_defaults():
    """Regression for the bug caught writing gas_closures.py: the corrugated
    functions' own Re_lo/Re_hi defaults (2000/4000) differ from
    CorrelationCoefficients' (2300/4000) -- using the wrong pair silently
    shifts the laminar/turbulent blend band. Pick a Re between the two
    defaults' laminar bounds where the two blends diverge."""
    Re = 2100.0  # below CorrelationCoefficients' Re_lo=2300 (laminar there),
    # but above the corrugated functions' own default Re_lo=2000 (turbulent there)
    ctx = _ctx(Re, 0.7, 1000.0, 800.0, extra={"corrugation_thickness_m": 3e-4, "corrugation_pitch_m": 4e-3})
    phi = (3e-4**2) / (4e-3 * D)

    got = tube_friction_closure("grooved").callable(ctx)
    using_corrcoeffs_thresholds = friction_corrugated_tube_vicente(
        Re, phi, Re_lo=CORR.Re_transition_lo, Re_hi=CORR.Re_transition_hi,
    )
    using_function_own_defaults = friction_corrugated_tube_vicente(Re, phi)

    assert got == using_corrcoeffs_thresholds
    assert got != using_function_own_defaults


def test_unknown_inside_tube_choice_raises():
    with pytest.raises(ValueError, match="inside_tube_choice"):
        tube_htc_closure("helical_insert")
    with pytest.raises(ValueError, match="inside_tube_choice"):
        tube_friction_closure("power_law")
