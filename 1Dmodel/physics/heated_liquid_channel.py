"""Deprecated location. Use
``hps_combustor.physics.liquid_flow.governing_equations`` instead. Kept as a
re-export shim during the Phase 1 restructure in
docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md.
"""

import warnings

from hps_combustor.physics.liquid_flow.governing_equations import (
    HeatedChannelCase,
    HeatedChannelDiagnostics,
    HeatedChannelProfileCase,
    HeatedChannelResult,
    HXGridHeatedChannelResult,
    heat_flux_profile,
    heated_channel_cell_fields,
    heated_channel_node_fields,
    inlet_enthalpy,
    solve_steady_heated_channel,
    solve_steady_heated_channel_on_hx_grid,
    solve_steady_heated_channel_profile,
    summarize_heated_channel_result,
    validate_heated_channel_case,
)

warnings.warn(
    "hps_combustor.physics.heated_liquid_channel has moved to "
    "hps_combustor.physics.liquid_flow.governing_equations; update imports.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "HeatedChannelCase",
    "HeatedChannelDiagnostics",
    "HeatedChannelProfileCase",
    "HeatedChannelResult",
    "HXGridHeatedChannelResult",
    "heat_flux_profile",
    "heated_channel_cell_fields",
    "heated_channel_node_fields",
    "inlet_enthalpy",
    "solve_steady_heated_channel",
    "solve_steady_heated_channel_on_hx_grid",
    "solve_steady_heated_channel_profile",
    "summarize_heated_channel_result",
    "validate_heated_channel_case",
]
