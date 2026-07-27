"""Deprecated location. Use ``hps_combustor.physics.liquid_flow.dispatch``
instead. Kept as a re-export shim during the Phase 1 restructure in
docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md.
"""

import warnings

from hps_combustor.physics.liquid_flow.dispatch import (
    CoolantClosureResult,
    CoolantState,
    coolant_inlet_state,
    coolant_state_from_Tp,
    coolant_state_from_ph,
    evaluate_coolant_closure,
)

warnings.warn(
    "hps_combustor.physics.coolant_models has moved to "
    "hps_combustor.physics.liquid_flow.dispatch; update imports.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = [
    "CoolantClosureResult",
    "CoolantState",
    "coolant_inlet_state",
    "coolant_state_from_Tp",
    "coolant_state_from_ph",
    "evaluate_coolant_closure",
]
