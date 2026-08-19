"""Registry adapter for `core` — Stage D, Slice 1
(docs/solver_design/FV_CORE_REWORK_PLAN.md section 3.3).

Currently covers shell-and-tube's tube-side gas closures, registered in
`physics/liquid_flow/gas_closures.py`. Selection there is FORCED by name
(`shellTubeProp.inside_tube_choice`), mirroring how the legacy adapter
branches -- not `select_supercritical`'s inferred ranking, which is for the
liquid/supercritical family where no single config field names "the"
correlation. `physics/liquid_flow/dispatch.py::evaluate_coolant_closure`
remains the entry point for that family; this module does not replace it.
"""

from __future__ import annotations

from hps_combustor.physics.liquid_flow import gas_closures  # noqa: F401 (registers on import)
from hps_combustor.physics.liquid_flow.registry import ClosureRecord, get_record

_TUBE_HTC_BY_CHOICE = {
    "smooth": "tube_straight_gnielinski_blended",
    "grooved": "tube_grooved_vicente",
}
_TUBE_FRICTION_BY_CHOICE = {
    "smooth": "tube_straight_friction_blended",
    "grooved": "tube_grooved_friction_vicente",
}


def tube_htc_closure(inside_tube_choice: str) -> ClosureRecord:
    """The registered Nu-type closure for shell-and-tube's tube-side hot gas,
    forced-selected by `shellTubeProp.inside_tube_choice`."""
    try:
        name = _TUBE_HTC_BY_CHOICE[inside_tube_choice]
    except KeyError:
        raise ValueError(
            f"no tube-side HTC closure registered for inside_tube_choice={inside_tube_choice!r} "
            f"-- expected one of {sorted(_TUBE_HTC_BY_CHOICE)}"
        ) from None
    return get_record(name)


def tube_friction_closure(inside_tube_choice: str) -> ClosureRecord:
    """The registered friction-type closure for shell-and-tube's tube-side
    hot gas, forced-selected by `shellTubeProp.inside_tube_choice`."""
    try:
        name = _TUBE_FRICTION_BY_CHOICE[inside_tube_choice]
    except KeyError:
        raise ValueError(
            f"no tube-side friction closure registered for inside_tube_choice={inside_tube_choice!r} "
            f"-- expected one of {sorted(_TUBE_FRICTION_BY_CHOICE)}"
        ) from None
    return get_record(name)
