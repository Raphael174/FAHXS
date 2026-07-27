"""Confirms the pre-Phase-1 physics module paths still work as re-export shims.

Phase 1 of docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md moved liquid
coolant physics into ``physics/liquid_flow/`` and split the ideal-gas
governing equations into ``physics/gas_flow/``. The old module paths are kept
as thin shims (deleted in Phase 4) so any external/legacy code importing them
keeps working. This test is the one place that deliberately imports the old
paths and expects a DeprecationWarning; every other importer in this repo has
been updated to the new canonical paths.
"""
import warnings


def test_liquid_coolant_shim_reexports_and_warns():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        from hps_combustor.physics import liquid_coolant as shim
        from hps_combustor.physics.liquid_flow import chf, correlations

    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    assert shim.saturation_state is correlations.saturation_state
    assert shim.groeneveld_2006_chf is chf.groeneveld_2006_chf


def test_coolant_models_shim_reexports_and_warns():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        from hps_combustor.physics import coolant_models as shim
        from hps_combustor.physics.liquid_flow import dispatch

    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    assert shim.evaluate_coolant_closure is dispatch.evaluate_coolant_closure


def test_heated_liquid_channel_shim_reexports_and_warns():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        from hps_combustor.physics import heated_liquid_channel as shim
        from hps_combustor.physics.liquid_flow import governing_equations

    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    assert shim.solve_steady_heated_channel is governing_equations.solve_steady_heated_channel


def test_liquid_hx_adapters_shim_reexports_and_warns():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        from hps_combustor.physics import liquid_hx_adapters as shim
        from hps_combustor.physics.liquid_flow import hx_adapters

    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    assert (
        shim.solve_helical_coil_liquid_from_duty
        is hx_adapters.solve_helical_coil_liquid_from_duty
    )


def test_gas_governing_equations_shim_reexports_and_warns():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        from hps_combustor.physics import governing_equations as shim
        from hps_combustor.physics.gas_flow import governing_equations as gas_ge

    assert any(issubclass(w.category, DeprecationWarning) for w in caught)
    assert shim.dT__dx_IdealGas is gas_ge.dT__dx_IdealGas
