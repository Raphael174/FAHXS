"""Stage D, Slice 3 of docs/solver_design/FV_CORE_REWORK_PLAN.md.

`core/hotgas/combustor.py` relocates the three gas-state providers
(`fpv_gas_state_provider`, `equilibrium_gas_state_provider`,
`oxygen_gas_state_provider`) from `transient_core/adapters_shelltube.py`
unchanged, confirmed geometry-independent before moving. Proves
bit-identical equivalence (same objects) to the shim's re-exports, the same
discipline every prior relocation this stage used.
"""
from __future__ import annotations

import pytest

from hps_combustor.core.hotgas.combustor import (
    GasStateProvider,
    ShellTubeGasState,
    _coerce_gas_state,
    equilibrium_gas_state_provider,
    fpv_gas_state_provider,
    oxygen_gas_state_provider,
)
from hps_combustor.transient_core import adapters_shelltube as legacy_shim


def test_adapters_shelltube_shim_reexports_are_identical_objects():
    assert legacy_shim.ShellTubeGasState is ShellTubeGasState
    assert legacy_shim.GasStateProvider is GasStateProvider
    assert legacy_shim.fpv_gas_state_provider is fpv_gas_state_provider
    assert legacy_shim.equilibrium_gas_state_provider is equilibrium_gas_state_provider
    assert legacy_shim.oxygen_gas_state_provider is oxygen_gas_state_provider
    assert legacy_shim._coerce_gas_state is _coerce_gas_state


def test_equilibrium_gas_state_provider_wraps_manifold():
    class FakeManifold:
        def at(self, h_removed):
            return (1000.0 - h_removed / 1000.0, 2.0, 3e-5, 0.05, 1200.0, 0.1, 0.2)

    provider, initial_progress = equilibrium_gas_state_provider(FakeManifold())
    assert initial_progress == 0.0
    state = provider(5000.0, 0.0, 0)
    assert isinstance(state, ShellTubeGasState)
    assert state.T == pytest.approx(995.0)
    assert state.progress_source == 0.0


def test_oxygen_gas_state_provider_cools_with_enthalpy_removed():
    provider, initial_progress = oxygen_gas_state_provider(T_inlet=120.0, pressure=2e5)
    assert initial_progress == 0.0
    state0 = provider(0.0, 0.0, 0)
    state_cooled = provider(5000.0, 0.0, 0)
    assert isinstance(state0, ShellTubeGasState)
    # removing enthalpy must not raise the temperature
    assert state_cooled.T <= state0.T


def test_coerce_gas_state_accepts_dict_and_validates():
    state = _coerce_gas_state({"T": 300.0, "rho": 1.0, "mu": 2e-5, "k": 0.03, "cp": 1000.0})
    assert isinstance(state, ShellTubeGasState)
    assert state.T == 300.0

    with pytest.raises(ValueError, match="must be positive"):
        _coerce_gas_state({"T": -1.0, "rho": 1.0, "mu": 2e-5, "k": 0.03, "cp": 1000.0})

    with pytest.raises(TypeError):
        _coerce_gas_state("not a state")
