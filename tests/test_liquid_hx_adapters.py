import numpy as np
import pytest
from types import SimpleNamespace

from hps_combustor.input_data import combustorProp, coolantProp, numericalProp, shellTubeProp
from hps_combustor.main_solve import main_solver
from hps_combustor.main_solve_shellntube import shellntube_solver
from hps_combustor.physics.liquid_flow.hx_adapters import (
    solve_helical_coil_liquid_from_data_master,
    solve_helical_coil_liquid_from_duty,
    solve_shelltube_shellside_liquid_from_duty,
    solve_shelltube_shellside_liquid_from_tube_result,
)
from hps_combustor.validation.liquid_solver_postprocess_audit import run_audit


def _water_coolant():
    return coolantProp(
        coolant="Water",
        coolant_model="equilibrium_liquid",
        mass_flow_c=0.12,
        T_in=420.0,
        p_in=1.0e6,
    )


def test_shelltube_shellside_liquid_adapter_coflow_maps_total_duty():
    stp = shellTubeProp()
    coolant = _water_coolant()
    dQ_per_tube = np.linspace(4.0, 7.0, 12)

    result = solve_shelltube_shellside_liquid_from_duty(
        coolant_prop=coolant,
        shelltube_prop=stp,
        dQ_profile_per_tube_W=dQ_per_tube,
        coolant_enters_at="z_min",
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    expected_heat = float(np.sum(dQ_per_tube) * stp.N_tubes)
    assert result.diagnostics.heat_rate_W == pytest.approx(expected_heat)
    assert result.diagnostics.energy_residual_ok is True
    assert result.diagnostics.pressure_drop_Pa > 0.0
    assert result.node_fields_hx_order["T_K"][-1] > result.node_fields_hx_order["T_K"][0]
    assert result.cell_fields_hx_order["heat_flux_W_m2"].shape == dQ_per_tube.shape


def test_shelltube_shellside_liquid_adapter_counterflow_keeps_hx_order():
    stp = shellTubeProp()
    coolant = _water_coolant()
    dQ_per_tube = np.linspace(4.0, 7.0, 12)

    result = solve_shelltube_shellside_liquid_from_duty(
        coolant_prop=coolant,
        shelltube_prop=stp,
        dQ_profile_per_tube_W=dQ_per_tube,
        coolant_enters_at="z_max",
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert result.coolant_enters_at == "z_max"
    assert result.node_fields_hx_order["z_m"][0] == pytest.approx(0.0)
    assert result.node_fields_hx_order["z_m"][-1] == pytest.approx(stp.L_tube)
    assert result.node_fields_hx_order["T_K"][-1] == pytest.approx(coolant.T_in)
    assert result.node_fields_hx_order["T_K"][0] > result.node_fields_hx_order["T_K"][-1]
    assert result.diagnostics.heat_rate_W == pytest.approx(float(np.sum(dQ_per_tube) * stp.N_tubes))


def test_shelltube_shellside_liquid_adapter_rejects_bad_duty_shape():
    with pytest.raises(ValueError, match="non-empty 1D"):
        solve_shelltube_shellside_liquid_from_duty(
            coolant_prop=_water_coolant(),
            shelltube_prop=shellTubeProp(),
            dQ_profile_per_tube_W=np.zeros((2, 2)),
            coolant_enters_at="z_min",
        )


def test_shelltube_shellside_liquid_adapter_accepts_solver_tube_result():
    stp = shellTubeProp()
    coolant = _water_coolant()
    dQ_per_tube = np.linspace(4.0, 7.0, 12)

    result = solve_shelltube_shellside_liquid_from_tube_result(
        coolant_prop=coolant,
        shelltube_prop=stp,
        tube_result={"dQ": dQ_per_tube},
        coolant_enters_at="z_min",
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert result.diagnostics.heat_rate_W == pytest.approx(float(np.sum(dQ_per_tube) * stp.N_tubes))
    with pytest.raises(KeyError, match="dQ"):
        solve_shelltube_shellside_liquid_from_tube_result(
            coolant_prop=coolant,
            shelltube_prop=stp,
            tube_result={},
            coolant_enters_at="z_min",
        )


def test_helical_coil_liquid_adapter_coflow_uses_total_parallel_area_and_duty():
    coolant = _water_coolant()
    combustor = combustorProp(flow_config="co", N_coils=2, Dh_coil=0.006)
    z_edges = np.linspace(0.0, 0.8, 9)
    dQ = np.linspace(80.0, 120.0, 8)

    result = solve_helical_coil_liquid_from_duty(
        coolant_prop=coolant,
        combustor_prop=combustor,
        numerical_prop=numericalProp(),
        dQ_profile_W=dQ,
        z_edges_m=z_edges,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert result.coolant_enters_at == "z_min"
    assert result.diagnostics.heat_rate_W == pytest.approx(float(np.sum(dQ)))
    assert result.diagnostics.energy_residual_ok is True
    assert result.node_fields_hx_order["T_K"][-1] > result.node_fields_hx_order["T_K"][0]
    assert result.diagnostics.pressure_drop_Pa > 0.0


def test_helical_coil_liquid_adapter_counterflow_defaults_from_combustor_flow_config():
    coolant = _water_coolant()
    combustor = combustorProp(flow_config="counter", N_coils=1, Dh_coil=0.006)
    z_edges = np.linspace(0.0, 0.8, 9)
    dQ = np.linspace(80.0, 120.0, 8)

    result = solve_helical_coil_liquid_from_duty(
        coolant_prop=coolant,
        combustor_prop=combustor,
        numerical_prop=numericalProp(),
        dQ_profile_W=dQ,
        z_edges_m=z_edges,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert result.coolant_enters_at == "z_max"
    assert result.node_fields_hx_order["z_m"][0] == pytest.approx(0.0)
    assert result.node_fields_hx_order["z_m"][-1] == pytest.approx(0.8)
    assert result.node_fields_hx_order["T_K"][-1] == pytest.approx(coolant.T_in)
    assert result.node_fields_hx_order["T_K"][0] > result.node_fields_hx_order["T_K"][-1]
    assert result.diagnostics.heat_rate_W == pytest.approx(float(np.sum(dQ)))


def test_helical_coil_liquid_adapter_rejects_bad_grid_length():
    with pytest.raises(ValueError, match="z_edges_m"):
        solve_helical_coil_liquid_from_duty(
            coolant_prop=_water_coolant(),
            combustor_prop=combustorProp(flow_config="co"),
            numerical_prop=numericalProp(),
            dQ_profile_W=np.ones(4),
            z_edges_m=np.linspace(0.0, 1.0, 4),
        )


def test_helical_coil_liquid_adapter_accepts_data_master_output():
    coolant = _water_coolant()
    combustor = combustorProp(flow_config="co", N_coils=1, Dh_coil=0.006)
    dQ = np.linspace(80.0, 120.0, 8)
    data_master = {
        "dQ": dQ,
        "L_ch": np.linspace(0.0, 0.7, dQ.size),
    }

    result = solve_helical_coil_liquid_from_data_master(
        coolant_prop=coolant,
        combustor_prop=combustor,
        numerical_prop=numericalProp(),
        data_master=data_master,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert result.diagnostics.heat_rate_W == pytest.approx(float(np.sum(dQ)))
    assert result.node_fields_hx_order["z_m"][0] == pytest.approx(0.0)
    assert result.node_fields_hx_order["z_m"][-1] == pytest.approx(0.8)
    with pytest.raises(KeyError, match="dQ"):
        solve_helical_coil_liquid_from_data_master(
            coolant_prop=coolant,
            combustor_prop=combustor,
            numerical_prop=numericalProp(),
            data_master={},
        )


def test_helical_solver_postprocess_method_uses_data_master_output():
    coolant = _water_coolant()
    combustor = combustorProp(flow_config="co", N_coils=1, Dh_coil=0.006)
    dQ = np.linspace(80.0, 120.0, 8)
    solver_like = SimpleNamespace(
        coolantProp=coolant,
        combustorProp=combustor,
        numericalProp=numericalProp(),
        data_master={"dQ": dQ, "L_ch": np.linspace(0.0, 0.7, dQ.size)},
    )

    result = main_solver.liquid_coolant_postprocess(
        solver_like,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert solver_like.liquid_coolant is result
    assert result.diagnostics.heat_rate_W == pytest.approx(float(np.sum(dQ)))
    assert result.diagnostics.energy_residual_ok is True


def test_shelltube_solver_postprocess_method_uses_tube_output():
    stp = shellTubeProp()
    coolant = _water_coolant()
    dQ_per_tube = np.linspace(4.0, 7.0, 12)
    solver_like = SimpleNamespace(
        coolantProp=coolant,
        stp=stp,
        tube={"dQ": dQ_per_tube},
        flow_config="counter",
    )

    result = shellntube_solver.liquid_coolant_postprocess(
        solver_like,
        lut_path="docs/reference/external/2006LUTdata.txt",
    )

    assert solver_like.liquid_coolant is result
    assert result.coolant_enters_at == "z_max"
    assert result.diagnostics.heat_rate_W == pytest.approx(float(np.sum(dQ_per_tube) * stp.N_tubes))
    assert result.diagnostics.energy_residual_ok is True


def test_liquid_solver_postprocess_audit_generates_checks(tmp_path):
    report = run_audit(output=tmp_path / "audit.json")

    assert report["checks"]["all_passed"] is True
    assert report["helical"]["heat_rate_W"] > 0.0
    assert report["shelltube"]["heat_rate_W"] > 0.0
    assert report["helical"]["energy_residual_ok"] is True
    assert report["shelltube"]["energy_residual_ok"] is True
    assert (tmp_path / "audit.json").exists()
