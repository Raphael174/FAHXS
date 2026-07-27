from types import SimpleNamespace

import pytest

from hps_combustor.schedule_inputs import apply_schedule_file


def make_transient():
    return SimpleNamespace(
        schedule_mass_flow_c=None,
        schedule_p_c_in=None,
        schedule_T_c_in=None,
        schedule_T_lox_in=None,
        schedule_mass_flow_lox=None,
        schedule_mass_flow_diesel=None,
        schedule_ignition_state=None,
        schedule_mass_flow_g=None,
        schedule_OF=None,
        ignition_time=0.0,
    )


def make_hotgas():
    return SimpleNamespace(
        mass_flow_g=0.1,
        mixing_ratio=2.5,
        T_inj_LOX=120.0,
    )


def make_coolant():
    return SimpleNamespace(
        mass_flow_c=0.15,
        p_in=8.0e6,
        T_in=90.0,
    )


def test_schedule_csv_accepts_decimal_commas_with_semicolon_delimiter(tmp_path):
    path = tmp_path / "schedule.csv"
    path.write_text(
        "\n".join(
            [
                "time_s;helium_m_dot_kg_s;helium_T_in_K;helium_p_in_Pa;lox_m_dot_kg_s;diesel_m_dot_kg_s;ignition",
                "0;0;90,5;8,0E6;0,03;0,01;0",
                "0,0025;1,88E-07;91,25;7,95E6;0,04;0,02;1",
            ]
        ),
        encoding="utf-8",
    )

    transient = make_transient()
    hotgas = make_hotgas()
    coolant = make_coolant()

    apply_schedule_file(transient, hotgas, path, coolant=coolant)

    assert transient.schedule_mass_flow_c == ((0.0, 0.0), (0.0025, 1.88e-7))
    assert transient.schedule_T_c_in == ((0.0, 90.5), (0.0025, 91.25))
    assert transient.schedule_p_c_in == ((0.0, 8.0e6), (0.0025, 7.95e6))
    assert transient.schedule_mass_flow_g == ((0.0, 0.04), (0.0025, 0.06))
    assert transient.schedule_OF == ((0.0, 3.0), (0.0025, 2.0))
    assert transient.ignition_time == 0.0025
    assert coolant.mass_flow_c == pytest.approx(1.88e-7)
    assert hotgas.mass_flow_g == pytest.approx(0.06)


def test_schedule_excel_rows_accept_decimal_comma_strings(monkeypatch, tmp_path):
    pd = pytest.importorskip("pandas")

    class FakeFrame:
        def notna(self):
            return self

        def where(self, _condition, _replacement):
            return self

        def to_dict(self, orient):
            assert orient == "records"
            return [
                {"time_s": "0,0", "m_dot_kg_s": "0,15", "T_in_K": "90,0"},
                {"time_s": "1,5", "m_dot_kg_s": "0,05", "T_in_K": "95,5"},
            ]

    monkeypatch.setattr(pd, "read_excel", lambda path, sheet_name=None: {"helium": FakeFrame()})
    path = tmp_path / "schedule.xlsx"
    path.write_bytes(b"not a real workbook; pandas is monkeypatched")

    transient = make_transient()
    hotgas = make_hotgas()

    apply_schedule_file(transient, hotgas, path)

    assert transient.schedule_mass_flow_c == ((0.0, 0.15), (1.5, 0.05))
    assert transient.schedule_T_c_in == ((0.0, 90.0), (1.5, 95.5))
