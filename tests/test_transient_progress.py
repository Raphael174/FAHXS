from hps_combustor.transient_core.progress import TransientProgressPrinter


def test_transient_progress_prints_requested_fields(capsys):
    progress = TransientProgressPrinter(
        total_steps=4,
        enabled=True,
        interval_steps=2,
    )

    progress.update(
        step=1,
        time_s=0.25,
        T_wall=[300.0, 310.0],
        T_coolant_outlet=95.0,
        p_coolant_outlet=8.0e6,
        T_gas_outlet=1200.0,
    )
    assert capsys.readouterr().out == ""

    progress.update(
        step=2,
        time_s=0.50,
        T_wall=[305.0, 315.0],
        T_coolant_outlet=100.0,
        p_coolant_outlet=7.8e6,
        T_gas_outlet=1180.0,
    )
    out = capsys.readouterr().out
    assert "step/total" in out
    assert "Tmat min/max" in out
    assert "2/4" in out
    assert "305.0/315.0" in out
    assert "100.0" in out
    assert "78.000" in out
    assert "1180.0" in out


def test_transient_progress_can_be_disabled(capsys):
    progress = TransientProgressPrinter(total_steps=1, enabled=False)
    progress.update(
        step=1,
        time_s=1.0,
        T_wall=[300.0],
        T_coolant_outlet=90.0,
        p_coolant_outlet=8.0e6,
        T_gas_outlet=1000.0,
    )
    assert capsys.readouterr().out == ""


def test_transient_progress_accepts_custom_pressure_label(capsys):
    progress = TransientProgressPrinter(
        total_steps=1,
        enabled=True,
        pressure_label="dpHe [bar]",
    )
    progress.update(
        step=1,
        time_s=1.0,
        T_wall=[300.0],
        T_coolant_outlet=90.0,
        p_coolant_outlet=7.5e5,
        T_gas_outlet=1000.0,
    )
    out = capsys.readouterr().out
    assert "dpHe [bar]" in out
    assert "7.500" in out
