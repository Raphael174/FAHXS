"""Result packaging for steady and transient solver runs.

The archive is meant to be cheap to store and sufficient to re-plot later:
JSON metadata/input presets plus compressed numeric arrays. Images are not
written by default.
"""

from __future__ import annotations

import csv
import io
import json
import math
import re
import shutil
import zipfile
from dataclasses import asdict, is_dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np


def package_steady_run(solver, inputs: dict[str, Any], summary: dict[str, Any] | None = None):
    run = inputs["run"]
    run_dir, metadata = _make_run_dir(run.run_name, run.output_root, "steady")

    _write_json(run_dir / "metadata.json", metadata)
    _write_inputs(run_dir, inputs)
    _write_json(run_dir / "summary.json", summary or {})
    _write_steady_data(run_dir, solver, save_csv=run.save_csv)
    _write_steady_figures(run_dir, solver)

    if run.save_input_snapshot:
        _copy_input_data_snapshot(run_dir)

    return _zip_if_requested(run_dir, run.make_archive)


def package_transient_run(solver, inputs: dict[str, Any], summary: dict[str, Any] | None = None):
    run = inputs["run"]
    run_dir, metadata = _make_run_dir(run.run_name, run.output_root, "transient")
    time_series = get_transient_time_series(solver)
    final_summary = summary or _transient_summary(time_series)

    _write_json(run_dir / "metadata.json", metadata)
    _write_inputs(run_dir, inputs)
    _write_json(run_dir / "summary.json", final_summary)
    _write_transient_data(run_dir, time_series, save_csv=run.save_csv)
    _write_hx_performance_report(run_dir, time_series, inputs, final_summary)

    if run.save_input_snapshot:
        _copy_input_data_snapshot(run_dir)

    return _zip_if_requested(run_dir, run.make_archive)


def load_transient_time_series(path: str | Path):
    """Load a packaged transient time-series from a .npz, run folder, or zip file."""
    path = Path(path)
    if path.is_dir():
        return _load_transient_npz(path / "data" / "transient_timeseries.npz")
    if path.suffix.lower() == ".npz":
        return _load_transient_npz(path)
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            member = _find_zip_member(archive, "data/transient_timeseries.npz")
            if member is None:
                raise FileNotFoundError(
                    f"No data/transient_timeseries.npz found inside archive: {path}"
                )
            with archive.open(member) as handle:
                return _load_transient_npz_bytes(handle.read())
    raise ValueError(f"Expected a transient .npz file, run folder, or .zip archive: {path}")


def get_transient_time_series(solver):
    """Return or build a dashboard-friendly transient time-series dict."""
    if hasattr(solver, "time_series"):
        return solver.time_series

    if not hasattr(solver, "sol"):
        raise ValueError("Transient solver has no time_series and no solve_ivp result.")

    times = np.asarray(solver.sol.t, dtype=float)
    n_time = len(times)
    n_nodes = int(solver.N)
    x = np.arange(n_nodes) * float(solver.dx)
    tube_count = float(getattr(getattr(solver, "stp", None), "N_tubes", 1))

    first = solver.fluid_pass(solver.sol.y[:, 0], solver._bc_at(times[0]))
    field_names = [k for k, v in first.items() if _is_node_array(v, n_nodes)]
    scalar_names = [k for k, v in first.items() if not _is_node_array(v, n_nodes)]

    fields = {"Tbar": np.zeros((n_time, n_nodes))}
    fields.update({name: np.zeros((n_time, n_nodes)) for name in field_names})
    scalars = {name: np.zeros(n_time) for name in scalar_names}
    scalars.update({
        "Q_hot_kW": np.zeros(n_time),
        "Q_cold_kW": np.zeros(n_time),
        "mdot_c": np.zeros(n_time),
        "mdot_g": np.zeros(n_time),
    })

    for i, t in enumerate(times):
        wall = solver.sol.y[:, i]
        bc = solver._bc_at(float(t))
        res = solver.fluid_pass(wall, bc)

        fields["Tbar"][i] = wall
        for name in field_names:
            fields[name][i] = np.asarray(res[name], dtype=float)
        for name in scalar_names:
            scalars[name][i] = float(res[name])

        if "dq_hot__dx" in res:
            scalars["Q_hot_kW"][i] = np.sum(res["dq_hot__dx"]) * solver.dx * tube_count / 1e3
        if "dq_cold__dx" in res:
            scalars["Q_cold_kW"][i] = np.sum(res["dq_cold__dx"]) * solver.dx * tube_count / 1e3
        scalars["mdot_c"][i] = float(bc.get("mdot_c", math.nan))
        scalars["mdot_g"][i] = float(bc.get("mdot_g", math.nan))

    solver.time_series = {"t": times, "x": x, "fields": fields, "scalars": scalars}
    solver.fields = fields
    solver.scalars = scalars
    return solver.time_series


def _load_transient_npz(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Transient data file not found: {path}")

    with np.load(path) as data:
        return _time_series_from_npz(data)


def _load_transient_npz_bytes(payload: bytes):
    with np.load(io.BytesIO(payload)) as data:
        return _time_series_from_npz(data)


def _time_series_from_npz(data):
    time_series = {
        "t": np.asarray(data["t"], dtype=float),
        "x": np.asarray(data["x"], dtype=float),
        "fields": {},
        "scalars": {},
    }
    for name in data.files:
        if name.startswith("field_"):
            time_series["fields"][name.removeprefix("field_")] = np.asarray(data[name], dtype=float)
        elif name.startswith("scalar_"):
            time_series["scalars"][name.removeprefix("scalar_")] = np.asarray(data[name], dtype=float)
    return time_series


def _find_zip_member(archive: zipfile.ZipFile, suffix: str):
    normalized_suffix = suffix.replace("\\", "/")
    for member in archive.namelist():
        if member.replace("\\", "/").endswith(normalized_suffix):
            return member
    return None


def _make_run_dir(run_name: str, output_root: str, mode: str):
    now = datetime.now()
    safe_name = _safe_name(run_name)
    safe_stamp = now.strftime("%d-%m-%Y_%Hh%M")
    display_date = now.strftime("%d/%m/%Y")
    display_time = now.strftime("%Hh%M")

    root = Path(output_root)
    run_dir = _unique_run_dir(root / f"{safe_name}_{mode}_{safe_stamp}")
    run_dir.mkdir(parents=True, exist_ok=False)

    metadata = {
        "run_name": run_name,
        "mode": mode,
        "created_date": display_date,
        "created_time_24h": display_time,
        "folder_timestamp": safe_stamp,
    }
    return run_dir, metadata


def _write_inputs(run_dir: Path, inputs: dict[str, Any]):
    inputs_dir = run_dir / "inputs"
    inputs_dir.mkdir(exist_ok=True)
    _write_json(inputs_dir / "input_preset.json", inputs)


def _unique_run_dir(base: Path):
    if not base.exists():
        return base
    index = 2
    while True:
        candidate = base.with_name(f"{base.name}_{index:02d}")
        if not candidate.exists():
            return candidate
        index += 1


def _write_steady_data(run_dir: Path, solver, save_csv: bool):
    data_dir = run_dir / "data"
    data_dir.mkdir(exist_ok=True)

    arrays = {}
    if hasattr(solver, "data_master"):
        arrays = _arrays_from_dict(solver.data_master)
    elif hasattr(solver, "tube"):
        arrays = _arrays_from_dict({f"tube_{k}": v for k, v in solver.tube.items()})
        arrays["x"] = np.arange(solver.N) * solver.dx
        arrays["T_shell"] = np.asarray(solver.T_shell, dtype=float)

    if not arrays:
        raise ValueError("No steady numeric data found on solver.")

    np.savez_compressed(data_dir / "steady_data.npz", **arrays)
    if save_csv:
        _write_csv_table(data_dir / "steady_table.csv", arrays)


def _write_steady_figures(run_dir: Path, solver):
    """Render the ``HXDashboard`` themed figures for a solved steady solver.

    Both maintained steady solvers are covered, and both go through the same
    dashboard code (model_data_process/data_plotting.py):

      * shell-and-tube (duck-typed via ``solver.tube``/``solver.stp``) —
        through the ``data_master`` adapter in data_plotting_shellntube.py;
      * helical (duck-typed via ``solver.data_master``) — directly.

    Both saver entry points render under a temporarily-forced Agg backend and
    restore the caller's backend afterwards, so this can never block on a GUI
    window nor leave an interactive session headless.

    Never raises: a plotting failure must not lose an otherwise-successful
    numeric archive. On failure, writes a short note instead of the figures.
    """
    fig_dir = run_dir / "figures"
    try:
        if hasattr(solver, "tube") and hasattr(solver, "stp"):
            from .model_data_process.data_plotting_shellntube import save_shelltube_dashboard
            save_shelltube_dashboard(solver, fig_dir)
        elif getattr(solver, "data_master", None):
            from .model_data_process.data_plotting import save_dashboard
            save_dashboard(solver.data_master, fig_dir,
                           coolant_name=solver.coolantProp.coolant)
        else:
            return
    except Exception as exc:  # pragma: no cover - defensive, see docstring
        fig_dir.mkdir(exist_ok=True)
        (fig_dir / "PLOTTING_FAILED.txt").write_text(
            f"Automatic figure generation failed: {exc!r}\n"
            "The rest of this run's numeric archive is unaffected.",
            encoding="utf-8",
        )


def _write_transient_data(run_dir: Path, time_series: dict[str, Any], save_csv: bool):
    data_dir = run_dir / "data"
    data_dir.mkdir(exist_ok=True)

    arrays = {
        "t": np.asarray(time_series["t"], dtype=float),
        "x": np.asarray(time_series["x"], dtype=float),
    }
    arrays.update({f"field_{k}": np.asarray(v, dtype=float)
                   for k, v in time_series["fields"].items()})
    arrays.update({f"scalar_{k}": np.asarray(v, dtype=float)
                   for k, v in time_series["scalars"].items()})
    np.savez_compressed(data_dir / "transient_timeseries.npz", **arrays)

    if save_csv:
        scalar_arrays = {"t": arrays["t"]}
        scalar_arrays.update({k: v for k, v in arrays.items() if k.startswith("scalar_")})
        _write_csv_table(data_dir / "transient_scalars.csv", scalar_arrays)


def _write_hx_performance_report(
    run_dir: Path,
    time_series: dict[str, Any],
    inputs: dict[str, Any],
    summary: dict[str, Any],
):
    report = _build_hx_performance_report(time_series, inputs, summary)
    (run_dir / "HX_performance_summary.txt").write_text(report, encoding="utf-8")


def _write_csv_table(path: Path, arrays: dict[str, np.ndarray]):
    one_d = {k: np.asarray(v) for k, v in arrays.items() if np.asarray(v).ndim == 1}
    if not one_d:
        return
    length = max(len(v) for v in one_d.values())
    columns = {k: v for k, v in one_d.items() if len(v) == length}
    if not columns:
        return

    names = list(columns)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for i in range(length):
            writer.writerow([_format_float(columns[name][i]) for name in names])


def _interp_time_array(t_source: np.ndarray, values: np.ndarray, t_export: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 1:
        return np.interp(t_export, t_source, arr)
    if arr.ndim == 2:
        out = np.empty((len(t_export), arr.shape[1]), dtype=float)
        for j in range(arr.shape[1]):
            out[:, j] = np.interp(t_export, t_source, arr[:, j])
        return out
    return arr


def _derived_hx_metrics(
    time_series: dict[str, Any],
    inputs: dict[str, Any],
    *,
    t_export: np.ndarray,
) -> dict[str, np.ndarray]:
    t_source = np.asarray(time_series["t"], dtype=float)
    fields = time_series.get("fields", {})
    scalars = time_series.get("scalars", {})
    n = len(t_export)
    nan = np.full(n, np.nan)

    def scalar(name):
        if name not in scalars:
            return nan.copy()
        return _interp_time_array(t_source, np.asarray(scalars[name], dtype=float), t_export)

    def field(name):
        if name not in fields:
            return None
        return _interp_time_array(t_source, np.asarray(fields[name], dtype=float), t_export)

    Q_hot_kW = scalar("Q_hot_kW")
    T_g = field("T_g")
    T_c = field("T_c")
    cp_c = field("cp_c")
    h_g = field("h_g")
    h_c = field("h_c")
    T_g_out = scalar("T_g_out")
    T_c_out = scalar("T_c_out")
    mdot_c = scalar("mdot_c")

    flow_config = getattr(inputs.get("combustor"), "flow_config", "co")
    coolant_inlet_index = -1 if flow_config == "counter" else 0

    T_h_in = T_g[:, 0] if T_g is not None else nan.copy()
    if T_c is not None:
        T_c_in = T_c[:, coolant_inlet_index]
    else:
        T_c_in = nan.copy()
    if cp_c is not None:
        cp_c_ref = cp_c[:, coolant_inlet_index]
    else:
        cp_c_ref = nan.copy()

    dT1 = T_h_in - T_c_out
    dT2 = T_g_out - T_c_in
    if flow_config != "counter":
        dT1 = T_h_in - T_c_in
        dT2 = T_g_out - T_c_out
    lmtd = _lmtd_array(dT1, dT2)
    heat_duty_W = Q_hot_kW * 1.0e3
    UA_W_K = _safe_divide(heat_duty_W, lmtd)
    C_cold_W_K = mdot_c * cp_c_ref
    NTU_cold_ref = _safe_divide(UA_W_K, C_cold_W_K)
    q_max_cold_ref_W = C_cold_W_K * (T_h_in - T_c_in)
    effectiveness_cold_ref = _safe_divide(heat_duty_W, q_max_cold_ref_W)

    metrics = {
        "heat_duty_W": heat_duty_W,
        "LMTD_K": lmtd,
        "UA_W_K": UA_W_K,
        "NTU_cold_capacity_ref": NTU_cold_ref,
        "effectiveness_cold_capacity_ref": effectiveness_cold_ref,
        "T_hot_in_K": T_h_in,
        "T_hot_out_K": T_g_out,
        "T_cold_in_state_K": T_c_in,
        "T_cold_out_K": T_c_out,
        "C_cold_W_K": C_cold_W_K,
    }

    if h_g is not None:
        metrics["h_g_mean_W_m2K"] = np.nanmean(h_g, axis=1)
    if h_c is not None:
        metrics["h_c_mean_W_m2K"] = np.nanmean(h_c, axis=1)
    return metrics


def _lmtd_array(dT1: np.ndarray, dT2: np.ndarray) -> np.ndarray:
    dT1 = np.asarray(dT1, dtype=float)
    dT2 = np.asarray(dT2, dtype=float)
    out = np.full_like(dT1, np.nan, dtype=float)
    valid = (dT1 > 0.0) & (dT2 > 0.0) & np.isfinite(dT1) & np.isfinite(dT2)
    close = valid & (np.abs(dT1 - dT2) < 1e-9)
    out[close] = 0.5 * (dT1[close] + dT2[close])
    normal = valid & ~close
    out[normal] = (dT1[normal] - dT2[normal]) / np.log(dT1[normal] / dT2[normal])
    return out


def _safe_divide(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full(np.broadcast_shapes(num.shape, den.shape), np.nan, dtype=float)
    return np.divide(num, den, out=out, where=np.isfinite(den) & (np.abs(den) > 1e-30))


def _build_hx_performance_report(
    time_series: dict[str, Any],
    inputs: dict[str, Any],
    summary: dict[str, Any],
) -> str:
    t = np.asarray(time_series["t"], dtype=float)
    derived = _derived_hx_metrics(time_series, inputs, t_export=t)
    scalars = time_series.get("scalars", {})
    fields = time_series.get("fields", {})

    sample_indices = _report_sample_indices(t)
    lines = [
        "Combustor-HX transient performance summary",
        "=" * 48,
        "",
        "Scope:",
        "- Generated from packaged transient numeric time-series.",
        "- Raw transient fields are saved once in data/transient_timeseries.npz.",
        "- Wall radial/axial transient temperatures are available there as T_wg, Tbar, and T_wc when produced by the solver.",
        "- Axial derivative fields from data_processing.py are intentionally omitted.",
        "",
        "Run:",
        f"- HX configuration: {getattr(inputs.get('combustor'), 'HX_config', 'unknown')}",
        f"- Flow configuration: {getattr(inputs.get('combustor'), 'flow_config', 'unknown')}",
        f"- Fluid model: {getattr(inputs.get('transient'), 'fluid_model', 'unknown')}",
        f"- Chemistry: {getattr(inputs.get('transient'), 'chemistry_transient', 'unknown')}",
        f"- Time span: {float(t[0]):.6g} to {float(t[-1]):.6g} s",
        f"- Stored time points: {len(t)}",
        f"- Axial nodes: {len(np.asarray(time_series['x']))}",
        "",
        "Global extrema / averages:",
    ]
    for label, values, unit in (
        ("Heat duty Q_hot", scalars.get("Q_hot_kW"), "kW"),
        ("Cold-side heat to helium", scalars.get("Q_cold_kW"), "kW"),
        ("HX effectiveness, cold-capacity reference", derived.get("effectiveness_cold_capacity_ref"), "-"),
        ("NTU, cold-capacity reference", derived.get("NTU_cold_capacity_ref"), "-"),
        ("LMTD", derived.get("LMTD_K"), "K"),
        ("UA", derived.get("UA_W_K"), "W/K"),
        ("Max wall temperature", scalars.get("T_wall_max"), "K"),
        ("Min wall temperature", scalars.get("T_wall_min"), "K"),
        ("Helium outlet temperature", scalars.get("T_c_out"), "K"),
        ("Hot outlet temperature", scalars.get("T_g_out"), "K"),
        ("Gas pressure drop", scalars.get("dp_g_total_Pa"), "Pa"),
        ("Shell pressure drop estimate", scalars.get("dp_shell_total_Pa"), "Pa"),
        ("Gas Reynolds max", scalars.get("Re_g_max"), "-"),
        ("Shell Reynolds max", scalars.get("Re_shell_max"), "-"),
        ("FPV progress outlet", scalars.get("progress_g_out"), "-"),
        ("Energy residual", scalars.get("energy_residual_J"), "J"),
    ):
        lines.append(_stat_line(label, values, unit))

    if "T_wg" in fields and "T_wc" in fields:
        wall_delta = np.asarray(fields["T_wg"], dtype=float) - np.asarray(fields["T_wc"], dtype=float)
        lines.append(_stat_line("Wall radial delta T, T_wg - T_wc", wall_delta, "K"))

    lines.extend(["", "Selected transient points:", _sample_header()])
    for idx in sample_indices:
        lines.append(_sample_line(t[idx], scalars, derived, idx))

    lines.extend([
        "",
        "Notes:",
        "- Effectiveness and NTU are coolant-capacity referenced because transient hot-side heat capacity is not fully reconstructed in the packaged data.",
        "- LMTD uses terminal temperatures and is reported only when terminal temperature differences are positive.",
        "- Mechanical strain/stress extrema are reported when corresponding scalar histories are present; transient_core currently does not compute them.",
    ])
    mechanical = _mechanical_lines(scalars)
    if mechanical:
        lines.extend(["", "Mechanical histories present:"])
        lines.extend(mechanical)
    return "\n".join(lines) + "\n"


def _report_sample_indices(t: np.ndarray) -> list[int]:
    if len(t) <= 6:
        return list(range(len(t)))
    targets = np.linspace(float(t[0]), float(t[-1]), 6)
    return sorted(set(int(np.argmin(np.abs(t - target))) for target in targets))


def _sample_header() -> str:
    return (
        "time_s | Q_hot_kW | eff_cold_ref | NTU_cold_ref | LMTD_K | "
        "T_wall_max_K | T_c_out_K | T_g_out_K"
    )


def _sample_line(t, scalars, derived, idx: int) -> str:
    def at(source, name):
        values = source.get(name)
        if values is None:
            return np.nan
        arr = np.asarray(values, dtype=float)
        return arr[idx] if idx < arr.size else np.nan

    return (
        f"{t:.6g} | "
        f"{at(scalars, 'Q_hot_kW'):.6g} | "
        f"{at(derived, 'effectiveness_cold_capacity_ref'):.6g} | "
        f"{at(derived, 'NTU_cold_capacity_ref'):.6g} | "
        f"{at(derived, 'LMTD_K'):.6g} | "
        f"{at(scalars, 'T_wall_max'):.6g} | "
        f"{at(scalars, 'T_c_out'):.6g} | "
        f"{at(scalars, 'T_g_out'):.6g}"
    )


def _stat_line(label: str, values, unit: str) -> str:
    if values is None:
        return f"- {label}: not available"
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return f"- {label}: not available"
    return (
        f"- {label}: final={arr[-1]:.6g} {unit}, "
        f"min={np.min(arr):.6g} {unit}, max={np.max(arr):.6g} {unit}, "
        f"mean={np.mean(arr):.6g} {unit}"
    )


def _mechanical_lines(scalars: dict[str, Any]) -> list[str]:
    names = [name for name in scalars if "strain" in name.lower() or "stress" in name.lower()]
    return [_stat_line(name, scalars[name], "") for name in sorted(names)]


def _copy_input_data_snapshot(run_dir: Path):
    source = Path(__file__).with_name("input_data.py")
    if source.exists():
        shutil.copy2(source, run_dir / "inputs" / "input_data_snapshot.py")


def _zip_if_requested(run_dir: Path, make_archive: bool):
    archive = None
    if make_archive:
        archive = shutil.make_archive(
            str(run_dir),
            "zip",
            root_dir=run_dir.parent,
            base_dir=run_dir.name,
        )
    return {"folder": str(run_dir), "archive": archive}


def _write_json(path: Path, payload: Any):
    path.write_text(json.dumps(_jsonable(payload), indent=2), encoding="utf-8")


def _arrays_from_dict(data: dict[str, Any]):
    arrays = {}
    for key, value in data.items():
        if value is None:
            continue
        try:
            arr = np.asarray(value, dtype=float)
        except (TypeError, ValueError):
            continue
        if arr.size:
            arrays[key] = arr
    return arrays


def _transient_summary(time_series: dict[str, Any]):
    scalars = time_series.get("scalars", {})
    summary = {}
    for name, values in scalars.items():
        arr = np.asarray(values, dtype=float)
        if arr.size:
            summary[f"{name}_final"] = float(arr[-1])
            summary[f"{name}_max"] = float(np.nanmax(arr))
            summary[f"{name}_min"] = float(np.nanmin(arr))
    return summary


def _is_node_array(value: Any, n_nodes: int):
    try:
        arr = np.asarray(value)
    except Exception:
        return False
    return arr.ndim == 1 and len(arr) == n_nodes


def _safe_name(name: str):
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip())
    return name.strip("_") or "run"


def _format_float(value: Any):
    try:
        return f"{float(value):.10g}"
    except (TypeError, ValueError):
        return value


def _jsonable(value: Any):
    if is_dataclass(value):
        return _jsonable(asdict(value))
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, float):
        return value if math.isfinite(value) else str(value)
    return value
