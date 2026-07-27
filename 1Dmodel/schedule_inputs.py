"""Load transient boundary-condition schedules from CSV or Excel files.

Excel may use two sheets:
  helium:      time_s, m_dot_kg_s, p_in_Pa, T_in_K
  propellants: time_s, m_dot_lox_kg_s, m_dot_diesel_kg_s

CSV or single-sheet Excel may also use descriptive column names:
  time_s, helium_m_dot_kg_s, helium_T_in_K, helium_p_in_Pa,
  lox_m_dot_kg_s, diesel_m_dot_kg_s, lox_T_in_K, ignition

Direct hot-gas columns `m_dot_g_kg_s` and `OF` may be used instead of
LOX/diesel columns.

Numeric cells may use decimal dots or decimal commas. CSV schedules with
decimal commas should use semicolon or tab field delimiters.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


def apply_schedule_file(transient, hotgas, path, coolant=None):
    if not path:
        return transient

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Schedule file not found: {path}")

    if path.suffix.lower() in (".xlsx", ".xls"):
        tables = _read_excel(path)
        if "helium" not in tables and "propellants" not in tables and tables:
            first_rows = next(iter(tables.values()))
            tables = {"helium": first_rows, "propellants": first_rows}
    else:
        rows = _read_csv(path)
        tables = {"helium": rows, "propellants": rows}

    helium = tables.get("helium", [])
    propellants = tables.get("propellants", [])

    transient.schedule_mass_flow_c = _schedule_any(
        helium, "time_s", ("m_dot_kg_s", "helium_m_dot_kg_s", "he_m_dot_kg_s")
    )
    transient.schedule_p_c_in = _schedule_any(
        helium, "time_s", ("p_in_Pa", "helium_p_in_Pa", "he_p_in_Pa")
    )
    transient.schedule_p_c_out = _schedule_any(
        helium, "time_s", ("p_out_Pa", "helium_p_out_Pa", "he_p_out_Pa")
    )
    transient.schedule_T_c_in = _schedule_any(
        helium, "time_s", ("T_in_K", "helium_T_in_K", "he_T_in_K")
    )

    gas_direct = _schedule_any(propellants, "time_s", ("m_dot_g_kg_s", "hotgas_m_dot_kg_s"))
    of_direct = _schedule_any(propellants, "time_s", ("OF", "of", "o_f"))
    transient.schedule_T_lox_in = _schedule_any(
        propellants, "time_s", ("lox_T_in_K", "T_lox_in_K", "lox_temperature_K")
    )
    transient.schedule_mass_flow_lox = _schedule_any(
        propellants, "time_s", ("m_dot_lox_kg_s", "lox_m_dot_kg_s")
    )
    transient.schedule_mass_flow_diesel = _schedule_any(
        propellants, "time_s", ("m_dot_diesel_kg_s", "diesel_m_dot_kg_s")
    )
    transient.schedule_ignition_state = _schedule_any(
        propellants, "time_s", ("ignition", "ignited", "combustion_on")
    )

    if gas_direct is not None:
        transient.schedule_mass_flow_g = gas_direct
    else:
        transient.schedule_mass_flow_g = _propellant_total_schedule(propellants)

    if of_direct is not None:
        transient.schedule_OF = of_direct
    else:
        transient.schedule_OF = _of_schedule(propellants)

    ignition_time = _first_enabled_time(transient.schedule_ignition_state)
    if ignition_time is not None:
        transient.ignition_time = ignition_time

    # Keep steady defaults aligned with the final scheduled operating point.
    # The transient schedule still controls startup/chilldown; these defaults are
    # mainly used for chemistry setup and steady-reference comparisons.
    if transient.schedule_mass_flow_g:
        hotgas.mass_flow_g = _last_positive_value(transient.schedule_mass_flow_g)
    if transient.schedule_OF:
        hotgas.mixing_ratio = transient.schedule_OF[-1][1]
    if transient.schedule_T_lox_in:
        hotgas.T_inj_LOX = transient.schedule_T_lox_in[-1][1]
    if coolant is not None:
        if transient.schedule_mass_flow_c:
            coolant.mass_flow_c = transient.schedule_mass_flow_c[-1][1]
        if transient.schedule_p_c_in:
            coolant.p_in = transient.schedule_p_c_in[-1][1]
        if transient.schedule_T_c_in:
            coolant.T_in = transient.schedule_T_c_in[-1][1]

    return transient


def _read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
        except csv.Error:
            dialect = csv.excel
        return list(csv.DictReader(handle, dialect=dialect))


def _read_excel(path: Path):
    try:
        import pandas as pd
    except ImportError as exc:
        raise RuntimeError(
            "Excel schedules require pandas and openpyxl. Use CSV, or install "
            "those packages in the project environment."
        ) from exc

    sheets = pd.read_excel(path, sheet_name=None)
    return {
        name.lower(): frame.where(frame.notna(), None).to_dict(orient="records")
        for name, frame in sheets.items()
    }


def _schedule(rows, time_key, value_key):
    pairs = []
    for row in rows:
        if value_key not in row or _is_blank(row.get(value_key)):
            continue
        time_value = row.get(time_key)
        if _is_blank(time_value):
            continue
        pairs.append((_to_float(time_value), _to_float(row[value_key])))
    return tuple(pairs) if pairs else None


def _schedule_any(rows, time_key, value_keys):
    for value_key in value_keys:
        schedule = _schedule(rows, time_key, value_key)
        if schedule is not None:
            return schedule
    return None


def _propellant_total_schedule(rows):
    pairs = []
    for row in rows:
        time_value = row.get("time_s")
        lox = _row_value(row, ("m_dot_lox_kg_s", "lox_m_dot_kg_s"))
        diesel = _row_value(row, ("m_dot_diesel_kg_s", "diesel_m_dot_kg_s"))
        if _is_blank(time_value) or _is_blank(lox) or _is_blank(diesel):
            continue
        pairs.append((_to_float(time_value), _to_float(lox) + _to_float(diesel)))
    return tuple(pairs) if pairs else None


def _of_schedule(rows):
    pairs = []
    for row in rows:
        time_value = row.get("time_s")
        lox = _row_value(row, ("m_dot_lox_kg_s", "lox_m_dot_kg_s"))
        diesel = _row_value(row, ("m_dot_diesel_kg_s", "diesel_m_dot_kg_s"))
        if _is_blank(time_value) or _is_blank(lox) or _is_blank(diesel):
            continue
        diesel = _to_float(diesel)
        if diesel <= 0:
            continue
        pairs.append((_to_float(time_value), _to_float(lox) / diesel))
    return tuple(pairs) if pairs else None


def _row_value(row, keys):
    for key in keys:
        value = row.get(key)
        if not _is_blank(value):
            return value
    return None


def _last_positive_value(schedule):
    for _, value in reversed(schedule):
        if value > 0:
            return value
    return schedule[-1][1]


def _first_enabled_time(schedule):
    if not schedule:
        return None
    for time_value, state in schedule:
        if state >= 0.5:
            return time_value
    return None


def _to_float(value):
    """Parse numeric schedule values, accepting decimal commas from Excel/CSV."""

    if isinstance(value, (int, float)):
        if isinstance(value, float) and math.isnan(value):
            raise ValueError("NaN numeric schedule value")
        return float(value)
    text = str(value).strip().replace("\u00a0", "").replace(" ", "")
    if not text:
        raise ValueError("empty numeric schedule value")
    if "," in text and "." in text:
        comma = text.rfind(",")
        dot = text.rfind(".")
        if comma > dot:
            # European thousands + decimal: 1.234,56
            text = text.replace(".", "").replace(",", ".")
        else:
            # US/UK thousands + decimal: 1,234.56
            text = text.replace(",", "")
    elif "," in text:
        # Common European decimal format: 0,0025 or 1,88E-07.
        text = text.replace(",", ".")
    return float(text)


def _is_blank(value) -> bool:
    if value in (None, ""):
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False
