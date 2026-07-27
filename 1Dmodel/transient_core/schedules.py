"""Schedule helpers shared by transient-core adapters."""

from __future__ import annotations

import numpy as np


def interp_schedule(schedule, t: float, default: float) -> float:
    """Linearly interpolate a project-style schedule with flat end holds.

    `schedule` is `None` or an iterable of `(time_s, value)` pairs. The behavior
    matches the legacy transient solvers' private `_interp_schedule` helpers.
    """

    if schedule is None:
        return float(default)
    pts = _schedule_array(schedule)
    if pts.size == 0:
        return float(default)

    time = float(t)
    if time <= pts[0, 0]:
        return float(pts[0, 1])
    if time >= pts[-1, 0]:
        return float(pts[-1, 1])
    return float(np.interp(time, pts[:, 0], pts[:, 1]))


def schedule_times(*schedules, t_min: float | None = None, t_max: float | None = None) -> np.ndarray:
    """Return sorted unique time points from one or more schedules."""

    times = []
    for schedule in schedules:
        if not schedule:
            continue
        pts = _schedule_array(schedule)
        if pts.size:
            times.append(pts[:, 0])
    if not times:
        return np.array([], dtype=float)

    out = np.concatenate(times)
    if t_min is not None:
        out = out[out >= float(t_min)]
    if t_max is not None:
        out = out[out <= float(t_max)]
    return np.unique(np.round(out[np.isfinite(out)], decimals=12))


def collect_transient_schedule_times(transient, names, *, t_min=0.0, t_max=None) -> np.ndarray:
    """Collect breakpoint times from named schedule attributes on a transient config."""

    schedules = [getattr(transient, name, None) for name in names]
    return schedule_times(*schedules, t_min=t_min, t_max=t_max)


def _schedule_array(schedule) -> np.ndarray:
    pts = np.asarray(list(schedule), dtype=float)
    if pts.size == 0:
        return np.empty((0, 2), dtype=float)
    if pts.ndim != 2 or pts.shape[1] < 2:
        raise ValueError("schedule must be an iterable of (time, value) pairs")
    pts = pts[:, :2]
    if not np.all(np.isfinite(pts)):
        raise ValueError("schedule contains non-finite values")
    order = np.argsort(pts[:, 0], kind="mergesort")
    return pts[order]
