"""Terminal progress reporting for transient simulations."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass
class TransientProgressPrinter:
    """Small throttled terminal reporter for long transient runs."""

    total_steps: int
    enabled: bool = True
    interval_steps: int | None = None
    interval_time_s: float | None = None
    pressure_label: str = "pHe,out [bar]"

    def __post_init__(self):
        self.total_steps = max(int(self.total_steps), 0)
        if self.interval_steps is None:
            self.interval_steps = max(1, math.ceil(max(self.total_steps, 1) / 20))
        else:
            self.interval_steps = max(1, int(self.interval_steps))
        if self.interval_time_s is not None and self.interval_time_s <= 0.0:
            self.interval_time_s = None
        self._next_step = self.interval_steps
        self._next_time = self.interval_time_s
        self._header_printed = False

    @classmethod
    def from_config(cls, config, *, total_steps: int):
        """Build from a transientProp-like object."""

        return cls(
            total_steps=total_steps,
            enabled=bool(getattr(config, "progress_print", True)),
            interval_steps=getattr(config, "progress_interval_steps", None),
            interval_time_s=getattr(config, "progress_interval_time_s", None),
        )

    def update(
        self,
        *,
        step: int,
        time_s: float,
        T_wall,
        T_coolant_outlet: float,
        p_coolant_outlet: float | None,
        T_gas_outlet: float | None,
    ) -> None:
        """Print one progress line if the step/time throttle says to."""

        if not self.enabled or self.total_steps <= 0:
            return

        step = int(step)
        time_s = float(time_s)
        due_by_step = step >= self._next_step
        due_by_time = self._next_time is not None and time_s >= self._next_time
        final_step = step >= self.total_steps
        if not (due_by_step or due_by_time or final_step):
            return

        if not self._header_printed:
            print(
                "  transient progress:"
                " step/total | time [s] | Tmat min/max [K] |"
                f" THe,out [K] | {self.pressure_label} | Tgas,out [K]"
            )
            self._header_printed = True

        wall = np.asarray(T_wall, dtype=float)
        wall_min = _finite_stat(wall, np.nanmin)
        wall_max = _finite_stat(wall, np.nanmax)
        p_bar = np.nan if p_coolant_outlet is None else float(p_coolant_outlet) / 1.0e5
        Tg = np.nan if T_gas_outlet is None else float(T_gas_outlet)
        print(
            f"  {step:5d}/{self.total_steps:<5d} |"
            f" {time_s:8.3f} |"
            f" {wall_min:8.1f}/{wall_max:<8.1f} |"
            f" {float(T_coolant_outlet):10.1f} |"
            f" {p_bar:12.3f} |"
            f" {Tg:11.1f}"
        )

        while self._next_step <= step:
            self._next_step += self.interval_steps
        if self._next_time is not None and self.interval_time_s is not None:
            while self._next_time <= time_s:
                self._next_time += self.interval_time_s


def _finite_stat(values: np.ndarray, func) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    return float(func(finite))
