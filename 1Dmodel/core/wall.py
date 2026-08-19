"""Wall conduction: radial resistance stack + analytic fin/rib — Stage C.

See docs/solver_design/FV_CORE_REWORK_PLAN.md section 3.4. Generalizes
``physics/heat_conduction.py``'s
``OneDimensionalSteadyConduction_ShellnHelicalTube`` while preserving the two
things that must not regress:

1. **``hot_side`` orientation.** ``"outer"`` = hot fluid outside the tube
   (helical coil config); ``"inner"`` = hot fluid inside the tube
   (shell-and-tube, where combustion gas is in the tubes). Getting this
   backwards silently swaps which perimeter each flux uses. CLAUDE.md flags
   the corrected mapping as a do-not-regress item.
2. **The quadratic quasi-static radial profile** (``a2 = s/2k``,
   ``a6 = s/6k``) and its closed-form 2x2 face solve, validated to <2 K
   against a fully resolved radial PDE (``docs/solver_design/
   check_wall_quasi_static_validity.py``). Reproduced here exactly.

What is NEW in this module is the **rib/fin path** for channel geometries
(rocket-nozzle regen channels, multi-start helical wraps): the land between
adjacent channels acts as a fin, so the cold side's effective heated
perimeter is not simply the channel's wetted perimeter. Without this, peak
hot-gas-wall temperature is under-predicted on ribbed high-aspect-ratio
channels -- the quantity that actually sizes a regen wall.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class FinResult:
    """Analytic fin solution for the land between adjacent channels."""

    efficiency: float          # eta_fin, adiabatic-tip
    m_1_m: float               # fin parameter [1/m]
    P_effective_m: float       # effective cold-side heated perimeter per channel


def rectangular_fin_efficiency(
    *, h_c: float, k_w: float, t_rib: float, h_channel: float
) -> tuple[float, float]:
    """Adiabatic-tip efficiency of a straight rectangular fin.

    ``m = sqrt(2*h_c/(k_w*t_rib))``, ``eta = tanh(m*L)/(m*L)``.

    The factor of **2** is deliberate and is the correct choice for a rib
    BETWEEN TWO CHANNELS: coolant wets both faces of the land, so the
    perimeter-to-area ratio of the fin cross-section is ``2/t_rib``, not
    ``1/t_rib``. ``physics/heat_conduction.py`` contains both conventions --
    ``compute_fin_efficiency_ch`` (factor 1, for a fin wetted on one side)
    and ``compute_eta_fin_rectangular`` (factor 2). This module uses the
    factor-2 form, matching ``compute_eta_fin_rectangular``; the design doc
    flags checking this convention as an explicit Stage C task.
    """
    if k_w <= 0.0 or t_rib <= 0.0:
        raise ValueError("k_w and t_rib must be positive")
    if h_c <= 0.0:
        # No convection -> a fin transfers nothing extra, but is not "0%
        # efficient" in the tanh sense; the limit eta->1 with zero driving
        # potential. Return 1.0 so P_effective degenerates gracefully.
        return 1.0, 0.0
    m = np.sqrt(2.0 * h_c / (k_w * t_rib))
    mL = m * h_channel
    if mL < 1e-12:
        return 1.0, float(m)
    return float(np.tanh(mL) / mL), float(m)


def channel_effective_perimeter(
    *,
    w_channel: float,
    h_channel: float,
    t_rib: float,
    h_c: float,
    k_w: float,
    include_closeout: bool = True,
) -> FinResult:
    """Effective cold-side heated perimeter for one rectangular channel.

        P_eff = w_ch + 2 * eta_fin * h_ch  [+ w_ch if the closeout conducts]

    The channel-root face (``w_ch``, adjacent to the hot wall) is at full
    effectiveness; the two side walls are the ribs/lands, derated by
    ``eta_fin``. ``include_closeout`` adds the outer (back) face at full
    effectiveness -- appropriate for a closed-out channel whose outer wall is
    conducting, conservative to disable if the closeout is a separate jacket
    with poor thermal contact.
    """
    eta, m = rectangular_fin_efficiency(
        h_c=h_c, k_w=k_w, t_rib=t_rib, h_channel=h_channel
    )
    P_eff = w_channel + 2.0 * eta * h_channel
    if include_closeout:
        P_eff += w_channel
    return FinResult(efficiency=eta, m_1_m=m, P_effective_m=float(P_eff))


@dataclass(frozen=True)
class WallFaces:
    """Reconstructed wall state at one axial cell."""

    T_wg: float          # [K] hot-gas-side face
    T_wc: float          # [K] coolant-side face (channel root for ribbed walls)
    q_hot_W_m: float     # [W/m] hot fluid -> wall, per unit axial length
    q_cold_W_m: float    # [W/m] wall -> cold fluid, per unit axial length
    k_w: float           # [W/m/K] conductivity at the mean wall temperature


class CylindricalWall:
    """Single-layer cylindrical wall with temperature-dependent conductivity.

    ``fluxes_at_Tbar`` is the transient/FV entry point: given the lumped
    thickness-mean temperature ``T_bar`` (the one integrated state per axial
    cell), reconstruct both face temperatures from the quadratic quasi-static
    profile and return the per-unit-length face fluxes. The driver integrates
    ``(rho*cp*A_wall) dT_bar/dt = q_hot - q_cold``.

    ``solve_steady`` is the steady counterpart: the 3-resistance network with
    ``k_w`` converged against the mean wall temperature.
    """

    def __init__(
        self,
        *,
        D_inner: float,
        thickness: float,
        k_of_T,
        hot_side: str = "outer",
        rib: dict | None = None,
    ):
        if D_inner <= 0.0 or thickness <= 0.0:
            raise ValueError("D_inner and thickness must be positive")
        if hot_side not in ("inner", "outer"):
            raise ValueError("hot_side must be 'inner' or 'outer'")
        self.D_inner = float(D_inner)
        self.thickness = float(thickness)
        self.k_of_T = k_of_T
        self.hot_side = hot_side
        # rib: dict(w_channel=, h_channel=, t_rib=, n_channels=, include_closeout=)
        self.rib = rib

    # -- geometry -------------------------------------------------------
    def perimeters(self, h_c: float | None = None, k_w: float | None = None) -> tuple[float, float]:
        """(hot-side, cold-side) heated perimeters [m] for this wall.

        For a plain tube these are pi*D_inner / pi*(D_inner+2s) assigned by
        ``hot_side``. For a ribbed channel wall, the cold-side perimeter is
        replaced by the fin-derated effective perimeter summed over all
        channels (requires ``h_c`` and ``k_w`` -- the fin efficiency depends
        on the film coefficient).
        """
        P_inner = np.pi * self.D_inner
        P_outer = np.pi * (self.D_inner + 2.0 * self.thickness)
        P_hot, P_cold = (P_inner, P_outer) if self.hot_side == "inner" else (P_outer, P_inner)

        if self.rib is not None:
            if h_c is None or k_w is None:
                raise ValueError("ribbed wall requires h_c and k_w to size the fin")
            fin = channel_effective_perimeter(
                w_channel=self.rib["w_channel"],
                h_channel=self.rib["h_channel"],
                t_rib=self.rib["t_rib"],
                h_c=h_c,
                k_w=k_w,
                include_closeout=self.rib.get("include_closeout", True),
            )
            P_cold = fin.P_effective_m * self.rib["n_channels"]
        return float(P_hot), float(P_cold)

    def _faces_from_hgeff(self, T_bar, h_g_eff, h_c, T_g, T_c, a2, a6):
        """Closed-form 2x2 solve — identical algebra to
        ``physics/heat_conduction.py::_faces_from_hgeff`` (validated <2 K vs a
        resolved radial PDE). Kept explicit rather than calling np.linalg.solve
        because this sits on the transient hot path."""
        m00 = 1.0 + h_g_eff * (a2 - a6)
        m01 = -h_c * a6
        m10 = -(1.0 + h_g_eff * a2)
        m11 = 1.0 + h_c * a2
        r0 = T_bar + h_g_eff * T_g * (a2 - a6) - h_c * T_c * a6
        r1 = h_c * T_c * a2 - h_g_eff * T_g * a2
        det = m00 * m11 - m01 * m10
        T_wg = (r0 * m11 - m01 * r1) / det
        T_wc = (m00 * r1 - r0 * m10) / det
        return T_wg, T_wc

    def fluxes_at_Tbar(
        self, *, T_bar: float, T_g: float, T_c: float, h_g: float, h_c: float,
        h_g_rad: float = 0.0,
    ) -> WallFaces:
        """Face temperatures and per-unit-length fluxes at a given mean wall T."""
        k_w = float(self.k_of_T(T_bar))
        a2 = self.thickness / (2.0 * k_w)
        a6 = self.thickness / (6.0 * k_w)
        h_g_eff = h_g + h_g_rad
        T_wg, T_wc = self._faces_from_hgeff(T_bar, h_g_eff, h_c, T_g, T_c, a2, a6)
        P_hot, P_cold = self.perimeters(h_c=h_c, k_w=k_w)
        return WallFaces(
            T_wg=float(T_wg),
            T_wc=float(T_wc),
            q_hot_W_m=float(h_g_eff * P_hot * (T_g - T_wg)),
            q_cold_W_m=float(h_c * P_cold * (T_wc - T_c)),
            k_w=k_w,
        )

    def solve_steady(
        self, *, T_g: float, T_c: float, h_g: float, h_c: float,
        h_g_rad: float = 0.0, max_iter: int = 50, tol: float = 1e-9,
    ) -> WallFaces:
        """Steady 3-resistance network, iterating k_w(T_mean) to convergence.

        At steady state q_hot == q_cold by construction (both equal the
        series-network heat rate), which is asserted in the Stage C tests --
        the same self-consistency guardrail the legacy module documents.
        """
        h_g_eff = h_g + h_g_rad
        r_in = 0.5 * self.D_inner
        r_out = r_in + self.thickness
        T_bar = 0.5 * (T_g + T_c)
        for _ in range(max_iter):
            k_w = float(self.k_of_T(T_bar))
            P_hot, P_cold = self.perimeters(h_c=h_c, k_w=k_w)
            R_wall = np.log(r_out / r_in) / (2.0 * np.pi * k_w)   # [K*m/W]
            R_hot = 1.0 / (h_g_eff * P_hot)
            R_cold = 1.0 / (h_c * P_cold)
            q_per_m = (T_g - T_c) / (R_hot + R_wall + R_cold)
            T_wg = T_g - q_per_m * R_hot
            T_wc = T_wg - q_per_m * R_wall
            T_bar_new = 0.5 * (T_wg + T_wc)
            if abs(T_bar_new - T_bar) < tol:
                T_bar = T_bar_new
                break
            T_bar = T_bar_new
        return WallFaces(
            T_wg=float(T_wg), T_wc=float(T_wc),
            q_hot_W_m=float(q_per_m), q_cold_W_m=float(q_per_m), k_w=float(k_w),
        )

    def heat_capacity_per_length(self, rho: float, cp_of_T, T_bar: float) -> float:
        """[J/m/K] wall thermal inertia per unit axial length."""
        r_in = 0.5 * self.D_inner
        r_out = r_in + self.thickness
        A_wall = np.pi * (r_out**2 - r_in**2)
        return float(rho * float(cp_of_T(T_bar)) * A_wall)
