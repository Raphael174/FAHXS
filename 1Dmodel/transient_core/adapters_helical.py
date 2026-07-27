"""Helical-coil adapter helpers for the transient core.

The functions here do not run the legacy helical solver. They translate helical
geometry and material/fluid properties into the geometry-neutral arrays used by
`transient_core`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..physics.friction_correlations import dispatch_friction_coil
from ..physics.heat_conduction import OneDimensionalSteadyConduction_ShellnHelicalTube
from ..physics.heat_transfer_correlations import dispatch_nu_coil
from .grid import AxialGrid


@dataclass(frozen=True)
class HelicalCoreGeometry:
    """Geometry bridge from shell-and-helical-tube inputs to `transient_core`."""

    grid: AxialGrid
    tube_inner_diameter: float
    tube_outer_diameter: float
    wall_thickness: float
    n_parallel: int

    @property
    def single_coolant_area(self) -> float:
        return float(np.pi * self.tube_inner_diameter**2 / 4.0)

    @property
    def single_wall_area(self) -> float:
        return float(
            np.pi * (self.tube_outer_diameter**2 - self.tube_inner_diameter**2) / 4.0
        )


@dataclass(frozen=True)
class HelicalCoolantFilm:
    """Coolant-side hydraulic and heat-transfer state for one time/grid state."""

    velocity_m_s: np.ndarray
    reynolds: np.ndarray
    prandtl: np.ndarray
    friction_factor: np.ndarray
    nusselt: np.ndarray
    h_W_m2K: np.ndarray
    conductance_W_K: np.ndarray


@dataclass(frozen=True)
class HelicalWallFlux:
    """Wall reconstruction and total per-cell heat rates for helical cells."""

    hot_heat_W: np.ndarray
    cold_heat_W: np.ndarray
    dq_hot_per_length_W_m: np.ndarray
    dq_cold_per_length_W_m: np.ndarray
    T_wg: np.ndarray
    T_wc: np.ndarray
    h_g_rad: np.ndarray
    q_w_rad: np.ndarray
    k_wall: np.ndarray


def flow_direction_from_config(flow_config: str) -> int:
    """Return transient-core coolant direction from the project flow string."""

    if flow_config == "co":
        return 1
    if flow_config == "counter":
        return -1
    raise ValueError("flow_config must be 'co' or 'counter'")


def build_helical_core_geometry(
    *,
    pipe_length: float,
    n_cells: int,
    tube_inner_diameter: float,
    wall_thickness: float,
    n_parallel: int = 1,
    flow_config: str = "co",
) -> HelicalCoreGeometry:
    """Build transient-core geometry arrays for a helical coil.

    `pipe_length` is the helical tube arc length, not combustor axial length.
    Areas and perimeters stored in the returned `AxialGrid` are totals across
    `n_parallel` identical coils.
    """

    if tube_inner_diameter <= 0.0:
        raise ValueError("tube_inner_diameter must be positive")
    if wall_thickness <= 0.0:
        raise ValueError("wall_thickness must be positive")
    if int(n_parallel) <= 0:
        raise ValueError("n_parallel must be positive")

    n = int(n_parallel)
    tube_outer_diameter = tube_inner_diameter + 2.0 * wall_thickness
    coolant_area = n * np.pi * tube_inner_diameter**2 / 4.0
    wall_area = n * np.pi * (tube_outer_diameter**2 - tube_inner_diameter**2) / 4.0
    hot_perimeter = n * np.pi * tube_outer_diameter
    coolant_perimeter = n * np.pi * tube_inner_diameter

    grid = AxialGrid.uniform(
        length=pipe_length,
        n_cells=n_cells,
        coolant_area=coolant_area,
        wall_area=wall_area,
        hot_perimeter=hot_perimeter,
        coolant_perimeter=coolant_perimeter,
        flow_direction=flow_direction_from_config(flow_config),
    )
    return HelicalCoreGeometry(
        grid=grid,
        tube_inner_diameter=float(tube_inner_diameter),
        tube_outer_diameter=float(tube_outer_diameter),
        wall_thickness=float(wall_thickness),
        n_parallel=n,
    )


def build_helical_core_geometry_from_solver(solver, *, n_cells: int | None = None):
    """Build core geometry from an initialized legacy helical solver object.

    The object must expose the geometry attributes initialized by
    `main_solve.main_solver` or `main_solve_transient.transient_solver`.
    This helper intentionally reads attributes only; it does not call
    `solver.solver()` or perform chemistry/radiation setup.
    """

    if n_cells is None:
        if not hasattr(solver, "N"):
            raise ValueError("n_cells must be provided when solver.N is unavailable")
        n_cells = int(solver.N)

    return build_helical_core_geometry(
        pipe_length=float(solver.L_ch_max),
        n_cells=int(n_cells),
        tube_inner_diameter=float(solver.Dh_ch),
        wall_thickness=float(solver.combustorProp.thickness_coil_wall),
        n_parallel=int(solver.N_ch),
        flow_config=str(solver.combustorProp.flow_config),
    )


def wall_heat_capacity_J_K(grid: AxialGrid, density: float, cp) -> np.ndarray:
    """Return per-cell total wall heat capacity."""

    cp_arr = _cell_property("cp", cp, grid.n_cells)
    if density <= 0.0:
        raise ValueError("density must be positive")
    if np.any(cp_arr <= 0.0):
        raise ValueError("cp must be positive")
    return float(density) * cp_arr * grid.wall_volume


def coolant_heat_capacity_J_K(grid: AxialGrid, density, cp) -> np.ndarray:
    """Return per-cell total coolant heat capacity."""

    rho_arr = _cell_property("density", density, grid.n_cells)
    cp_arr = _cell_property("cp", cp, grid.n_cells)
    if np.any(rho_arr <= 0.0) or np.any(cp_arr <= 0.0):
        raise ValueError("density and cp must be positive")
    return rho_arr * cp_arr * grid.coolant_volume


def conductance_from_h(grid: AxialGrid, h_coolant) -> np.ndarray:
    """Return total wall-to-coolant conductance per cell from `h_c`."""

    h = _cell_property("h_coolant", h_coolant, grid.n_cells)
    if np.any(h < 0.0):
        raise ValueError("h_coolant must be non-negative")
    return h * grid.coolant_perimeter * grid.dx


def helical_coolant_film(
    geometry: HelicalCoreGeometry,
    *,
    mdot_total: float,
    rho,
    mu,
    k,
    cp,
    friction_selector: str,
    nusselt_selector: str,
    roughness: float,
    coil_radius: float,
    corrCoeffs,
    x_for_developing=None,
    friction_error_factor: float = 1.0,
    nusselt_error_factor: float = 1.0,
) -> HelicalCoolantFilm:
    """Compute helical coolant-side film coefficients from supplied properties.

    The formulas mirror the legacy helical solver call convention:

    ```text
    U  = mdot_total / (rho * A_coolant,total)
    Re = rho * U * D_i / mu
    Pr = cp * mu / k
    f  = dispatch_friction_coil(...)
    Nu = dispatch_nu_coil(...)
    h  = Nu * k / D_i
    G  = h * P_coolant,total * dx
    ```

    `rho`, `mu`, `k`, and `cp` may be scalars or per-cell arrays. This helper
    does not call CoolProp; the future production adapter should supply
    properties from the selected thermodynamic backend.
    """

    grid = geometry.grid
    rho_arr = _cell_property("rho", rho, grid.n_cells)
    mu_arr = _cell_property("mu", mu, grid.n_cells)
    k_arr = _cell_property("k", k, grid.n_cells)
    cp_arr = _cell_property("cp", cp, grid.n_cells)
    if np.any(rho_arr <= 0.0) or np.any(mu_arr <= 0.0):
        raise ValueError("rho and mu must be positive")
    if np.any(k_arr <= 0.0) or np.any(cp_arr <= 0.0):
        raise ValueError("k and cp must be positive")
    if coil_radius <= 0.0:
        raise ValueError("coil_radius must be positive")

    if x_for_developing is None:
        x_arr = np.full(grid.n_cells, 10e10)
    else:
        x_arr = _cell_property("x_for_developing", x_for_developing, grid.n_cells)

    velocity = abs(float(mdot_total)) / (rho_arr * grid.coolant_area)
    reynolds = rho_arr * velocity * geometry.tube_inner_diameter / mu_arr
    prandtl = cp_arr * mu_arr / k_arr

    friction = np.empty(grid.n_cells, dtype=float)
    nusselt = np.empty(grid.n_cells, dtype=float)
    for i in range(grid.n_cells):
        friction[i] = dispatch_friction_coil(
            friction_selector,
            Re=float(reynolds[i]),
            Dh=geometry.tube_inner_diameter,
            Rc=coil_radius,
            roughness=roughness,
            x=float(x_arr[i]),
            error_factor=friction_error_factor,
            corrCoeffs=corrCoeffs,
        )
        nusselt[i] = dispatch_nu_coil(
            nusselt_selector,
            Re=float(reynolds[i]),
            Pr=float(prandtl[i]),
            d=geometry.tube_inner_diameter,
            R=coil_radius,
            f_fd=float(friction[i]),
            x=float(x_arr[i]),
            error_factor=nusselt_error_factor,
            corrCoeffs=corrCoeffs,
        )

    h = nusselt * k_arr / geometry.tube_inner_diameter
    return HelicalCoolantFilm(
        velocity_m_s=velocity,
        reynolds=reynolds,
        prandtl=prandtl,
        friction_factor=friction,
        nusselt=nusselt,
        h_W_m2K=h,
        conductance_W_K=conductance_from_h(grid, h),
    )


def helical_wall_flux(
    geometry: HelicalCoreGeometry,
    *,
    Tbar_wall,
    T_coolant,
    T_gas,
    h_gas,
    h_coolant,
    wall_conductivity_at_T,
    h_g_rad=None,
) -> HelicalWallFlux:
    """Reconstruct helical wall faces and total per-cell heat rates.

    The existing conduction object returns `dq_*__dx` on a single-tube,
    per-length basis. This helper converts those values to the transient-core
    convention:

    ```text
    Q_cell,total = dq__dx_single_tube * dx * n_parallel
    ```

    Radiation should be supplied as a precomputed `h_g_rad` array for the fast
    path. Internal radiation-backend iteration is intentionally not used here.
    """

    grid = geometry.grid
    Tbar = _cell_property("Tbar_wall", Tbar_wall, grid.n_cells)
    Tc = _cell_property("T_coolant", T_coolant, grid.n_cells)
    Tg = _cell_property("T_gas", T_gas, grid.n_cells)
    hg = _cell_property("h_gas", h_gas, grid.n_cells)
    hc = _cell_property("h_coolant", h_coolant, grid.n_cells)
    hrad = np.zeros(grid.n_cells) if h_g_rad is None else _cell_property(
        "h_g_rad", h_g_rad, grid.n_cells
    )
    if np.any(hg < 0.0) or np.any(hc < 0.0) or np.any(hrad < 0.0):
        raise ValueError("heat-transfer coefficients must be non-negative")

    dq_hot = np.empty(grid.n_cells, dtype=float)
    dq_cold = np.empty(grid.n_cells, dtype=float)
    T_wg = np.empty(grid.n_cells, dtype=float)
    T_wc = np.empty(grid.n_cells, dtype=float)
    q_w_rad = np.empty(grid.n_cells, dtype=float)
    k_wall = np.empty(grid.n_cells, dtype=float)

    for i in range(grid.n_cells):
        node = OneDimensionalSteadyConduction_ShellnHelicalTube(
            h_g=float(hg[i]),
            h_c=float(hc[i]),
            T_c=float(Tc[i]),
            T_g=float(Tg[i]),
            s_w=geometry.wall_thickness,
            Dh_ch=geometry.tube_inner_diameter,
            f_kw_at_T=wall_conductivity_at_T,
            T_wg_0=float(Tbar[i]),
            T_wc_0=float(Tbar[i]),
            T_c_check_0=float(Tc[i]),
            dx=float(grid.dx[i]),
            rad_enabled=False,
            hot_side="outer",
        )
        result = node.fluxes_at_Tbar(float(Tbar[i]), h_g_rad=float(hrad[i]))
        dq_hot[i] = result["dq_hot__dx"]
        dq_cold[i] = result["dq_cold__dx"]
        T_wg[i] = result["T_wg"]
        T_wc[i] = result["T_wc"]
        q_w_rad[i] = result["q_w_rad"]
        k_wall[i] = result["k_w"]

    multiplier = geometry.n_parallel * grid.dx
    return HelicalWallFlux(
        hot_heat_W=dq_hot * multiplier,
        cold_heat_W=dq_cold * multiplier,
        dq_hot_per_length_W_m=dq_hot,
        dq_cold_per_length_W_m=dq_cold,
        T_wg=T_wg,
        T_wc=T_wc,
        h_g_rad=hrad,
        q_w_rad=q_w_rad,
        k_wall=k_wall,
    )


def _cell_property(name: str, value, n_cells: int) -> np.ndarray:
    if np.isscalar(value):
        arr = np.full(n_cells, float(value))
    else:
        arr = np.asarray(value, dtype=float)
    if arr.ndim != 1 or arr.size != n_cells:
        raise ValueError(f"{name} must be scalar or shape ({n_cells},)")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr
