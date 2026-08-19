"""Shell-and-tube adapter helpers for the transient core."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
from CoolProp.CoolProp import PropsSI

from ..physics.bell_delaware import bell_delaware_shell
from ..physics.friction_correlations import (
    dispatch_friction_tube_straight,
    friction_corrugated_tube_vicente,
)
from ..physics.heat_conduction import OneDimensionalSteadyConduction_ShellnHelicalTube
from ..physics.heat_transfer_correlations import (
    dispatch_nu_tube_straight,
    nu_corrugated_tube_vicente,
)
from .grid import AxialGrid
from .integrator import (
    WallCoolantIntegrationResult,
    WallCoolantStepInputs,
    fixed_time_grid,
    integrate_wall_coolant_fixed_step,
)
from .progress import TransientProgressPrinter
from .schedules import interp_schedule
from .compressible_coolant import (
    CoolantThermodynamicState,
    coolprop_state_from_mass_energy,
    enforce_density_bounds,
    enforce_internal_energy_floor,
    initial_mass_energy_from_TP,
)
from .wall_compressible_coolant import (
    WallCompressibleCoolantStepResult,
    semi_implicit_wall_compressible_coolant_step,
)
from .wall_coolant import implicit_wall_coolant_step


@dataclass(frozen=True)
class ShellTubeCoreGeometry:
    """Geometry bridge from baffled shell-and-tube inputs to `transient_core`."""

    grid: AxialGrid
    tube_inner_diameter: float
    tube_outer_diameter: float
    wall_thickness: float
    n_tubes: int
    shell_inner_diameter: float

    @property
    def single_tube_inner_area(self) -> float:
        return float(np.pi * self.tube_inner_diameter**2 / 4.0)

    @property
    def single_tube_wall_area(self) -> float:
        return float(
            np.pi * (self.tube_outer_diameter**2 - self.tube_inner_diameter**2) / 4.0
        )


@dataclass(frozen=True)
class ShellTubeShellFilm:
    """Shell-side hydraulic and heat-transfer state for one time/grid state."""

    mass_flux_kg_m2s: np.ndarray
    reynolds: np.ndarray
    prandtl: np.ndarray
    h_W_m2K: np.ndarray
    conductance_W_K: np.ndarray
    dp_shell_Pa: np.ndarray


@dataclass(frozen=True)
class ShellTubeTubeGasFilm:
    """Tube-side hot-gas hydraulic and heat-transfer state."""

    velocity_m_s: np.ndarray
    reynolds: np.ndarray
    prandtl: np.ndarray
    friction_factor: np.ndarray
    nusselt: np.ndarray
    h_W_m2K: np.ndarray
    hot_conductance_W_K: np.ndarray
    dp_per_length_Pa_m: np.ndarray


@dataclass(frozen=True)
class ShellTubeWallFlux:
    """Wall reconstruction and total per-cell heat rates for shell-and-tube cells."""

    hot_heat_W: np.ndarray
    cold_heat_W: np.ndarray
    dq_hot_per_length_W_m: np.ndarray
    dq_cold_per_length_W_m: np.ndarray
    T_wg: np.ndarray
    T_wc: np.ndarray
    h_g_rad: np.ndarray
    q_w_rad: np.ndarray
    k_wall: np.ndarray


@dataclass(frozen=True)
class ShellTubeGasState:
    """Thermophysical state returned by a shell-and-tube hot-gas provider."""

    T: float
    rho: float
    mu: float
    k: float
    cp: float
    progress_source: float = 0.0


@dataclass(frozen=True)
class ShellTubeHotGasMarch:
    """Sequential representative-tube hot-gas march result."""

    T_gas: np.ndarray
    h_gas_W_m2K: np.ndarray
    gas_velocity_m_s: np.ndarray
    reynolds: np.ndarray
    prandtl: np.ndarray
    friction_factor: np.ndarray
    nusselt: np.ndarray
    dp_per_length_Pa_m: np.ndarray
    wall_flux: ShellTubeWallFlux
    enthalpy_removed_J_kg: np.ndarray
    progress_variable: np.ndarray
    T_gas_outlet: float
    enthalpy_removed_outlet_J_kg: float
    progress_outlet: float


@dataclass(frozen=True)
class ShellTubeFluidProperties:
    """Cellwise thermophysical properties used by shell-and-tube adapters."""

    rho: np.ndarray
    mu: np.ndarray
    k: np.ndarray
    cp: np.ndarray


@dataclass(frozen=True)
class ShellTubeStepInputDiagnostics:
    """Diagnostics returned with one assembled wall/coolant step input."""

    wall_coolant_inputs: WallCoolantStepInputs
    hot_gas_march: ShellTubeHotGasMarch
    shell_film: ShellTubeShellFilm
    coolant_properties: ShellTubeFluidProperties
    wall_heat_capacity_J_K: np.ndarray


@dataclass(frozen=True)
class ShellTubeTransientCoreResult:
    """Result from a shell-and-tube transient-core wall/coolant run."""

    integration: WallCoolantIntegrationResult
    step_diagnostics: tuple[ShellTubeStepInputDiagnostics, ...]


@dataclass(frozen=True)
class ShellTubeCompressibleIntegrationResult:
    """Time history from wall + compressible coolant mass/energy integration."""

    t: np.ndarray
    T_wall: np.ndarray
    T_coolant: np.ndarray
    T_coolant_outlet: np.ndarray
    coolant_mass_kg: np.ndarray
    coolant_internal_energy_J: np.ndarray
    coolant_pressure_Pa: np.ndarray
    coolant_density_kg_m3: np.ndarray
    coolant_specific_enthalpy_J_kg: np.ndarray
    face_mdot_kg_s: np.ndarray
    hot_heat_added_J: np.ndarray
    advective_energy_in_J: np.ndarray
    advective_energy_out_J: np.ndarray
    energy_residual_J: np.ndarray
    mass_residual_kg: np.ndarray
    heat_wall_to_coolant_W: np.ndarray
    last_step: WallCompressibleCoolantStepResult | None


GasStateProvider = Callable[[float, float, int], ShellTubeGasState | dict]
GasProviderAtTime = Callable[[float], tuple[GasStateProvider, float]]
ScalarAtTime = Callable[[float], float]
CoolantPropertyProvider = Callable[[np.ndarray, float], ShellTubeFluidProperties | dict]


def shelltube_flow_direction(flow_config: str) -> int:
    """Return shell-side coolant direction from the project flow string."""

    if flow_config == "co":
        return 1
    if flow_config == "counter":
        return -1
    raise ValueError("flow_config must be 'co' or 'counter'")


def build_shelltube_core_geometry(
    *,
    tube_length: float,
    n_cells: int,
    shell_inner_diameter: float,
    tube_outer_diameter: float,
    wall_thickness: float,
    n_tubes: int,
    flow_config: str = "co",
) -> ShellTubeCoreGeometry:
    """Build transient-core geometry arrays for the shell-and-tube config.

    The coolant inventory area is approximated as shell bore area minus total
    tube outside area. Wall and tube-side perimeters are totals across all tubes.
    """

    if tube_length <= 0.0:
        raise ValueError("tube_length must be positive")
    if shell_inner_diameter <= 0.0:
        raise ValueError("shell_inner_diameter must be positive")
    if tube_outer_diameter <= 0.0:
        raise ValueError("tube_outer_diameter must be positive")
    if wall_thickness <= 0.0:
        raise ValueError("wall_thickness must be positive")
    if int(n_tubes) <= 0:
        raise ValueError("n_tubes must be positive")

    n = int(n_tubes)
    tube_inner_diameter = tube_outer_diameter - 2.0 * wall_thickness
    if tube_inner_diameter <= 0.0:
        raise ValueError("tube inner diameter must be positive")

    shell_area = np.pi * shell_inner_diameter**2 / 4.0
    displaced_tube_area = n * np.pi * tube_outer_diameter**2 / 4.0
    coolant_area = shell_area - displaced_tube_area
    if coolant_area <= 0.0:
        raise ValueError("shell-side coolant area must be positive")

    wall_area = n * np.pi * (tube_outer_diameter**2 - tube_inner_diameter**2) / 4.0
    hot_perimeter = n * np.pi * tube_inner_diameter
    coolant_perimeter = n * np.pi * tube_outer_diameter

    grid = AxialGrid.uniform(
        length=tube_length,
        n_cells=n_cells,
        coolant_area=coolant_area,
        wall_area=wall_area,
        hot_perimeter=hot_perimeter,
        coolant_perimeter=coolant_perimeter,
        flow_direction=shelltube_flow_direction(flow_config),
    )
    return ShellTubeCoreGeometry(
        grid=grid,
        tube_inner_diameter=float(tube_inner_diameter),
        tube_outer_diameter=float(tube_outer_diameter),
        wall_thickness=float(wall_thickness),
        n_tubes=n,
        shell_inner_diameter=float(shell_inner_diameter),
    )


def build_shelltube_core_geometry_from_solver(solver, *, n_cells: int | None = None):
    """Build core geometry from an initialized legacy shell-and-tube solver."""

    if n_cells is None:
        if not hasattr(solver, "N"):
            raise ValueError("n_cells must be provided when solver.N is unavailable")
        n_cells = int(solver.N)

    return build_shelltube_core_geometry(
        tube_length=float(solver.stp.L_tube),
        n_cells=int(n_cells),
        shell_inner_diameter=float(solver.stp.D_shell_inner),
        tube_outer_diameter=float(solver.stp.D_tube_outer),
        wall_thickness=float(solver.stp.thickness_tube_wall),
        n_tubes=int(solver.stp.N_tubes),
        flow_config=str(getattr(solver, "flow_config", "co")),
    )


def shelltube_wall_heat_capacity_J_K(grid: AxialGrid, density: float, cp) -> np.ndarray:
    """Return per-cell total tube-wall heat capacity."""

    cp_arr = _cell_property("cp", cp, grid.n_cells)
    if density <= 0.0:
        raise ValueError("density must be positive")
    if np.any(cp_arr <= 0.0):
        raise ValueError("cp must be positive")
    return float(density) * cp_arr * grid.wall_volume


def shelltube_coolant_heat_capacity_J_K(grid: AxialGrid, density, cp) -> np.ndarray:
    """Return per-cell shell-side coolant heat capacity."""

    rho_arr = _cell_property("density", density, grid.n_cells)
    cp_arr = _cell_property("cp", cp, grid.n_cells)
    if np.any(rho_arr <= 0.0) or np.any(cp_arr <= 0.0):
        raise ValueError("density and cp must be positive")
    return rho_arr * cp_arr * grid.coolant_volume


def shelltube_conductance_from_h(grid: AxialGrid, h_shell) -> np.ndarray:
    """Return total shell-side wall-to-coolant conductance per cell."""

    h = _cell_property("h_shell", h_shell, grid.n_cells)
    if np.any(h < 0.0):
        raise ValueError("h_shell must be non-negative")
    return h * grid.coolant_perimeter * grid.dx


def shelltube_shell_film(
    geometry: ShellTubeCoreGeometry,
    bell_geometry: dict,
    *,
    mdot_shell: float,
    rho,
    mu,
    k,
    cp,
    corrCoeffs=None,
    mu_ratio=1.0,
    mdot_floor: float = 1.0e-12,
) -> ShellTubeShellFilm:
    """Compute shell-side Bell-Delaware film coefficients from supplied properties.

    This mirrors the maintained shell-and-tube transient call convention:

    ```text
    G_s = mdot_shell / S_m
    Re_s = D_o * G_s / mu
    Pr_s = cp * mu / k
    h_shell = bell_delaware_shell(...)[h_shell]
    G_i = h_shell_i * P_outer,total * dx_i
    ```

    The thermophysical properties may be scalar or per-cell arrays. This helper
    does not call CoolProp.
    """

    grid = geometry.grid
    rho_arr = _cell_property("rho", rho, grid.n_cells)
    mu_arr = _cell_property("mu", mu, grid.n_cells)
    k_arr = _cell_property("k", k, grid.n_cells)
    cp_arr = _cell_property("cp", cp, grid.n_cells)
    mu_ratio_arr = _cell_property("mu_ratio", mu_ratio, grid.n_cells)
    if np.any(rho_arr <= 0.0) or np.any(mu_arr <= 0.0):
        raise ValueError("rho and mu must be positive")
    if np.any(k_arr <= 0.0) or np.any(cp_arr <= 0.0):
        raise ValueError("k and cp must be positive")
    if "S_m" not in bell_geometry or bell_geometry["S_m"] <= 0.0:
        raise ValueError("bell_geometry must contain positive S_m")

    mdot_abs = abs(float(mdot_shell))
    G_s = mdot_abs / float(bell_geometry["S_m"])
    mass_flux = np.full(grid.n_cells, G_s, dtype=float)
    reynolds = geometry.tube_outer_diameter * mass_flux / mu_arr
    prandtl = cp_arr * mu_arr / k_arr
    h = np.empty(grid.n_cells, dtype=float)
    dp = np.empty(grid.n_cells, dtype=float)

    if mdot_abs <= float(mdot_floor):
        h = k_arr / geometry.tube_outer_diameter
        return ShellTubeShellFilm(
            mass_flux_kg_m2s=mass_flux,
            reynolds=reynolds,
            prandtl=prandtl,
            h_W_m2K=h,
            conductance_W_K=shelltube_conductance_from_h(grid, h),
            dp_shell_Pa=np.zeros(grid.n_cells, dtype=float),
        )

    for i in range(grid.n_cells):
        geom = dict(bell_geometry)
        geom["rho_s"] = float(rho_arr[i])
        result = bell_delaware_shell(
            geom,
            Re_s=float(reynolds[i]),
            Pr_s=float(prandtl[i]),
            k_s=float(k_arr[i]),
            cp_s=float(cp_arr[i]),
            mu_s=float(mu_arr[i]),
            mdot_s=mdot_abs,
            mu_ratio=float(mu_ratio_arr[i]),
            corrCoeffs=corrCoeffs,
        )
        h[i] = result["h_shell"]
        dp[i] = result.get("dp_shell", 0.0)

    return ShellTubeShellFilm(
        mass_flux_kg_m2s=mass_flux,
        reynolds=reynolds,
        prandtl=prandtl,
        h_W_m2K=h,
        conductance_W_K=shelltube_conductance_from_h(grid, h),
        dp_shell_Pa=dp,
    )


def shelltube_tube_gas_film(
    geometry: ShellTubeCoreGeometry,
    *,
    mdot_hot_total: float,
    rho,
    mu,
    k,
    cp,
    inside_tube_choice: str,
    nusselt_selector: str,
    roughness: float,
    corrCoeffs=None,
    corrugation_thickness: float = 0.0,
    corrugation_pitch: float = 1.0,
    x_for_developing=None,
    T_bulk=None,
    T_wall=None,
) -> ShellTubeTubeGasFilm:
    """Compute representative tube-side hot-gas film coefficients.

    Properties are supplied as scalar or per-cell arrays. The mass flow is the
    total bundle hot-side flow; the adapter divides it by `N_tubes`, matching
    the maintained shell-and-tube solver convention.
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
    if roughness < 0.0:
        raise ValueError("roughness must be non-negative")

    if x_for_developing is None:
        x_arr = grid.x_centers
    else:
        x_arr = _cell_property("x_for_developing", x_for_developing, grid.n_cells)
    Tb_arr = None if T_bulk is None else _cell_property("T_bulk", T_bulk, grid.n_cells)
    Tw_arr = None if T_wall is None else _cell_property("T_wall", T_wall, grid.n_cells)

    mdot_tube = abs(float(mdot_hot_total)) / geometry.n_tubes
    velocity = mdot_tube / (rho_arr * geometry.single_tube_inner_area)
    reynolds = rho_arr * velocity * geometry.tube_inner_diameter / mu_arr
    prandtl = cp_arr * mu_arr / k_arr

    nu_factor, f_factor = _tube_surface_factors(inside_tube_choice, corrCoeffs)
    phi = _corrugation_severity(
        corrugation_thickness,
        corrugation_pitch,
        geometry.tube_inner_diameter,
    )
    re_lo = getattr(corrCoeffs, "Re_transition_lo", 2300.0)
    re_hi = getattr(corrCoeffs, "Re_transition_hi", 4000.0)

    friction = np.empty(grid.n_cells, dtype=float)
    nusselt = np.empty(grid.n_cells, dtype=float)
    for i in range(grid.n_cells):
        if inside_tube_choice == "grooved":
            friction[i] = friction_corrugated_tube_vicente(
                float(reynolds[i]),
                phi,
                Re_lo=re_lo,
                Re_hi=re_hi,
            )
            nusselt[i] = nu_corrugated_tube_vicente(
                float(reynolds[i]),
                float(prandtl[i]),
                phi,
                D_i=geometry.tube_inner_diameter,
                x=float(x_arr[i]),
                Re_lo=re_lo,
                Re_hi=re_hi,
            )
        else:
            friction[i] = dispatch_friction_tube_straight(
                float(reynolds[i]),
                roughness,
                geometry.tube_inner_diameter,
                x=float(x_arr[i]),
                Re_lo=re_lo,
                Re_hi=re_hi,
            )
            nusselt[i] = dispatch_nu_tube_straight(
                nusselt_selector,
                Re=float(reynolds[i]),
                Pr=float(prandtl[i]),
                d=geometry.tube_inner_diameter,
                x=float(x_arr[i]),
                f_fd=float(friction[i]),
                T_bulk=None if Tb_arr is None else float(Tb_arr[i]),
                T_wall=None if Tw_arr is None else float(Tw_arr[i]),
                error_factor=1.0,
                corrCoeffs=corrCoeffs,
            )
    friction *= f_factor
    nusselt *= nu_factor

    h = nusselt * k_arr / geometry.tube_inner_diameter
    dpdx = friction * rho_arr * velocity**2 / (2.0 * geometry.tube_inner_diameter)
    return ShellTubeTubeGasFilm(
        velocity_m_s=velocity,
        reynolds=reynolds,
        prandtl=prandtl,
        friction_factor=friction,
        nusselt=nusselt,
        h_W_m2K=h,
        hot_conductance_W_K=h * grid.hot_perimeter * grid.dx,
        dp_per_length_Pa_m=dpdx,
    )


def shelltube_wall_flux(
    geometry: ShellTubeCoreGeometry,
    *,
    Tbar_wall,
    T_coolant,
    T_gas,
    h_gas,
    h_shell,
    wall_conductivity_at_T,
    h_g_rad=None,
) -> ShellTubeWallFlux:
    """Reconstruct shell-and-tube wall faces and total per-cell heat rates.

    The maintained conduction object returns `dq_*__dx` for one representative
    tube. This helper converts those values to the transient-core total-cell
    convention:

    ```text
    Q_cell,total = dq__dx_single_tube * dx * N_tubes
    ```

    The shell-and-tube orientation is `hot_side="inner"`: hot gas is inside the
    tubes, while helium/coolant is outside on the shell side.
    """

    grid = geometry.grid
    Tbar = _cell_property("Tbar_wall", Tbar_wall, grid.n_cells)
    Tc = _cell_property("T_coolant", T_coolant, grid.n_cells)
    Tg = _cell_property("T_gas", T_gas, grid.n_cells)
    hg = _cell_property("h_gas", h_gas, grid.n_cells)
    hc = _cell_property("h_shell", h_shell, grid.n_cells)
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
            hot_side="inner",
        )
        result = node.fluxes_at_Tbar(float(Tbar[i]), h_g_rad=float(hrad[i]))
        dq_hot[i] = result["dq_hot__dx"]
        dq_cold[i] = result["dq_cold__dx"]
        T_wg[i] = result["T_wg"]
        T_wc[i] = result["T_wc"]
        q_w_rad[i] = result["q_w_rad"]
        k_wall[i] = result["k_w"]

    multiplier = geometry.n_tubes * grid.dx
    return ShellTubeWallFlux(
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


def shelltube_hot_gas_march(
    geometry: ShellTubeCoreGeometry,
    *,
    Tbar_wall,
    T_coolant,
    h_shell,
    mdot_hot_total: float,
    gas_state_at: GasStateProvider,
    wall_conductivity_at_T,
    inside_tube_choice: str,
    nusselt_selector: str,
    roughness: float,
    corrCoeffs=None,
    corrugation_thickness: float = 0.0,
    corrugation_pitch: float = 1.0,
    progress_initial: float = 0.0,
    h_g_rad=None,
) -> ShellTubeHotGasMarch:
    """March hot gas through one representative tube against current wall/coolant.

    `gas_state_at(h_removed, progress, i)` supplies thermochemistry and transport
    properties. The march updates:

    ```text
    h_removed_{i+1} = h_removed_i + dq_hot__dx_i * dx_i / mdot_tube
    progress_{i+1}  = progress_i + progress_source_i * dx_i / U_i
    ```

    Heat rates returned in `wall_flux` are total bundle-cell watts, while
    `enthalpy_removed` is per unit mass in the representative tube.
    """

    grid = geometry.grid
    if mdot_hot_total <= 0.0:
        raise ValueError("mdot_hot_total must be positive for hot-gas marching")
    Tbar = _cell_property("Tbar_wall", Tbar_wall, grid.n_cells)
    Tc = _cell_property("T_coolant", T_coolant, grid.n_cells)
    hs = _cell_property("h_shell", h_shell, grid.n_cells)
    hrad = np.zeros(grid.n_cells) if h_g_rad is None else _cell_property(
        "h_g_rad", h_g_rad, grid.n_cells
    )
    if np.any(hs < 0.0) or np.any(hrad < 0.0):
        raise ValueError("heat-transfer coefficients must be non-negative")

    n = grid.n_cells
    T_g = np.empty(n, dtype=float)
    rho = np.empty(n, dtype=float)
    mu = np.empty(n, dtype=float)
    k = np.empty(n, dtype=float)
    cp = np.empty(n, dtype=float)
    progress_source = np.empty(n, dtype=float)
    h_g = np.empty(n, dtype=float)
    velocity = np.empty(n, dtype=float)
    reynolds = np.empty(n, dtype=float)
    prandtl = np.empty(n, dtype=float)
    friction = np.empty(n, dtype=float)
    nusselt = np.empty(n, dtype=float)
    dpdx = np.empty(n, dtype=float)
    dq_hot = np.empty(n, dtype=float)
    dq_cold = np.empty(n, dtype=float)
    T_wg = np.empty(n, dtype=float)
    T_wc = np.empty(n, dtype=float)
    q_w_rad = np.empty(n, dtype=float)
    k_wall = np.empty(n, dtype=float)
    h_removed_history = np.empty(n, dtype=float)
    progress_history = np.empty(n, dtype=float)

    mdot_tube = abs(float(mdot_hot_total)) / geometry.n_tubes
    h_removed = 0.0
    progress = float(progress_initial)

    for i in range(n):
        h_removed_history[i] = h_removed
        progress_history[i] = progress
        state = _coerce_gas_state(gas_state_at(h_removed, progress, i))
        T_g[i] = state.T
        rho[i] = state.rho
        mu[i] = state.mu
        k[i] = state.k
        cp[i] = state.cp
        progress_source[i] = state.progress_source

        film = shelltube_tube_gas_film(
            _single_cell_shelltube_geometry(geometry, i),
            mdot_hot_total=float(mdot_hot_total),
            rho=state.rho,
            mu=state.mu,
            k=state.k,
            cp=state.cp,
            inside_tube_choice=inside_tube_choice,
            nusselt_selector=nusselt_selector,
            roughness=roughness,
            corrCoeffs=corrCoeffs,
            corrugation_thickness=corrugation_thickness,
            corrugation_pitch=corrugation_pitch,
            x_for_developing=float(grid.x_centers[i]),
            T_bulk=state.T,
            T_wall=float(Tbar[i]),
        )
        h_g[i] = film.h_W_m2K[0]
        velocity[i] = film.velocity_m_s[0]
        reynolds[i] = film.reynolds[0]
        prandtl[i] = film.prandtl[0]
        friction[i] = film.friction_factor[0]
        nusselt[i] = film.nusselt[0]
        dpdx[i] = film.dp_per_length_Pa_m[0]

        node = OneDimensionalSteadyConduction_ShellnHelicalTube(
            h_g=float(h_g[i]),
            h_c=float(hs[i]),
            T_c=float(Tc[i]),
            T_g=float(T_g[i]),
            s_w=geometry.wall_thickness,
            Dh_ch=geometry.tube_inner_diameter,
            f_kw_at_T=wall_conductivity_at_T,
            T_wg_0=float(Tbar[i]),
            T_wc_0=float(Tbar[i]),
            T_c_check_0=float(Tc[i]),
            dx=float(grid.dx[i]),
            rad_enabled=False,
            hot_side="inner",
        )
        result = node.fluxes_at_Tbar(float(Tbar[i]), h_g_rad=float(hrad[i]))
        dq_hot[i] = result["dq_hot__dx"]
        dq_cold[i] = result["dq_cold__dx"]
        T_wg[i] = result["T_wg"]
        T_wc[i] = result["T_wc"]
        q_w_rad[i] = result["q_w_rad"]
        k_wall[i] = result["k_w"]

        h_removed += dq_hot[i] * grid.dx[i] / mdot_tube
        if velocity[i] > 0.0:
            progress += progress_source[i] * grid.dx[i] / velocity[i]

    wall_flux = ShellTubeWallFlux(
        hot_heat_W=dq_hot * geometry.n_tubes * grid.dx,
        cold_heat_W=dq_cold * geometry.n_tubes * grid.dx,
        dq_hot_per_length_W_m=dq_hot,
        dq_cold_per_length_W_m=dq_cold,
        T_wg=T_wg,
        T_wc=T_wc,
        h_g_rad=hrad,
        q_w_rad=q_w_rad,
        k_wall=k_wall,
    )
    outlet_state = _coerce_gas_state(gas_state_at(h_removed, progress, n))
    return ShellTubeHotGasMarch(
        T_gas=T_g,
        h_gas_W_m2K=h_g,
        gas_velocity_m_s=velocity,
        reynolds=reynolds,
        prandtl=prandtl,
        friction_factor=friction,
        nusselt=nusselt,
        dp_per_length_Pa_m=dpdx,
        wall_flux=wall_flux,
        enthalpy_removed_J_kg=h_removed_history,
        progress_variable=progress_history,
        T_gas_outlet=float(outlet_state.T),
        enthalpy_removed_outlet_J_kg=float(h_removed),
        progress_outlet=float(progress),
    )


def shelltube_step_inputs(
    geometry: ShellTubeCoreGeometry,
    bell_geometry: dict,
    *,
    Tbar_wall,
    T_coolant,
    mdot_coolant: float,
    T_coolant_inlet: float,
    p_coolant: float,
    mdot_hot_total: float,
    gas_state_at: GasStateProvider,
    coolant_properties_at: CoolantPropertyProvider,
    wall_density: float,
    wall_cp,
    wall_conductivity_at_T,
    inside_tube_choice: str,
    nusselt_selector: str,
    tube_roughness: float,
    corrCoeffs=None,
    corrugation_thickness: float = 0.0,
    corrugation_pitch: float = 1.0,
    progress_initial: float = 0.0,
    flow_direction: int | None = None,
    h_g_rad=None,
) -> ShellTubeStepInputDiagnostics:
    """Assemble one transient-core wall/coolant step for shell-and-tube.

    This is the production-dispatch seam between geometry/thermochemistry and
    the geometry-neutral wall/coolant integrator. It evaluates coolant
    properties, Bell-Delaware shell film, quasi-steady representative-tube
    hot-gas march, wall heat capacity, and returns `WallCoolantStepInputs`.
    """

    grid = geometry.grid
    if p_coolant <= 0.0:
        raise ValueError("p_coolant must be positive")
    flow = grid.flow_direction if flow_direction is None else int(flow_direction)
    if flow not in (-1, 1):
        raise ValueError("flow_direction must be +1 or -1")

    Tbar = _cell_property("Tbar_wall", Tbar_wall, grid.n_cells)
    Tc = _cell_property("T_coolant", T_coolant, grid.n_cells)
    coolant_props = _coerce_fluid_properties(
        coolant_properties_at(Tc, float(p_coolant)),
        grid.n_cells,
    )
    shell_film = shelltube_shell_film(
        geometry,
        bell_geometry,
        mdot_shell=mdot_coolant,
        rho=coolant_props.rho,
        mu=coolant_props.mu,
        k=coolant_props.k,
        cp=coolant_props.cp,
        corrCoeffs=corrCoeffs,
    )
    if mdot_hot_total > 0.0:
        hot_march = shelltube_hot_gas_march(
            geometry,
            Tbar_wall=Tbar,
            T_coolant=Tc,
            h_shell=shell_film.h_W_m2K,
            mdot_hot_total=mdot_hot_total,
            gas_state_at=gas_state_at,
            wall_conductivity_at_T=wall_conductivity_at_T,
            inside_tube_choice=inside_tube_choice,
            nusselt_selector=nusselt_selector,
            roughness=tube_roughness,
            corrCoeffs=corrCoeffs,
            corrugation_thickness=corrugation_thickness,
            corrugation_pitch=corrugation_pitch,
            progress_initial=progress_initial,
            h_g_rad=h_g_rad,
        )
    else:
        zeros = np.zeros(grid.n_cells, dtype=float)
        if callable(wall_conductivity_at_T):
            k_wall = np.array([wall_conductivity_at_T(float(Ti)) for Ti in Tbar], dtype=float)
        else:
            k_wall = _cell_property("k_wall", wall_conductivity_at_T, grid.n_cells)
        hot_march = ShellTubeHotGasMarch(
            T_gas=np.full(grid.n_cells, np.nan),
            h_gas_W_m2K=zeros.copy(),
            gas_velocity_m_s=zeros.copy(),
            reynolds=zeros.copy(),
            prandtl=zeros.copy(),
            friction_factor=zeros.copy(),
            nusselt=zeros.copy(),
            dp_per_length_Pa_m=zeros.copy(),
            wall_flux=ShellTubeWallFlux(
                hot_heat_W=zeros.copy(),
                cold_heat_W=zeros.copy(),
                dq_hot_per_length_W_m=zeros.copy(),
                dq_cold_per_length_W_m=zeros.copy(),
                T_wg=Tbar.copy(),
                T_wc=Tbar.copy(),
                h_g_rad=zeros.copy(),
                q_w_rad=zeros.copy(),
                k_wall=k_wall,
            ),
            enthalpy_removed_J_kg=zeros.copy(),
            progress_variable=zeros.copy(),
            T_gas_outlet=np.nan,
            enthalpy_removed_outlet_J_kg=0.0,
            progress_outlet=float(progress_initial),
        )
    wall_cp_value = wall_cp(Tbar) if callable(wall_cp) else wall_cp
    wall_capacity = shelltube_wall_heat_capacity_J_K(grid, wall_density, wall_cp_value)
    coolant_capacity = shelltube_coolant_heat_capacity_J_K(
        grid,
        coolant_props.rho,
        coolant_props.cp,
    )
    inputs = WallCoolantStepInputs(
        wall_heat_capacity=wall_capacity,
        coolant_heat_capacity=coolant_capacity,
        coolant_cp=coolant_props.cp,
        mdot_coolant=float(mdot_coolant),
        T_coolant_inlet=float(T_coolant_inlet),
        hot_heat_W=hot_march.wall_flux.hot_heat_W,
        wall_to_coolant_conductance_W_per_K=shell_film.conductance_W_K,
        flow_direction=flow,
    )
    return ShellTubeStepInputDiagnostics(
        wall_coolant_inputs=inputs,
        hot_gas_march=hot_march,
        shell_film=shell_film,
        coolant_properties=coolant_props,
        wall_heat_capacity_J_K=wall_capacity,
    )


def _shelltube_initial_pressure_profile(
    geometry: ShellTubeCoreGeometry,
    *,
    inlet_pressure: float,
    pressure_drop: float,
    flow_direction: int,
) -> np.ndarray:
    """Return a linear pressure profile consistent with nominal flow direction."""

    n = geometry.grid.n_cells
    if n == 1:
        return np.array([float(inlet_pressure) - 0.5 * float(pressure_drop)])
    frac = np.linspace(0.0, 1.0, n)
    if flow_direction == 1:
        return float(inlet_pressure) - float(pressure_drop) * frac
    return float(inlet_pressure) - float(pressure_drop) * frac[::-1]


def _shelltube_boundary_pressure_profile(
    geometry: ShellTubeCoreGeometry,
    *,
    inlet_pressure: float,
    outlet_pressure: float,
    flow_direction: int,
) -> np.ndarray:
    """Return a cell pressure profile interpolated between hydraulic boundaries."""

    n = geometry.grid.n_cells
    if n == 1:
        return np.array([0.5 * (float(inlet_pressure) + float(outlet_pressure))])
    frac = np.linspace(0.0, 1.0, n)
    if flow_direction == 1:
        left = float(inlet_pressure)
        right = float(outlet_pressure)
    else:
        left = float(outlet_pressure)
        right = float(inlet_pressure)
    return left + (right - left) * frac


def _shelltube_nominal_pressure_drop(
    geometry: ShellTubeCoreGeometry,
    bell_geometry: dict,
    *,
    mdot_shell: float,
    T_coolant: np.ndarray,
    p_coolant: float,
    coolant_properties_at: CoolantPropertyProvider,
    corrCoeffs,
    mdot_floor: float,
) -> float:
    """Estimate a positive whole-exchanger shell pressure drop for resistance calibration."""

    props = _coerce_fluid_properties(
        coolant_properties_at(T_coolant, p_coolant),
        geometry.grid.n_cells,
    )
    film = shelltube_shell_film(
        geometry,
        bell_geometry,
        mdot_shell=max(abs(float(mdot_shell)), float(mdot_floor)),
        rho=props.rho,
        mu=props.mu,
        k=props.k,
        cp=props.cp,
        corrCoeffs=corrCoeffs,
        mdot_floor=mdot_floor,
    )
    finite_dp = np.asarray(film.dp_shell_Pa, dtype=float)
    finite_dp = finite_dp[np.isfinite(finite_dp)]
    if finite_dp.size == 0:
        return 1.0
    # Bell-Delaware returns a whole-shell estimate for each local property set.
    return max(float(np.nanmean(np.abs(finite_dp))), 1.0)


def _shelltube_face_resistance(
    pressure_drop: float,
    density: np.ndarray,
    mdot_reference: float,
    *,
    n_faces: int,
) -> np.ndarray:
    """Distribute a whole-shell quadratic resistance over internal/outlet faces."""

    rho_mean = max(float(np.nanmean(density)), 1.0e-9)
    mdot_ref = max(abs(float(mdot_reference)), 1.0e-9)
    dp_face = max(float(pressure_drop), 1.0) / max(int(n_faces), 1)
    return np.full(max(int(n_faces), 1), rho_mean * dp_face / mdot_ref**2)


def _orifice_mdot_from_resistance(dp: float, rho: float, resistance: float) -> float:
    if dp == 0.0:
        return 0.0
    return float(np.sign(dp) * np.sqrt(max(float(rho) * abs(float(dp)) / float(resistance), 0.0)))


def _shelltube_face_inertance(geometry: ShellTubeCoreGeometry) -> np.ndarray:
    """Return low-Mach face inertance coefficients for mass-flow states.

    For a 1D low-Mach momentum balance written in mass-flow form,

    ```text
    I_f d(mdot_f)/dt = p_left - p_right - K_f mdot_f |mdot_f|
    I_f ~= L_f / A_f
    ```

    This keeps pressure/inertance dynamics without resolving acoustic waves.
    """

    grid = geometry.grid
    area = np.asarray(grid.coolant_area, dtype=float)
    dx = np.asarray(grid.dx, dtype=float)
    n = grid.n_cells
    inertance = np.empty(n + 1, dtype=float)
    inertance[0] = 0.5 * dx[0] / max(area[0], 1.0e-30)
    for i in range(n - 1):
        length = 0.5 * (dx[i] + dx[i + 1])
        face_area = 0.5 * (area[i] + area[i + 1])
        inertance[i + 1] = length / max(face_area, 1.0e-30)
    inertance[-1] = 0.5 * dx[-1] / max(area[-1], 1.0e-30)
    return inertance


def _shelltube_low_mach_momentum_faces(
    face_old: np.ndarray,
    pressure: np.ndarray,
    density: np.ndarray,
    resistance: np.ndarray,
    inertance: np.ndarray,
    *,
    dt: float,
    inlet_pressure: float,
    outlet_pressure: float,
    flow_direction: int,
    relaxation_time_s: float = 5.0e-3,
) -> np.ndarray:
    """Advance shell-side face mass flows with low-Mach inertance/friction."""

    p = np.asarray(pressure, dtype=float)
    rho = np.asarray(density, dtype=float)
    old = np.asarray(face_old, dtype=float)
    R = np.asarray(resistance, dtype=float)
    I = np.asarray(inertance, dtype=float)
    n = p.size
    if old.size != n + 1 or I.size != n + 1:
        raise ValueError("face_old and inertance must have length n_cells + 1")
    if R.size != n:
        raise ValueError("resistance must have length n_cells")
    if dt <= 0.0:
        return old.copy()

    left_boundary = float(inlet_pressure) if flow_direction == 1 else float(outlet_pressure)
    right_boundary = float(outlet_pressure) if flow_direction == 1 else float(inlet_pressure)
    dp = np.empty(n + 1, dtype=float)
    rho_face = np.empty(n + 1, dtype=float)
    resistance_face = np.empty(n + 1, dtype=float)

    dp[0] = left_boundary - p[0]
    rho_face[0] = rho[0]
    resistance_face[0] = R[0]
    for i in range(n - 1):
        dp[i + 1] = p[i] - p[i + 1]
        rho_face[i + 1] = 0.5 * (rho[i] + rho[i + 1])
        resistance_face[i + 1] = 0.5 * (R[i] + R[i + 1])
    dp[-1] = p[-1] - right_boundary
    rho_face[-1] = rho[-1]
    resistance_face[-1] = R[-1]

    new = np.empty(n + 1, dtype=float)
    for i in range(n + 1):
        K = max(float(resistance_face[i]), 1.0e-30) / max(float(rho_face[i]), 1.0e-30)
        new[i] = _implicit_quadratic_momentum_update(
            old_mdot=float(old[i]),
            pressure_drive=float(dp[i]),
            inertance=max(float(I[i]), 1.0e-30),
            resistance_over_density=K,
            dt=dt,
        )
    boundary_drive = float(inlet_pressure) - float(outlet_pressure)
    if boundary_drive != 0.0:
        throughflow_sign = float(flow_direction) if boundary_drive > 0.0 else -float(flow_direction)
        throughflow_component = np.maximum(throughflow_sign * new, 0.0)
        new = throughflow_sign * throughflow_component
    tau = max(float(relaxation_time_s), 0.0)
    if tau > 0.0:
        theta = min(float(dt) / (tau + float(dt)), 1.0)
        new = old + theta * (new - old)
    return new


def _shelltube_low_mach_lumped_faces(
    geometry: ShellTubeCoreGeometry,
    bell_geometry: dict,
    face_old: np.ndarray,
    temperature: np.ndarray,
    pressure_profile: np.ndarray,
    *,
    dt: float,
    inlet_pressure: float,
    outlet_pressure: float,
    flow_direction: int,
    mdot_reference: float,
    corrCoeffs,
    mdot_floor: float,
) -> np.ndarray:
    """Advance one pressure-driven shell-side through-flow momentum state."""

    old = np.asarray(face_old, dtype=float)
    n_faces = geometry.grid.n_cells + 1
    if old.size != n_faces:
        raise ValueError("face_old must have length n_cells + 1")
    if dt <= 0.0:
        return old.copy()

    flow = 1 if int(flow_direction) >= 0 else -1
    old_through = max(float(flow) * float(np.nanmean(old)), 0.0)
    drive = max(float(inlet_pressure) - float(outlet_pressure), 0.0)
    resistance = _shelltube_lumped_resistance(
        geometry,
        bell_geometry,
        temperature,
        pressure_profile,
        mdot_shell=max(old_through, min(float(mdot_reference), 0.2), float(mdot_floor)),
        corrCoeffs=corrCoeffs,
        mdot_floor=mdot_floor,
    )
    mdot_new = _implicit_quadratic_momentum_update(
        old_mdot=old_through,
        pressure_drive=drive,
        inertance=_shelltube_lumped_inertance(geometry, bell_geometry),
        resistance_over_density=resistance,
        dt=dt,
    )
    return np.full(n_faces, float(flow) * max(float(mdot_new), 0.0), dtype=float)


def _shelltube_lumped_resistance(
    geometry: ShellTubeCoreGeometry,
    bell_geometry: dict,
    temperature: np.ndarray,
    pressure_profile: np.ndarray,
    *,
    mdot_shell: float,
    corrCoeffs,
    mdot_floor: float,
) -> float:
    """Return Bell-Delaware shell resistance coefficient `DeltaP/mdot^2`."""

    mdot = max(abs(float(mdot_shell)), float(mdot_floor))
    props = _coolprop_fluid_properties_at_profile(temperature, pressure_profile, "Helium")
    film = shelltube_shell_film(
        geometry,
        bell_geometry,
        mdot_shell=mdot,
        rho=props.rho,
        mu=props.mu,
        k=props.k,
        cp=props.cp,
        corrCoeffs=corrCoeffs,
        mdot_floor=mdot_floor,
    )
    dp = max(float(np.nanmean(np.abs(film.dp_shell_Pa))), 1.0)
    return dp / max(mdot * mdot, 1.0e-30)


def _shelltube_lumped_inertance(geometry: ShellTubeCoreGeometry, bell_geometry: dict) -> float:
    """Return a conservative shell-side through-flow inertance scale."""

    flow_area = max(float(bell_geometry.get("S_m", np.nan)), 1.0e-12)
    crossings = max(float(bell_geometry.get("N_baffles", 0)) + 1.0, 1.0)
    path_length = crossings * max(float(geometry.shell_inner_diameter), float(geometry.grid.length))
    return path_length / flow_area


def _implicit_quadratic_momentum_update(
    *,
    old_mdot: float,
    pressure_drive: float,
    inertance: float,
    resistance_over_density: float,
    dt: float,
) -> float:
    """Solve `I(m-new-old)/dt = dp - K*m_new*abs(m_new)`."""

    b = float(inertance) / float(dt)
    c = float(pressure_drive) + b * float(old_mdot)
    if c == 0.0:
        return 0.0
    sign = 1.0 if c > 0.0 else -1.0
    a = max(float(resistance_over_density), 0.0)
    if a <= 0.0:
        return c / b
    magnitude = (-b + np.sqrt(b * b + 4.0 * a * abs(c))) / (2.0 * a)
    return float(sign * magnitude)


def _shelltube_quasi_steady_faces(
    pressure: np.ndarray,
    density: np.ndarray,
    resistance: np.ndarray,
    *,
    mdot_inlet: float,
    outlet_pressure: float,
    flow_direction: int,
    mdot_floor: float,
) -> np.ndarray:
    """Build face mass flows with scheduled inlet and pressure-driven residual flow."""

    p = np.asarray(pressure, dtype=float)
    rho = np.asarray(density, dtype=float)
    n = p.size
    if resistance.size != n:
        raise ValueError("resistance must contain n_cells values")
    face = np.zeros(n + 1, dtype=float)
    mdot_cmd = max(float(mdot_inlet), 0.0)
    if mdot_cmd > mdot_floor:
        return np.full(n + 1, mdot_cmd if flow_direction == 1 else -mdot_cmd, dtype=float)

    if flow_direction == 1:
        face[-1] = max(_orifice_mdot_from_resistance(
            p[-1] - float(outlet_pressure),
            rho[-1],
            resistance[-1],
        ), 0.0)
    else:
        face[0] = min(_orifice_mdot_from_resistance(
            float(outlet_pressure) - p[0],
            rho[0],
            resistance[-1],
        ), 0.0)

    return face


def _outlet_pressure_at(schedule, t: float, *, default: float) -> float:
    return float(interp_schedule(schedule, t, float(default)))


def _refine_time_grid_max_step(t_grid: np.ndarray, max_dt: float) -> np.ndarray:
    """Insert internal points so no interval exceeds `max_dt`."""

    t = np.asarray(t_grid, dtype=float)
    if t.size < 2 or max_dt <= 0.0:
        return t
    pieces = [np.array([t[0]])]
    for left, right in zip(t[:-1], t[1:]):
        dt = float(right - left)
        if dt <= 0.0:
            continue
        n_sub = max(1, int(np.ceil(dt / max_dt)))
        pieces.append(np.linspace(left, right, n_sub + 1)[1:])
    return np.unique(np.round(np.concatenate(pieces), decimals=12))


def _schedule_max_abs(schedule, default: float) -> float:
    values = [abs(float(default))]
    if schedule:
        for row in schedule:
            if row is not None and len(row) >= 2:
                values.append(abs(float(row[1])))
    return max(values)


def _cfl_stable_substep_count(
    mass: np.ndarray,
    face_mdot: np.ndarray,
    dt: float,
    *,
    safety: float = 0.25,
) -> int:
    """Number of equal substeps needed to keep the explicit coolant mass/energy
    advection (`conservative_mass_energy_step`, forward Euler in the conserved
    variables) within its stability limit for this macro step.

    This scheme is unconditionally unstable once a step advances a cell by more
    than roughly its residence time `mass / mdot` -- confirmed empirically
    2026-08-18 on the shell-and-tube bang-bang validation case: stable at
    dt/tau <= ~0.2, a fast-growing single-cell spike at dt/tau ~0.4, and a
    `FloatingPointError` from `enforce_internal_energy_bounds` within a handful
    of steps at dt/tau > ~1. `_limit_face_mdot_for_inventory` guards a DIFFERENT
    failure mode (a macro step draining a cell's mass net-zero) and is only
    engaged near valve closure; it does not help here, where mass stays exactly
    conserved throughout and the instability is in the per-cell TEMPERATURE
    field advecting faster than the step can resolve.

    Uses the worst case across cells (max |face flow|, min cell mass) rather
    than `mdot_effective`'s mean, since a single under-resolved cell is enough
    to trip the energy-bounds guard.
    """
    if dt <= 0.0:
        return 1
    m = np.asarray(mass, dtype=float)
    if m.size == 0:
        return 1
    mdot_max = float(np.max(np.abs(np.asarray(face_mdot, dtype=float))))
    if mdot_max <= 0.0:
        return 1
    tau_min = float(np.min(m)) / mdot_max
    if tau_min <= 0.0:
        return 1
    return max(int(np.ceil(dt / (safety * tau_min))), 1)


def _limit_face_mdot_for_inventory(
    mass: np.ndarray,
    face_mdot: np.ndarray,
    dt: float,
    *,
    internal_energy_J: np.ndarray | None = None,
    specific_enthalpy_J_kg: np.ndarray | None = None,
    keep_fraction: float = 0.05,
) -> np.ndarray:
    """Limit face flows so one macro step cannot empty a cell.

    This is a bounded-cost substitute for resolving every cell residence time.
    It preserves the sign pattern and still allows residual discharge after
    valve closure, but prevents explicit conserved-variable updates from leaving
    the thermodynamic domain when `max_step` is much larger than cell residence.
    """

    if dt <= 0.0:
        return np.asarray(face_mdot, dtype=float).copy()
    limited = np.asarray(face_mdot, dtype=float).copy()
    m = np.asarray(mass, dtype=float)
    U = None if internal_energy_J is None else np.asarray(internal_energy_J, dtype=float)
    h = None if specific_enthalpy_J_kg is None else np.asarray(specific_enthalpy_J_kg, dtype=float)
    min_remaining = np.maximum(float(keep_fraction), 0.0) * m
    removable = np.maximum(m - min_remaining, 0.0)
    removable_U = None if U is None else np.maximum(U - np.maximum(float(keep_fraction), 0.0) * U, 0.0)
    n = m.size
    for _ in range(4):
        changed = False
        for i in range(n):
            outgoing = max(-limited[i], 0.0) + max(limited[i + 1], 0.0)
            max_outgoing = removable[i] / dt
            scale = 1.0
            if outgoing > max_outgoing > 0.0:
                scale = min(scale, max_outgoing / outgoing)
            if U is not None and h is not None:
                outgoing_energy_rate = outgoing * max(float(h[i]), 1.0e-30)
                max_energy_rate = removable_U[i] / dt
                if outgoing_energy_rate > max_energy_rate > 0.0:
                    scale = min(scale, max_energy_rate / outgoing_energy_rate)
                elif outgoing_energy_rate > 0.0 and max_energy_rate <= 0.0:
                    scale = 0.0
            if scale < 1.0:
                if limited[i] < 0.0:
                    limited[i] *= scale
                if limited[i + 1] > 0.0:
                    limited[i + 1] *= scale
                changed = True
            elif outgoing > 0.0 and max_outgoing <= 0.0:
                if limited[i] < 0.0:
                    limited[i] = 0.0
                if limited[i + 1] > 0.0:
                    limited[i + 1] = 0.0
                changed = True
        if not changed:
            break
    return limited


def _limit_closed_valve_outlet_discharge_to_backpressure(
    geometry: ShellTubeCoreGeometry,
    state: CoolantThermodynamicState,
    mass: np.ndarray,
    face_mdot: np.ndarray,
    *,
    dt: float,
    outlet_pressure: float,
    flow_direction: int,
) -> np.ndarray:
    """Limit closed-inlet residual discharge by the downstream pressure state."""

    if dt <= 0.0:
        return np.asarray(face_mdot, dtype=float).copy()
    limited = np.asarray(face_mdot, dtype=float).copy()
    outlet_cell = geometry.grid.outlet_index
    outlet_face = -1 if int(flow_direction) == 1 else 0
    mdot = float(limited[outlet_face])
    outgoing = mdot > 0.0 if outlet_face == limited.size - 1 else mdot < 0.0
    if not outgoing:
        return limited

    try:
        rho_floor = PropsSI(
            "D",
            "T",
            float(state.temperature[outlet_cell]),
            "P",
            max(float(outlet_pressure), 1.0e3),
            "Helium",
        )
    except Exception:
        return limited
    if not np.isfinite(rho_floor) or rho_floor <= 0.0:
        return limited

    target_mass = float(rho_floor) * float(geometry.grid.coolant_volume[outlet_cell])
    removable = max(float(mass[outlet_cell]) - target_mass, 0.0)
    max_abs_mdot = removable / float(dt)
    if abs(mdot) > max_abs_mdot:
        limited[outlet_face] = np.sign(mdot) * max_abs_mdot
    return limited


def run_shelltube_transient_core(
    geometry: ShellTubeCoreGeometry,
    bell_geometry: dict,
    *,
    T_wall_initial,
    T_coolant_initial,
    t_end: float,
    max_step: float,
    coolant_properties_at: CoolantPropertyProvider,
    wall_density: float,
    wall_cp,
    wall_conductivity_at_T,
    inside_tube_choice: str,
    nusselt_selector: str,
    tube_roughness: float,
    mdot_coolant_default: float,
    T_coolant_inlet_default: float,
    p_coolant_default: float,
    mdot_hot_total_default: float,
    mdot_coolant_schedule=None,
    T_coolant_inlet_schedule=None,
    p_coolant_schedule=None,
    mdot_hot_total_schedule=None,
    gas_state_at: GasStateProvider | None = None,
    gas_provider_at_time: GasProviderAtTime | None = None,
    mdot_hot_total_at_time: ScalarAtTime | None = None,
    hot_side_schedules=(),
    t_eval=None,
    n_save: int | None = None,
    corrCoeffs=None,
    corrugation_thickness: float = 0.0,
    corrugation_pitch: float = 1.0,
    progress_initial: float = 0.0,
    flow_direction: int | None = None,
    h_g_rad=None,
    mdot_floor: float = 1e-12,
    coolant_state_model: str = "temperature",
    progress_config=None,
    progress_enabled: bool = True,
) -> ShellTubeTransientCoreResult:
    """Run the shell-and-tube transient-core wall/coolant model.

    This is an adapter-level runner, not yet the user-facing `main_transient.py`
    dispatch. Schedules use the project `(time_s, value)` convention and are
    inserted into the fixed-step grid.

    Use `gas_state_at` for a single combustion model over the full run. Use
    `gas_provider_at_time(t)` and `mdot_hot_total_at_time(t)` when the hot-side
    model changes with time, for example GOX chilldown before ignition followed
    by finite-rate combustion.
    """

    if gas_provider_at_time is None and gas_state_at is None:
        raise ValueError("Provide gas_state_at or gas_provider_at_time")

    if n_save is not None:
        if int(n_save) < 2:
            raise ValueError("n_save must be at least 2")
        save_times = np.linspace(0.0, float(t_end), int(n_save))
        if t_eval is None:
            t_eval = save_times
        else:
            t_eval = np.unique(np.concatenate([np.asarray(t_eval, dtype=float), save_times]))

    schedule_breakpoints = (
        (
            mdot_coolant_schedule,
            T_coolant_inlet_schedule,
            p_coolant_schedule,
            getattr(progress_config, "schedule_p_c_out", None),
            mdot_hot_total_schedule,
            *tuple(hot_side_schedules or ()),
        )
        if getattr(progress_config, "insert_schedule_breakpoints", True) else
        ()
    )
    t_grid = fixed_time_grid(
        t_end=float(t_end),
        max_step=float(max_step),
        t_eval=t_eval,
        schedules=schedule_breakpoints,
    )
    diagnostics: list[ShellTubeStepInputDiagnostics] = []
    progress = TransientProgressPrinter.from_config(
        progress_config,
        total_steps=max(len(t_grid) - 1, 0),
    ) if progress_config is not None and progress_enabled else None
    if progress is not None and coolant_state_model == "mass_energy":
        progress.pressure_label = "dpHe [bar]"

    if coolant_state_model in ("mass_energy", "low_mach_momentum"):
        return _run_shelltube_transient_core_mass_energy(
            geometry,
            bell_geometry,
            T_wall_initial=T_wall_initial,
            T_coolant_initial=T_coolant_initial,
            t_grid=t_grid,
            coolant_properties_at=coolant_properties_at,
            wall_density=wall_density,
            wall_cp=wall_cp,
            wall_conductivity_at_T=wall_conductivity_at_T,
            inside_tube_choice=inside_tube_choice,
            nusselt_selector=nusselt_selector,
            tube_roughness=tube_roughness,
            mdot_coolant_default=mdot_coolant_default,
            T_coolant_inlet_default=T_coolant_inlet_default,
            p_coolant_default=p_coolant_default,
            mdot_hot_total_default=mdot_hot_total_default,
            mdot_coolant_schedule=mdot_coolant_schedule,
            T_coolant_inlet_schedule=T_coolant_inlet_schedule,
            p_coolant_schedule=p_coolant_schedule,
            mdot_hot_total_schedule=mdot_hot_total_schedule,
            gas_state_at=gas_state_at,
            gas_provider_at_time=gas_provider_at_time,
            mdot_hot_total_at_time=mdot_hot_total_at_time,
            corrCoeffs=corrCoeffs,
            corrugation_thickness=corrugation_thickness,
            corrugation_pitch=corrugation_pitch,
            progress_initial=progress_initial,
            flow_direction=flow_direction,
            h_g_rad=h_g_rad,
            mdot_floor=mdot_floor,
            progress=progress,
            momentum_model=(
                "low_mach"
                if coolant_state_model == "low_mach_momentum" else
                "quasi_steady"
            ),
            outlet_pressure_schedule=getattr(progress_config, "schedule_p_c_out", None),
            outlet_pressure_fixed=getattr(
                progress_config,
                "transient_coolant_outlet_pressure",
                None,
            ),
        )
    if coolant_state_model != "temperature":
        raise ValueError(
            "coolant_state_model must be 'temperature', 'mass_energy', "
            "or 'low_mach_momentum'"
        )

    def step_builder(t: float, T_wall: np.ndarray, T_coolant: np.ndarray) -> WallCoolantStepInputs:
        provider_t = gas_state_at
        progress_initial_t = progress_initial
        if gas_provider_at_time is not None:
            provider_t, progress_initial_t = gas_provider_at_time(float(t))
        if provider_t is None:
            raise ValueError("gas_provider_at_time returned no gas-state provider")
        mdot_hot_t = (
            float(mdot_hot_total_at_time(float(t)))
            if mdot_hot_total_at_time is not None else
            interp_schedule(mdot_hot_total_schedule, t, mdot_hot_total_default)
        )

        assembled = shelltube_step_inputs(
            geometry,
            bell_geometry,
            Tbar_wall=T_wall,
            T_coolant=T_coolant,
            mdot_coolant=interp_schedule(mdot_coolant_schedule, t, mdot_coolant_default),
            T_coolant_inlet=interp_schedule(
                T_coolant_inlet_schedule,
                t,
                T_coolant_inlet_default,
            ),
            p_coolant=interp_schedule(p_coolant_schedule, t, p_coolant_default),
            mdot_hot_total=mdot_hot_t,
            gas_state_at=provider_t,
            coolant_properties_at=coolant_properties_at,
            wall_density=wall_density,
            wall_cp=wall_cp,
            wall_conductivity_at_T=wall_conductivity_at_T,
            inside_tube_choice=inside_tube_choice,
            nusselt_selector=nusselt_selector,
            tube_roughness=tube_roughness,
            corrCoeffs=corrCoeffs,
            corrugation_thickness=corrugation_thickness,
            corrugation_pitch=corrugation_pitch,
            progress_initial=progress_initial_t,
            flow_direction=flow_direction,
            h_g_rad=h_g_rad,
        )
        diagnostics.append(assembled)
        return assembled.wall_coolant_inputs

    integration = integrate_wall_coolant_fixed_step(
        T_wall_initial=T_wall_initial,
        T_coolant_initial=T_coolant_initial,
        t_eval=t_grid,
        step_inputs=step_builder,
        mdot_floor=mdot_floor,
    )
    return ShellTubeTransientCoreResult(
        integration=integration,
        step_diagnostics=tuple(diagnostics),
    )


def _run_shelltube_transient_core_mass_energy(
    geometry: ShellTubeCoreGeometry,
    bell_geometry: dict,
    *,
    T_wall_initial,
    T_coolant_initial,
    t_grid: np.ndarray,
    coolant_properties_at: CoolantPropertyProvider,
    wall_density: float,
    wall_cp,
    wall_conductivity_at_T,
    inside_tube_choice: str,
    nusselt_selector: str,
    tube_roughness: float,
    mdot_coolant_default: float,
    T_coolant_inlet_default: float,
    p_coolant_default: float,
    mdot_hot_total_default: float,
    mdot_coolant_schedule=None,
    T_coolant_inlet_schedule=None,
    p_coolant_schedule=None,
    mdot_hot_total_schedule=None,
    gas_state_at: GasStateProvider | None = None,
    gas_provider_at_time: GasProviderAtTime | None = None,
    mdot_hot_total_at_time: ScalarAtTime | None = None,
    corrCoeffs=None,
    corrugation_thickness: float = 0.0,
    corrugation_pitch: float = 1.0,
    progress_initial: float = 0.0,
    flow_direction: int | None = None,
    h_g_rad=None,
    mdot_floor: float = 1e-12,
    progress: TransientProgressPrinter | None = None,
    momentum_model: str = "quasi_steady",
    outlet_pressure_schedule=None,
    outlet_pressure_fixed: float | None = None,
) -> ShellTubeTransientCoreResult:
    """Run wall + compressible coolant mass/energy with quasi-steady momentum."""

    grid = geometry.grid
    flow = grid.flow_direction if flow_direction is None else int(flow_direction)
    if flow not in (-1, 1):
        raise ValueError("flow_direction must be +1 or -1")

    T_wall0 = _cell_property("T_wall_initial", T_wall_initial, grid.n_cells)
    T_coolant0 = _cell_property("T_coolant_initial", T_coolant_initial, grid.n_cells)
    p0_inlet = float(interp_schedule(p_coolant_schedule, 0.0, p_coolant_default))
    momentum_model = str(momentum_model)
    if momentum_model not in ("quasi_steady", "low_mach"):
        raise ValueError("momentum_model must be 'quasi_steady' or 'low_mach'")
    mdot0 = max(float(interp_schedule(mdot_coolant_schedule, 0.0, mdot_coolant_default)), 0.0)
    mdot_reference = max(
        _schedule_max_abs(mdot_coolant_schedule, mdot_coolant_default),
        mdot0,
        float(mdot_coolant_default),
        mdot_floor,
    )
    dp_nominal = _shelltube_nominal_pressure_drop(
        geometry,
        bell_geometry,
        mdot_shell=mdot_reference,
        T_coolant=T_coolant0,
        p_coolant=p0_inlet,
        coolant_properties_at=coolant_properties_at,
        corrCoeffs=corrCoeffs,
        mdot_floor=mdot_floor,
    )
    p_outlet0 = _outlet_pressure_at(
        outlet_pressure_schedule,
        0.0,
        default=(
            outlet_pressure_fixed
            if outlet_pressure_fixed is not None else
            max(p0_inlet - dp_nominal, 1.0e3)
        ),
    )
    if momentum_model == "low_mach":
        p_initial = _shelltube_boundary_pressure_profile(
            geometry,
            inlet_pressure=p0_inlet,
            outlet_pressure=p_outlet0,
            flow_direction=flow,
        )
    else:
        p_initial = _shelltube_initial_pressure_profile(
            geometry,
            inlet_pressure=p0_inlet,
            pressure_drop=dp_nominal,
            flow_direction=flow,
        )
    p_outlet = max(float(p_outlet0), 1.0e3)
    mass, internal_energy = initial_mass_energy_from_TP(
        T_coolant0,
        p_initial,
        grid.coolant_volume,
        "Helium",
    )
    state0 = coolprop_state_from_mass_energy(mass, internal_energy, grid.coolant_volume, "Helium")
    resistance = _shelltube_face_resistance(
        dp_nominal,
        state0.density,
        mdot_reference,
        n_faces=grid.n_cells,
    )
    inertance = _shelltube_face_inertance(geometry)

    t = np.asarray(t_grid, dtype=float)
    n_time = t.size
    n_cells = grid.n_cells
    T_wall = np.zeros((n_time, n_cells), dtype=float)
    T_coolant = np.zeros((n_time, n_cells), dtype=float)
    T_out = np.zeros(n_time, dtype=float)
    coolant_mass = np.zeros((n_time, n_cells), dtype=float)
    coolant_U = np.zeros((n_time, n_cells), dtype=float)
    pressure = np.zeros((n_time, n_cells), dtype=float)
    density = np.zeros((n_time, n_cells), dtype=float)
    enthalpy = np.zeros((n_time, n_cells), dtype=float)
    face_mdot = np.zeros((n_time, n_cells + 1), dtype=float)
    hot_heat_added_J = np.zeros(n_time, dtype=float)
    adv_in_J = np.zeros(n_time, dtype=float)
    adv_out_J = np.zeros(n_time, dtype=float)
    energy_residual_J = np.zeros(n_time, dtype=float)
    mass_residual_kg = np.zeros(n_time, dtype=float)
    heat_wall_to_coolant_W = np.zeros((n_time, n_cells), dtype=float)

    T_wall[0] = T_wall0
    coolant_mass[0] = mass
    coolant_U[0] = internal_energy
    T_coolant[0] = state0.temperature
    pressure[0] = state0.pressure
    density[0] = state0.density
    enthalpy[0] = state0.specific_enthalpy_J_kg
    T_out[0] = T_coolant[0, grid.outlet_index]

    diagnostics: list[ShellTubeStepInputDiagnostics] = []
    last_step = None

    for j in range(n_time - 1):
        tj = float(t[j])
        dt = float(t[j + 1] - t[j])
        state = coolprop_state_from_mass_energy(
            coolant_mass[j],
            coolant_U[j],
            grid.coolant_volume,
            "Helium",
        )
        mdot_cmd = max(float(interp_schedule(mdot_coolant_schedule, tj, mdot_coolant_default)), 0.0)
        p_inlet = float(interp_schedule(p_coolant_schedule, tj, p_coolant_default))
        p_outlet = max(
            _outlet_pressure_at(
                outlet_pressure_schedule,
                tj,
                default=(
                    outlet_pressure_fixed
                    if outlet_pressure_fixed is not None else
                    p_outlet
                ),
            ),
            1.0e3,
        )
        if momentum_model == "low_mach":
            p_transport = _shelltube_boundary_pressure_profile(
                geometry,
                inlet_pressure=p_inlet,
                outlet_pressure=p_outlet,
                flow_direction=flow,
            )
            faces = _shelltube_low_mach_lumped_faces(
                geometry,
                bell_geometry,
                face_mdot[j],
                state.temperature,
                p_transport,
                dt=dt,
                inlet_pressure=p_inlet,
                outlet_pressure=p_outlet,
                flow_direction=flow,
                mdot_reference=mdot_reference,
                corrCoeffs=corrCoeffs,
                mdot_floor=mdot_floor,
            )
        else:
            faces = _shelltube_quasi_steady_faces(
                state.pressure,
                state.density,
                resistance,
                mdot_inlet=mdot_cmd,
                outlet_pressure=p_outlet,
                flow_direction=flow,
                mdot_floor=mdot_floor,
            )
        mdot_scale = _schedule_max_abs(mdot_coolant_schedule, mdot_coolant_default)
        if momentum_model == "low_mach":
            mdot_cap = max(mdot_cmd, mdot_floor)
        else:
            mdot_cap = max(2.0 * mdot_scale, mdot_floor)
        faces = np.clip(faces, -mdot_cap, mdot_cap)
        if momentum_model != "low_mach" and mdot_cmd <= mdot_floor and j > 0:
            outlet_face = -1 if flow == 1 else 0
            outlet_cell = grid.outlet_index
            outlet_drive = float(state.pressure[outlet_cell]) - float(p_outlet)
            if flow == -1:
                outlet_drive = float(p_outlet) - float(state.pressure[outlet_cell])
            if outlet_drive > 0.0 and abs(faces[outlet_face]) <= mdot_floor:
                faces[outlet_face] = 0.5 * face_mdot[j, outlet_face]
            faces = _limit_closed_valve_outlet_discharge_to_backpressure(
                geometry,
                state,
                coolant_mass[j],
                faces,
                dt=dt,
                outlet_pressure=p_outlet,
                flow_direction=flow,
            )
        if momentum_model != "low_mach" and mdot_cmd <= mdot_floor:
            faces = _limit_face_mdot_for_inventory(
                coolant_mass[j],
                faces,
                dt,
                internal_energy_J=coolant_U[j],
                specific_enthalpy_J_kg=state.specific_enthalpy_J_kg,
            )
        mdot_effective = max(float(np.mean(np.abs(faces))), mdot_floor)
        if momentum_model == "low_mach":
            coolant_props_for_film = _coolprop_fluid_properties_at_profile(
                state.temperature,
                p_transport,
                "Helium",
            )
        else:
            p_transport = _shelltube_boundary_pressure_profile(
                geometry,
                inlet_pressure=p_inlet,
                outlet_pressure=p_outlet,
                flow_direction=flow,
            )
            coolant_props_for_film = _coolprop_fluid_properties_at_profile(
                state.temperature,
                p_transport,
                "Helium",
            )

        provider_t = gas_state_at
        progress_initial_t = progress_initial
        if gas_provider_at_time is not None:
            provider_t, progress_initial_t = gas_provider_at_time(tj)
        if provider_t is None:
            raise ValueError("gas_provider_at_time returned no gas-state provider")
        mdot_hot_t = (
            float(mdot_hot_total_at_time(tj))
            if mdot_hot_total_at_time is not None else
            interp_schedule(mdot_hot_total_schedule, tj, mdot_hot_total_default)
        )
        T_inlet = float(interp_schedule(T_coolant_inlet_schedule, tj, T_coolant_inlet_default))

        assembled = shelltube_step_inputs(
            geometry,
            bell_geometry,
            Tbar_wall=T_wall[j],
            T_coolant=state.temperature,
            mdot_coolant=mdot_effective,
            T_coolant_inlet=T_inlet,
            p_coolant=p_inlet,
            mdot_hot_total=mdot_hot_t,
            gas_state_at=provider_t,
            coolant_properties_at=lambda _T, _p, props=coolant_props_for_film: props,
            wall_density=wall_density,
            wall_cp=wall_cp,
            wall_conductivity_at_T=wall_conductivity_at_T,
            inside_tube_choice=inside_tube_choice,
            nusselt_selector=nusselt_selector,
            tube_roughness=tube_roughness,
            corrCoeffs=corrCoeffs,
            corrugation_thickness=corrugation_thickness,
            corrugation_pitch=corrugation_pitch,
            progress_initial=progress_initial_t,
            flow_direction=flow,
            h_g_rad=h_g_rad,
        )
        diagnostics.append(assembled)

        h_inlet = float(PropsSI("H", "T", T_inlet, "P", max(p_inlet, 1.0e3), "Helium"))

        # Subcycle the explicit coolant mass/energy advection when one macro
        # step would exceed its CFL-stable limit (see _cfl_stable_substep_count
        # docstring: forward-Euler advection, unconditionally unstable past
        # roughly one cell residence time). faces/hot_heat_W/conductance stay
        # frozen across substeps -- the same "quasi-steady over the fixed step"
        # assumption already used for the hot-gas march and momentum, just
        # applied at a finer grain for the part of the update that actually
        # needs it.
        n_sub = _cfl_stable_substep_count(coolant_mass[j], faces, dt)
        sub_dt = dt / n_sub
        Tw_cur = T_wall[j]
        m_cur = coolant_mass[j]
        U_cur = coolant_U[j]
        T_cur = state.temperature
        h_cur = state.specific_enthalpy_J_kg
        hot_heat_acc = 0.0
        adv_in_acc = 0.0
        adv_out_acc = 0.0
        energy_residual_acc = 0.0
        mass_residual_acc = 0.0
        thermal_step = None
        for _sub in range(n_sub):
            thermal_step = semi_implicit_wall_compressible_coolant_step(
                Tw_cur,
                assembled.wall_heat_capacity_J_K,
                m_cur,
                U_cur,
                T_cur,
                h_cur,
                faces,
                assembled.wall_coolant_inputs.hot_heat_W,
                assembled.wall_coolant_inputs.wall_to_coolant_conductance_W_per_K,
                sub_dt,
                inlet_enthalpy_J_kg=h_inlet,
                outlet_backflow_enthalpy_J_kg=h_inlet,
                mass_floor=1.0e-12,
            )
            m_candidate, U_candidate = enforce_density_bounds(
                np.maximum(thermal_step.coolant.mass_new, 1.0e-12),
                thermal_step.coolant.internal_energy_new_J,
                grid.coolant_volume,
            )
            if momentum_model == "low_mach":
                provisional_U = enforce_internal_energy_floor(
                    m_candidate,
                    U_candidate,
                    grid.coolant_volume,
                    "Helium",
                    clip=False,
                )
                provisional_state = coolprop_state_from_mass_energy(
                    m_candidate,
                    provisional_U,
                    grid.coolant_volume,
                    "Helium",
                )
                p_projected = _shelltube_boundary_pressure_profile(
                    geometry,
                    inlet_pressure=p_inlet,
                    outlet_pressure=p_outlet,
                    flow_direction=flow,
                )
                m_sub_new, U_sub_new = _coolant_mass_energy_from_TP_profile(
                    provisional_state.temperature,
                    p_projected,
                    grid.coolant_volume,
                    "Helium",
                )
            else:
                m_sub_new = m_candidate
                U_sub_new = enforce_internal_energy_floor(
                    m_sub_new,
                    U_candidate,
                    grid.coolant_volume,
                    "Helium",
                    clip=False,
                )
            hot_heat_acc += thermal_step.hot_heat_added_J
            adv_in_acc += thermal_step.coolant.advective_energy_in_J
            adv_out_acc += thermal_step.coolant.advective_energy_out_J
            energy_residual_acc += thermal_step.total_energy_residual_J
            mass_residual_acc += thermal_step.coolant.mass_residual_kg
            Tw_cur = thermal_step.T_wall_new
            m_cur = m_sub_new
            U_cur = U_sub_new
            if _sub < n_sub - 1:
                sub_state = coolprop_state_from_mass_energy(
                    m_cur, U_cur, grid.coolant_volume, "Helium"
                )
                T_cur = sub_state.temperature
                h_cur = sub_state.specific_enthalpy_J_kg

        m_new, U_new = m_cur, U_cur
        new_state = coolprop_state_from_mass_energy(
            m_new,
            U_new,
            grid.coolant_volume,
            "Helium",
        )
        T_wall[j + 1] = Tw_cur
        coolant_mass[j + 1] = m_new
        coolant_U[j + 1] = U_new
        T_coolant[j + 1] = new_state.temperature
        pressure[j + 1] = new_state.pressure
        density[j + 1] = new_state.density
        enthalpy[j + 1] = new_state.specific_enthalpy_J_kg
        face_mdot[j + 1] = faces
        T_out[j + 1] = new_state.temperature[grid.outlet_index]
        hot_heat_added_J[j + 1] = hot_heat_acc
        adv_in_J[j + 1] = adv_in_acc
        adv_out_J[j + 1] = adv_out_acc
        energy_residual_J[j + 1] = energy_residual_acc
        mass_residual_kg[j + 1] = mass_residual_acc
        heat_wall_to_coolant_W[j + 1] = thermal_step.heat_wall_to_coolant_W
        last_step = thermal_step
        if progress is not None:
            pressure_report = (
                float(np.nanmean(assembled.shell_film.dp_shell_Pa))
                if momentum_model != "low_mach" else
                float(new_state.pressure[grid.outlet_index])
            )
            progress.update(
                step=j + 1,
                time_s=float(t[j + 1]),
                T_wall=T_wall[j + 1],
                T_coolant_outlet=T_out[j + 1],
                p_coolant_outlet=pressure_report,
                T_gas_outlet=assembled.hot_gas_march.T_gas_outlet,
            )

    integration = ShellTubeCompressibleIntegrationResult(
        t=t,
        T_wall=T_wall,
        T_coolant=T_coolant,
        T_coolant_outlet=T_out,
        coolant_mass_kg=coolant_mass,
        coolant_internal_energy_J=coolant_U,
        coolant_pressure_Pa=pressure,
        coolant_density_kg_m3=density,
        coolant_specific_enthalpy_J_kg=enthalpy,
        face_mdot_kg_s=face_mdot,
        hot_heat_added_J=hot_heat_added_J,
        advective_energy_in_J=adv_in_J,
        advective_energy_out_J=adv_out_J,
        energy_residual_J=energy_residual_J,
        mass_residual_kg=mass_residual_kg,
        heat_wall_to_coolant_W=heat_wall_to_coolant_W,
        last_step=last_step,
    )
    return ShellTubeTransientCoreResult(
        integration=integration,
        step_diagnostics=tuple(diagnostics),
    )


def coolprop_fluid_properties(fluid: str) -> CoolantPropertyProvider:
    """Build a cellwise property provider backed by CoolProp."""

    def provider(T, pressure: float) -> ShellTubeFluidProperties:
        T_arr = np.asarray(T, dtype=float)
        return ShellTubeFluidProperties(
            rho=np.array([PropsSI("D", "T", float(Ti), "P", pressure, fluid) for Ti in T_arr]),
            mu=np.array([PropsSI("V", "T", float(Ti), "P", pressure, fluid) for Ti in T_arr]),
            k=np.array([PropsSI("L", "T", float(Ti), "P", pressure, fluid) for Ti in T_arr]),
            cp=np.array([PropsSI("C", "T", float(Ti), "P", pressure, fluid) for Ti in T_arr]),
        )

    return provider


def _coolprop_fluid_properties_at_profile(T, pressure, fluid: str) -> ShellTubeFluidProperties:
    """Return CoolProp properties from cellwise temperature and pressure arrays."""

    T_arr = np.asarray(T, dtype=float)
    p_arr = np.asarray(pressure, dtype=float)
    if T_arr.ndim != 1 or p_arr.ndim != 1 or T_arr.size != p_arr.size:
        raise ValueError("T and pressure must be one-dimensional arrays of equal length")
    return ShellTubeFluidProperties(
        rho=np.array([
            PropsSI("D", "T", float(Ti), "P", max(float(pi), 1.0e3), fluid)
            for Ti, pi in zip(T_arr, p_arr)
        ]),
        mu=np.array([
            PropsSI("V", "T", float(Ti), "P", max(float(pi), 1.0e3), fluid)
            for Ti, pi in zip(T_arr, p_arr)
        ]),
        k=np.array([
            PropsSI("L", "T", float(Ti), "P", max(float(pi), 1.0e3), fluid)
            for Ti, pi in zip(T_arr, p_arr)
        ]),
        cp=np.array([
            PropsSI("C", "T", float(Ti), "P", max(float(pi), 1.0e3), fluid)
            for Ti, pi in zip(T_arr, p_arr)
        ]),
    )


def _coolant_mass_energy_from_TP_profile(T, pressure, volume, fluid: str) -> tuple[np.ndarray, np.ndarray]:
    """Return cell mass and internal energy from cellwise T/p/volume arrays."""

    T_arr = np.asarray(T, dtype=float)
    p_arr = np.asarray(pressure, dtype=float)
    V_arr = np.asarray(volume, dtype=float)
    if T_arr.ndim != 1 or p_arr.ndim != 1 or V_arr.ndim != 1:
        raise ValueError("T, pressure, and volume must be one-dimensional arrays")
    if T_arr.size != p_arr.size or T_arr.size != V_arr.size:
        raise ValueError("T, pressure, and volume must have the same length")
    rho = np.array([
        PropsSI("D", "T", float(Ti), "P", max(float(pi), 1.0e3), fluid)
        for Ti, pi in zip(T_arr, p_arr)
    ])
    u = np.array([
        PropsSI("U", "T", float(Ti), "P", max(float(pi), 1.0e3), fluid)
        for Ti, pi in zip(T_arr, p_arr)
    ])
    mass = rho * V_arr
    return mass, mass * u


def fpv_gas_state_provider(fpv) -> tuple[GasStateProvider, float]:
    """Return a gas-state provider backed by an `FPVManifold`.

    The returned initial progress is `fpv.Yc_inlet()`, matching the maintained
    finite-rate transient convention.
    """

    progress_initial = float(fpv.Yc_inlet())

    def provider(h_removed: float, progress: float, _i: int) -> ShellTubeGasState:
        T, rho, mu, k, cp, _xH2O, _xCO2, omega = fpv.state(h_removed, progress)
        return ShellTubeGasState(T=T, rho=rho, mu=mu, k=k, cp=cp, progress_source=omega)

    return provider, progress_initial


def equilibrium_gas_state_provider(manifold) -> tuple[GasStateProvider, float]:
    """Return a gas-state provider backed by an equilibrium/frozen manifold."""

    def provider(h_removed: float, _progress: float, _i: int) -> ShellTubeGasState:
        T, rho, mu, k, cp, _xH2O, _xCO2 = manifold.at(h_removed)
        return ShellTubeGasState(T=T, rho=rho, mu=mu, k=k, cp=cp, progress_source=0.0)

    return provider, 0.0


def oxygen_gas_state_provider(
    *,
    T_inlet: float,
    pressure: float,
    fluid: str = "Oxygen",
    T_min: float = 95.0,
    T_max: float = 1200.0,
) -> tuple[GasStateProvider, float]:
    """Return a pre-ignition oxygen sensible-cooling gas-state provider.

    `h_removed` is interpreted as specific enthalpy removed from the incoming
    oxygen stream. Temperature is recovered from `(H, P)` through CoolProp and
    clipped to the requested bounds before transport properties are evaluated.
    """

    if T_inlet <= 0.0 or pressure <= 0.0:
        raise ValueError("T_inlet and pressure must be positive")
    if T_min <= 0.0 or T_max <= T_min:
        raise ValueError("temperature bounds are invalid")

    T0 = float(np.clip(T_inlet, T_min, T_max))
    h0 = float(PropsSI("H", "T", T0, "P", pressure, fluid))

    def provider(h_removed: float, _progress: float, _i: int) -> ShellTubeGasState:
        h = h0 - max(float(h_removed), 0.0)
        T = float(PropsSI("T", "H", h, "P", pressure, fluid))
        T = float(np.clip(T, T_min, T_max))
        return ShellTubeGasState(
            T=T,
            rho=float(PropsSI("D", "T", T, "P", pressure, fluid)),
            mu=float(PropsSI("V", "T", T, "P", pressure, fluid)),
            k=float(PropsSI("L", "T", T, "P", pressure, fluid)),
            cp=float(PropsSI("C", "T", T, "P", pressure, fluid)),
            progress_source=0.0,
        )

    return provider, 0.0


def _tube_surface_factors(choice: str, corrCoeffs) -> tuple[float, float]:
    if choice == "smooth":
        return 1.0, 1.0
    if choice == "grooved":
        return (
            getattr(corrCoeffs, "tube_grooved_Nu_factor", 1.0),
            getattr(corrCoeffs, "tube_grooved_f_factor", 1.0),
        )
    if choice == "intensification_factor":
        factor = getattr(corrCoeffs, "tube_intensification_factor", 1.0)
        return factor, factor
    if choice in ("helical_insert", "power_law"):
        raise NotImplementedError(
            f"inside_tube_choice={choice!r} is exposed but not implemented yet"
        )
    raise ValueError(f"Unsupported inside_tube_choice: {choice!r}")


def _corrugation_severity(thickness: float, pitch: float, tube_inner_diameter: float) -> float:
    e = max(float(thickness), 0.0)
    p = max(float(pitch), 1e-12)
    return e**2 / (p * max(float(tube_inner_diameter), 1e-12))


def _single_cell_shelltube_geometry(
    geometry: ShellTubeCoreGeometry,
    cell_index: int,
) -> ShellTubeCoreGeometry:
    grid = geometry.grid
    i = int(cell_index)
    cell_grid = AxialGrid(
        x_edges=np.array([0.0, float(grid.dx[i])]),
        coolant_area=np.array([float(grid.coolant_area[i])]),
        wall_area=np.array([float(grid.wall_area[i])]),
        hot_perimeter=np.array([float(grid.hot_perimeter[i])]),
        coolant_perimeter=np.array([float(grid.coolant_perimeter[i])]),
        flow_direction=grid.flow_direction,
    )
    return ShellTubeCoreGeometry(
        grid=cell_grid,
        tube_inner_diameter=geometry.tube_inner_diameter,
        tube_outer_diameter=geometry.tube_outer_diameter,
        wall_thickness=geometry.wall_thickness,
        n_tubes=geometry.n_tubes,
        shell_inner_diameter=geometry.shell_inner_diameter,
    )


def _coerce_gas_state(value) -> ShellTubeGasState:
    if isinstance(value, ShellTubeGasState):
        state = value
    elif isinstance(value, dict):
        state = ShellTubeGasState(
            T=value["T"],
            rho=value["rho"],
            mu=value["mu"],
            k=value["k"],
            cp=value["cp"],
            progress_source=value.get("progress_source", 0.0),
        )
    else:
        raise TypeError("gas_state_at must return ShellTubeGasState or dict")
    vals = np.array([state.T, state.rho, state.mu, state.k, state.cp], dtype=float)
    if not np.all(np.isfinite(vals)):
        raise ValueError("gas state contains non-finite values")
    if state.T <= 0.0 or state.rho <= 0.0 or state.mu <= 0.0:
        raise ValueError("gas T, rho, and mu must be positive")
    if state.k <= 0.0 or state.cp <= 0.0:
        raise ValueError("gas k and cp must be positive")
    return state


def _coerce_fluid_properties(value, n_cells: int) -> ShellTubeFluidProperties:
    if isinstance(value, ShellTubeFluidProperties):
        props = value
    elif isinstance(value, dict):
        props = ShellTubeFluidProperties(
            rho=value["rho"],
            mu=value["mu"],
            k=value["k"],
            cp=value["cp"],
        )
    else:
        raise TypeError("fluid property provider must return ShellTubeFluidProperties or dict")
    rho = _cell_property("rho", props.rho, n_cells)
    mu = _cell_property("mu", props.mu, n_cells)
    k = _cell_property("k", props.k, n_cells)
    cp = _cell_property("cp", props.cp, n_cells)
    if np.any(rho <= 0.0) or np.any(mu <= 0.0):
        raise ValueError("rho and mu must be positive")
    if np.any(k <= 0.0) or np.any(cp <= 0.0):
        raise ValueError("k and cp must be positive")
    return ShellTubeFluidProperties(rho=rho, mu=mu, k=k, cp=cp)


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
