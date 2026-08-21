"""
@ author : Raphaël Aubry

Shell-and-tube adapter for the HXDashboard themed figures.

The steady shell-and-tube solver (``main_solve_shellntube.shellntube_solver``)
does not build a ``data_master`` dict the way the helical solver does — it keeps
its results in ``solver.tube`` / ``solver.T_shell`` / ``solver.shell_liquid``.
This module maps those onto the ``data_master`` keys
``data_plotting.HXDashboard`` reads, so a shell-and-tube run is plotted by
exactly the same dashboard code as a helical run, rather than by a parallel set
of figures that would drift out of step with it.

Wired into ``result_package.package_steady_run()``, so every archived
shell-and-tube run gets a ``figures/`` folder with no manual plotting call.

What this configuration genuinely does not model, and is therefore absent (not
zero-filled) from the adapted dict:

  * radiation (``emissivity_g``/``absorptivity_g``/``X_CO2``/``X_H2O``/
    ``q_w_rad``/``h_g_rad``) — there is no participating-media radiation model
    on the shell-and-tube path at all, so ``radiation()`` skips itself and
    ``thermal()`` shows total wall flux only.
  * the coolant compressibility factor ``Z`` — not tracked on the shell side.

Shell-side properties are re-derived post-hoc from the SOLVED profile using the
same property calls the solver uses internally (real-fluid flash in liquid mode,
``IdealGasBackend`` in gas mode) — a re-evaluation of the converged state, not
an approximation of it.

Usage:
    from .data_plotting_shellntube import save_shelltube_dashboard
    save_shelltube_dashboard(solver, out_dir)
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import CoolProp.CoolProp as CP

from .data_plotting import HXDashboard, _agg_backend, _style
from ..physics.liquid_flow.coolprop_state_cache import coolprop_fluid_string, get_cached_state
from ..physics.liquid_flow.correlations import two_phase_sound_speed
from ..physics.liquid_flow.regime import (
    PSEUDO_CRITICAL_BAND_FRACTION,
    pseudo_critical_temperature,
)
from ..mechanical.loads import stress_external_pressure_tube, stress_thermal_tube

# Characterized material-property data range (°C -> K), i.e. the domain of the
# interp1d tables in mechanical/material_specs/material_temperature_strength.py.
# Outside this range the property functions clamp flat rather than extrapolate a
# real material law. Keep in sync with init_material_temperature_properties's
# supported material set.
_MATERIAL_DATA_RANGE_K = {
    "INCO718": (-240.0 + 273.15, 760.0 + 273.15),
    "ST316L": (27.0 + 273.15, 816.0 + 273.15),
}


class _ShellSideState:
    """Per-node shell-side thermophysical properties, derived post-hoc from the
    solved (T or p,h) profile using the SAME property calls the solver itself
    uses internally — not an approximation, a re-evaluation of the converged
    state.

    Gas mode (single_phase_coolprop) never marches shell pressure (see
    main_solve_shellntube.py's _shell_h_at docstring), so p is held at the
    nominal inlet value throughout — that is the solver's own behaviour, not an
    artifact introduced here.
    """

    def __init__(self, solver):
        s = solver
        N = s.N
        cool = s.coolantProp.coolant
        self.T = np.asarray(s.T_shell, dtype=float)
        rho = np.zeros(N); mu = np.zeros(N); k = np.zeros(N); cp = np.zeros(N)
        a = np.zeros(N); p = np.zeros(N)
        quality = np.full(N, np.nan)
        void = np.full(N, np.nan)

        if s._liquid_mode and s.shell_liquid is not None:
            cool_cp = coolprop_fluid_string(cool, s._liquid_backend)
            sl = s.shell_liquid
            p[:] = np.asarray(sl["p"], dtype=float)
            quality[:] = np.asarray(sl["quality"], dtype=float)
            void[:] = np.asarray(sl["void"], dtype=float)
            # Density comes straight from the solver's own converged march
            # (never re-flashed here), so it can never be blanked out by an
            # unrelated transport-property or sound-speed failure below.
            rho[:] = np.asarray(sl["rho"], dtype=float)
            h_arr = np.asarray(sl["h"], dtype=float)
            self.h = h_arr
            for i in range(N):
                # mu/k/cp: their own try/except, independent of sound speed.
                try:
                    flashed = get_cached_state(cool_cp).flash_ph(float(p[i]), float(h_arr[i]))
                    mu[i] = flashed.viscosity(); k[i] = flashed.conductivity()
                    cp[i] = flashed.cpmass()
                except (ValueError, RuntimeError):
                    mu[i] = k[i] = cp[i] = np.nan
                # Speed of sound is undefined for a single-EOS two-phase state
                # (it depends on the phase distribution, not a lone equilibrium
                # value) -- CoolProp raises there by design, not by error. Use
                # this codebase's own homogeneous-mixture closure (Wood's
                # equation, correlations.two_phase_sound_speed) inside the dome
                # instead of discarding the node; only fall back to the direct
                # EOS value where the state is genuinely single-phase.
                if 0.0 <= quality[i] <= 1.0:
                    try:
                        a[i] = two_phase_sound_speed(
                            p_Pa=float(p[i]), void_fraction=float(void[i]),
                            rho_mix_kg_m3=float(rho[i]), fluid=cool_cp)
                    except (ValueError, ZeroDivisionError):
                        a[i] = np.nan
                else:
                    try:
                        a[i] = get_cached_state(cool_cp).flash_ph(
                            float(p[i]), float(h_arr[i])).speed_sound()
                    except (ValueError, RuntimeError):
                        a[i] = np.nan
        else:
            self.h = np.full(N, np.nan)
            p[:] = s.coolantProp.p_in
            for i in range(N):
                T_i = float(self.T[i])
                rho[i] = s._thermo.density(cool, T_i, s.coolantProp.p_in)
                mu[i] = s._thermo.viscosity(cool, T_i, s.coolantProp.p_in)
                k[i] = s._thermo.conductivity(cool, T_i, s.coolantProp.p_in)
                cp[i] = s._thermo.cp(cool, T_i, s.coolantProp.p_in)
                try:
                    a[i] = CP.PropsSI("SPEED_OF_SOUND", "T", T_i, "P", s.coolantProp.p_in, cool)
                except ValueError:
                    a[i] = np.nan

        self.rho, self.mu, self.k, self.cp, self.a, self.p = rho, mu, k, cp, a, p
        self.quality, self.void = quality, void

        G_s = s.coolantProp.mass_flow_c / s.geom["S_m"]
        self.G_s = G_s
        with np.errstate(invalid="ignore", divide="ignore"):
            self.Re = s.stp.D_tube_outer * G_s / mu
            self.Pr = cp * mu / k
            self.U = G_s / rho
            # Whole-bundle volumetric flow (contrast with the tube side, which
            # this codebase books per representative tube).
            self.Vdot_Ls = s.coolantProp.mass_flow_c / rho * 1000.0
            # Not a rigorous 1D compressible Mach -- the shell side has no
            # single streamwise coordinate. This is a cross-flow velocity-scale
            # / local sound-speed ratio, an order-of-magnitude compressibility
            # gauge only, never a momentum-equation validity check the way the
            # tube-side Mach is.
            self.Mach_like = self.U / a


def shelltube_data_master(solver):
    """Adapt a SOLVED ``shellntube_solver`` onto the ``data_master`` keys read by
    ``data_plotting.HXDashboard``. Returns ``(data_master, shell_state)``.

    Requires ``solver.solve()`` to have run. Quantities this configuration does
    not model are left out of the dict entirely rather than zero-filled, so the
    dashboard's own guards blank those panels instead of drawing a physically
    meaningless flat line.
    """
    s = solver
    d = s.tube
    N = s.N
    shell = _ShellSideState(s)

    Di, Do = s.D_tube_i, s.stp.D_tube_outer
    sw, dx = s.stp.thickness_tube_wall, s.dx

    # --- resistance network (per node) --------------------------------------
    A_hot = np.pi * Di * dx          # hot side is INSIDE the tubes here
    A_cold = np.pi * Do * dx
    T_wg = np.asarray(d["T_wg"], dtype=float)
    T_wc = np.asarray(d["T_wc"], dtype=float)
    T_avg_C = (T_wg + T_wc) / 2.0 - 273.15
    k_w = np.array([s.k_t(t) for t in T_avg_C])
    R_w = np.log((Di / 2 + sw) / (Di / 2)) / (2 * np.pi * dx * k_w)
    h_g = np.asarray(d["h_g"], dtype=float)
    h_c = np.asarray(d["h_c"], dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        R_g = 1.0 / (h_g * A_hot)
        R_c = 1.0 / (h_c * A_cold)
        UA = 1.0 / (R_g + R_w + R_c)
        Biot_g = R_w / R_g
        Biot_c = R_w / R_c
        Nu_g = h_g * Di / np.asarray(d["k_g"], dtype=float)
        Nu_c = h_c * Do / shell.k

    # --- mechanical (per node) ----------------------------------------------
    # Hot gas is INSIDE the tubes (hot_side="inner"), so T_wg is the inner face
    # here -- the opposite argument order from the helical coil, where the
    # coolant is the inner fluid.
    CTE = np.array([s.CTE_t(t) for t in T_avg_C])
    Modulus = np.array([s.E_t(t) for t in T_avg_C])
    Yield = np.array([s.Yield_t(t) for t in (T_wg - 273.15)])   # hot face governs
    st_in = np.zeros(N); st_out = np.zeros(N)
    for i in range(N):
        st_in[i], st_out[i] = stress_thermal_tube(
            T_inner=T_wg[i] - 273.15, T_outer=T_wc[i] - 273.15,
            CTE=CTE[i], E=Modulus[i], poisson=s.poisson_t)
    # Node-local external-pressure stress from the ACTUAL local shell/tube
    # pressures, rather than the single global (p_in, p0) pair that
    # solver.compute_stress() uses -- it refines that scalar into a profile.
    dP_ext = shell.p - np.asarray(d["p_g"], dtype=float)
    stress_pressure = np.array([
        stress_external_pressure_tube(dP_ext[i], sw, Di) for i in range(N)])

    dm = {
        "L_HX": np.arange(N) * dx,
        # temperatures / heat transfer
        "T_g": np.asarray(d["T_g"], dtype=float),
        "T_wg": T_wg,
        "T_wc": T_wc,
        "T_c": shell.T,
        "q_w": np.asarray(d["q_w_shell"], dtype=float),
        "h_g_conv": h_g,
        "h_c": h_c,
        "Res_g": R_g, "Res_w": R_w, "Res_c": R_c,
        "UA": UA, "Biot_g": Biot_g, "Biot_c": Biot_c,
        "k_w": k_w,
        # coolant (shell side, whole bundle)
        "p_c": shell.p, "rho_c": shell.rho, "U_c": shell.U,
        "Mach_c": shell.Mach_like, "c_c": shell.a,
        "cp_c": shell.cp, "mu_c": shell.mu, "k_c": shell.k,
        "Re_c": shell.Re, "Pr_c": shell.Pr, "Nu_c": Nu_c,
        # hot gas (tube side, PER REPRESENTATIVE TUBE)
        "p_g": np.asarray(d["p_g"], dtype=float),
        "rho_g": np.asarray(d["rho_g"], dtype=float),
        "U_g": np.asarray(d["U_g"], dtype=float),
        "Mach_g": np.asarray(d["mach_g"], dtype=float),
        "cp_g": np.asarray(d["cp_g"], dtype=float),
        "mu_g": np.asarray(d["mu_g"], dtype=float),
        "k_g": np.asarray(d["k_g"], dtype=float),
        "Re_g": np.asarray(d["Re_g"], dtype=float),
        "Pr_g": np.asarray(d["Pr_g"], dtype=float),
        "Nu_g": Nu_g,
        # mechanical
        "CTE": CTE, "Modulus": Modulus, "Yield": Yield,
        "stress_pressure": stress_pressure,
        "stress_thermal_inner": st_in, "stress_thermal_outer": st_out,
        "stress_inner": st_in + stress_pressure,
        "stress_outer": st_out + stress_pressure,
    }

    # --- two-phase fields: only when the march actually enters the dome -----
    # Above the critical pressure the state closure still returns an *extended*
    # quality, (h - h_f)/h_fg carried on past x = 1 (the solver's own summary
    # prints e.g. "quality: 4.8 -> 20.1" for a supercritical N2 run). Feeding
    # that to boiling()/phase_change() would draw a dryout line at x = 1 under
    # a curve that never had a dome to leave. So the two-phase themes are
    # offered only when some node is genuinely inside 0 <= x <= 1; otherwise
    # they skip themselves and extras()' pseudo-critical panel carries the
    # diagnosis instead.
    in_dome = np.isfinite(shell.quality) & (shell.quality >= 0.0) & (shell.quality <= 1.0)
    if s._liquid_mode and s.shell_liquid is not None and in_dome.any():
        dp_fric = np.asarray(d["dp_shell"], dtype=float)
        dp_acc = np.asarray(s.shell_liquid["dp_accel"], dtype=float)
        dm.update({
            "enthalpy_c": shell.h,
            "quality_c": shell.quality,
            "void_c": shell.void,
            "chf_margin_c": np.asarray(d["chf_margin"], dtype=float),
            # The solver stores positive per-node DROPS; the dashboard plots a
            # gradient, so both become negative per-metre values here.
            "dp_c__dx": -(dp_fric + dp_acc) / dx,
            "dp_c__dx_accel": -dp_acc / dx,
        })

    return dm, shell


class ShellTubeDashboard(HXDashboard):
    """``HXDashboard`` fed by a shell-and-tube solver, plus one extra figure for
    the shell-and-tube-specific diagnostics the generic themes have no panel for
    (volumetric flow, thermal-stress margin against the temperature-dependent
    yield, wall temperature against the material's characterized range, and
    supercritical proximity to the pseudo-critical line).
    """

    def __init__(self, solver, save_dir=None, dpi=200):
        dm, shell = shelltube_data_master(solver)
        cool = solver.coolantProp.coolant
        name = (coolprop_fluid_string(cool, solver._liquid_backend)
                if solver._liquid_mode else cool)
        super().__init__(dm, coolant_name=name, save_dir=save_dir, dpi=dpi)
        self.solver = solver
        self.shell = shell

    # ------------------------------------------------------------------
    def extras(self):
        """Shell-and-tube specifics with no panel in the generic themes."""
        s, d, x, shell = self.solver, self.d, self.x, self.shell
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
        fig.suptitle("Shell-and-tube specifics", fontweight='bold')

        # (0,0) volumetric flow rates -- note the two different conventions
        ax = axes[0, 0]
        Vdot_g = s.mdot_tube / d["rho_g"] * 1000.0
        ax.plot(x, shell.Vdot_Ls, color="cornflowerblue",
                label=r"$\dot{V}_c$ (shell, whole bundle)")
        ax.plot(x, Vdot_g, color="red", label=r"$\dot{V}_g$ (per representative tube)")
        ax.set_yscale('log')
        ax.set_ylabel("Volumetric flow [L/s]")
        _style(ax, x=x)

        # (0,1) thermal-stress safety margin against the local yield
        ax = axes[0, 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            margin = d["Yield"] / np.abs(d["stress_thermal_inner"])
        ax.semilogy(x, margin, color="darkorange", label="Yield / thermal stress")
        ax.axhline(1.0, color='red', linestyle='--', linewidth=1, label="margin = 1")
        ax.set_ylabel("Safety margin [-]")
        ax.set_title("Thermal-stress margin"
                     + (f"  |  collapse margin = {s.collapse_margin:.4f}"
                        if hasattr(s, "collapse_margin") else ""), fontsize=9)
        _style(ax, x=x)

        # (1,0) wall temperature against the material's characterized range
        ax = axes[1, 0]
        ax.plot(x, d["T_wg"] - 273.15, color="darkorange", label=r"$T_{wg}$")
        ax.plot(x, d["T_wc"] - 273.15, color="mediumblue", label=r"$T_{wc}$")
        rng = _MATERIAL_DATA_RANGE_K.get(s.stp.material_tube)
        if rng is not None:
            ceiling_C = rng[1] - 273.15
            ax.axhline(ceiling_C, color='crimson', linestyle=':', linewidth=1.2,
                       label=f"{s.stp.material_tube} data ceiling ({ceiling_C:.0f} °C)")
            beyond = d["T_wg"] > rng[1]
            if beyond.any():
                ax.axvspan(x[beyond].min(), x[beyond].max(), color="crimson",
                           alpha=0.10, zorder=0)
        ax.set_ylabel("Wall temperature [°C]")
        ax.set_title("Beyond the ceiling the property tables clamp flat, "
                     "they do not extrapolate", fontsize=8)
        _style(ax, x=x)

        # (1,1) supercritical proximity to the pseudo-critical line, if it applies
        ax = axes[1, 1]
        T_pc = np.full(len(x), np.nan)
        if s._liquid_mode:
            cool_cp = coolprop_fluid_string(s.coolantProp.coolant, s._liquid_backend)
            for i in range(len(x)):
                try:
                    T_pc[i] = pseudo_critical_temperature(cool_cp, float(shell.p[i]))
                except (ValueError, RuntimeError):
                    T_pc[i] = np.nan
        fin = np.isfinite(T_pc)
        if fin.any():
            band = PSEUDO_CRITICAL_BAND_FRACTION * T_pc
            ax.plot(x, d["T_c"] - 273.15, color="cornflowerblue", label=r"$T_c$")
            ax.plot(x[fin], T_pc[fin] - 273.15, color="purple", linestyle='--',
                    label=r"$T_{pc}(p)$")
            ax.fill_between(x[fin], (T_pc - band)[fin] - 273.15,
                            (T_pc + band)[fin] - 273.15,
                            color="purple", alpha=0.15, label="pseudo-critical band")
            ax.set_ylabel("Temperature [°C]")
            ax.set_title("Supercritical proximity  |  HTD-risk nodes: "
                         f"{getattr(s, '_sc_htd_nodes', 0)}", fontsize=9)
            _style(ax, x=x)
        else:
            ax.text(0.5, 0.5, "no supercritical coolant states on this run",
                    ha='center', va='center', transform=ax.transAxes,
                    fontsize=9, color='gray')
            ax.set_xticks([]); ax.set_yticks([])

        self._finish(fig, "shelltube_extras.png")

    # ------------------------------------------------------------------
    def all(self):
        super().all()
        self.extras()
        return list(self.written)


def save_shelltube_dashboard(solver, out_dir) -> list[str]:
    """Render the full dashboard for a SOLVED ``shellntube_solver`` into
    ``out_dir`` (created if missing). Returns the list of written paths.

    Forces the Agg backend for the duration and restores the previous one, so
    calling this mid-session never leaves the user's interactive plots headless.
    """
    out_dir = str(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    with _agg_backend():
        return ShellTubeDashboard(solver, save_dir=out_dir).all()
