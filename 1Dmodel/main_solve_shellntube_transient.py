"""
Dynamic (transient) solver — baffled shell-and-tube config (WP2/3).

@ author : Raphaël Aubry

Same architecture as main_solve_transient.py (lumped wall-energy ODE per axial
node, quasi-steady fluid `fluid_pass` per RHS evaluation, solve_ivp/BDF
integration) — see DESIGN_PLAN_shellntube_transient.md section 4 and 4.5
("both configs ship transient — hard requirement").

Counter-flow is handled by a warm-started coolant-profile relaxation inside each
fluid pass. The transient wall integrator uses a bounded-cost linearly-implicit
fixed step by default so long counter-flow runs do not pay BDF's Jacobian-probing
cost.

No tube-side radiation (optically thin at this Ø, per the steady solver's design
note) — the Nu_g/T_wall fixed point from the helical transient is retained only
in case a user sets a non-zero n_tube_gas (Kays-Crawford correction knob).
"""
# Current implementation supports both co-flow and counter-flow. Counter-flow
# uses the same profile-relaxation idiom as the helical transient solver.
if __name__ == "__main__" and __package__ is None:
    import sys as _sys, os as _os, importlib.util as _ilu, runpy as _rp
    _pkg_dir = _os.path.dirname(_os.path.abspath(__file__))
    _parent = _os.path.dirname(_pkg_dir)
    _alias = "_hps"
    if _alias not in _sys.modules:
        _spec = _ilu.spec_from_file_location(
            _alias, _os.path.join(_pkg_dir, "__init__.py"),
            submodule_search_locations=[_pkg_dir])
        _pkg = _ilu.module_from_spec(_spec)
        _pkg.__path__ = [_pkg_dir]
        _pkg.__package__ = _alias
        _sys.modules[_alias] = _pkg
        _spec.loader.exec_module(_pkg)
    if _parent not in _sys.path:
        _sys.path.insert(0, _parent)
    _rp.run_module(f"{_alias}.main_solve_shellntube_transient", run_name="__main__", alter_sys=True)
    raise SystemExit(0)

import numpy as np
from types import SimpleNamespace
from scipy.integrate import solve_ivp

from .main_solve_shellntube import shellntube_solver
from .physics.heat_conduction import OneDimensionalSteadyConduction_ShellnHelicalTube
from .physics.bell_delaware import bell_delaware_shell
from .physics.combustion_chemistry.gas_manifold import build_equilibrium_manifold, EquilibriumManifold
from .physics.combustion_chemistry.fpv_manifold import build_fpv_manifold, FPVManifold
from .transient_core import (
    build_shelltube_core_geometry_from_solver,
    coolprop_fluid_properties,
    equilibrium_gas_state_provider,
    fpv_gas_state_provider,
    oxygen_gas_state_provider,
    run_shelltube_transient_core,
)


def _interp_schedule(schedule, t, default):
    if schedule is None:
        return default
    pts = list(schedule)
    if t <= pts[0][0]:
        return pts[0][1]
    if t >= pts[-1][0]:
        return pts[-1][1]
    for (t0, v0), (t1, v1) in zip(pts[:-1], pts[1:]):
        if t0 <= t <= t1:
            w = (t - t0) / (t1 - t0)
            return v0 + w * (v1 - v0)
    return pts[-1][1]


class shellntube_transient_solver(shellntube_solver):
    def __init__(self, coolantProp, hotgasProp, shellTubeProp, numericalProp,
                system_requirements, transientProp, corrCoeffs=None, N_axial=80,
                flow_config="co"):
        super().__init__(coolantProp, hotgasProp, shellTubeProp, numericalProp,
                         system_requirements, corrCoeffs=corrCoeffs, N_axial=N_axial,
                         flow_config=flow_config)
        if flow_config not in ("co", "counter"):
            raise ValueError(f"flow_config must be 'co' or 'counter', got {flow_config!r}")
        self.transientProp = transientProp
        self._counter_Tc_profile = None
        self._counter_fallback_mdot = 0.02 * max(self.coolantProp.mass_flow_c, 1e-9)

        # wall thermal mass (tube annulus)
        Do, Di = self.stp.D_tube_outer, self.D_tube_i
        self.A_wall = np.pi * (Do ** 2 - Di ** 2) / 4.0
        self.rho_w = self.rho_t

        self._build_chem_manifold()

    # ------------------------------------------------------------------
    def _build_chem_manifold(self):
        mode = self.transientProp.chemistry_transient
        self._chem_mode = mode
        Tg0, pg0, Yg0 = self._gas_inlet_TPY
        if mode == "finite_rate":
            iCO = self.gas_phase.species_index("CO")
            self._fpv = FPVManifold(build_fpv_manifold(
                self.gas_phase, Y_inlet=Yg0, T_inlet=Tg0, p=pg0,
                species_index={"CO2": self.gas_phase.species_index("CO2"),
                              "H2O": self.gas_phase.species_index("H2O"), "CO": iCO},
                n_h=self.numericalProp.fpv_n_h,
                n_c=self.numericalProp.fpv_n_c,
                t_relax=self.numericalProp.fpv_t_relax,
                n_t=self.numericalProp.fpv_n_t,
                cache_dir=self.numericalProp.fpv_cache_dir))
        else:
            self._eqm = EquilibriumManifold(build_equilibrium_manifold(
                self.gas_phase, T_inlet=Tg0, p_inlet=pg0, Y_inlet=Yg0, mode=mode, n_h=200))
        self.gas_phase.TPY = Tg0, pg0, Yg0  # restore inlet

    # ------------------------------------------------------------------
    def _bc_at(self, t):
        tp = self.transientProp
        return dict(
            mdot_c=_interp_schedule(tp.schedule_mass_flow_c, t, self.coolantProp.mass_flow_c),
            mdot_g=_interp_schedule(tp.schedule_mass_flow_g, t, self.hotgasProp.mass_flow_g),
            mdot_lox=_interp_schedule(getattr(tp, "schedule_mass_flow_lox", None), t, 0.0),
            mdot_diesel=_interp_schedule(getattr(tp, "schedule_mass_flow_diesel", None), t, 0.0),
            T_c_in=_interp_schedule(tp.schedule_T_c_in, t, self.coolantProp.T_in),
            p_c_in=_interp_schedule(tp.schedule_p_c_in, t, self.coolantProp.p_in),
            T_lox_in=_interp_schedule(tp.schedule_T_lox_in, t, self.hotgasProp.T_inj_LOX),
            ignited=(
                _interp_schedule(getattr(tp, "schedule_ignition_state", None), t, 1.0 if t >= tp.ignition_time else 0.0) >= 0.5
            ),
        )

    # ------------------------------------------------------------------
    def fluid_pass(self, Tbar_vec, bc):
        """Dispatch quasi-steady fluid pass for co-flow or counter-flow."""
        if self.flow_config == "co":
            return self._march_fluids(Tbar_vec, bc)
        return self._relax_counter_flow(Tbar_vec, bc)

    def _relax_counter_flow(self, Tbar_vec, bc, max_iter=None, tol=None, omega0=0.5):
        """Counter-flow shell-side profile relaxation at fixed wall field."""
        if bc["mdot_c"] < self._counter_fallback_mdot:
            return self._march_fluids(Tbar_vec, bc)

        N = self.N
        cool = self.coolantProp.coolant
        warm_started = self._counter_Tc_profile is not None and len(self._counter_Tc_profile) == N
        if warm_started:
            T_c_profile = self._counter_Tc_profile.copy()
        else:
            T_c_profile = np.linspace(0.5 * (bc["T_c_in"] + float(Tbar_vec[0])),
                                      bc["T_c_in"], N)
        if max_iter is None:
            max_iter = (
                self.transientProp.counterflow_warm_relax_iter
                if warm_started else
                self.transientProp.counterflow_initial_relax_iter
            )
        if tol is None:
            tol = self.transientProp.counterflow_relax_tol_K

        res = None
        omega = omega0
        delta_prev = None
        for _it in range(max_iter):
            res = self._march_fluids(Tbar_vec, bc, T_c_profile=T_c_profile)
            dq_cold = res["dq_cold__dx"]
            T_c_new = np.empty(N)
            T_c_new[N - 1] = bc["T_c_in"]
            for j in range(N - 1, 0, -1):
                cp_c = self._thermo.cp(cool, T_c_new[j], bc["p_c_in"])
                T_c_new[j - 1] = T_c_new[j] + (dq_cold[j] * self.stp.N_tubes) * self.dx / \
                    (bc["mdot_c"] * cp_c)
                T_c_new[j - 1] = float(np.clip(T_c_new[j - 1], 100.0, 1500.0))

            delta = float(np.max(np.abs(T_c_new - T_c_profile)))
            if delta_prev is not None and delta > delta_prev:
                omega = max(omega * 0.5, 0.02)
            T_c_profile = (1.0 - omega) * T_c_profile + omega * T_c_new
            delta_prev = delta
            if delta < tol:
                break

        self._counter_Tc_profile = T_c_profile.copy()
        res["T_c"] = T_c_profile
        res["T_c_out"] = float(T_c_profile[0])
        return res

    def _march_fluids(self, Tbar_vec, bc, T_c_profile=None):
        """Forward gas march with either co-flow coolant marching or a prescribed
        counter-flow coolant profile."""
        N, dx = self.N, self.dx
        cool = self.coolantProp.coolant
        fr = (self._chem_mode == "finite_rate")

        prescribed_Tc = T_c_profile is not None
        T_c = float(T_c_profile[0]) if prescribed_Tc else bc["T_c_in"]
        p_c = bc["p_c_in"]
        mdot_c = bc["mdot_c"]
        combustion_on = bool(bc["ignited"])
        mdot_hot_total = bc["mdot_g"] if combustion_on else bc.get("mdot_lox", 0.0)
        mdot_tube = mdot_hot_total / self.stp.N_tubes   # per-tube flow; totals are schedule-level
        h_removed = 0.0
        Yc = self._fpv.Yc_inlet() if fr else 0.0
        T_g_gox = float(np.clip(bc.get("T_lox_in", self.hotgasProp.T_inj_LOX), 95.0, 1200.0))

        dq_hot = np.zeros(N); dq_cold = np.zeros(N)
        T_wg_a = np.zeros(N); T_wc_a = np.zeros(N)
        T_g_a = np.zeros(N); T_c_a = np.zeros(N)
        h_g_a = np.zeros(N); h_c_a = np.zeros(N)

        for i in range(N):
            if prescribed_Tc:
                T_c = float(T_c_profile[i])
            if combustion_on and fr:
                T_g, rho_g, mu_g, k_g, cp_g, X_H2O, X_CO2, omega_Yc = self._fpv.state(h_removed, Yc)
            elif combustion_on:
                T_g, rho_g, mu_g, k_g, cp_g, X_H2O, X_CO2 = self._eqm.at(h_removed)
            else:
                T_g = T_g_gox
                p_gox = max(self.hotgasProp.p0, 1e4)
                rho_g = self._thermo.density("Oxygen", T_g, p_gox)
                mu_g = self._thermo.viscosity("Oxygen", T_g, p_gox)
                k_g = self._thermo.conductivity("Oxygen", T_g, p_gox)
                cp_g = self._thermo.cp("Oxygen", T_g, p_gox)
                omega_Yc = 0.0
            hot_flowing = bool(mdot_tube > 0)
            U_g = mdot_tube / (rho_g * self.A_tube_i) if hot_flowing else 0.0
            Re_g = rho_g * U_g * self.D_tube_i / mu_g if hot_flowing else 1.0
            Pr_g = cp_g * mu_g / k_g

            rho_c = self._thermo.density(cool, T_c, p_c)
            mu_c = self._thermo.viscosity(cool, T_c, p_c)
            k_c = self._thermo.conductivity(cool, T_c, p_c)
            cp_c = self._thermo.cp(cool, T_c, p_c)
            Pr_c = cp_c * mu_c / k_c
            G_s = mdot_c / self.geom["S_m"]
            Re_c = self.stp.D_tube_outer * G_s / mu_c
            geom = dict(self.geom); geom["rho_s"] = rho_c
            bell = bell_delaware_shell(geom, Re_s=Re_c, Pr_s=Pr_c, k_s=k_c, cp_s=cp_c,
                                       mu_s=mu_c, mdot_s=mdot_c, corrCoeffs=self.corrCoeffs)
            h_c = bell["h_shell"]

            x_local = (i + 0.5) * dx
            T_wg_est = float(Tbar_vec[i])
            for _fp in range(4):
                if hot_flowing:
                    f_g, Nu_g = self._tube_side_hydraulics(
                        Re_g, Pr_g, x_local, T_g, T_wg_est)
                    h_g = Nu_g * k_g / self.D_tube_i
                else:
                    h_g = 1.0
                node = OneDimensionalSteadyConduction_ShellnHelicalTube(
                    h_g=h_g, h_c=h_c, T_c=T_c, T_g=T_g,
                    s_w=self.stp.thickness_tube_wall, Dh_ch=self.D_tube_i,
                    f_kw_at_T=self.k_t, T_wg_0=Tbar_vec[i], T_wc_0=Tbar_vec[i],
                    T_c_check_0=T_c, dx=dx, rad_enabled=False, hot_side="inner")
                r = node.fluxes_at_Tbar(Tbar_vec[i], h_g_rad=0.0)
                if abs(r["T_wg"] - T_wg_est) < 1e-4:
                    break
                T_wg_est = r["T_wg"]

            dq_hot[i] = r["dq_hot__dx"]; dq_cold[i] = r["dq_cold__dx"]
            T_wg_a[i] = r["T_wg"]; T_wc_a[i] = r["T_wc"]; T_g_a[i] = T_g; T_c_a[i] = T_c
            h_g_a[i] = h_g; h_c_a[i] = h_c

            # advance fluids (per-tube dq_hot; shell side sees N_tubes tubes -> total flux)
            if combustion_on and hot_flowing:
                h_removed += dq_hot[i] * dx / mdot_tube
                if fr and U_g > 0:
                    Yc += omega_Yc * dx / U_g
            elif hot_flowing:
                T_g_gox -= dq_hot[i] * dx / (mdot_tube * cp_g)
                T_g_gox = float(np.clip(T_g_gox, 95.0, 1200.0))
            if mdot_c > 0 and not prescribed_Tc:
                T_c += (dq_cold[i] * self.stp.N_tubes) * dx / (mdot_c * cp_c)
                # Guard against solve_ivp's internal Newton/Jacobian probing calling
                # fluid_pass with a wall field far from physical (e.g. near the cold-start
                # T_wall_initial): an explicit per-node march is not unconditionally bounded,
                # and PropsSI hard-crashes on unphysical T. Clamp to a generous physical
                # range (CoolProp Helium is safely evaluable throughout) instead of letting
                # a transient excursion propagate into a solver crash; the integrator's own
                # step-size control will reject/shrink steps that hit this bound.
                T_c = float(np.clip(T_c, 100.0, 1500.0))
            p_c -= 0.0  # shell dp not tracked per-node in the transient march (see steady solver)

        return dict(dq_hot__dx=dq_hot, dq_cold__dx=dq_cold, T_wg=T_wg_a, T_wc=T_wc_a,
                    T_g=T_g_a, T_c=T_c_a, h_g=h_g_a, h_c=h_c_a, T_c_out=T_c,
                    T_g_out=(
                        self._fpv.state(h_removed, Yc)[0]
                        if combustion_on and fr else
                        self._eqm.at(h_removed)[0]
                        if combustion_on else
                        T_g_gox
                    ))

    # ------------------------------------------------------------------
    def _wall_rhs(self, t, Tbar_vec):
        return self._wall_rate(t, Tbar_vec)[0]

    def _wall_rate(self, t, Tbar_vec):
        bc = self._bc_at(t)
        res = self.fluid_pass(Tbar_vec, bc)
        cp_w = np.array([self.cp_t(T - 273.15) for T in Tbar_vec])
        denom = self.rho_w * cp_w * self.A_wall
        dTbar_dt = (res["dq_hot__dx"] - res["dq_cold__dx"]) / denom
        P_h = np.pi * self.D_tube_i
        P_c = np.pi * self.stp.D_tube_outer
        lam = (np.asarray(res["h_g"]) * P_h + np.asarray(res["h_c"]) * P_c) / denom
        return dTbar_dt, lam

    # ------------------------------------------------------------------
    def solve_transient(self, verbose=True):
        tp = self.transientProp
        Tbar0 = np.full(self.N, tp.T_wall_initial)
        t_eval = np.linspace(0.0, tp.t_end, tp.n_save)
        if verbose:
            print(f"Integrating shell-and-tube wall ODE: {self.N} nodes, t_end={tp.t_end}s, "
                  f"method={tp.solver_method}, max_step={tp.max_step}s ...")
        if tp.solver_method.lower() in ("fixed_step", "euler", "explicit_euler"):
            sol = self._solve_fixed_step(Tbar0, t_eval)
            self.sol = sol
            if verbose:
                self._print_final_state(sol)
            return sol
        ivp_options = {}
        if tp.use_sparse_jacobian and tp.solver_method in ("BDF", "Radau"):
            ivp_options["jac_sparsity"] = np.eye(self.N, dtype=bool)
        sol = solve_ivp(self._wall_rhs, (0.0, tp.t_end), Tbar0, method=tp.solver_method,
                        t_eval=t_eval, max_step=tp.max_step, rtol=1e-5, atol=1e-2,
                        **ivp_options)
        if not sol.success:
            print(f"  WARNING: solve_ivp did not fully succeed: {sol.message}")
        self.sol = sol
        if verbose:
            self._print_final_state(sol)
        return sol

    def solve_transient_core(self, verbose=True):
        """Run the finite-volume wall + helium transient core.

        Hot gas remains quasi-steady over each fixed step. The hot-side property
        provider can switch between GOX chilldown and combustion at schedule
        breakpoints, while helium temperature is an integrated state.
        """
        tp = self.transientProp
        geometry = build_shelltube_core_geometry_from_solver(self)
        bc0 = self._bc_at(0.0)
        T_coolant_initial = (
            bc0["T_c_in"]
            if getattr(tp, "T_coolant_initial", None) is None else
            float(tp.T_coolant_initial)
        )
        T_wall0 = np.full(self.N, float(tp.T_wall_initial))
        T_coolant0 = np.full(self.N, float(T_coolant_initial))

        if self._chem_mode == "finite_rate":
            combustion_provider, combustion_progress0 = fpv_gas_state_provider(self._fpv)
        else:
            combustion_provider, combustion_progress0 = equilibrium_gas_state_provider(self._eqm)

        def gas_provider_at_time(t):
            bc = self._bc_at(float(t))
            if bc["ignited"]:
                return combustion_provider, combustion_progress0
            return oxygen_gas_state_provider(
                T_inlet=bc["T_lox_in"],
                pressure=max(float(self.hotgasProp.p0), 1.0e4),
            )

        def mdot_hot_total_at_time(t):
            bc = self._bc_at(float(t))
            return bc["mdot_g"] if bc["ignited"] else bc.get("mdot_lox", 0.0)

        if verbose:
            print(
                "Integrating shell-and-tube transient core: "
                f"{self.N} nodes, t_end={tp.t_end}s, max_step={tp.max_step}s, "
                f"chemistry={self._chem_mode}, flow={self.flow_config}"
            )

        result = run_shelltube_transient_core(
            geometry,
            self.geom,
            T_wall_initial=T_wall0,
            T_coolant_initial=T_coolant0,
            t_end=float(tp.t_end),
            max_step=float(tp.max_step),
            coolant_properties_at=coolprop_fluid_properties(self.coolantProp.coolant),
            wall_density=float(self.rho_t),
            wall_cp=lambda T: np.array([self.cp_t(float(Ti) - 273.15) for Ti in np.asarray(T)]),
            wall_conductivity_at_T=self.k_t,
            inside_tube_choice=self.stp.inside_tube_choice,
            nusselt_selector=self.stp.Nusselt_tube,
            tube_roughness=self.stp.tube_roughness,
            mdot_coolant_default=float(self.coolantProp.mass_flow_c),
            T_coolant_inlet_default=float(self.coolantProp.T_in),
            p_coolant_default=float(self.coolantProp.p_in),
            mdot_hot_total_default=float(self.hotgasProp.mass_flow_g),
            mdot_coolant_schedule=tp.schedule_mass_flow_c,
            T_coolant_inlet_schedule=tp.schedule_T_c_in,
            p_coolant_schedule=tp.schedule_p_c_in,
            mdot_hot_total_schedule=tp.schedule_mass_flow_g,
            gas_provider_at_time=gas_provider_at_time,
            mdot_hot_total_at_time=mdot_hot_total_at_time,
            hot_side_schedules=(
                tp.schedule_mass_flow_g,
                getattr(tp, "schedule_mass_flow_lox", None),
                getattr(tp, "schedule_T_lox_in", None),
                getattr(tp, "schedule_ignition_state", None),
                getattr(tp, "schedule_p_c_out", None),
            ),
            n_save=int(tp.n_save),
            corrCoeffs=self.corrCoeffs,
            corrugation_thickness=self.stp.corrugation_thickness,
            corrugation_pitch=self.stp.corrugation_pitch,
            flow_direction=geometry.grid.flow_direction,
            mdot_floor=getattr(tp, "transient_coolant_mdot_floor", 1e-9),
            coolant_state_model=(
                "low_mach_momentum"
                if getattr(tp, "coolant_momentum_model", "quasi_steady") == "low_mach" else
                "mass_energy"
            ),
            progress_config=tp,
            progress_enabled=verbose,
        )

        self.core_result = result
        self.sol = SimpleNamespace(
            t=result.integration.t,
            y=result.integration.T_wall.T,
            success=True,
            message="transient_core completed",
            nfev=len(result.step_diagnostics),
        )
        self.time_series = self._time_series_from_core_result(result, geometry)
        if verbose:
            self._print_core_final_state(result)
        return result

    def _time_series_from_core_result(self, result, geometry):
        integ = result.integration
        t = np.asarray(integ.t, dtype=float)
        n_time, n_nodes = integ.T_wall.shape
        fields = {
            "Tbar": integ.T_wall.copy(),
            "T_c": integ.T_coolant.copy(),
            "heat_wall_to_coolant_W": integ.heat_wall_to_coolant_W.copy(),
            "T_g": np.full((n_time, n_nodes), np.nan),
            "T_wg": np.full((n_time, n_nodes), np.nan),
            "T_wc": np.full((n_time, n_nodes), np.nan),
            "h_g": np.full((n_time, n_nodes), np.nan),
            "h_c": np.full((n_time, n_nodes), np.nan),
            "U_g": np.full((n_time, n_nodes), np.nan),
            "Re_g": np.full((n_time, n_nodes), np.nan),
            "Pr_g": np.full((n_time, n_nodes), np.nan),
            "Nu_g": np.full((n_time, n_nodes), np.nan),
            "f_g": np.full((n_time, n_nodes), np.nan),
            "dp_g_per_m": np.full((n_time, n_nodes), np.nan),
            "h_removed_g": np.full((n_time, n_nodes), np.nan),
            "progress_g": np.full((n_time, n_nodes), np.nan),
            "G_shell": np.full((n_time, n_nodes), np.nan),
            "Re_shell": np.full((n_time, n_nodes), np.nan),
            "Pr_shell": np.full((n_time, n_nodes), np.nan),
            "dp_shell": np.full((n_time, n_nodes), np.nan),
            "rho_c": np.full((n_time, n_nodes), np.nan),
            "mu_c": np.full((n_time, n_nodes), np.nan),
            "k_c": np.full((n_time, n_nodes), np.nan),
            "cp_c": np.full((n_time, n_nodes), np.nan),
        }
        optional_field_map = {
            "coolant_mass_kg": "coolant_mass_kg",
            "coolant_internal_energy_J": "coolant_internal_energy_J",
            "p_c": "coolant_pressure_Pa",
            "rho_c_state": "coolant_density_kg_m3",
            "h_c_state": "coolant_specific_enthalpy_J_kg",
            "face_mdot_c": "face_mdot_kg_s",
        }
        for output_name, attr_name in optional_field_map.items():
            if hasattr(integ, attr_name):
                fields[output_name] = np.asarray(getattr(integ, attr_name), dtype=float).copy()
        scalars = {
            "Q_hot_kW": np.zeros(n_time),
            "Q_cold_kW": np.sum(integ.heat_wall_to_coolant_W, axis=1) / 1.0e3,
            "T_c_out": integ.T_coolant_outlet.copy(),
            "T_g_out": np.full(n_time, np.nan),
            "T_wall_max": np.max(integ.T_wall, axis=1),
            "T_wall_min": np.min(integ.T_wall, axis=1),
            "T_c_min": np.min(integ.T_coolant, axis=1),
            "T_c_max": np.max(integ.T_coolant, axis=1),
            "mdot_c": np.zeros(n_time),
            "mdot_c_effective": np.full(n_time, np.nan),
            "mdot_c_inlet_face": np.full(n_time, np.nan),
            "mdot_c_outlet_face": np.full(n_time, np.nan),
            "mdot_g": np.zeros(n_time),
            "energy_residual_J": integ.energy_residual_J.copy(),
            "h_removed_g_out_J_kg": np.full(n_time, np.nan),
            "progress_g_out": np.full(n_time, np.nan),
            "dp_g_total_Pa": np.full(n_time, np.nan),
            "dp_shell_total_Pa": np.full(n_time, np.nan),
            "Re_g_max": np.full(n_time, np.nan),
            "Re_shell_max": np.full(n_time, np.nan),
        }
        if hasattr(integ, "mass_residual_kg"):
            scalars["coolant_mass_residual_kg"] = np.asarray(
                integ.mass_residual_kg,
                dtype=float,
            ).copy()
        diagnostics = result.step_diagnostics
        for i in range(n_time):
            if not diagnostics:
                bc = self._bc_at(float(t[i]))
                scalars["mdot_c"][i] = bc["mdot_c"]
                scalars["mdot_g"][i] = bc["mdot_g"] if bc["ignited"] else bc.get("mdot_lox", 0.0)
                continue
            diag = diagnostics[min(i, len(diagnostics) - 1)]
            march = diag.hot_gas_march
            fields["T_g"][i] = march.T_gas
            fields["T_wg"][i] = march.wall_flux.T_wg
            fields["T_wc"][i] = march.wall_flux.T_wc
            fields["h_g"][i] = march.h_gas_W_m2K
            fields["h_c"][i] = diag.shell_film.h_W_m2K
            fields["U_g"][i] = march.gas_velocity_m_s
            fields["Re_g"][i] = march.reynolds
            fields["Pr_g"][i] = march.prandtl
            fields["Nu_g"][i] = march.nusselt
            fields["f_g"][i] = march.friction_factor
            fields["dp_g_per_m"][i] = march.dp_per_length_Pa_m
            fields["h_removed_g"][i] = march.enthalpy_removed_J_kg
            fields["progress_g"][i] = march.progress_variable
            fields["G_shell"][i] = diag.shell_film.mass_flux_kg_m2s
            fields["Re_shell"][i] = diag.shell_film.reynolds
            fields["Pr_shell"][i] = diag.shell_film.prandtl
            fields["dp_shell"][i] = diag.shell_film.dp_shell_Pa
            fields["rho_c"][i] = diag.coolant_properties.rho
            fields["mu_c"][i] = diag.coolant_properties.mu
            fields["k_c"][i] = diag.coolant_properties.k
            fields["cp_c"][i] = diag.coolant_properties.cp
            scalars["Q_hot_kW"][i] = np.sum(march.wall_flux.hot_heat_W) / 1.0e3
            scalars["T_g_out"][i] = march.T_gas_outlet
            scalars["h_removed_g_out_J_kg"][i] = march.enthalpy_removed_outlet_J_kg
            scalars["progress_g_out"][i] = march.progress_outlet
            scalars["dp_g_total_Pa"][i] = np.sum(march.dp_per_length_Pa_m * geometry.grid.dx)
            scalars["dp_shell_total_Pa"][i] = _finite_mean(diag.shell_film.dp_shell_Pa)
            scalars["Re_g_max"][i] = _finite_max(march.reynolds)
            scalars["Re_shell_max"][i] = _finite_max(diag.shell_film.reynolds)
            bc = self._bc_at(float(t[i]))
            scalars["mdot_c"][i] = bc["mdot_c"]
            scalars["mdot_c_effective"][i] = diag.wall_coolant_inputs.mdot_coolant
            if "face_mdot_c" in fields:
                faces = fields["face_mdot_c"][i]
                if geometry.grid.flow_direction == 1:
                    scalars["mdot_c_inlet_face"][i] = faces[0]
                    scalars["mdot_c_outlet_face"][i] = faces[-1]
                else:
                    scalars["mdot_c_inlet_face"][i] = -faces[-1]
                    scalars["mdot_c_outlet_face"][i] = -faces[0]
            scalars["mdot_g"][i] = bc["mdot_g"] if bc["ignited"] else bc.get("mdot_lox", 0.0)
        return {
            "t": t,
            "x": geometry.grid.x_centers.copy(),
            "fields": fields,
            "scalars": scalars,
        }

    def _print_core_final_state(self, result):
        ts = self.time_series
        scalars = ts["scalars"]
        print(
            f"  final: T_g_out={scalars['T_g_out'][-1]:.1f} K  "
            f"T_c_out={scalars['T_c_out'][-1]:.1f} K  "
            f"Q_hot={scalars['Q_hot_kW'][-1]:.1f} kW  "
            f"T_wall_max={scalars['T_wall_max'][-1]:.1f} K"
        )

    def _solve_fixed_step(self, Tbar0, t_eval):
        """Bounded-cost linearly-implicit wall integrator.

        This is the fast production path for shell-and-tube counter-flow
        transients: one fluid pass per time step, with schedule breakpoints
        inserted into the grid so ignition/flow jumps are not stepped across.
        The local film stiffness is treated implicitly per wall node:

            Tbar_i^{n+1} = Tbar_i^n + dt * R_i(Tbar^n) / (1 + dt * lambda_i)

        This preserves the one-fluid-pass-per-step cost while preventing the
        helium film time constant from destabilizing practical time steps.
        """
        grid = self._fixed_time_grid(t_eval)
        y = np.zeros((self.N, len(grid)))
        y[:, 0] = Tbar0
        nfev = 0
        for j in range(len(grid) - 1):
            t = float(grid[j])
            dt = float(grid[j + 1] - grid[j])
            rhs, lam = self._wall_rate(t, y[:, j])
            nfev += 1
            y[:, j + 1] = np.clip(y[:, j] + dt * rhs / (1.0 + dt * lam), 80.0, 4000.0)
        return SimpleNamespace(
            t=grid,
            y=y,
            success=True,
            message="fixed_step (linearly-implicit) completed",
            nfev=nfev,
        )

    def _fixed_time_grid(self, t_eval):
        tp = self.transientProp
        base = np.arange(0.0, tp.t_end + 0.5 * tp.max_step, tp.max_step)
        points = [base, np.asarray(t_eval, dtype=float), np.array([0.0, tp.t_end])]
        for name in (
            "schedule_mass_flow_c", "schedule_mass_flow_g", "schedule_mass_flow_lox",
            "schedule_mass_flow_diesel", "schedule_T_c_in", "schedule_p_c_in",
            "schedule_T_lox_in", "schedule_ignition_state", "schedule_OF",
        ):
            schedule = getattr(tp, name, None)
            if schedule:
                points.append(np.array([float(t) for t, _ in schedule]))
        grid = np.unique(np.clip(np.concatenate(points), 0.0, tp.t_end))
        return grid

    def _print_final_state(self, sol):
        bc = self._bc_at(float(sol.t[-1]))
        res = self.fluid_pass(sol.y[:, -1], bc)
        Q = np.sum(res["dq_hot__dx"]) * self.dx * self.stp.N_tubes / 1e3
        print(f"  final: T_g_out={res['T_g_out']:.1f} K  T_c_out={res['T_c_out']:.1f} K  "
              f"Q={Q:.1f} kW  wall dT max={np.max(res['T_wg']-res['T_wc']):.1f} K")


def _finite_max(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.max(arr)) if arr.size else np.nan


def _finite_mean(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.mean(arr)) if arr.size else np.nan


if __name__ == "__main__":
    from .input_data import (coolantProp, hotgasProp, shellTubeProp, numericalProp,
                             system_requirements, transientProp)
    tp = transientProp()
    s = shellntube_transient_solver(coolantProp=coolantProp(), hotgasProp=hotgasProp(),
                                    shellTubeProp=shellTubeProp(), numericalProp=numericalProp(),
                                    system_requirements=system_requirements(), transientProp=tp,
                                    N_axial=80)
    s.solve_transient(verbose=True)
