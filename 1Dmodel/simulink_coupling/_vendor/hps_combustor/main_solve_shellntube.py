"""
Steady solver — baffled shell-and-tube heat exchanger (WP2).

@ author : Raphaël Aubry

Hot combustion gas flows axially inside straight tubes; coolant flows on the
shell side in zig-zag cross-flow around segmental baffles (the EchTherm-style
geometry, see doc/DESIGN_PLAN_shellntube_transient.md sections 2-3).

Architecture (section 1.3 / WP2):
  - Tube side (representative single tube): quasi-1D march, reusing the same
    combustion chemistry (Cantera) and blended tube-side Nu/friction
    correlations as the rest of the codebase.
  - Shell side: 1D energy march (axial), 0D-per-node Bell-Delaware h — NOT a
    momentum ODE (shell flow has no single streamwise coordinate;
    Bell-Delaware is the correct altitude, section 1.2). Gas mode holds
    pressure at the nominal inlet value (h unaffected either way). Liquid
    mode marches an actual (lagged) pressure profile — the lumped
    whole-bundle Bell-Delaware dp apportioned per node outside the two-phase
    dome, the Muller-Steinhagen-Heck friction gradient inside it — since
    local Tsat(p) matters once the coolant boils (see _shell_h_at).
  - Wall: reuses OneDimensionalSteadyConduction_ShellnHelicalTube with
    hot_side="inner" (hot fluid inside the tube — the opposite radial
    arrangement from the helical coil).
  - Both fluid inlets are prescribed (a genuine two-point BVP, unlike the
    helical solver's prescribed-outlet shortcut) → solved by the predictive
    sweep iteration of section 1.3: march the tube gas against the current
    shell-temperature field, march the shell energy against the resulting
    duty field, under-relax, repeat.

Run directly: python main_solve_shellntube.py
"""
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
    _rp.run_module(f"{_alias}.main_solve_shellntube", run_name="__main__", alter_sys=True)
    raise SystemExit(0)

import numpy as np
from CoolProp.CoolProp import PropsSI

from .physics.combustion_chemistry.combustion_gas import combustion_gas_solve, choose_fuel
from .physics.combustion_chemistry.fpv_manifold import build_fpv_manifold, FPVManifold
from .physics.heat_transfer_correlations import (
    dispatch_nu_tube_straight,
    nu_corrugated_tube_vicente,
)
from .physics.friction_correlations import (
    dispatch_friction_tube_straight,
    friction_corrugated_tube_vicente,
    getFrictionColebrook1939,
)
from .physics.heat_conduction import OneDimensionalSteadyConduction_ShellnHelicalTube
from .physics.bell_delaware import bell_delaware_shell
from .physics.liquid_flow.hx_adapters import solve_shelltube_shellside_liquid_from_tube_result
from .physics.liquid_flow.dispatch import evaluate_coolant_closure
from .physics.liquid_flow.correlations import (
    bergles_rohsenow_onb_wall_superheat,
    saturation_state,
)
from .physics.liquid_flow.regime import real_fluid_state_ph
from .physics.liquid_flow.coolprop_state_cache import coolprop_fluid_string, get_cached_state
from .physics.liquid_flow.chf import chf_regime
from .mechanical.geometry.shelltube_geometry import compute_bell_delaware_geometry
from .mechanical.material_specs.material_temperature_strength import init_material_temperature_properties
from .mechanical.loads import (stress_external_pressure_tube, collapse_pressure_thin_tube,
                               stress_thermal_tube)
from .input_data import CorrelationCoefficients


class shellntube_solver:
    def __init__(self, coolantProp, hotgasProp, shellTubeProp, numericalProp,
                system_requirements, corrCoeffs=None, N_axial=200,
                flow_config="counter"):
        self.coolantProp = coolantProp
        self.hotgasProp = hotgasProp
        self.stp = shellTubeProp
        self.numericalProp = numericalProp
        self.system_requirements = system_requirements
        self.corrCoeffs = corrCoeffs if corrCoeffs is not None else CorrelationCoefficients()
        self.N = int(N_axial)
        self.flow_config = flow_config   # "co" | "counter" (tube-gas always forward)
        self.dx = self.stp.L_tube / self.N

        # --- geometry ---
        self.geom = compute_bell_delaware_geometry(
            D_shell_inner=self.stp.D_shell_inner, D_tube_outer=self.stp.D_tube_outer,
            pitch_ratio=self.stp.pitch_ratio, layout=self.stp.layout,
            N_tubes=self.stp.N_tubes, N_baffles=self.stp.N_baffles,
            baffle_cut=self.stp.baffle_cut, L_tube=self.stp.L_tube,
            clearance_tube_baffle=self.stp.clearance_tube_baffle,
            clearance_baffle_shell=self.stp.clearance_baffle_shell,
            clearance_bundle_shell=self.stp.clearance_bundle_shell,
            N_sealing_strip_pairs=self.stp.N_sealing_strip_pairs,
            baffle_spacing=self.stp.baffle_spacing,
            L_inlet_spacing=self.stp.L_inlet_spacing,
            L_outlet_spacing=self.stp.L_outlet_spacing)

        self.D_tube_i = self.stp.D_tube_outer - 2 * self.stp.thickness_tube_wall
        self.A_tube_i = np.pi * self.D_tube_i ** 2 / 4.0
        self.P_tube_i = np.pi * self.D_tube_i

        # --- materials ---
        (self.CTE_t, self.E_t, self.Yield_t, self.k_t, self.rho_t, self.poisson_t, self.cp_t
         ) = init_material_temperature_properties(self.stp.material_tube)

        # --- combustion inlet state (tube side) ---
        self.chem_mech_path, self.Y_fuel, self.Hv_fuel = choose_fuel(self.hotgasProp.fuel)
        self.combustion_node = combustion_gas_solve(
            fuel=self.hotgasProp.fuel, oxidizer=self.hotgasProp.oxidizer,
            OF=self.hotgasProp.mixing_ratio, T_inj_LOX=self.hotgasProp.T_inj_LOX,
            T_g_init=self.hotgasProp.T_g_init, p0=self.hotgasProp.p0,
            chem_mech_path=self.chem_mech_path, Hv_fuel=self.Hv_fuel, Y_fuel=self.Y_fuel)
        self.combustion_node.solve()
        self.gas_phase = self.combustion_node.phase
        # cache the inlet combustion state ONCE — _tube_side_march must reset to this
        # every sweep, since gas_phase is a persistent object that the previous sweep's
        # march left in a cooled-down state.
        self._gas_inlet_TPY = (float(self.gas_phase.T), float(self.hotgasProp.p0),
                               np.array(self.gas_phase.Y, dtype=float))
        self._fpv = None
        if self.numericalProp.chemistry_model == "finite_rate":
            self._setup_fpv_manifold()

        # per-tube mass flow
        self.mdot_tube = self.hotgasProp.mass_flow_g / self.stp.N_tubes
        self._tube_Nu_factor, self._tube_f_factor = self._tube_surface_factors()
        self._corrugation_phi = self._tube_corrugation_severity()

        # Shell-side liquid coupling (Phase 3 of
        # docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md): when
        # coolant_model == "equilibrium_liquid", _shell_side_march advances
        # (p, h) instead of T (no phase-change/latent-heat accounting in a
        # plain T-based march would otherwise let "boiling" coolant heat past
        # 1000+ K, since cp*dT alone has nowhere to put the latent heat), and
        # _shell_h_at switches from Bell-Delaware (single-phase only, no
        # boiling model) to the same validated Gungor-Winterton/CHF closure
        # used by the helical solver once inside the two-phase dome. Shell
        # pressure is now a real (lagged) running profile fed by both the
        # two-phase Muller-Steinhagen-Heck friction gradient and the
        # single-phase Bell-Delaware total dp, apportioned per axial node
        # (see _shell_h_at/_shell_side_march) — no longer held at the nominal
        # coolantProp.p_in everywhere, though this remains a gas-mode
        # simplification (gas mode is untouched). Cross-flow enhancement of
        # boiling HTC (a real, shell-specific effect Gungor-Winterton does not
        # capture, since it is a tube-flow correlation) is still NOT modeled —
        # see docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md Phase 3.
        self._liquid_mode = self.coolantProp.coolant_model == "equilibrium_liquid"
        # Lagged shell-side wall heat flux for the boiling HTC's heat-flux
        # (Bo) term - one sweep-iteration behind, seeded at zero, same
        # first-order-accurate lagged-closure pattern used by the helical
        # solver's per-node boiling HTC (see main_solve.py).
        self._shell_qw_lagged = np.zeros(self.N)
        # Lagged per-node shell-side pressure drop [Pa], same one-sweep-behind
        # pattern as _shell_qw_lagged - see _shell_h_at's dp_node return value.
        self._shell_dp_lagged = np.zeros(self.N)
        # Lagged per-node shell-side (coolant-side) wall temperature [K], same
        # one-sweep-behind pattern - feeds the supercritical property-ratio
        # closures' (T_s/T_b) correction (see _shell_h_at). Seeded at the
        # physical inlet temperature (a neutral "wall ~ bulk" start, same
        # cold-start convention as q_w/dp being seeded at zero).
        self._shell_Tw_lagged = np.full(self.N, float(self.coolantProp.T_in))
        # Worst (largest) Bergles-Rohsenow ONB wall-superheat margin seen on
        # the shell side across all sweeps - same diagnostic as main_solve.py.
        self._shell_onb_max_margin = float("-inf")
        # Supercritical-regime diagnostics (populated only when the coolant is
        # above its critical pressure; see _shell_h_at's supercritical branch).
        self._sc_closure_name = None
        self._sc_regimes = set()
        self._sc_htd_nodes = 0
        self._sc_extrapolated = False
        self._sc_extrap_message = None
        # Opt-in tabulated CoolProp property backend for the liquid march
        # (default "HEOS" = exact, unchanged behavior). Validated once here
        # (fail fast on a typo) rather than per march node.
        self._liquid_backend = getattr(self.coolantProp, "liquid_property_backend", "HEOS")
        if self._liquid_mode:
            coolprop_fluid_string(self.coolantProp.coolant, self._liquid_backend)

    def _tube_surface_factors(self):
        """Return EchTherm-style tube-side Nu/friction multipliers.

        Smooth tubes use the published smooth-tube correlations directly.
        Grooved tubes use the Vicente/Cruz corrugated-tube path plus optional
        calibration multipliers. Intensification-factor remains a direct
        multiplier hook.
        """
        choice = self.stp.inside_tube_choice
        if choice == "smooth":
            return 1.0, 1.0
        if choice == "grooved":
            return self.corrCoeffs.tube_grooved_Nu_factor, self.corrCoeffs.tube_grooved_f_factor
        if choice == "intensification_factor":
            return self.corrCoeffs.tube_intensification_factor, self.corrCoeffs.tube_intensification_factor
        if choice in ("helical_insert", "power_law"):
            raise NotImplementedError(
                f"inside_tube_choice={choice!r} is exposed from EchTherm but not implemented yet."
            )
        raise ValueError(f"Unsupported inside_tube_choice: {choice!r}")

    def _tube_corrugation_severity(self):
        """Vicente/Cruz helical corrugation severity index phi = e^2/(p*D_i)."""
        e = max(float(self.stp.corrugation_thickness), 0.0)
        p = max(float(self.stp.corrugation_pitch), 1e-12)
        return e ** 2 / (p * max(self.D_tube_i, 1e-12))

    def _tube_side_hydraulics(self, Re_g, Pr_g, x_local, T_bulk, T_wall):
        """Return tube-side Darcy friction factor and Nusselt number."""
        choice = self.stp.inside_tube_choice
        if choice == "grooved":
            f_g = friction_corrugated_tube_vicente(
                Re_g, self._corrugation_phi,
                Re_lo=self.corrCoeffs.Re_transition_lo,
                Re_hi=self.corrCoeffs.Re_transition_hi,
            )
            Nu_g = nu_corrugated_tube_vicente(
                Re_g, Pr_g, self._corrugation_phi,
                D_i=self.D_tube_i, x=x_local,
                Re_lo=self.corrCoeffs.Re_transition_lo,
                Re_hi=self.corrCoeffs.Re_transition_hi,
            )
        else:
            f_g = dispatch_friction_tube_straight(
                Re_g, self.stp.tube_roughness, self.D_tube_i, x=x_local,
                Re_lo=self.corrCoeffs.Re_transition_lo,
                Re_hi=self.corrCoeffs.Re_transition_hi,
            )
            Nu_g = dispatch_nu_tube_straight(
                self.stp.Nusselt_tube, Re=Re_g, Pr=Pr_g, d=self.D_tube_i, x=x_local,
                f_fd=f_g, T_bulk=T_bulk, T_wall=T_wall, error_factor=1.0,
                corrCoeffs=self.corrCoeffs,
            )
        return f_g * self._tube_f_factor, Nu_g * self._tube_Nu_factor

    def _setup_fpv_manifold(self):
        """Build the finite-rate FPV cooling manifold for the tube-side gas."""
        Tg0, pg0, Yg0 = self._gas_inlet_TPY
        gas = self.gas_phase
        self._fpv = FPVManifold(build_fpv_manifold(
            gas,
            Y_inlet=Yg0,
            T_inlet=Tg0,
            p=pg0,
            species_index={
                "CO2": gas.species_index("CO2"),
                "H2O": gas.species_index("H2O"),
                "CO": gas.species_index("CO"),
            },
            n_h=self.numericalProp.fpv_n_h,
            n_c=self.numericalProp.fpv_n_c,
            t_relax=self.numericalProp.fpv_t_relax,
            n_t=self.numericalProp.fpv_n_t,
            cache_dir=self.numericalProp.fpv_cache_dir,
        ))

    # ------------------------------------------------------------------
    def _tube_side_march(self, T_shell_profile, liquid_state=None):
        """March the representative tube from x=0 (gas inlet) to L, given the CURRENT
        shell-temperature field T_shell_profile[i] (K, same axial grid). Returns
        per-node dQ [W, per tube], T_g, T_wg, T_wc, p_g arrays, and gas outlet state.

        ``liquid_state``, when not None, is the enthalpy/quality field from the
        LAST sweep's ``_shell_side_march`` (dict with ``h`` and ``quality``
        arrays, same axial index convention as T_shell_profile) - only used
        when ``coolantProp.coolant_model == "equilibrium_liquid"``, to pick
        the correct coolant-side HTC branch in ``_shell_h_at``.
        """
        N, dx = self.N, self.dx
        Tg0, pg0, Yg0 = self._gas_inlet_TPY
        self.gas_phase.TPY = Tg0, pg0, Yg0
        T_g, p_g = Tg0, pg0
        h_removed = 0.0
        Yc = self._fpv.Yc_inlet() if self._fpv is not None else 0.0
        W_g = self.gas_phase.mean_molecular_weight
        P_tube_o = np.pi * self.stp.D_tube_outer

        dQ = np.zeros(N); T_g_a = np.zeros(N); T_wg_a = np.zeros(N); T_wc_a = np.zeros(N)
        h_g_a = np.zeros(N); q_w_shell_a = np.zeros(N); chf_margin_a = np.full(N, np.nan)
        p_g_a = np.zeros(N); h_c_a = np.zeros(N); dp_shell_a = np.zeros(N)
        for i in range(N):
            if self._fpv is not None:
                T_g, rho_g, mu_g, k_g, cp_g, _xH2O, _xCO2, omega_Yc = self._fpv.state(
                    h_removed,
                    Yc,
                )
            else:
                g = self.gas_phase
                cp_g = g.cp; mu_g = g.viscosity; k_g = g.thermal_conductivity; rho_g = g.density
            U_g = self.mdot_tube / (rho_g * self.A_tube_i)
            Re_g = rho_g * U_g * self.D_tube_i / mu_g
            Pr_g = cp_g * mu_g / k_g
            x_local = (i + 0.5) * dx
            T_s_local = T_shell_profile[i]
            f_g, Nu_g = self._tube_side_hydraulics(Re_g, Pr_g, x_local, T_g, T_s_local)
            h_g = Nu_g * k_g / self.D_tube_i

            # Shell-side coolant flow length from ITS OWN entrance (co-flow:
            # index 0; counter-flow: index N-1 - same convention as the
            # h/quality/p axial index everywhere else in this file), over the
            # tube outer diameter - Taylor's supercritical entrance correction.
            x_local_shell = (i * dx) if self.flow_config == "co" else ((N - 1 - i) * dx)
            x_over_D_shell = max(x_local_shell, 0.5 * dx) / self.stp.D_tube_outer

            # shell-side h at this node (recomputed locally — see module docstring)
            h_c, dp_c_node = self._shell_h_at(
                T_s_local,
                h_c_enthalpy=(liquid_state["h"][i] if liquid_state is not None else None),
                quality_local=(liquid_state["quality"][i] if liquid_state is not None else None),
                q_w_lagged=self._shell_qw_lagged[i],
                p_local=(liquid_state["p"][i] if (liquid_state is not None and "p" in liquid_state) else None),
                wall_temp_K=self._shell_Tw_lagged[i],
                x_over_D=x_over_D_shell,
            )
            chf_margin_a[i] = self._last_shell_chf_margin if self._last_shell_chf_margin is not None else np.nan
            h_c_a[i] = h_c
            dp_shell_a[i] = dp_c_node
            if self._last_shell_sc_info is not None:
                info = self._last_shell_sc_info
                self._sc_closure_name = info["closure_name"]
                self._sc_regimes.add(info["regime"])
                if info["htd_risk"]:
                    self._sc_htd_nodes += 1
                rep = info["extrapolation_report"]
                if rep is not None and not rep.in_range:
                    self._sc_extrapolated = True
                    if self._sc_extrap_message is None:
                        self._sc_extrap_message = rep.message()

            node = OneDimensionalSteadyConduction_ShellnHelicalTube(
                h_g=h_g, h_c=h_c, T_c=T_s_local, T_g=T_g,
                s_w=self.stp.thickness_tube_wall, Dh_ch=self.D_tube_i,
                f_kw_at_T=self.k_t, T_wg_0=T_g, T_wc_0=T_s_local, T_c_check_0=T_s_local,
                dx=dx, rad_enabled=False, hot_side="inner")
            node.Solve1Dconduction()

            dQ[i] = node.dQ
            T_g_a[i] = T_g; T_wg_a[i] = node.T_wg; T_wc_a[i] = node.T_wc; h_g_a[i] = h_g
            p_g_a[i] = p_g
            # Wall heat flux into the shell-side coolant at this node (per unit
            # outer tube area) - stored for NEXT sweep's lagged boiling Bo term.
            q_w_shell_a[i] = node.dQ / (P_tube_o * dx)

            dh_g = node.dQ / self.mdot_tube
            if self._fpv is not None:
                h_removed += dh_g
                if U_g > 0:
                    Yc += omega_Yc / U_g * dx
            else:
                if self.numericalProp.chemistry_model not in ("equilibrium", "frozen"):
                    raise ValueError(
                        f"Unsupported steady chemistry_model: {self.numericalProp.chemistry_model!r}"
                    )
                self.combustion_node.remove_energy(
                    dh=dh_g,
                    updated_pressure=p_g,
                    equilibrium_dh_gas_ON=(self.numericalProp.chemistry_model == "equilibrium"),
                )
                T_g = self.gas_phase.T
            p_g -= f_g * rho_g * U_g ** 2 / (2 * self.D_tube_i) * dx

        self._shell_qw_lagged = q_w_shell_a
        self._shell_dp_lagged = dp_shell_a
        self._shell_Tw_lagged = T_wc_a
        return dict(dQ=dQ, T_g=T_g_a, T_wg=T_wg_a, T_wc=T_wc_a, h_g=h_g_a, h_c=h_c_a, p_g=p_g_a,
                   q_w_shell=q_w_shell_a, chf_margin=chf_margin_a, dp_shell=dp_shell_a,
                   T_g_out=T_g, p_g_out=p_g)

    def _shell_h_at(self, T_shell_local, h_c_enthalpy=None, quality_local=None, q_w_lagged=0.0,
                    p_local=None, wall_temp_K=None, x_over_D=None):
        """Shell-side coolant-side heat transfer coefficient at one axial node.

        Gas mode (unchanged): Bell-Delaware baffled cross-flow HTC, properties
        from a (T, p) CoolProp flash, nominal pressure (no gas-side shell dp
        march - out of scope, see module docstring).

        Liquid mode: outside the two-phase dome AND outside the supercritical
        range (subcooled liquid or superheated vapor), properties come from a
        (p, h) flash instead of (T, p) - required for phase-change consistency,
        since T alone cannot distinguish liquid/vapor branches near saturation
        - and Bell-Delaware is still used (it is a valid single-phase
        correlation either side of the dome). Inside the dome (0<=quality<=1)
        OR at supercritical pressure (quality is NaN there - no dome exists),
        Bell-Delaware has no boiling/supercritical physics at all, so this
        switches to evaluate_coolant_closure, which internally dispatches to
        the validated Gungor-Winterton/CHF closure (subcritical dome) or the
        McCarthy-Wolf/Taylor property-ratio registry (supercritical) the
        helical solver uses, with the tube outer diameter and shell-side
        cross-flow mass flux (mdot_c / S_m) as the characteristic scale — a
        documented simplification either way (shell-side cross-flow boiling
        AND supercritical enhancement are real, shell-specific effects that
        neither Gungor-Winterton nor McCarthy-Wolf/Taylor, all tube-flow
        correlations, capture; see __init__'s Phase 3 note and
        docs/solver_design/FLUID_AGNOSTIC_CLOSURES_AND_SUPERCRITICAL_PLAN.md).
        ``p_local`` is the running (lagged) shell pressure at this node in
        liquid mode - falls back to the nominal coolantProp.p_in when None
        (gas mode, or the first sweep before a liquid_state pressure field
        exists). ``wall_temp_K``/``x_over_D`` are the lagged shell-side wall
        temperature and flow-length/diameter ratio, forwarded to the
        supercritical property-ratio closures only (ignored subcritically).

        Returns (h_shell [W/m2K], dp_node_Pa): the latter is this node's
        share of shell-side pressure drop, used by _shell_side_march to march
        an actual (lagged) pressure profile instead of holding pressure flat.
        Two-phase dome / supercritical: dp_node = evaluate_coolant_closure's
        friction gradient (already computed for the HTC call above, previously
        discarded) times dx - a straight-pipe correlation, the same documented
        tube-flow-vs-cross-flow simplification as the boiling/supercritical HTC
        itself. Single-phase subcritical branch: the lumped whole-bundle
        Bell-Delaware dp_shell (computed with the LOCAL density/viscosity at
        this node - important once density swings through the boiling
        transition), apportioned evenly over the N axial nodes - an
        approximation (assumes each axial slice carries an equal share of the
        total drop; a true baffle-by-baffle discretization would need
        re-deriving Bell-Delaware's per-crossing terms) but at least uses the
        correct cross-flow-bundle friction model and the correct local
        density, instead of the previous constant-p_in simplification. Gas
        mode returns dp_node=0.0 (unchanged; no gas-side pressure march).
        """
        cool = self.coolantProp.coolant
        cool_cp = coolprop_fluid_string(cool, self._liquid_backend) if self._liquid_mode else cool
        p_c = p_local if p_local is not None else self.coolantProp.p_in

        self._last_shell_sc_info = None
        # quality_local is NaN at supercritical pressure (no dome, undefined
        # quality - see regime.supercritical_state_ph) - route it through
        # evaluate_coolant_closure exactly like the subcritical two-phase dome,
        # since that function already dispatches correctly on both internally.
        if self._liquid_mode and quality_local is not None and (
            (0.0 <= quality_local <= 1.0) or np.isnan(quality_local)
        ):
            mass_flux = self.coolantProp.mass_flow_c / self.geom["S_m"]
            closure = evaluate_coolant_closure(
                coolant_prop=self.coolantProp,
                p_Pa=p_c,
                h_J_kg=h_c_enthalpy,
                mass_flux_kg_m2_s=mass_flux,
                hydraulic_diameter_m=self.stp.D_tube_outer,
                heat_flux_W_m2=max(q_w_lagged, 0.0),
                lut_path=self.coolantProp.liquid_chf_lut_path,
                wall_temp_K=wall_temp_K,
                # Shell-side cross-flow has no vertical/horizontal tube
                # orientation of its own; "horizontal" is an arbitrary but
                # harmless default (no registered supercritical closure hard-
                # filters on orientation, only regime does - see registry.py).
                orientation="horizontal",
                geometry="shell_crossflow",
                x_over_D=x_over_D,
            )
            self._last_shell_chf_margin = closure.chf_margin
            self._last_bell = None
            if closure.regime is not None:
                self._last_shell_sc_info = dict(
                    regime=closure.regime, closure_name=closure.closure_name,
                    htd_risk=closure.htd_risk, extrapolation_report=closure.extrapolation_report,
                )
            return closure.htc_W_m2_K, closure.dpdz_friction_Pa_m * self.dx

        self._last_shell_chf_margin = None
        if self._liquid_mode:
            flashed = get_cached_state(cool_cp).flash_ph(p_c, h_c_enthalpy)
            rho = flashed.rhomass()
            mu = flashed.viscosity()
            k = flashed.conductivity()
            cp = flashed.cpmass()
        else:
            rho = PropsSI('D', 'T', T_shell_local, 'P', p_c, cool)
            mu = PropsSI('V', 'T', T_shell_local, 'P', p_c, cool)
            k = PropsSI('L', 'T', T_shell_local, 'P', p_c, cool)
            cp = PropsSI('C', 'T', T_shell_local, 'P', p_c, cool)
        Pr = cp * mu / k
        G_s = self.coolantProp.mass_flow_c / self.geom["S_m"]
        Re_s = self.stp.D_tube_outer * G_s / mu
        geom = dict(self.geom); geom["rho_s"] = rho
        r = bell_delaware_shell(geom, Re_s=Re_s, Pr_s=Pr, k_s=k, cp_s=cp, mu_s=mu,
                                mdot_s=self.coolantProp.mass_flow_c, corrCoeffs=self.corrCoeffs)
        self._last_bell = r

        if self._liquid_mode and q_w_lagged > 0.0 and quality_local is not None and quality_local < 0.0:
            # Bergles-Rohsenow ONB check (see main_solve.py for the same
            # criterion on the helical solver): estimate the local wall
            # temperature from this node's convective h_shell and lagged wall
            # flux, and flag if nucleate boiling is physically expected even
            # though bulk quality is still subcooled.
            T_wall_est = T_shell_local + q_w_lagged / max(r["h_shell"], 1e-9)
            sat_local = saturation_state(cool_cp, p_c)
            dT_wall_superheat = T_wall_est - sat_local.T_sat_K
            dT_onb = bergles_rohsenow_onb_wall_superheat(p_Pa=p_c, heat_flux_W_m2=q_w_lagged)
            self._shell_onb_max_margin = max(self._shell_onb_max_margin, dT_wall_superheat - dT_onb)

        dp_node = (r["dp_shell"] / self.N) if self._liquid_mode else 0.0
        return r["h_shell"], dp_node

    def _shell_side_march(self, dQ_profile):
        """March shell-side coolant energy from ITS inlet using the duty field
        dQ_profile [W, per tube] -> total per node = dQ_profile * N_tubes.
        Direction depends on flow_config (co: same end as gas inlet; counter:
        opposite end). Returns T_shell array on the SAME axial index convention
        as dQ_profile (index 0 = gas inlet end, regardless of coolant direction),
        the coolant outlet temperature, and a liquid_state dict (h, quality,
        void arrays in the same index convention) when
        coolantProp.coolant_model == "equilibrium_liquid", else None.

        Liquid mode marches ENTHALPY (h += dQ/mdot), not T (T_cur += dQ/(mdot*cp)):
        a plain cp*dT update has no latent-heat sink, so "boiling" coolant would
        just keep heating past 1000+ K instead of pinning near T_sat - exactly
        the failure this branch exists to avoid. Shell pressure now marches
        alongside enthalpy using the PREVIOUS sweep's per-node dp
        (self._shell_dp_lagged, populated by _shell_h_at during
        _tube_side_march) - the same one-sweep-lagged pattern as the boiling
        HTC's heat-flux term, replacing the old flat coolantProp.p_in
        assumption (gas mode is unaffected)."""
        N, dx = self.N, self.dx
        cool = self.coolantProp.coolant
        dQ_total = dQ_profile * self.stp.N_tubes

        T = np.zeros(N)

        if not self._liquid_mode:
            if self.flow_config == "co":
                T_cur = self.coolantProp.T_in
                for i in range(N):
                    T[i] = T_cur
                    cp_c = PropsSI('C', 'T', T_cur, 'P', self.coolantProp.p_in, cool)
                    T_cur += dQ_total[i] / (self.coolantProp.mass_flow_c * cp_c)
                T_out = T_cur
            else:  # counter: coolant enters at the gas-OUTLET end (index N-1) and marches to 0
                T_cur = self.coolantProp.T_in
                for i in range(N - 1, -1, -1):
                    T[i] = T_cur
                    cp_c = PropsSI('C', 'T', T_cur, 'P', self.coolantProp.p_in, cool)
                    T_cur += dQ_total[i] / (self.coolantProp.mass_flow_c * cp_c)
                T_out = T_cur  # exits at index 0 (gas inlet end)
            return T, T_out, None

        cool_cp = coolprop_fluid_string(cool, self._liquid_backend)
        p_in = self.coolantProp.p_in
        h_in = PropsSI('H', 'T', self.coolantProp.T_in, 'P', p_in, cool)
        quality = np.zeros(N); void = np.zeros(N); h_arr = np.zeros(N); p_arr = np.zeros(N)
        # Very low coolant flow relative to duty can superheat the coolant
        # past the fluid's EOS validity ceiling (e.g. CoolProp's Water backend
        # is only valid up to 3000 K) - a real physical edge case (the duty
        # available vastly exceeds what this little coolant flow can absorb
        # even as superheated vapor), not a numerical bug. Hold the last
        # EOS-valid (state, h, p) for the remainder of the march instead of
        # crashing.
        last_valid = {"state": None, "h": None, "p": None}

        def _advance(h_cur, p_cur, idx):
            try:
                # real_fluid_state_ph: dome-based below p_crit (bit-identical
                # to the former equilibrium_state_ph call), single-phase
                # real-EOS above it (supercritical coolant does not crash on
                # the ValueError this still needs to catch for a genuine
                # EOS-ceiling overshoot, e.g. very high superheat).
                state = real_fluid_state_ph(cool_cp, p_cur, h_cur)
                last_valid["state"] = state
                last_valid["h"] = h_cur
                last_valid["p"] = p_cur
            except ValueError:
                # Freeze h AND p at the last EOS-valid values - otherwise they
                # keep drifting node to node, and the frozen-state/drifting
                # combination would re-trigger the same ValueError on the
                # NEXT sweep's _shell_h_at call (which reads them back out of
                # liquid_state["h"]/["p"]).
                state = last_valid["state"]
                h_cur = last_valid["h"]
                p_cur = last_valid["p"]
            T[idx] = state.T_K
            h_arr[idx] = h_cur
            p_arr[idx] = p_cur
            quality[idx] = state.quality
            void[idx] = state.void_fraction
            h_next = h_cur + dQ_total[idx] / self.coolantProp.mass_flow_c
            p_next = max(p_cur - self._shell_dp_lagged[idx], 1.0)
            return h_next, p_next

        h_cur, p_cur = h_in, p_in
        if self.flow_config == "co":
            for i in range(N):
                h_cur, p_cur = _advance(h_cur, p_cur, i)
        else:  # counter: coolant enters at the gas-OUTLET end (index N-1) and marches to 0
            for i in range(N - 1, -1, -1):
                h_cur, p_cur = _advance(h_cur, p_cur, i)
        try:
            T_out = real_fluid_state_ph(cool_cp, p_cur, h_cur).T_K
        except ValueError:
            T_out = last_valid["state"].T_K
        return T, T_out, dict(h=h_arr, quality=quality, void=void, p=p_arr, h_out=h_cur, p_out=p_cur)

    # ------------------------------------------------------------------
    def solve(self, max_sweeps=25, omega=0.5, tol=0.05, verbose=True):
        """Predictive sweep iteration (section 1.3): converge the shell temperature
        field against the tube-side duty field it produces."""
        N = self.N
        liquid_state = None
        if self._liquid_mode:
            # coolantProp.T_out is not meaningful for
            # coolant_model="equilibrium_liquid" (same reason the helical
            # solver's counter-flow shooting reference exists - a single (T,p)
            # pair cannot seed a genuine two-phase state), so seed uniformly
            # at the physical subcooled inlet instead of a T_in/T_out guess.
            T_shell = np.full(N, self.coolantProp.T_in)
            p_c = self.coolantProp.p_in
            h_in = PropsSI('H', 'T', self.coolantProp.T_in, 'P', p_c, self.coolantProp.coolant)
            eq0 = real_fluid_state_ph(self.coolantProp.coolant, p_c, h_in)
            liquid_state = dict(h=np.full(N, h_in), quality=np.full(N, eq0.quality),
                                void=np.full(N, eq0.void_fraction), p=np.full(N, p_c),
                                h_out=h_in)
        else:
            # initial guess: linear profile between the two prescribed inlets
            T_guess_out = 0.5 * (self.coolantProp.T_in + self.coolantProp.T_out)
            T_shell = np.linspace(self.coolantProp.T_in, T_guess_out, N)

        for sweep in range(max_sweeps):
            tube = self._tube_side_march(T_shell, liquid_state=liquid_state)
            T_shell_new, T_c_out, liquid_state = self._shell_side_march(tube["dQ"])
            delta = np.max(np.abs(T_shell_new - T_shell))
            T_shell = (1 - omega) * T_shell + omega * T_shell_new
            if verbose:
                print(f"  sweep {sweep:2d}: max|dT_shell|={delta:7.3f} K  "
                      f"T_g_out={tube['T_g_out']:.1f}  T_c_out={T_c_out:.1f}")
            if delta < tol:
                break
        else:
            print(f"  WARNING: sweep did not converge to {tol} K in {max_sweeps} sweeps "
                  f"(last delta={delta:.3f} K)")

        self.tube = tube
        self.T_shell = T_shell
        self.T_c_out = T_c_out
        # Liquid-mode (p,h) shell-side result (h, quality, void arrays, HX
        # axial order) from the final converged sweep - None in gas mode.
        # tube["chf_margin"] is the matching per-node CHF margin (also NaN
        # outside the two-phase dome or without coolantProp.liquid_chf_lut_path).
        self.shell_liquid = liquid_state
        self.Q_tot = np.sum(tube["dQ"]) * self.stp.N_tubes
        self.n_sweeps = sweep + 1
        return self

    # ------------------------------------------------------------------
    def compute_stress(self):
        """External-pressure stress/collapse check for the tubes (90 bar shell
        outside, low tube-side pressure — opposite loading from the helical coil)."""
        T_wg_max = np.max(self.tube["T_wg"])
        dP_ext = self.coolantProp.p_in - self.hotgasProp.p0  # shell minus tube pressure
        sigma_ext = stress_external_pressure_tube(dP_ext, self.stp.thickness_tube_wall,
                                                  self.D_tube_i)
        P_cr = collapse_pressure_thin_tube(self.E_t((T_wg_max + np.min(self.tube["T_wc"])) / 2 - 273),
                                           self.stp.thickness_tube_wall, self.D_tube_i, self.poisson_t)
        self.sigma_ext = sigma_ext
        self.P_cr = P_cr
        self.collapse_margin = abs(dP_ext) / P_cr
        return dict(sigma_ext=sigma_ext, P_cr=P_cr, collapse_margin=self.collapse_margin)

    def liquid_coolant_postprocess(self, lut_path=None, min_pressure_Pa=1.0):
        """Map completed shell-and-tube duty into the liquid coolant solver.

        This opt-in bridge consumes ``self.tube["dQ"]`` from the converged steady
        shell-and-tube solve and returns liquid p-h fields/diagnostics in the
        solver's axial order. It does not alter the existing shell-side
        temperature-only march.
        """
        result = solve_shelltube_shellside_liquid_from_tube_result(
            coolant_prop=self.coolantProp,
            shelltube_prop=self.stp,
            tube_result=self.tube,
            coolant_enters_at="z_min" if self.flow_config == "co" else "z_max",
            lut_path=lut_path,
            min_pressure_Pa=min_pressure_Pa,
        )
        self.liquid_coolant = result
        return result

    def print_summary(self):
        d = self.tube
        print("=" * 55)
        print("SHELL-AND-TUBE RESULTS")
        print("=" * 55)
        print(f"  N_tubes={self.stp.N_tubes}  N_baffles={self.stp.N_baffles}  "
              f"sweeps={self.n_sweeps}")
        print(f"  Q_tot = {self.Q_tot/1e3:.2f} kW")
        print(f"  T_g:  {d['T_g'][0]:.1f} -> {d['T_g_out']:.1f} K")
        print(f"  T_c:  {self.coolantProp.T_in:.1f} -> {self.T_c_out:.1f} K")
        print(f"  T_wg max = {np.max(d['T_wg']):.1f} K  T_wc max = {np.max(d['T_wc']):.1f} K")
        if self._last_bell is not None:
            print(f"  Bell-Delaware (last node): Jc={self._last_bell['Jc']:.3f} "
                  f"Jl={self._last_bell['Jl']:.3f} Jb={self._last_bell['Jb']:.3f}")
        else:
            print("  Bell-Delaware: n/a (last node was in the two-phase boiling closure)")
        if self.shell_liquid is not None:
            q = np.asarray(self.shell_liquid["quality"], dtype=float)
            v = np.asarray(self.shell_liquid["void"], dtype=float)
            p_arr = np.asarray(self.shell_liquid.get("p", []), dtype=float)
            cm = np.asarray(d["chf_margin"], dtype=float)
            finite_mask = np.isfinite(cm)
            finite_cm = cm[finite_mask]
            if finite_cm.size:
                worst_idx = np.flatnonzero(finite_mask)[np.argmin(finite_cm)]
                regime = chf_regime(float(q[worst_idx]))
                cm_str = f"{finite_cm.min():.2f} ({regime})"
            else:
                cm_str = "n/a"
            finite_q = q[np.isfinite(q)]
            if finite_q.size:
                print(f"  quality: {finite_q.min():.3f} -> {finite_q.max():.3f}  |  "
                      f"void max = {v.max():.3f}  |  min CHF margin = {cm_str}")
            if p_arr.size:
                print(f"  shell pressure: {p_arr.max()/1e5:.3f} -> {p_arr.min()/1e5:.3f} bar "
                      f"(inlet -> outlet, lagged march)")
            if self._shell_onb_max_margin > 0.0:
                print(f"  WARNING: Bergles-Rohsenow ONB criterion exceeded on the shell side "
                      f"(max wall-superheat margin = {self._shell_onb_max_margin:.2f} K) -- "
                      f"subcooled nucleate boiling physically expected before the numerical "
                      f"blend window")
            if self._sc_closure_name is not None:  # supercritical coolant somewhere in the march
                regimes = ", ".join(sorted(self._sc_regimes))
                print(
                    f"  Supercritical coolant: closure={self._sc_closure_name} | "
                    f"regimes traversed: {regimes} | HTD-risk nodes={self._sc_htd_nodes}"
                )
                print(
                    "  WARNING: shell-side cross-flow has no validated supercritical closure "
                    "(McCarthy-Wolf/Taylor are straight-tube fits) -- applied here as an "
                    "unvalidated geometry extrapolation, same caveat as the shell-side boiling gap."
                )
                if self._sc_extrapolated:
                    print(f"  WARNING: {self._sc_extrap_message} (validity extrapolation)")
                if self._sc_htd_nodes > 0:
                    print(
                        f"  WARNING: heat-transfer deterioration (HTD) risk flagged at "
                        f"{self._sc_htd_nodes} node(s) (McEligot-Jackson buoyancy/flow-"
                        f"acceleration criteria) -- deteriorated-HTC magnitude is NOT modeled; "
                        f"treat these nodes as outside the validated envelope."
                    )
        if hasattr(self, "collapse_margin"):
            print(f"  external-pressure collapse margin = {self.collapse_margin:.3f} "
                  f"(dP/P_cr; <1 = safe)")
        print("=" * 55)


if __name__ == "__main__":
    from .input_data import coolantProp, hotgasProp, shellTubeProp, numericalProp, system_requirements

    solver = shellntube_solver(
        coolantProp=coolantProp(), hotgasProp=hotgasProp(), shellTubeProp=shellTubeProp(),
        numericalProp=numericalProp(), system_requirements=system_requirements(),
        N_axial=150, flow_config="counter")
    solver.solve(verbose=True)
    solver.compute_stress()
    solver.print_summary()
