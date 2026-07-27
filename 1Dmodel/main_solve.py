"""
@ author : Raphaël Aubry
"""

# Allow running this file directly (play button, `python main_solve.py`).
# When executed as a script __package__ is None, so relative imports fail.
# This block registers the package under a stable alias (the folder name may
# start with a digit and therefore cannot be a Python identifier), then
# re-runs the module via runpy with the package context set up.
if __name__ == "__main__" and __package__ is None:
    import sys as _sys, os as _os, importlib.util as _ilu, runpy as _rp
    _pkg_dir = _os.path.dirname(_os.path.abspath(__file__))
    _parent  = _os.path.dirname(_pkg_dir)
    _alias   = "_hps"  # stable alias regardless of the on-disk folder name
    if _alias not in _sys.modules:
        _spec = _ilu.spec_from_file_location(
            _alias,
            _os.path.join(_pkg_dir, "__init__.py"),
            submodule_search_locations=[_pkg_dir],
        )
        _pkg = _ilu.module_from_spec(_spec)
        _pkg.__path__ = [_pkg_dir]
        _pkg.__package__ = _alias
        _sys.modules[_alias] = _pkg
        _spec.loader.exec_module(_pkg)
    if _parent not in _sys.path:
        _sys.path.insert(0, _parent)
    _rp.run_module(f"{_alias}.main_solve", run_name="__main__", alter_sys=True)
    raise SystemExit(0)

# Config / dataclasses
from .input_data import CorrelationCoefficients
# Physics
from .physics.combustion_chemistry.combustion_gas import combustion_gas_solve, choose_fuel
from .physics.combustion_chemistry.fpv_manifold import build_fpv_manifold, FPVManifold
from .physics.friction_correlations import getFrictionColebrook1939, getFrictionDeveloping, dispatch_friction_coil
from .physics.heat_transfer_correlations import dispatch_nu_coil, dispatch_nu_shell
from .physics.heat_conduction import OneDimensionalSteadyConduction_ShellnHelicalTube
from .physics.gas_flow.governing_equations import *
from .physics.liquid_flow.hx_adapters import solve_helical_coil_liquid_from_data_master
from .physics.liquid_flow.dispatch import BOILING_ONSET_BLEND_HALF_WIDTH, evaluate_coolant_closure
from .physics.liquid_flow.correlations import equilibrium_state_ph
from .physics.liquid_flow.regime import real_fluid_state_ph
from .physics.liquid_flow.sanity_checks import check_liquid_march
from .physics.radiation_model.radiation_build import make_ehlme_backend
from .physics.combustion_chemistry.other import print_mass_of_N_most_abundant_species
# Mechanical
from .mechanical.material_specs.material_temperature_strength import init_material_temperature_properties
from .mechanical.geometry.helix_geometry import *
from .mechanical.loads import *
# Data
from .model_data_process.data_processing import make_solver_data
from .model_data_process.data_plotting import HXDashboard
# Libraries
from CoolProp.CoolProp import PropsSI
from CoolProp import AbstractState
import numpy as np
from scipy.constants import gas_constant
from scipy.integrate import simpson
from pathlib import Path
import copy


class main_solver :

    def __init__ (  self,
                    coolantProp,
                    hotgasProp,
                    combustorProp,
                    numericalProp,
                    system_requirements,
                    corrCoeffs=None,
                    _liquid_enthalpy_hot_end_override=None):

        # main_solver is the shell-and-helical-tube backend and its axial-length
        # bookkeeping (_advance_state's L_HX update) branches on this exact
        # field: HX_config == "shellnHelicalTube" uses the true helical
        # arc-length-to-axial-position mapping, anything else silently falls
        # back to a naive linear dx accumulation that is wrong for a wound
        # coil. "shellntube" is shellntube_solver's own label, not a valid
        # main_solver configuration — constructing main_solver with it produces
        # silently-incorrect L_HX-based loop termination and mass/length
        # results, with no error. Fail fast instead.
        if combustorProp.HX_config != "shellnHelicalTube":
            raise ValueError(
                "main_solver requires combustorProp.HX_config == 'shellnHelicalTube' "
                f"for correct axial-length accounting; got {combustorProp.HX_config!r}. "
                "Use shellntube_solver (main_solve_shellntube.py) for 'shellntube'."
            )

        # extract dataclasses
        self.coolantProp = coolantProp
        self.hotgasProp = hotgasProp
        self.combustorProp = combustorProp
        self.numericalProp = numericalProp
        self.system_requirements = system_requirements
        self.corrCoeffs = corrCoeffs if corrCoeffs is not None else CorrelationCoefficients()

        # Coolant molar mass (gas-mode sound-speed/Mach calc only): looked up
        # from CoolProp for whatever fluid is configured, not a hardcoded
        # per-fluid constant — coolantProp.coolant may be any CoolProp fluid.
        self._coolant_molar_mass_g_mol = PropsSI('MOLAR_MASS', self.coolantProp.coolant) * 1000.0

        # fuel set up
        self.chem_mech_path, self.Y_fuel, self.Hv_fuel = choose_fuel(self.hotgasProp.fuel)

        """ 
        DEFINE HELICAL COIL GEOMETRY
        """
        self.N_ch = self.combustorProp.N_coils
        self.Dh_ch = self.combustorProp.Dh_coil 
        self.A_ch = np.pi*self.combustorProp.Dh_coil **2/4
        self.coil_pitch = self.Dh_ch + 2*self.combustorProp.thickness_coil_wall + self.combustorProp.coil_gap
        self.D_coil = self.combustorProp.inner_diameter - 2*self.combustorProp.gap_shell2coil - self.combustorProp.Dh_coil - 2*self.combustorProp.thickness_coil_wall
        self.D_inner_coil_passage = self.D_coil - self.combustorProp.Dh_coil - 2*self.combustorProp.thickness_coil_wall
        self.D_tube = self.Dh_ch + 2*self.combustorProp.thickness_coil_wall

        self.Rc = self.D_coil/2*(1 + (self.coil_pitch/(np.pi*self.D_coil))**2)  # coil curvature
        if self.D_coil < 0 : 
            raise Exception("Coil center-to-center diameter is negative - check your geometry")

        # coil axial length from pipe centers
        self.L_coil = (self.numericalProp.L_HX_max-self.combustorProp.mixing_length) - 2*self.combustorProp.length_2_coil - (self.Dh_ch+2*self.combustorProp.thickness_coil_wall)
        #* get functions to translate tube length to axial combustor length and get total tube length
        self.func_s_to_x, self.func_s_to_theta, self.L_ch_max = HelixGeometryRadiusCST(   coil_pitch=self.coil_pitch,
                                                                                D_coil=self.D_coil, 
                                                                                L_coil=self.L_coil)
        self.Dh_cc = compute_Dh_shell(D_coil=self.D_coil, 
                                d_coil_outer=self.Dh_ch+2*self.combustorProp.thickness_coil_wall, 
                                shell_diameter=self.combustorProp.inner_diameter, 
                                coil_pitch=self.coil_pitch)
        self.Ap_cc = np.pi*self.Dh_cc**2/4

        self.area_g_square_PT = self.combustorProp.inner_diameter**2 - np.pi*( (self.combustorProp.inner_diameter - 2*self.combustorProp.gap_shell2coil)**2 - (self.D_coil-self.combustorProp.Dh_coil - 2*self.combustorProp.thickness_coil_wall)**2)/4
        self.area_g_round_PT = np.pi*self.combustorProp.inner_diameter**2/4 - np.pi*( (self.combustorProp.inner_diameter - 2*self.combustorProp.gap_shell2coil)**2 - (self.D_coil-self.combustorProp.Dh_coil - 2*self.combustorProp.thickness_coil_wall)**2)/4
        self.relative_area_dif = (self.area_g_square_PT-self.area_g_round_PT)/self.area_g_round_PT
        print(f'relative diff square vs round PT={self.relative_area_dif}')

        #* HX pipe wall port area for mass computation
        self.Ap_HX = np.pi*((self.Dh_ch+2*self.combustorProp.thickness_coil_wall)**2 - self.Dh_ch**2)/4

        self.P_wc_unit = np.pi*self.combustorProp.Dh_coil 
        # Arc step = one full coil turn divided by N_arc_steps_per_turn.
        # N=1 (default): one full-turn per march step — dQ uses T_g_start for the whole turn.
        # N>1: finer arc discretisation; T_g updates N times per turn, reducing first-order
        #      Euler bias (~3% over-prediction of Q at N=1 for 27-turn HX at design point).
        # frac = 1.0 always correct in the 1D model — see doc/CALIBRATION_METHODOLOGY.md.
        self.numericalProp.dx = np.pi * self.D_coil / self.numericalProp.N_arc_steps_per_turn

        #! remove restriction on Helium inlet properties
        #* ensures respected geometry of HX

        """
        MODEL ---------------------------------- RADIATION BUILD
        """
        if self.numericalProp.radiation_ON == True : 
            #! RADIATIVE PARAMETERS FOR CO2-H2O MIXTURE 
            self.radiation_backend = make_ehlme_backend(Path(__file__).parent/"physics/radiation_model/ehlme2025_mixture.json")  # if only mixture is present
            #! mean beam length radiative heat
            # taking the approximate effective gas volume around each coil
            V_tot_1_turn = np.pi*self.combustorProp.inner_diameter**2/4 * self.numericalProp.L_HX_max 
            V_pipe = np.pi*(self.Dh_ch+2*self.combustorProp.thickness_coil_wall)**2/4*self.L_ch_max
            A_wg = np.pi*(self.Dh_ch+2*self.combustorProp.thickness_coil_wall)*self.L_ch_max
            self.Le = self.corrCoeffs.mbl_factor*(V_tot_1_turn-V_pipe)/A_wg
            #!----------------------------------
            print("Le=", self.Le, "V_tot_1_turn-V_pipe=", V_tot_1_turn-V_pipe)
            print(f"Pipe max length = {self.L_ch_max} m")

        else:
            self.radiation_backend = None
        
        
        #! using the shell hydraulic diameter to compute the gas bulk velocity

        #* material properties divided between conductive zone (HX) and burner (CC)
        self.func_CTE_HX, self.func_E_HX, self.func_Yield_HX, self.func_conductivity_HX, self.density_HX, self.poisson_HX, self.func_cp_HX = init_material_temperature_properties(self.combustorProp.material_HX)
        _, _, _, _, self.density_CC, _, _ = init_material_temperature_properties(self.combustorProp.material_CC)


        # initialize coolant — start from the end of the coil that coincides with the gas INLET
        # counter-flow: He exits at the gas-inlet end  → start at T_He_out, p_He_out; march sign = -1
        # co-flow:      He enters at the gas-inlet end → start at T_He_in,  p_He_in;  march sign = +1
        if self.combustorProp.flow_config == "co":
            self.T_c, self.p_c = self.coolantProp.T_in, self.coolantProp.p_in
            self._flow_sign = 1
        else:  # "counter" (default)
            self.T_c, self.p_c = self.coolantProp.T_out, self.coolantProp.p_out
            self._flow_sign = -1
        self.rho_c = PropsSI('D', 'T', self.T_c, 'P', self.p_c, self.coolantProp.coolant)
        self.U_c = self.coolantProp.mass_flow_c / (self.rho_c * self.A_ch * self.N_ch)

        # Liquid (boiling) coolant path: state is (p,h), never (T,p) — see
        # docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md Design Decision 1.
        # This is a first coupled-march wiring (Phase 2); horizontal-orientation
        # boiling-correction and geometry-specific (helical) enhancement are not
        # yet implemented — see physics/liquid_flow/correlations.py.
        self._liquid_mode = self.coolantProp.coolant_model == "equilibrium_liquid"
        if self._liquid_mode:
            self._liquid_min_pressure_Pa = 1.0
            # Worst (largest) Bergles-Rohsenow ONB wall-superheat margin seen
            # across the march (see evaluate_coolant_closure). Positive means
            # subcooled nucleate boiling was physically expected at some node
            # despite bulk quality < 0 - a real gap in the numerical
            # quality-only blend window, surfaced as a warning in
            # _check_global() rather than silently ignored.
            self._onb_max_margin = float("-inf")
            # Supercritical-regime diagnostics (populated only when the coolant
            # is above its critical pressure; see the supercritical branch of
            # evaluate_coolant_closure). Surfaced in _check_global().
            self._sc_closure_name = None       # which registry closure was used
            self._sc_regimes = set()           # regime labels seen along the march
            self._sc_htd_nodes = 0             # nodes flagged for HTD onset
            self._sc_extrapolated = False      # any node outside closure validity
            self._sc_extrap_message = None     # a representative extrapolation msg
            # Lower enthalpy safety floor, mirroring the pressure floor above.
            # Without it, a coarse arc-length step (or a shooting search probe
            # that starts too close to the physical inlet enthalpy) can drive
            # (p, h) past the fluid's valid CoolProp range in a single step,
            # crashing PropsSI deep inside the march instead of being caught
            # by the loop guard. The floor is generously below the physical
            # cold inlet enthalpy (never legitimately reached — hitting it
            # means the march has run away nonphysically) and is verified
            # evaluable at __init__ time, backing off geometrically if not.
            h_in_reference = PropsSI(
                'H', 'T', float(self.coolantProp.T_in), 'P', float(self.coolantProp.p_in),
                self.coolantProp.coolant,
            )
            margin = 5.0e5
            while margin > 1.0e3:
                candidate_floor = h_in_reference - margin
                try:
                    # real_fluid_state_ph so this succeeds for a supercritical
                    # inlet too (equilibrium_state_ph would raise, forcing the
                    # else-floor); subcritically identical to the old call.
                    real_fluid_state_ph(self.coolantProp.coolant, float(self.coolantProp.p_in), candidate_floor)
                    break
                except ValueError:
                    margin *= 0.5
            else:
                candidate_floor = h_in_reference - 1.0e3
            self._liquid_min_enthalpy_J_kg = candidate_floor
            if _liquid_enthalpy_hot_end_override is not None:
                # Internal hook for solve_counterflow_liquid_reference(): the
                # legacy T_out/p_out pair cannot express a two-phase starting
                # state (only single-phase (T,P) is well-posed), so the
                # shooting helper supplies the starting enthalpy directly
                # instead of going through PropsSI(T,p). Not a public
                # coolantProp field — this is solver march mechanics, not a
                # physical input.
                #
                # Pressure: coolantProp.p_out is NOT used here (it is
                # unrelated to this shooting mode and defaults to a value that
                # has nothing to do with p_in). Liquid friction pressure drop
                # over the coil is typically small relative to line pressure
                # (order 0.01-1 bar in the validated water cases), so the
                # hot-end starting pressure is approximated as p_in; the
                # actual cold-end pressure will differ from p_in by that
                # small friction drop. This is a documented approximation,
                # not a two-variable (h, p) shooting solve.
                self.enthalpy_c = float(_liquid_enthalpy_hot_end_override)
                self.p_c = float(self.coolantProp.p_in)
            else:
                self.enthalpy_c = PropsSI('H', 'T', self.T_c, 'P', self.p_c, self.coolantProp.coolant)
            self.q_w = 0.0  # seed for the lagged heat-flux term in the first node's boiling HTC
            initial_closure = evaluate_coolant_closure(
                coolant_prop=self.coolantProp,
                p_Pa=self.p_c,
                h_J_kg=self.enthalpy_c,
                mass_flux_kg_m2_s=self.coolantProp.mass_flow_c / (self.A_ch * self.N_ch),
                hydraulic_diameter_m=self.Dh_ch,
                heat_flux_W_m2=0.0,
                lut_path=self.coolantProp.liquid_chf_lut_path,
            )
            self.quality_c = initial_closure.state.quality
            self.void_c = initial_closure.state.void_fraction
            if initial_closure.onb_wall_superheat_margin_K is not None:
                self._onb_max_margin = max(self._onb_max_margin, initial_closure.onb_wall_superheat_margin_K)
            # self.rho_c (set generically above from T_c/p_c before the
            # equilibrium_liquid branch) is stale here for the
            # _liquid_enthalpy_hot_end_override path (counter-flow shooting):
            # it was flashed from the pre-override T_out/p_out, not the
            # actual starting (p, h) state. Refresh it so the first march
            # node's accelerational-pressure-drop term (which reads this as
            # its "previous node" density reference) starts from the correct
            # value instead of an unrelated one.
            self.rho_c = initial_closure.state.rho_kg_m3

            # Fail fast and clearly if the starting state is already past
            # complete vaporization / at the pressure or enthalpy floor —
            # otherwise the march loop guard (_coolant_flow_continues) simply
            # never executes a single node, data_master stays empty, and the
            # first downstream np.max() call crashes with an opaque
            # "zero-size array" error far from the real cause. This is the
            # single most common liquid-mode misconfiguration: T_out/p_out
            # (only meaningful for the single-phase-gas legacy march) left at
            # values that correspond to superheated vapor at the requested
            # pressure once coolant_model="equilibrium_liquid" is set.
            if not self._coolant_flow_continues():
                fluid = self.coolantProp.coolant
                if self.combustorProp.flow_config == "counter" and _liquid_enthalpy_hot_end_override is None:
                    which = "coolantProp.T_out / p_out"
                else:
                    which = "coolantProp.T_in / p_in"
                raise ValueError(
                    f"Liquid coolant march cannot start: the initial state "
                    f"(fluid={fluid!r}, T={self.T_c:.1f} K, p={self.p_c/1e5:.2f} bar) "
                    f"has quality={self.quality_c:.3f} (>=1.0 means already fully "
                    f"vaporized), p={self.p_c:.1f} Pa, h={self.enthalpy_c:.1f} J/kg. "
                    f"Check {which}: for coolant_model='equilibrium_liquid' these "
                    f"must correspond to a subcooled or saturated liquid state at "
                    f"the requested pressure, not a value tuned for a single-phase "
                    f"gas coolant. For counter-flow, consider "
                    f"solve_counterflow_liquid_reference() instead of guessing "
                    f"T_out/p_out - see docs/USER_GUIDE.md Section 0.5/0.6."
                )

        # initialize hot gas
        self.p_g = np.copy(self.hotgasProp.p0)
        self.combustion_node = combustion_gas_solve(fuel=self.hotgasProp.fuel, oxidizer=self.hotgasProp.oxidizer, OF=self.hotgasProp.mixing_ratio,
                                T_inj_LOX=self.hotgasProp.T_inj_LOX, T_g_init=self.hotgasProp.T_g_init, 
                                p0=self.p_g, 
                                chem_mech_path=self.chem_mech_path, Hv_fuel=self.Hv_fuel, Y_fuel=self.Y_fuel)
        self.combustion_node.solve() # combustion gas properties extractable using cantera
        self.gas_phase = self.combustion_node.phase
        self.T_g  = self.gas_phase.T
        self.rho_g = self.gas_phase.density
        self.W_g = self.gas_phase.mean_molecular_weight
        self._fpv = None
        self._fpv_h_removed = 0.0
        self._fpv_Yc = 0.0
        self._fpv_omega_Yc = 0.0
        if self.numericalProp.chemistry_model == "finite_rate":
            self._setup_fpv_manifold()
            self._update_fpv_gas_state()
        self.U_g = self.hotgasProp.mass_flow_g/(self.rho_g*self.Ap_cc)

        """ 
        COMBUSTION DATA
        """
        print(self.T_g, self.rho_g, self.gas_phase.cp, self.gas_phase.cp/self.gas_phase.cv, self.gas_phase.mean_molecular_weight, self.gas_phase.viscosity, self.gas_phase.thermal_conductivity)
        print("U_g = ", self.U_g)
        print_mass_of_N_most_abundant_species(self.gas_phase, 10)
        # initialize numerical parameters 
        self.L_HX = 0
        if self.combustorProp.HX_config=="shellnHelicalTube" :
            self.L_HX+=self.Dh_ch+2*self.combustorProp.thickness_coil_wall
        self.L_ch = 0
        self.T_wg, self.T_wc, self.T_c_check = np.copy(self.T_g), np.copy(self.T_c), np.copy(self.T_c) # initialize temperatures for 1D conduction

        # fresh data dict per solver instance — never shared across runs
        self.data_master = make_solver_data()

        #! species index for radiation, molar fractions
        self.index_H2O = self.gas_phase.species_index("H2O")
        self.index_CO2 = self.gas_phase.species_index("CO2")

        AS = AbstractState('HEOS',self.coolantProp.coolant)
        self.T_c_min_coolprop = AS.Tmin() # K

    def _setup_fpv_manifold(self):
        """Build the steady finite-rate FPV cooling manifold from the combustor inlet."""
        gas = self.gas_phase
        Tg0 = float(gas.T)
        pg0 = float(self.hotgasProp.p0)
        Yg0 = np.array(gas.Y, dtype=float)
        iCO = gas.species_index("CO")
        self._fpv = FPVManifold(build_fpv_manifold(
            gas,
            Y_inlet=Yg0,
            T_inlet=Tg0,
            p=pg0,
            species_index={
                "CO2": gas.species_index("CO2"),
                "H2O": gas.species_index("H2O"),
                "CO": iCO,
            },
            n_h=self.numericalProp.fpv_n_h,
            n_c=self.numericalProp.fpv_n_c,
            t_relax=self.numericalProp.fpv_t_relax,
            n_t=self.numericalProp.fpv_n_t,
            cache_dir=self.numericalProp.fpv_cache_dir,
        ))
        self._fpv_h_removed = 0.0
        self._fpv_Yc = self._fpv.Yc_inlet()

    def _update_fpv_gas_state(self):
        """Refresh gas properties from the FPV table at the current h/Yc state."""
        T, rho, mu, k, cp, xH2O, xCO2, omega = self._fpv.state(
            self._fpv_h_removed,
            self._fpv_Yc,
        )
        self.T_g = T
        self.rho_g = rho
        self.mu_g = mu
        self.k_g = k
        self.cp_g = cp
        r_specific = gas_constant * 1e3 / self.W_g
        self.cv_g = max(cp - r_specific, 1.0)
        self.gamma_g = self.cp_g / self.cv_g
        self.X_H2O = xH2O
        self.X_CO2 = xCO2
        self._fpv_omega_Yc = omega
    # ------------------------------------------------------------------
    # SANITY CHECKS
    # ------------------------------------------------------------------

    def _check_node(self):
        """Per-node physical sanity checks — only called when debug_verbose=True."""
        p = self.numericalProp
        tag = f"[x={self.L_HX:.3f}m]"
        if p.check_temperature_ordering:
            if self.T_c > self.T_wc + 0.5:
                print(f"WARNING {tag} T_c ({self.T_c:.1f} K) > T_wc ({self.T_wc:.1f} K) — wall colder than coolant")
            if self.T_wg > self.T_g + 0.5:
                print(f"WARNING {tag} T_wg ({self.T_wg:.1f} K) > T_g ({self.T_g:.1f} K) — wall hotter than gas")
        if p.check_mach_limits:
            if self.Mach_c > 0.5:
                print(f"WARNING {tag} Mach_c = {self.Mach_c:.3f} > 0.5 — compressibility significant in coolant")
            if self.Mach_g > 0.5:
                print(f"WARNING {tag} Mach_g = {self.Mach_g:.3f} > 0.5 — compressibility significant in hot gas")
        if p.check_Re_regime:
            if self.Re_c < 4000:
                print(f"WARNING {tag} Re_c = {self.Re_c:.0f} < 4000 — turbulent Nu correlation in laminar regime")
        if p.check_Z_deviation:
            if abs(self.Z - 1) > p.Z_tolerance:
                print(f"WARNING {tag} Z = {self.Z:.4f} (|Z-1| > {p.Z_tolerance}) — real-gas deviation; ideal-gas equations in use")
        if p.check_stress_limits:
            ratio = max(abs(self.stress_inner), abs(self.stress_outer)) / max(self.Yield, 1.0)
            if ratio > 0.8:
                print(f"WARNING {tag} stress/yield = {ratio:.2f} > 0.8")

    def _check_global(self):
        """Run-end energy balance and limit checks — called from HX_sizing_brief."""
        import numpy as np
        from CoolProp.CoolProp import PropsSI
        p = self.numericalProp

        if self._liquid_mode:
            # Gas-path checks below call PropsSI('H','T',...,'P',...), which is
            # ill-posed inside the two-phase dome; run the liquid-specific gates
            # instead (see docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md
            # Design Decision 6) and skip the gas-only checks entirely.
            self.liquid_sanity_report = check_liquid_march(
                self.data_master,
                fluid=self.coolantProp.coolant,
                mass_flow_c=self.coolantProp.mass_flow_c,
                energy_balance_tol=p.energy_balance_tol,
            )
            report = self.liquid_sanity_report
            print(
                f"Liquid coolant sanity gates: PASSED={report.passed} | "
                f"energy imbalance={report.energy_balance_rel_error*100:.2f}% | "
                f"min CHF margin={report.min_chf_margin:.2f} "
                f"({report.chf_regime_at_min_margin}) | "
                f"dryout_risk={report.dryout_risk} | "
                f"Mach_c_max={report.mach_c_max:.3f}"
            )
            for message in report.messages:
                print(f"  WARNING: {message}")
            if self._onb_max_margin > 0.0:
                print(
                    f"  WARNING: Bergles-Rohsenow ONB criterion exceeded in the "
                    f"subcooled region (max wall-superheat margin = "
                    f"{self._onb_max_margin:.2f} K) -- subcooled nucleate boiling "
                    f"is physically expected before the numerical blend window "
                    f"(bulk quality > -{BOILING_ONSET_BLEND_HALF_WIDTH:.3f}); wall "
                    f"temperature in that span may be under-predicted by the pure "
                    f"single-phase correlation."
                )
            if self._sc_closure_name is not None:  # supercritical run
                regimes = ", ".join(sorted(self._sc_regimes))
                print(
                    f"Supercritical coolant: closure={self._sc_closure_name} | "
                    f"regimes traversed: {regimes} | HTD-flagged nodes="
                    f"{self._sc_htd_nodes}"
                )
                if self._sc_extrapolated:
                    print(f"  WARNING: {self._sc_extrap_message} (validity extrapolation)")
                if self._sc_htd_nodes > 0:
                    print(
                        f"  WARNING: heat-transfer deterioration (HTD) risk flagged at "
                        f"{self._sc_htd_nodes} node(s) (McEligot-Jackson buoyancy/flow-"
                        f"acceleration criteria, Urbano & Nasuti 2013) -- the forced-"
                        f"convection property-ratio closure's negligible-buoyancy "
                        f"assumption is violated there; deteriorated-HTC magnitude is NOT "
                        f"modeled, treat these nodes as outside the validated envelope."
                    )
            return

        if p.check_energy_balance:
            Q_nodes = sum(self.data_master["dQ"])
            # Coolant enthalpy rise: integrate cp dT along the run
            T_c_arr = np.array(self.data_master["T_c"])
            p_c_arr = np.array(self.data_master["p_c"])
            h_c_out = PropsSI('H', 'T', T_c_arr[0],  'P', p_c_arr[0],  self.coolantProp.coolant)
            h_c_in  = PropsSI('H', 'T', T_c_arr[-1], 'P', p_c_arr[-1], self.coolantProp.coolant)
            Q_coolant = self.coolantProp.mass_flow_c * abs(h_c_out - h_c_in)
            imbalance = abs(Q_nodes - Q_coolant) / max(Q_nodes, 1.0)
            print(f"Energy balance: Q_nodes={Q_nodes/1e3:.2f} kW | Q_coolant={Q_coolant/1e3:.2f} kW | imbalance={imbalance*100:.1f}%")
            if imbalance > p.energy_balance_tol:
                print(f"  WARNING: energy imbalance {imbalance*100:.1f}% > {p.energy_balance_tol*100:.0f}% threshold")

        if p.check_temperature_ordering:
            T_wg = np.array(self.data_master["T_wg"])
            T_wc = np.array(self.data_master["T_wc"])
            T_c  = np.array(self.data_master["T_c"])
            T_g  = np.array(self.data_master["T_g"])
            n_bad_wc = np.sum(T_c > T_wc + 0.5)
            n_bad_wg = np.sum(T_wg > T_g + 0.5)
            if n_bad_wc > 0:
                print(f"WARNING: T_c > T_wc at {n_bad_wc} nodes")
            if n_bad_wg > 0:
                print(f"WARNING: T_wg > T_g at {n_bad_wg} nodes")

        if p.check_mach_limits:
            Mach_c_max = max(self.data_master["Mach_c"])
            Mach_g_max = max(self.data_master["Mach_g"])
            if Mach_c_max > 0.5:
                print(f"WARNING: Mach_c_max = {Mach_c_max:.3f} > 0.5")
            if Mach_g_max > 0.5:
                print(f"WARNING: Mach_g_max = {Mach_g_max:.3f} > 0.5")

        if p.check_stress_limits:
            ratio_max = np.max([
                np.max(np.abs(self.data_master["stress_inner"]) / np.array(self.data_master["Yield"])),
                np.max(np.abs(self.data_master["stress_outer"]) / np.array(self.data_master["Yield"]))
            ])
            if ratio_max > 0.8:
                print(f"WARNING: max stress/yield = {ratio_max:.2f} > 0.8 — approaching yield")

    def _coolant_flow_continues(self) -> bool:
        """Loop guard for the coolant march.

        Gas mode: stop if the coolant drops below the CoolProp safety floor.
        Liquid mode: stop on pressure floor or enthalpy floor (runaway
        nonphysical cooling) only. Complete vaporization (quality >= 1) is NOT
        a stop condition: evaluate_coolant_closure()'s single-phase branch
        already gives a physically correct closure past that point (real
        CoolProp vapor properties, a phase-agnostic Nusselt correlation), so
        the march continues into superheated vapor rather than silently
        ending at the dryout boundary. The coil-length limit in solver()'s
        while-loop condition remains the outer bound either way.
        """
        if self._liquid_mode:
            return (
                self.p_c > self._liquid_min_pressure_Pa
                and self.enthalpy_c > self._liquid_min_enthalpy_J_kg
            )
        return self.T_c > self.T_c_min_coolprop * self.numericalProp.T_c_safety_factor

    def solver (self) :

        # iterate through HX using numerical parameters
        #* stay in loop as long as :
            #! HX below max length
            #! coolant flow can still continue (see _coolant_flow_continues)
        while   self.L_HX  <= (self.numericalProp.L_HX_max - self.combustorProp.mixing_length - 2*self.combustorProp.length_2_coil) \
                and self._coolant_flow_continues():

            """ 
            Resolve heat transfer 
            """

            if self._liquid_mode:
                # (p,h)-state liquid/boiling closure — see Design Decision 1 in
                # docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md. The
                # boiling HTC's heat-flux dependence (Bo term) uses the PREVIOUS
                # node's wall flux (self.q_w before this node's conduction solve)
                # since the wall flux at THIS node is only known after solving
                # conduction with h_c as an input — a one-node-lagged closure,
                # first-order accurate and consistent as dx -> 0.
                # self.rho_c still holds the PREVIOUS node's density here (about
                # to be overwritten below) — captured for the accelerational
                # pressure-drop term, same one-node-lagged pattern.
                rho_c_prev = self.rho_c
                mass_flux_c = self.coolantProp.mass_flow_c / (self.A_ch * self.N_ch)
                closure = evaluate_coolant_closure(
                    coolant_prop=self.coolantProp,
                    p_Pa=self.p_c,
                    h_J_kg=self.enthalpy_c,
                    mass_flux_kg_m2_s=mass_flux_c,
                    hydraulic_diameter_m=self.Dh_ch,
                    heat_flux_W_m2=max(self.q_w, 0.0),
                    lut_path=self.coolantProp.liquid_chf_lut_path,
                    # self.T_wc holds the PREVIOUS node's coolant-side wall
                    # temperature (updated post-conduction each node; seeded at
                    # T_c) -- the lagged wall temp the supercritical property-
                    # ratio corrections need, same one-node-lag pattern as q_w.
                    # Ignored on the subcritical/boiling path (bit-identical).
                    # geometry: the helical coil is applied as a mild
                    # extrapolation of the straight-tube supercritical fits (see
                    # supercritical.py); orientation defaults to vertical (Cheng's
                    # basis) -- a documented approximation for a coil.
                    wall_temp_K=self.T_wc,
                    geometry="helical_coil",
                    # Flow length / diameter from the heated inlet, for Taylor's
                    # entrance-effect exponent (self.L_ch is arc length BEFORE
                    # this node, accumulated at the end of _advance_state()).
                    x_over_D=self.L_ch / self.Dh_ch if self.Dh_ch > 0.0 else None,
                )
                state = closure.state
                self.T_c = state.T_K
                self.quality_c = state.quality
                self.void_c = state.void_fraction
                self.rho_c = state.rho_kg_m3
                self.cp_c = state.cp_J_kg_K
                self.mu_c = state.mu_Pa_s
                self.k_c = state.k_W_m_K
                self.Pr_c = state.Pr
                self.h_c = closure.htc_W_m2_K
                self.c_c = closure.sound_speed_m_s if closure.sound_speed_m_s is not None else float('nan')
                self.chf_margin_c = closure.chf_margin if closure.chf_margin is not None else float('nan')
                self._dpdz_friction_c = closure.dpdz_friction_Pa_m
                if closure.onb_wall_superheat_margin_K is not None:
                    self._onb_max_margin = max(self._onb_max_margin, closure.onb_wall_superheat_margin_K)
                if closure.regime is not None:  # supercritical node
                    self._sc_closure_name = closure.closure_name
                    self._sc_regimes.add(closure.regime)
                    if closure.htd_risk:
                        self._sc_htd_nodes += 1
                    rep = closure.extrapolation_report
                    if rep is not None and not rep.in_range:
                        self._sc_extrapolated = True
                        if self._sc_extrap_message is None:
                            self._sc_extrap_message = rep.message()
                self.U_c = mass_flux_c / self.rho_c
                self.Re_c = mass_flux_c * self.Dh_ch / self.mu_c
                # Accelerational pressure drop (HEM): -dP/dz = friction + accel,
                # accel = G^2 * d(1/rho)/dz evaluated in the FLOW direction. Uses
                # the same one-node-lagged pattern as the boiling HTC above
                # (rho_c_prev is the previous node's density, dx is constant per
                # march step) rather than an implicit/iterative derivative -
                # first-order accurate and consistent as dx -> 0. See "Pressure
                # drop budget shifts" in
                # docs/solver_design/water_coolant_conversion_plan.md section 4:
                # this term was previously omitted as "small vs friction
                # pre-CHF" but can dominate as density collapses through the
                # two-phase dome and into superheated vapor.
                inv_rho_grad_index_order = (
                    1.0 / self.rho_c - 1.0 / rho_c_prev
                ) / self.numericalProp.dx
                inv_rho_grad_flow_order = self._flow_sign * inv_rho_grad_index_order
                self._dpdz_accel_c = mass_flux_c**2 * inv_rho_grad_flow_order
                # Gas-only diagnostics do not apply to the liquid path.
                self.Z = float('nan')
                self.gamma_c = float('nan')
                self.cv_c = float('nan')
                self.f_c = float('nan')
                self.Nu_c = float('nan')
                self.De = float('nan')
                self.He = float('nan')
            else:
                #* check for compressibility factor
                self.Z = PropsSI('Z','T',self.T_c,'P',self.p_c,self.coolantProp.coolant)
                # TODO(Z-correction): governing_equations.py uses ideal-gas formulations.
                # At He supercritical conditions (90 bar / 120 K) Z ~ 1.04-1.06: small but
                # non-zero. To include: multiply p = Z*rho*R*T in dU/dx and dp/dx derivatives.
                # start with initialized (U, p, T, rho)_c and T_g

                # Coolant thermodynamics + convection coef
                self.cp_c = PropsSI('C','T', self.T_c,'P',self.p_c,self.coolantProp.coolant)
                self.cv_c = PropsSI('CVMASS','T', self.T_c,'P',self.p_c,self.coolantProp.coolant)
                self.gamma_c = self.cp_c/self.cv_c
                self.mu_c = PropsSI('V','T',self.T_c,'P',self.p_c,self.coolantProp.coolant)
                self.k_c = PropsSI('L','T',self.T_c,'P',self.p_c,self.coolantProp.coolant)
                # flow characteristics
                self.Re_c=  self.rho_c*self.U_c*self.Dh_ch/self.mu_c
                self.De = self.Re_c*np.sqrt(self.Dh_ch/self.D_coil) # Dean number
                self.He = self.Re_c*np.sqrt(self.Dh_ch/(2*self.Rc)) # Helical number
                # Coil-side friction and Nusselt
                self.Pr_c = self.cp_c * self.mu_c / self.k_c
                self.f_c = dispatch_friction_coil(
                    self.combustorProp.friction_coil,
                    Re=self.Re_c, Dh=self.Dh_ch, Rc=self.D_coil / 2,
                    roughness=self.combustorProp.channel_roughness, x=10e10,
                    error_factor=1 + self.numericalProp.artificial_error_friction_cold,
                    corrCoeffs=self.corrCoeffs,
                )
                self.Nu_c = dispatch_nu_coil(
                    self.combustorProp.Nusselt_coil,
                    Re=self.Re_c, Pr=self.Pr_c, d=self.Dh_ch, R=self.D_coil / 2,
                    f_fd=self.f_c, x=10e10,
                    error_factor=1 + self.numericalProp.artificial_error_Nu_cold,
                    corrCoeffs=self.corrCoeffs,
                )
                self.h_c = self.Nu_c * self.k_c / self.Dh_ch

            # hot gas thermodynamics + convective effect
            if self._fpv is not None:
                self._update_fpv_gas_state()
            else:
                self.cp_g, self.cv_g = self.gas_phase.cp, self.gas_phase.cv
                self.gamma_g = self.cp_g/self.cv_g
                self.mu_g = self.gas_phase.viscosity 
                self.k_g = self.gas_phase.thermal_conductivity
                self.rho_g = self.gas_phase.density
                self.W_g = self.gas_phase.mean_molecular_weight
            # flow characteristics
            self.U_g = self.hotgasProp.mass_flow_g/(self.rho_g*self.Ap_cc)
            self.Re_g= self.rho_g*self.U_g*self.Dh_cc/self.mu_g
            self.f_fd_g = getFrictionColebrook1939(self.Re_g, self.combustorProp.combustor_roughness/self.Dh_cc)
            # check for fully rough region 
            self.f_g = getFrictionDeveloping(self.f_fd_g, self.Dh_cc, 10e10)
            # coolant heat transfer properties
            self.Pr_g = self.cp_g*self.mu_g/self.k_g

            #! radiation inputs
            if self._fpv is None:
                self.X_H2O = self.combustion_node.phase.X[self.index_H2O]
                self.X_CO2 = self.combustion_node.phase.X[self.index_CO2]


            self.Re_sh = self.rho_g*self.U_g*(self.Dh_ch+2*self.combustorProp.thickness_coil_wall)/self.mu_g

            # Shell-side Nusselt
            self.Nu_g, self.h_g = dispatch_nu_shell(
                self.combustorProp.Nusselt_shell,
                Re_sh=self.Re_sh, Re_g=self.Re_g, Pr_g=self.Pr_g, k_g=self.k_g,
                U_g=self.U_g, rho_g=self.rho_g, mu_g=self.mu_g,
                coil_pitch=self.coil_pitch, Dh_cc=self.Dh_cc,
                Dh_ch=self.Dh_ch, D_coil=self.D_coil,
                thickness_wall=self.combustorProp.thickness_coil_wall,
                T_bulk=float(self.T_g), T_wall=float(self.T_wg),
                Nusselt_correction=self.combustorProp.Nusselt_correction,
                error_factor=1 + self.numericalProp.artificial_error_Nu_hot,
                corrCoeffs=self.corrCoeffs,
            )
            


            if self._liquid_mode:
                # self.c_c was already set from the real-EOS/Wood's-equation
                # closure earlier this node (property-evaluation block) -
                # a real fluid sound speed, not an ideal-gas relation (which
                # would not apply to water/liquids).
                self.Mach_c = self.U_c / self.c_c
            else:
                self.c_c = np.sqrt(gas_constant*1e3/self._coolant_molar_mass_g_mol*self.T_c*self.cp_c/self.cv_c)
                self.Mach_c = self.U_c/self.c_c

            self.c_g = np.sqrt(gas_constant*1e3/self.W_g*self.T_g*self.cp_g/self.cv_g)
            self.Mach_g = self.U_g/self.c_g
            # 1D heat conduction 

            rad_state = None
            if self.numericalProp.radiation_ON:
                rad_state = {
                    'p': self.p_g,
                    'yH2O': self.X_H2O,
                    'yCO2': self.X_CO2,
                    'Le': self.Le,
                }

            self.heat_transfer_node =  OneDimensionalSteadyConduction_ShellnHelicalTube(
                                        h_g=self.h_g, h_c=self.h_c,
                                        T_c=self.T_c, T_g=self.T_g, 
                                        s_w=self.combustorProp.thickness_coil_wall,
                                        Dh_ch=self.Dh_ch,
                                        f_kw_at_T=self.func_conductivity_HX,
                                        T_wg_0=self.T_wg, T_wc_0=self.T_wc, T_c_check_0=self.T_c_check,
                                        dx=self.numericalProp.dx,
                                    # radiation:
                                        rad_enabled=self.numericalProp.radiation_ON,
                                        eps_s=self.corrCoeffs.emissivity_wall,
                                        rad_backend=self.radiation_backend,
                                        rad_state=rad_state)
            

            self.heat_transfer_node.Solve1Dconduction()
            self._unpack_heat_transfer_node()
            
            """ 
            STRESSES
            """
            #!extract functions
            #* properties at mean wall temperatures
            self.CTE = self.func_CTE_HX((self.T_wg+self.T_wc)/2-273)
            self.Modulus = self.func_E_HX((self.T_wg+self.T_wc)/2-273)
            # yield_at_hot_wall=True (default): evaluate yield at T_wg — the hot face is the
            # weakest point and governs failure; using mean wall temp is slightly non-conservative.
            if self.numericalProp.yield_at_hot_wall:
                self.Yield = self.func_Yield_HX(self.T_wg - 273)
            else:
                self.Yield = self.func_Yield_HX((self.T_wg + self.T_wc)/2 - 273)
            #* stresses
            self.stress_pressure = stress_pressure_tube(P=self.p_c, thickness_pipe=self.combustorProp.thickness_coil_wall, Dh_pipe=self.Dh_ch)
            self.stress_thermal_inner,  self.stress_thermal_outer= stress_thermal_tube(T_inner=self.T_wc, T_outer=self.T_wg,
                                                                                        CTE=self.CTE,
                                                                                        E=self.Modulus,
                                                                                        poisson=self.poisson_HX)
            self.stress_inner, self.stress_outer = [self.stress_thermal_inner+self.stress_pressure,  self.stress_thermal_outer+self.stress_pressure]

            # compute derivatives of coolant
            if self._liquid_mode:
                # (p,h) governing equations: dh/dx = dQ/(mdot);
                # -dp/dx = friction + acceleration (HEM). self._dpdz_accel_c
                # was computed earlier this node (property-evaluation block,
                # alongside the boiling closure) since it needs the previous
                # node's density, not this node's derivative outputs.
                self.dh_c__dx = (self.dq__dx / self.N_ch) / (self.coolantProp.mass_flow_c / self.N_ch)
                self.dp_c__dx_accel = -self._dpdz_accel_c
                self.dp_c__dx = -self._dpdz_friction_c + self.dp_c__dx_accel
                self.dU_c__dx = float('nan')
                self.drho_c__dx = float('nan')
                self.dT_c__dx = float('nan')
            else:
                #* ignore impact of pressure on temperature for now
                #! using power gradient dq__dx in order to conserve the effective perimeters (i.e. heat exchange widths)
                self.dU_c__dx = dU__dx_IdealGas(U=self.U_c, A=self.A_ch, p=self.p_c, P_w=1, q_w=self.dq__dx/self.N_ch, T=self.T_c, cp=self.cp_c, m_dot=self.coolantProp.mass_flow_c/self.N_ch, dA__dx=0, f=self.f_c, Dh=self.Dh_ch)
                self.dT_c__dx = dT__dx_IdealGas(P_w=1, q_w=self.dq__dx/self.N_ch, m_dot=self.coolantProp.mass_flow_c/self.N_ch, U=self.U_c, dU__dx=self.dU_c__dx, cp=self.cp_c)
                self.drho_c__dx = drho__dx_IdealGas_logical(rho=self.rho_c, U=self.U_c, dU__dx=self.dU_c__dx, A=self.A_ch, dA__dx=0)
                self.dp_c__dx = dp__dx_IdealGas_logical(p=self.p_c, T=self.T_c, dT__dx=self.dT_c__dx, rho=self.rho_c, drho__dx=self.drho_c__dx)
            # compute dp hot gas
            self.dp_g__dx = -self.f_g*self.rho_g*self.U_g**2/(2*self.Dh_cc)

            # per-node sanity checks (only when debug_verbose=True)
            if self.numericalProp.debug_verbose:
                self._check_node()

            # record data
            for key in self.data_master:
                if hasattr(self, key):
                    self.data_master[key].append(getattr(self, key))

            # advance all state variables by one dx step
            self._advance_state()
            


    
    def _unpack_heat_transfer_node(self) -> None:
        """Copy all outputs from the heat-conduction node onto self after Solve1Dconduction()."""
        n = self.heat_transfer_node
        self.UP = n.UP;  self.UA = n.UA
        self.dQ = n.dQ;  self.q_w = n.q_w;  self.dq__dx = n.dq__dx
        self.Res_g = n.Res_g;  self.Res_c = n.Res_c;  self.Res_w = n.Res_w
        self.Biot_g = n.Res_w / n.Res_g;  self.Biot_c = n.Res_w / n.Res_c
        self.q_w_rad = n.q_w_rad;  self.emissivity_g = n.eps_emit
        self.absorptivity_g = n.eps_abs;  self.h_g_rad = n.h_g_rad;  self.h_g_conv = n.h_g
        self.T_wg = n.T_wg_new;  self.T_wc = n.T_wc_new
        self.k_w = n.k_w;  self.T_c_check = n.T_c_check_f
        self.dh_g = self.dQ / self.hotgasProp.mass_flow_g

    def _advance_state(self) -> None:
        """Advance all state variables by one dx step after recording node data."""
        dx = self.numericalProp.dx
        # lengths
        self.L_ch += dx
        if self.combustorProp.HX_config in ("shellnHelicalTube", "coolingcoil"):
            self.L_HX = self.func_s_to_x(self.L_ch) + self.Dh_ch + 2 * self.combustorProp.thickness_coil_wall
        else:
            self.L_HX += dx
        # coolant
        sign = self._flow_sign
        if self._liquid_mode:
            self.enthalpy_c = max(
                self.enthalpy_c + self.dh_c__dx * dx * sign, self._liquid_min_enthalpy_J_kg
            )
            self.p_c = max(self.p_c + self.dp_c__dx * dx * sign, self._liquid_min_pressure_Pa)
            # Keep T_c/quality_c/void_c in sync with the just-advanced (p,h)
            # state. The next loop iteration's per-node closure overwrites
            # these anyway, but after the LAST iteration (loop exit) this is
            # the only place these get refreshed — without it, self.T_c stays
            # one node stale relative to self.enthalpy_c/p_c post-solve,
            # which is wrong for anything reading solver.T_c after solve()
            # returns (print_summary, the counterflow liquid shooting
            # residual, etc).
            # real_fluid_state_ph: dome-based below p_crit (identical to the
            # former equilibrium_state_ph call), single-phase real-EOS above it
            # (supercritical coolant does not crash; quality/void are NaN/0).
            eq_state = real_fluid_state_ph(self.coolantProp.coolant, self.p_c, self.enthalpy_c)
            self.T_c = eq_state.T_K
            self.quality_c = eq_state.quality
            self.void_c = eq_state.void_fraction
        else:
            self.T_c   += self.dT_c__dx   * dx * sign
            self.p_c   += self.dp_c__dx   * dx * sign
            self.rho_c += self.drho_c__dx * dx * sign
            self.U_c   += self.dU_c__dx   * dx * sign
        # gas
        self.p_g += self.dp_g__dx * dx
        if self._fpv is not None:
            self._fpv_h_removed += self.dh_g
            if self.U_g > 0:
                self._fpv_Yc += self._fpv_omega_Yc / self.U_g * dx
            self._update_fpv_gas_state()
        else:
            eq_on = self.numericalProp.chemistry_model == "equilibrium"
            if self.numericalProp.chemistry_model not in ("equilibrium", "frozen"):
                raise ValueError(
                    f"Unsupported steady chemistry_model: {self.numericalProp.chemistry_model!r}"
                )
            self.combustion_node.remove_energy(dh=self.dh_g, updated_pressure=self.p_g,
                                               equilibrium_dh_gas_ON=eq_on)
            g = self.combustion_node.phase
            self.T_g = g.T
            self.cp_g = g.cp;  self.cv_g = g.cv
            self.mu_g = g.viscosity;  self.k_g = g.thermal_conductivity;  self.rho_g = g.density

    def compute_performance(self):
        """Compute all post-run scalar metrics. Stores results as attributes and returns a dict."""
        d = self.data_master

        # HX effectiveness
        self.C_g_avg = np.average(np.array(d["cp_g"]) * self.hotgasProp.mass_flow_g)
        self.C_c_avg = np.average(np.array(d["cp_c"]) * self.coolantProp.mass_flow_c)
        self.T_wg_max = np.max(d["T_wg"])
        self.T_wc_max = np.max(d["T_wc"])
        self.dT_max   = np.max(np.array(d["T_wg"]) - np.array(d["T_wc"]))
        # He cold inlet is at data[-1] (counter-flow) or data[0] (co-flow)
        T_c_cold_inlet = d["T_c"][0] if self._flow_sign == 1 else d["T_c"][-1]
        self.eta_HX   = (self.C_g_avg * (d["T_g"][0] - d["T_g"][-1])
                         / (min(self.C_c_avg, self.C_g_avg) * (d["T_g"][0] - T_c_cold_inlet)))
        self.Mach_g_max = np.max(d["Mach_g"])
        self.Q_tot    = sum(d["dQ"]) * 1e-3   # kW
        if self._liquid_mode:
            # simpson(cp, T) is invalid through a phase change (T plateaus at
            # Tsat while enthalpy keeps rising); use the enthalpy difference
            # directly instead.
            self.Q_He = np.abs(
                self.coolantProp.mass_flow_c * (d["enthalpy_c"][-1] - d["enthalpy_c"][0])
            ) * 1e-3  # kW
        else:
            self.Q_He     = np.abs(self.coolantProp.mass_flow_c
                                   * simpson(np.array(d["cp_c"]), x=np.array(d["T_c"]))) * 1e-3  # kW
        self.T_fin_max = np.max(d["T_g"])
        # p_He_inlet - p_He_outlet > 0 regardless of flow direction
        # counter: data[0]=He_outlet, data[-1]=He_inlet  → [-1]-[0] > 0
        # co-flow: data[0]=He_inlet,  data[-1]=He_outlet → [0]-[-1] > 0
        self.dp_c_tot = abs(d["p_c"][-1] - d["p_c"][0])

        # Propulsion (Humble 1995, Space Propulsion Analysis and Design)
        self.R_g   = gas_constant * 1e3 / self.W_g
        self.Gam_g = self.gamma_g * (2 / (self.gamma_g + 1)) ** ((self.gamma_g + 1) / (2 * (self.gamma_g - 1)))
        self.c0    = np.sqrt(self.gamma_g * self.R_g * self.T_g)
        self.exit_area = np.pi * self.combustorProp.exhaust_diameter ** 2 / 4
        self.p0    = self.c0 * self.hotgasProp.mass_flow_g / (self.Gam_g * self.exit_area)
        self.T_exit = self.T_g / (1 + (self.gamma_g - 1) / 2)
        self.p_exit = self.p0 / (1 + (self.gamma_g - 1) / 2) ** (1 / (self.gamma_g - 1))
        self.v_exit = np.sqrt(2 * self.gamma_g * self.R_g * self.T_g / (self.gamma_g - 1)
                              * (1 - (self.p_exit / self.p0) ** ((self.gamma_g - 1) / self.gamma_g)))
        self.Thrust_SL  = self.v_exit * self.hotgasProp.mass_flow_g + (self.p_exit - 101325) * self.exit_area
        self.Thrust_VAC = self.v_exit * self.hotgasProp.mass_flow_g + self.p_exit * self.exit_area

        # Mass budget
        self.F = 1 / (self.hotgasProp.mixing_ratio + 1)
        self.O = 1 - self.F
        self.mass_kerosene = self.system_requirements.burn_time * self.hotgasProp.mass_flow_g * self.F
        self.mass_LOx      = self.system_requirements.burn_time * self.hotgasProp.mass_flow_g * self.O
        self.mass_injector_plate = (self.density_CC
                                    * np.pi * self.combustorProp.inner_diameter ** 2 / 4
                                    * self.combustorProp.wall_thickness_inj)
        self.mass_gas_eject = (self.density_CC
                               * 1.1 * np.pi * self.combustorProp.inner_diameter ** 2 / 4
                               * self.combustorProp.wall_thickness_cc)
        self.mass_mixing_zone = (self.density_CC * self.combustorProp.mixing_length
                                 * np.pi * ((self.combustorProp.inner_diameter + 2 * self.combustorProp.wall_thickness_cc) ** 2
                                            - self.combustorProp.inner_diameter ** 2) / 4)
        self.N_turns = self.L_coil / self.coil_pitch
        self.mass_HX = (self.L_ch * self.density_HX
                        * np.pi * ((self.Dh_ch + 2 * self.combustorProp.thickness_coil_wall) ** 2
                                   - self.Dh_ch ** 2) / 4)
        self.mass_shell_walls = (self.density_CC * (self.L_HX - self.combustorProp.mixing_length)
                                 * np.pi * ((self.combustorProp.inner_diameter + 2 * self.combustorProp.wall_thickness_cc) ** 2
                                            - self.combustorProp.inner_diameter ** 2) / 4)
        self.mass_combustor = (self.mass_injector_plate + self.mass_gas_eject
                               + self.mass_mixing_zone + self.mass_shell_walls)
        self.mass_tot = self.mass_kerosene + self.mass_LOx + self.mass_HX + self.mass_combustor

        # Stress
        self.max_stress__yield = np.max([
            np.max(np.abs(d["stress_inner"]) / np.array(d["Yield"])),
            np.max(np.abs(d["stress_outer"]) / np.array(d["Yield"]))
        ])

        return dict(
            eta_HX=self.eta_HX, Q_tot_kW=self.Q_tot, Q_He_kW=self.Q_He,
            T_wg_max=self.T_wg_max, T_wc_max=self.T_wc_max, dT_wall_max=self.dT_max,
            Mach_g_max=self.Mach_g_max, dp_c_bar=self.dp_c_tot / 1e5,
            N_turns=self.N_turns, L_HX=self.L_HX, L_pipe=self.L_ch,
            Thrust_SL_N=self.Thrust_SL, Thrust_VAC_N=self.Thrust_VAC,
            p0_bar=self.p0 / 1e5, v_exit=self.v_exit,
            mass_HX_kg=self.mass_HX, mass_combustor_kg=self.mass_combustor,
            mass_tot_kg=self.mass_tot, max_stress__yield=self.max_stress__yield,
        )

    def print_summary(self):
        """Print formatted summary of all performance metrics. Call compute_performance() first."""
        d = self.data_master
        print()
        print("=" * 55)
        print("RESULTS")
        print("=" * 55)
        print(f"  Config:       {self.combustorProp.HX_config}  |  N_coils = {self.combustorProp.N_coils}")
        print(f"  L_pipe = {self.L_ch:.3f} m  |  L_HX = {self.L_HX:.3f} m")
        print(f"  D_coil = {self.D_coil*1e3:.1f} mm  |  N_turns = {self.N_turns:.1f}")
        print(f"  D_coil/Dh = {self.D_coil/self.Dh_ch:.2f}  |  Rc = {self.Rc*1e3:.1f} mm")
        print()
        print("  Thermal performance")
        print(f"    eta_HX   = {self.eta_HX:.3f}")
        print(f"    Q_tot    = {self.Q_tot:.2f} kW  |  Q_He (Simpson) = {self.Q_He:.2f} kW")
        print(f"    T_wg_max = {self.T_wg_max:.1f} K  |  T_wc_max = {self.T_wc_max:.1f} K")
        print(f"    dT_wall_max = {self.dT_max:.1f} K")
        print(f"    T_g: {np.max(d['T_g']):.1f} K -> {np.min(d['T_g']):.1f} K")
        print(f"    T_c_in = {self.T_c:.1f} K  |  p_c_in = {self.p_c/1e5:.2f} bar")
        print()
        print("  Flow")
        print(f"    Re_c_max = {np.max(d['Re_c'])/1e6:.2f}e6")
        print(f"    De_max = {np.max(d['De']):.0f}  |  De_min = {np.min(d['De']):.0f}")
        print(f"    He_max = {np.max(d['He']):.0f}  |  He_min = {np.min(d['He']):.0f}")
        print(f"    Nu_c: {np.max(d['Nu_c']):.0f} -> {np.min(d['Nu_c']):.0f}")
        print(f"    Nu_g: {np.max(d['Nu_g']):.0f} -> {np.min(d['Nu_g']):.0f}")
        print(f"    Biot_g: {np.max(d['Biot_g']):.3f} -> {np.min(d['Biot_g']):.3f}")
        print(f"    Biot_c: {np.max(d['Biot_c']):.3f} -> {np.min(d['Biot_c']):.3f}")
        print(f"    Mach_g_max = {self.Mach_g_max:.4f}")
        print(f"    dP_c = {self.dp_c_tot/1e5:.3f} bar")
        print()
        print("  Geometry")
        print(f"    D_inner_coil_passage = {self.D_inner_coil_passage*1e3:.2f} mm")
        print(f"    D_tube               = {self.D_tube*1e3:.2f} mm")
        print(f"    D_inner/D_tube       = {self.D_inner_coil_passage/self.D_tube:.3f}")
        print()
        print("  Mass budget")
        print(f"    HX pipe      = {self.mass_HX:.3f} kg")
        print(f"    Combustor    = {self.mass_combustor:.3f} kg")
        print(f"    Fuel         = {self.mass_kerosene:.2f} kg")
        print(f"    LOX          = {self.mass_LOx:.2f} kg")
        print(f"    TOTAL        = {self.mass_tot:.2f} kg")
        print()
        print("  Propulsion")
        print(f"    p0 = {self.p0/1e5:.2f} bar  |  v_exit = {self.v_exit:.1f} m/s")
        print(f"    T_exit = {self.T_exit:.1f} K  |  p_exit = {self.p_exit/1e5:.3f} bar")
        print(f"    Thrust SL = {self.Thrust_SL:.2f} N  |  Thrust VAC = {self.Thrust_VAC:.2f} N")
        print()
        print("  Stresses")
        print(f"    max |sigma_inner| = {np.max(np.abs(d['stress_inner']))/1e6:.1f} MPa")
        print(f"    max |sigma_outer| = {np.max(np.abs(d['stress_outer']))/1e6:.1f} MPa")
        print(f"    max sigma/yield   = {self.max_stress__yield:.3f}")
        print("=" * 55)

    def HX_sizing_brief(self, plotON=True, printON=True):
        """Convenience wrapper: compute metrics, run checks, print, and plot."""
        self._check_global()
        self.compute_performance()
        if printON:
            self.print_summary()
        if plotON:
            HXDashboard(self.data_master).all()

    def liquid_coolant_postprocess(self, lut_path=None, min_pressure_Pa=1.0):
        """Map the completed helical steady duty profile into the liquid coolant solver.

        This is an opt-in integration bridge. It does not alter the maintained
        helium march; it consumes the already-recorded ``data_master["dQ"]`` and
        returns liquid p-h fields/diagnostics in the helical solver's axial order.
        """
        result = solve_helical_coil_liquid_from_data_master(
            coolant_prop=self.coolantProp,
            combustor_prop=self.combustorProp,
            numerical_prop=self.numericalProp,
            data_master=self.data_master,
            lut_path=lut_path,
            min_pressure_Pa=min_pressure_Pa,
        )
        self.liquid_coolant = result
        return result


def solve_counterflow_physical_reference(
    coolantProp,
    hotgasProp,
    combustorProp,
    numericalProp,
    system_requirements,
    corrCoeffs=None,
):
    """Shoot the helical counter-flow outlet temperature to match T_cold,in.

    The legacy counter-flow steady march starts at the gas-inlet end, where the
    helium exits. For a physical comparison with transient counter-flow, the
    known boundary is the cold helium inlet at the gas-outlet end. This helper
    repeatedly runs the existing march and adjusts `coolantProp.T_out`.
    """
    if combustorProp.flow_config != "counter":
        solver = main_solver(coolantProp, hotgasProp, combustorProp, numericalProp,
                             system_requirements, corrCoeffs=corrCoeffs)
        solver.solver()
        return solver

    target = float(coolantProp.T_in)
    max_iter = int(getattr(numericalProp, "counterflow_reference_max_iter", 4))
    tol = float(getattr(numericalProp, "counterflow_reference_tol_K", 1.0))

    def _run(T_hot_end):
        coolant = copy.copy(coolantProp)
        coolant.T_out = float(T_hot_end)
        solver = main_solver(coolant, hotgasProp, combustorProp, numericalProp,
                             system_requirements, corrCoeffs=corrCoeffs)
        solver.solver()
        residual = float(solver.data_master["T_c"][-1]) - target
        return solver, residual

    T_user = float(getattr(coolantProp, "T_out", target + 100.0))
    guesses = [max(target + 1.0, T_user), target + 0.5 * max(T_user - target, 20.0)]
    history = []
    best_solver = None
    best_residual = float("inf")

    for guess in guesses:
        solver, residual = _run(guess)
        history.append((float(guess), float(residual)))
        if abs(residual) < abs(best_residual):
            best_solver, best_residual = solver, residual
        if abs(residual) <= tol:
            best_solver.counterflow_reference_residual_K = residual
            best_solver.counterflow_reference_T_out_guess_K = guess
            return best_solver

    for _ in range(max(0, max_iter - len(history))):
        (x0, f0), (x1, f1) = history[-2], history[-1]
        if abs(f1 - f0) < 1e-9:
            x2 = x1 + 25.0
        else:
            x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x2 = float(np.clip(x2, target, 2000.0))
        solver, residual = _run(x2)
        history.append((x2, float(residual)))
        if abs(residual) < abs(best_residual):
            best_solver, best_residual = solver, residual
        if abs(residual) <= tol:
            best_solver.counterflow_reference_residual_K = residual
            best_solver.counterflow_reference_T_out_guess_K = x2
            return best_solver

    best_solver.counterflow_reference_residual_K = best_residual
    best_solver.counterflow_reference_T_out_guess_K = min(history, key=lambda pair: abs(pair[1]))[0]
    return best_solver


def solve_counterflow_liquid_reference(
    coolantProp,
    hotgasProp,
    combustorProp,
    numericalProp,
    system_requirements,
    corrCoeffs=None,
):
    """Shoot the helical counter-flow liquid march's hot-end starting enthalpy
    to match the user's physical cold coolant inlet (T_in, p_in).

    The plain liquid counter-flow march (``main_solver`` with
    ``coolant_model == "equilibrium_liquid"``) starts from
    ``coolantProp.T_out``/``p_out`` as a single-phase (T,P) guess — the same
    legacy shortcut the gas march uses. For liquids this is a worse interface
    than for gas: a (T,P) pair cannot express a genuine two-phase starting
    state at all (T and P are not independent inside the dome), so the user
    would have to guess a single-phase hot-end state by trial and error. This
    helper removes that guess: it shoots on the hot-end starting ENTHALPY
    directly (never temperature — invalid inside the dome, see Design
    Decision 1 in docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md) so
    that the march's cold end matches the physically supplied
    ``coolantProp.T_in``/``p_in`` inlet. ``coolantProp.T_out``/``p_out`` are
    not used by this path at all.

    Pressure caveat: only enthalpy is shot; the hot-end starting pressure is
    approximated as ``p_in`` (liquid friction pressure drop over the coil is
    typically small relative to line pressure in the validated water cases —
    order 0.01-1 bar). The converged cold-end pressure therefore differs from
    ``p_in`` by that friction drop, which can shift the converged cold-end
    temperature by more than the raw energy tolerance would suggest,
    especially near saturation where ``T = Tsat(p)`` is pressure-sensitive.
    This is a documented one-variable-shooting approximation, not a coupled
    (h, p) 2D root-find; `solver.counterflow_reference_residual_J_kg` reports
    the true (enthalpy-space) convergence achieved.
    """
    if combustorProp.flow_config != "counter" or coolantProp.coolant_model != "equilibrium_liquid":
        solver = main_solver(coolantProp, hotgasProp, combustorProp, numericalProp,
                             system_requirements, corrCoeffs=corrCoeffs)
        solver.solver()
        return solver

    fluid = coolantProp.coolant
    target_h = PropsSI('H', 'T', float(coolantProp.T_in), 'P', float(coolantProp.p_in), fluid)
    # Bisection halves the bracket each iteration (linear convergence, unlike
    # the superlinear secant method used for the gas reference), so it needs
    # more iterations for the same tolerance: with rough_span ~ cp*150K and
    # tol_h ~2000 J/kg, about 8-10 bisections are needed after the initial
    # bracket. Default overridden here (not shared with the gas path's
    # counterflow_reference_max_iter=6, which is tuned for secant).
    max_iter = int(getattr(numericalProp, "counterflow_liquid_reference_max_iter", 20))
    tol_h = float(getattr(numericalProp, "counterflow_reference_tol_J_kg", 2.0e3))

    def _run(h_hot_end):
        solver = main_solver(coolantProp, hotgasProp, combustorProp, numericalProp,
                             system_requirements, corrCoeffs=corrCoeffs,
                             _liquid_enthalpy_hot_end_override=h_hot_end)
        solver.solver()
        # solver.enthalpy_c (the live attribute) holds the true post-march
        # state after the last _advance_state() call. data_master["enthalpy_c"]
        # only records PRE-advance states (one node "behind" the true march
        # endpoint — see the energy-balance note in
        # physics/liquid_flow/sanity_checks.py), so it must not be used here.
        residual = float(solver.enthalpy_c) - target_h
        return solver, residual

    # Bracketed bisection, not secant/Newton: h_cold_end(h_hot_end) can be
    # non-smooth across boiling onset/CHF (see "Known Risks" in
    # docs/solver_design/LIQUID_COOLANT_INTEGRATION_PLAN.md), and a secant
    # step can overshoot into a nonphysical enthalpy guess there (observed:
    # a wild low-enthalpy guess crashed CoolProp mid-march). Bisection keeps
    # every guess strictly inside a bracket with confirmed opposite-sign
    # residuals, so it cannot overshoot. residual(h_hot_end) was confirmed
    # monotonic increasing empirically (tests/test_liquid_counterflow_reference.py).
    #
    # The required margin above target_h varies over orders of magnitude with
    # mass flow (a high-flow case may need only ~2 K worth of margin; a
    # low-flow boiling case may need >150 K worth) — a single fixed-size
    # bracket cannot cover both, so the bracket itself is found adaptively:
    # start from a small margin and geometrically expand (if still short of
    # target) or shrink (if already past it) until the sign of the residual
    # flips, then bisect between the last two opposite-sign points.
    cp_ref = PropsSI('C', 'T', float(coolantProp.T_in), 'P', float(coolantProp.p_in), fluid)
    initial_margin = max(cp_ref * 2.0, 500.0)  # ~2 K worth of margin, a physically small starting probe
    max_bracket_attempts = 20

    lo = target_h + initial_margin
    solver_lo, f_lo = _run(lo)
    best_solver, best_residual = solver_lo, f_lo
    if abs(f_lo) <= tol_h:
        solver_lo.counterflow_reference_residual_J_kg = f_lo
        solver_lo.counterflow_reference_h_hot_end_guess_J_kg = lo
        return solver_lo

    hi, f_hi = lo, f_lo
    bracketed = False
    if f_lo < 0.0:
        # Not enough margin yet — expand hi upward until the residual turns positive.
        offset = initial_margin
        for _ in range(max_bracket_attempts):
            offset *= 2.0
            hi = target_h + offset
            solver_hi, f_hi = _run(hi)
            if abs(f_hi) < abs(best_residual):
                best_solver, best_residual = solver_hi, f_hi
            if abs(f_hi) <= tol_h:
                solver_hi.counterflow_reference_residual_J_kg = f_hi
                solver_hi.counterflow_reference_h_hot_end_guess_J_kg = hi
                return solver_hi
            if f_hi > 0.0:
                bracketed = True
                break
    else:
        # Already past target even at the smallest sane margin — shrink lo
        # toward target_h until the residual turns negative.
        hi, f_hi = lo, f_lo
        offset = initial_margin
        for _ in range(max_bracket_attempts):
            offset *= 0.5
            lo = target_h + offset
            solver_lo, f_lo = _run(lo)
            if abs(f_lo) < abs(best_residual):
                best_solver, best_residual = solver_lo, f_lo
            if abs(f_lo) <= tol_h:
                solver_lo.counterflow_reference_residual_J_kg = f_lo
                solver_lo.counterflow_reference_h_hot_end_guess_J_kg = lo
                return solver_lo
            if f_lo < 0.0:
                bracketed = True
                break

    if not bracketed:
        # Could not find a sign change within the attempt budget; return the
        # closest evaluated guess rather than bisecting a non-bracket.
        best_solver.counterflow_reference_residual_J_kg = best_residual
        best_solver.counterflow_reference_h_hot_end_guess_J_kg = (
            lo if abs(f_lo) < abs(f_hi) else hi
        )
        return best_solver

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        solver_mid, f_mid = _run(mid)
        if abs(f_mid) < abs(best_residual):
            best_solver, best_residual = solver_mid, f_mid
        if abs(f_mid) <= tol_h:
            solver_mid.counterflow_reference_residual_J_kg = f_mid
            solver_mid.counterflow_reference_h_hot_end_guess_J_kg = mid
            return solver_mid
        if (f_mid > 0.0) == (f_hi > 0.0):
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid

    best_solver.counterflow_reference_residual_J_kg = best_residual
    best_solver.counterflow_reference_h_hot_end_guess_J_kg = 0.5 * (lo + hi)
    return best_solver


if __name__ == "__main__":
    from .input_data import coolantProp, hotgasProp, combustorProp, numericalProp, system_requirements

    combustor = main_solver(hotgasProp=hotgasProp(),
                            coolantProp=coolantProp(),
                            combustorProp=combustorProp(),
                            system_requirements=system_requirements(),
                            numericalProp=numericalProp())

    combustor.solver()
    combustor.compute_performance()
    combustor.print_summary()
    combustor.HX_sizing_brief()
