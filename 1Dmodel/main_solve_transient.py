"""
Dynamic (transient) heat-exchanger solver — start-up / shut-down.

@ author : Raphaël Aubry  (transient extension)

Architecture (see doc/DESIGN_PLAN_shellntube_transient.md section 4):
  - Legacy `fluid_model = "quasi_steady"` integrates only the lumped
    thickness-mean wall temperature per axial node, T_bar_i(t).
  - Production `fluid_model = "transient_coolant"` integrates wall temperature
    plus conserved helium mass/internal energy in a 1D finite-volume coolant
    model. The coolant face fluxes can be prescribed/quasi-steady or advanced
    with a low-Mach pressure-driven momentum closure.
  - Hot gas remains quasi-steady over its sub-ms residence time: one spatial
    march from inlet to outlet at the current wall/coolant state gives the
    per-node hot-side heat fluxes. Chemistry defaults to finite-rate FPV.
  - Wall ODE:  (rho_w cp_w A_wall) dT_bar_i/dt = dq_hot__dx_i - dq_cold__dx_i,
    integrated with scipy.solve_ivp (adaptive).

This module reuses `main_solver` wholesale for geometry, radiation, combustion
and material-property setup; it only replaces the steady wall `fsolve` with the
transient `fluxes_at_Tbar` reconstruction and adds the time integration.

Config-agnostic by construction (section 4.5): `fluid_pass` is implemented here
for the shell-and-helical-tube config; the shell-and-tube config provides its
own `fluid_pass` with the identical signature (WP2).

Run directly:  python main_solve_transient.py
"""

# --- self-contained bootstrap (folder name starts with a digit), mirrors main_solve.py ---
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
    _rp.run_module(f"{_alias}.main_solve_transient", run_name="__main__", alter_sys=True)
    raise SystemExit(0)

import numpy as np
from types import SimpleNamespace
from scipy.integrate import solve_ivp

from .main_solve import main_solver, solve_counterflow_physical_reference
from .transient_core.compressible_coolant import (
    coolprop_state_from_mass_energy,
    enforce_density_bounds,
    enforce_internal_energy_floor,
    initial_mass_energy_from_TP,
)
from .transient_core.integrator import fixed_time_grid
from .transient_core.progress import TransientProgressPrinter
from .transient_core.adapters_shelltube import (
    _coolant_mass_energy_from_TP_profile,
    _implicit_quadratic_momentum_update,
    _limit_face_mdot_for_inventory,
)
from .transient_core.wall_compressible_coolant import (
    semi_implicit_wall_compressible_coolant_step,
)
from .physics.heat_conduction import OneDimensionalSteadyConduction_ShellnHelicalTube
from .physics.friction_correlations import getFrictionColebrook1939, getFrictionDeveloping, dispatch_friction_coil
from .physics.heat_transfer_correlations import dispatch_nu_coil, dispatch_nu_shell
from .physics.radiation_model.radiation_equations import qrad_net_mbl, hrad_from_q
from .physics.combustion_chemistry.fpv_manifold import build_fpv_manifold, FPVManifold


def _interp_schedule(schedule, t, default):
    """Linear interpolation of a list[(t, value)] schedule, flat-held outside range.
    `schedule is None` -> return `default` (hold the steady input)."""
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


class transient_solver(main_solver):
    """Dynamic solver for the shell-and-helical-tube config.

    Composes `main_solver` (all steady setup) and adds the lumped-wall transient.
    Supports both `combustorProp.flow_config in ("co", "counter")`.

    Co-flow: both fluids enter at x=0 — a well-posed initial-value problem,
    solved by a single forward march (`_march_fluids`).

    Counter-flow: gas enters at x=0, coolant enters at x=L. Even though the wall
    field `Tbar_vec` is fixed (given, not solved for, within one RHS evaluation),
    the two fluid temperatures are still mutually dependent along the whole
    length (the per-node wall flux needs both simultaneously) — a genuine
    two-point BVP in the fluid temperatures. Resolved by a 1-D shooting method
    (`_march_fluids` + `brentq` on the coolant's x=0 guess vs. its true inlet at
    x=L) — cheaper than the steady solver's full sweep because the wall does
    NOT need co-solving here, only the two fluid profiles. Warm-started from the
    previous RHS call's converged guess since `Tbar` changes little between
    adjacent `solve_ivp` evaluations.
    """

    def __init__(self, coolantProp, hotgasProp, combustorProp, numericalProp,
                 system_requirements, transientProp, corrCoeffs=None):
        super().__init__(coolantProp, hotgasProp, combustorProp, numericalProp,
                         system_requirements, corrCoeffs=corrCoeffs)
        self.transientProp = transientProp
        if self.combustorProp.flow_config not in ("co", "counter"):
            raise ValueError(f"flow_config must be 'co' or 'counter', got "
                             f"{self.combustorProp.flow_config!r}")
        self._counter_Tc_profile = None  # warm-start cache for the counter-flow relaxation
        # Below this coolant flow, skip the counter-flow relaxation and march co-flow
        # style: the flow is so low that quasi-steady is already invalid (He residence
        # >> transient timescale — section 4.6 flags the outlet unreliable there anyway),
        # AND the cold-start wall field makes the relaxation ill-behaved. NOT an NTU
        # criterion (NTU ~ mdot^-0.2 here, nearly flow-independent, so it can't
        # discriminate the low-flow regime — an earlier NTU-threshold attempt wrongly
        # disabled counter-flow at the OPERATING point where co/counter differ by ~160 K).
        self._counter_fallback_mdot = 0.02 * max(self.coolantProp.mass_flow_c, 1e-9)

        # fixed axial grid + per-node geometry, captured from one steady march
        self._build_axial_grid()

        # wall thermal-mass geometry (uniform annulus, per node)
        Dh, s = self.Dh_ch, self.combustorProp.thickness_coil_wall
        self.A_wall = np.pi * ((Dh + 2 * s) ** 2 - Dh ** 2) / 4.0   # [m^2] tube cross-section
        self.rho_w = self.density_HX

        # cache the inlet burnt-gas state
        self._gas_inlet_TPY = (float(self.gas_phase.T), float(self.p_g),
                               np.array(self.gas_phase.Y, dtype=float))

        # Build the chemistry manifold (state vs enthalpy-removed) + radiation table.
        # This is the C=C_eq(h) equilibrium slice of the FGM (see design doc 4b/C0):
        # the gas thermochemical path is precomputed ONCE so the march never calls
        # Cantera — making the (physically required) equilibrium mode as fast as frozen.
        # The manifold is per-unit-mass (h in J/kg) => independent of mass flow, so one
        # build covers the whole mdot ramp; only O/F or p changes force a rebuild.
        self._build_chem_tables()

    # ------------------------------------------------------------------
    def _build_chem_tables(self, n_h=250, T_floor=340.0):
        """Chemistry manifold: tabulate gas state (+ radiation emissivity) so the march
        makes zero Cantera calls. Mode from transientProp.chemistry_transient:
          "equilibrium" (default) : C=C_eq(h) 1-D manifold vs enthalpy-removed (fast, correct)
          "frozen"                : composition fixed (validation only)
          "finite_rate"           : FPV (h, Yc) manifold + progress-variable transport
        """
        mode = self.transientProp.chemistry_transient
        self._chem_mode = mode
        g = self.gas_phase
        Tg0 = float(g.T); pg0 = float(self.p_g); Yg0 = np.array(g.Y, dtype=float)
        h0 = float(g.enthalpy_mass)
        self._gas_h0 = h0
        self._gas_p0 = pg0

        if mode == "finite_rate":
            self._chem_eq = True  # radiation built from the manifold's equilibrium edge
            iCO = self.gas_phase.species_index("CO")
            self._fpv = FPVManifold(build_fpv_manifold(
                g, Y_inlet=Yg0, T_inlet=Tg0, p=pg0,
                species_index={"CO2": self.index_CO2, "H2O": self.index_H2O, "CO": iCO},
                n_h=self.numericalProp.fpv_n_h,
                n_c=self.numericalProp.fpv_n_c,
                T_floor=T_floor,
                t_relax=self.numericalProp.fpv_t_relax,
                n_t=self.numericalProp.fpv_n_t,
                cache_dir=self.numericalProp.fpv_cache_dir))
            g.TPY = Tg0, pg0, Yg0  # restore inlet
            # radiation table over (T_eval, h) using the equilibrium-edge composition (c=1)
            m = self._fpv.m
            if self.numericalProp.radiation_ON and self.radiation_backend is not None:
                self._rad_Teval = np.linspace(280.0, 3400.0, 120)
                self._rad_hgrid = m["h_grid"]
                eps2d = np.zeros((len(self._rad_Teval), len(self._rad_hgrid)))
                for j in range(len(self._rad_hgrid)):
                    yH2O = float(m["xH2O"][j, -1]); yCO2 = float(m["xCO2"][j, -1])
                    for iT, Te in enumerate(self._rad_Teval):
                        eps2d[iT, j] = self.radiation_backend(
                            T_eval=float(Te), p=pg0, yH2O=yH2O, yCO2=yCO2, Le=self.Le)
                self._rad_eps2d = eps2d
            else:
                self._rad_Teval = None
            return

        eq = (mode == "equilibrium")
        self._chem_eq = eq

        hr, T, rho, mu, k, cp, xh2o, xco2 = ([] for _ in range(8))

        def _rec(h_removed):
            hr.append(h_removed); T.append(float(g.T)); rho.append(float(g.density))
            mu.append(float(g.viscosity)); k.append(float(g.thermal_conductivity))
            cp.append(float(g.cp)); xh2o.append(float(g.X[self.index_H2O]))
            xco2.append(float(g.X[self.index_CO2]))

        _rec(0.0)
        dh = 2200.0 * (Tg0 - T_floor) / n_h   # ~cp*dT spread over n_h steps
        cur = h0
        for _ in range(3 * n_h):
            cur -= dh
            g.HP = cur, pg0
            if eq:
                g.equilibrate('HP')
            _rec(h0 - cur)
            if g.T < T_floor:
                break
        # restore inlet state (TPY freezes composition back to the inlet burnt mix)
        g.TPY = Tg0, pg0, Yg0

        self._gas_hgrid = np.array(hr)
        self._gas_T = np.array(T); self._gas_rho = np.array(rho)
        self._gas_mu = np.array(mu); self._gas_k = np.array(k); self._gas_cp = np.array(cp)
        self._gas_xh2o = np.array(xh2o); self._gas_xco2 = np.array(xco2)

        # --- radiation emissivity table ---
        if not (self.numericalProp.radiation_ON and self.radiation_backend is not None):
            self._rad_Teval = None
            return
        self._rad_Teval = np.linspace(280.0, 3400.0, 120)
        if eq:
            # 2-D: composition varies along the path -> eps(T_eval, h_removed)
            self._rad_hgrid = np.linspace(hr[0], hr[-1], 60)
            eps2d = np.zeros((len(self._rad_Teval), len(self._rad_hgrid)))
            for j, hh in enumerate(self._rad_hgrid):
                yH2O = float(np.interp(hh, self._gas_hgrid, self._gas_xh2o))
                yCO2 = float(np.interp(hh, self._gas_hgrid, self._gas_xco2))
                for iT, Te in enumerate(self._rad_Teval):
                    eps2d[iT, j] = self.radiation_backend(
                        T_eval=float(Te), p=pg0, yH2O=yH2O, yCO2=yCO2, Le=self.Le)
            self._rad_eps2d = eps2d
        else:
            # frozen: composition fixed -> 1-D eps(T_eval)
            yH2O = float(self._gas_xh2o[0]); yCO2 = float(self._gas_xco2[0])
            self._rad_eps1d = np.array([self.radiation_backend(
                T_eval=float(Te), p=pg0, yH2O=yH2O, yCO2=yCO2, Le=self.Le)
                for Te in self._rad_Teval])

    def _gas_at(self, h_removed):
        """Interpolated gas state at a given specific enthalpy removed [J/kg]."""
        hg = self._gas_hgrid
        return (float(np.interp(h_removed, hg, self._gas_T)),
                float(np.interp(h_removed, hg, self._gas_rho)),
                float(np.interp(h_removed, hg, self._gas_mu)),
                float(np.interp(h_removed, hg, self._gas_k)),
                float(np.interp(h_removed, hg, self._gas_cp)),
                float(np.interp(h_removed, hg, self._gas_xh2o)),
                float(np.interp(h_removed, hg, self._gas_xco2)))

    def _eps_at(self, T_eval, h_removed):
        """Interpolated gas emissivity (2-D for equilibrium, 1-D for frozen)."""
        if self._chem_eq:
            # bilinear over (T_eval, h_removed)
            j = np.interp(h_removed, self._rad_hgrid, np.arange(len(self._rad_hgrid)))
            j0 = int(np.clip(np.floor(j), 0, len(self._rad_hgrid) - 2)); wj = j - j0
            col = self._rad_eps2d[:, j0] * (1 - wj) + self._rad_eps2d[:, j0 + 1] * wj
            return float(np.interp(T_eval, self._rad_Teval, col))
        return float(np.interp(T_eval, self._rad_Teval, self._rad_eps1d))

    def _h_g_rad(self, T_g, T_wg, h_removed):
        """Radiation coefficient [W/m^2K] from tabulated emissivity."""
        if self._rad_Teval is None:
            return 0.0
        eps_emit = self._eps_at(T_g, h_removed)
        eps_abs = self._eps_at(T_wg, h_removed)
        q = qrad_net_mbl(T_g, T_wg, eps_emit, eps_abs, self.corrCoeffs.emissivity_wall)
        return hrad_from_q(T_g, T_wg, q)

    # ------------------------------------------------------------------
    def _build_axial_grid(self):
        """Run one steady march to fix the total arc length, then coarsen to a
        transient grid of `transientProp.n_axial` uniform nodes (the wall field is
        smooth; a coarse grid keeps each per-RHS fluid march cheap)."""
        if getattr(self.transientProp, "skip_steady_reference_probe", False):
            self.N = int(self.transientProp.n_axial)
            self.dx = float(self.L_ch_max) / self.N
            self.Tbar_steady = np.full(self.N, np.nan)
            self._probe_steady = None
            return

        # Run the probe with chemistry matching the transient so the
        # steady reference is consistent AND setup is fast (frozen skips per-node equilibrate).
        import copy
        np_probe = copy.copy(self.numericalProp)
        np_probe.chemistry_model = self.transientProp.chemistry_transient
        np_probe.equilibrium_dh_gas_ON = (self.transientProp.chemistry_transient == "equilibrium")
        if (
            self.combustorProp.flow_config == "counter"
            and getattr(np_probe, "counterflow_physical_steady_reference", False)
        ):
            probe = solve_counterflow_physical_reference(
                self.coolantProp, self.hotgasProp, self.combustorProp,
                np_probe, self.system_requirements,
                corrCoeffs=self.corrCoeffs,
            )
        else:
            probe = main_solver(self.coolantProp, self.hotgasProp, self.combustorProp,
                                np_probe, self.system_requirements,
                                corrCoeffs=self.corrCoeffs)
            probe.solver()
        N_steady = len(probe.data_master["T_g"])
        L_total = N_steady * self.numericalProp.dx        # total arc length marched

        self.N = int(self.transientProp.n_axial)
        self.dx = L_total / self.N
        # steady wall field, resampled onto the coarse grid (settle-check reference)
        Tbar_fine = 0.5 * (np.array(probe.data_master["T_wg"])
                           + np.array(probe.data_master["T_wc"]))
        x_fine = np.arange(N_steady) * self.numericalProp.dx
        x_coarse = (np.arange(self.N) + 0.5) * self.dx
        self.Tbar_steady = np.interp(x_coarse, x_fine, Tbar_fine)
        self._probe_steady = probe

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
            OF=_interp_schedule(tp.schedule_OF, t, self.hotgasProp.mixing_ratio),
            T_lox_in=_interp_schedule(tp.schedule_T_lox_in, t, self.hotgasProp.T_inj_LOX),
            ignited=(
                _interp_schedule(getattr(tp, "schedule_ignition_state", None), t, 1.0 if t >= tp.ignition_time else 0.0) >= 0.5
            ),
        )

    # ------------------------------------------------------------------
    def fluid_pass(self, Tbar_vec, bc, record=False):
        """Public entry point: dispatch on flow_config.

        Co-flow: single forward march, coolant starts at its true inlet (x=0) —
        a well-posed IVP, no iteration.
        Counter-flow: profile relaxation (successive substitution). A naive
        single-scalar SHOOTING method (guess T_c(0), march forward, root-find
        against the true T_c(L)) was tried first and found numerically
        ill-conditioned at high NTU — the outlet becomes nearly insensitive to
        the inlet guess (a standard heat-exchanger-effectiveness fact), so the
        shooting residual goes flat and bracketing fails. Iterating on the full
        T_c(x) PROFILE instead (march gas forward using the current guess, march
        coolant backward using the resulting flux, under-relax, repeat) does not
        have this conditioning problem — same idiom as the steady solver's own
        sweep, but only 2 fields iterate here (Tbar is fixed, not co-solved).
        See class docstring.
        """
        if self.combustorProp.flow_config == "co":
            return self._march_fluids(bc["T_c_in"], bc["p_c_in"], Tbar_vec, bc)
        return self._relax_counter_flow(Tbar_vec, bc)

    # ------------------------------------------------------------------
    def _relax_counter_flow(self, Tbar_vec, bc, max_iter=None, tol=None, omega0=0.5):
        """Counter-flow: converge the coolant's T_c(x) profile by successive
        substitution with warm-starting. Gas enters at x=0 (i=0), coolant enters
        at x=L (i=N-1). Counter-flow genuinely differs from co-flow here by
        ~120-185 K in coolant outlet (verified), so this MUST actually run at
        operating conditions — an earlier NTU-threshold fallback wrongly skipped
        it there. With a realistic wall field it converges in a few iterations
        (~0.5 s cold-start, faster warm-started); the low-flow fallback (below)
        handles the ill-behaved cold-start/near-stagnant regime.
        """
        if bc["mdot_c"] < self._counter_fallback_mdot:
            return self._march_fluids(bc["T_c_in"], bc["p_c_in"], Tbar_vec, bc)

        N = self.N
        T_c_in_true = bc["T_c_in"]
        cool = self.coolantProp.coolant

        warm_started = self._counter_Tc_profile is not None and len(self._counter_Tc_profile) == N
        if warm_started:
            T_c_profile = self._counter_Tc_profile.copy()
        else:
            # initial guess: linear between the true inlet (i=N-1) and a rough
            # outlet estimate (midway between inlet and local wall temperature)
            T_c_profile = np.linspace(0.5 * (T_c_in_true + float(Tbar_vec[0])),
                                      T_c_in_true, N)
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
            res = self._march_fluids(None, bc["p_c_in"], Tbar_vec, bc, T_c_profile=T_c_profile)
            dq_cold = res["dq_cold__dx"]
            # backward march: coolant's true downstream is i=N-1 -> i=0
            T_c_new = np.empty(N)
            T_c_new[N - 1] = T_c_in_true
            for j in range(N - 1, 0, -1):
                cp_c = self._thermo.cp(cool, T_c_new[j], bc["p_c_in"])
                T_c_new[j - 1] = T_c_new[j] + (dq_cold[j] / self.N_ch) * self.dx / \
                    (bc["mdot_c"] / self.N_ch * cp_c)
                # clip inside the loop (not just after) — an intermediate value can go
                # unphysical before the backward march even finishes, crashing the next
                # step's PropsSI call. Same rationale as the co-flow per-node clip.
                T_c_new[j - 1] = float(np.clip(T_c_new[j - 1], 100.0, 4000.0))
            delta = float(np.max(np.abs(T_c_new - T_c_profile)))
            # adaptive damping: if the residual grew (a sign of oscillation, not
            # convergence), shrink the relaxation factor — but still APPLY an update
            # (at the smaller step), don't just re-probe the same unchanged state.
            if delta_prev is not None and delta > delta_prev:
                omega = max(omega * 0.5, 0.02)
            T_c_profile = (1 - omega) * T_c_profile + omega * T_c_new
            delta_prev = delta
            if delta < tol:
                break

        self._counter_Tc_profile = T_c_profile.copy()
        res["T_c"] = T_c_profile
        res["T_c_out"] = float(T_c_profile[0])   # coolant's true exit is at x=0
        return res

    # ------------------------------------------------------------------
    def _march_fluids(self, Tc0, pc0, Tbar_vec, bc, T_c_profile=None, p_c_profile=None):
        """One quasi-steady forward spatial march (x=0 -> x=L) at frozen wall field
        Tbar_vec and BCs `bc`.

        If `T_c_profile` is None (co-flow): the coolant's x=0 state is (Tc0, pc0)
        and it is MARCHED forward alongside the gas (a well-posed IVP).
        If `T_c_profile` is given (counter-flow relaxation): T_c at each node is
        READ from the profile (not marched) at a uniform pressure `pc0` — used by
        `_relax_counter_flow` to get a self-consistent dq_cold(x) for its own
        backward coolant march; `Tc0` is unused in this mode.

        Returns dict of per-node arrays: dq_hot__dx, dq_cold__dx [W/m], T_wg, T_wc,
        T_g, T_c, h_g, h_c, plus T_c_out/T_g_out (fluid states at x=L for co-flow;
        for counter-flow's profile mode T_c_out is not meaningful here — see
        `_relax_counter_flow`, which overrides it with the true x=0 exit).
        """
        N = self.N
        dx = self.dx
        cool = self.coolantProp.coolant
        s_w = self.combustorProp.thickness_coil_wall
        prescribed_Tc = (T_c_profile is not None)
        p_profile = None
        if p_c_profile is not None:
            p_profile = np.asarray(p_c_profile, dtype=float)
            if p_profile.shape != (N,):
                raise ValueError("p_c_profile must have shape (N,)")

        # --- initialise He at x=0 ---
        if prescribed_Tc:
            T_c, p_c = float(T_c_profile[0]), pc0
            if p_profile is not None:
                p_c = float(p_profile[0])
        else:
            T_c, p_c = Tc0, pc0
        rho_c = self._thermo.density(cool, T_c, p_c)
        U_c = bc["mdot_c"] / (rho_c * self.A_ch * self.N_ch)

        # --- initialise gas at inlet on the chemistry manifold (no Cantera in the march) ---
        pg0 = self._gas_p0
        h_removed = 0.0                       # cumulative specific enthalpy removed [J/kg]
        fr = (self._chem_mode == "finite_rate")
        Yc = self._fpv.Yc_inlet() if fr else 0.0   # progress variable (inlet = local eq)
        if fr:
            T_g, rho_g, mu_g, k_g, cp_g, X_H2O, X_CO2, omega_Yc = self._fpv.state(0.0, Yc)
        else:
            T_g, rho_g, mu_g, k_g, cp_g, X_H2O, X_CO2 = self._gas_at(0.0)
        p_g = pg0
        combustion_on = bool(bc["ignited"])
        mdot_g = bc["mdot_g"] if combustion_on else bc.get("mdot_lox", 0.0)
        T_g_gox = float(np.clip(bc.get("T_lox_in", self.hotgasProp.T_inj_LOX), 95.0, 1200.0))

        dq_hot = np.zeros(N); dq_cold = np.zeros(N)
        T_wg_a = np.zeros(N); T_wc_a = np.zeros(N)
        T_g_a = np.zeros(N); T_c_a = np.zeros(N)
        h_g_a = np.zeros(N); h_c_a = np.zeros(N)

        for i in range(N):
            if prescribed_Tc:
                T_c = float(T_c_profile[i])   # read-only in this mode; not marched below
                if p_profile is not None:
                    p_c = float(p_profile[i])
                rho_c = self._thermo.density(cool, T_c, p_c)
                U_c = bc["mdot_c"] / (rho_c * self.A_ch * self.N_ch)
            # ---- He (coil-side) properties + convection ----
            cp_c = self._thermo.cp(cool, T_c, p_c)
            mu_c = self._thermo.viscosity(cool, T_c, p_c)
            k_c = self._thermo.conductivity(cool, T_c, p_c)
            Re_c = rho_c * U_c * self.Dh_ch / mu_c
            Pr_c = cp_c * mu_c / k_c
            f_c = dispatch_friction_coil(self.combustorProp.friction_coil, Re=Re_c,
                                         Dh=self.Dh_ch, Rc=self.D_coil / 2,
                                         roughness=self.combustorProp.channel_roughness,
                                         x=10e10, error_factor=1.0, corrCoeffs=self.corrCoeffs)
            Nu_c = dispatch_nu_coil(self.combustorProp.Nusselt_coil, Re=Re_c, Pr=Pr_c,
                                    d=self.Dh_ch, R=self.D_coil / 2, f_fd=f_c, x=10e10,
                                    error_factor=1.0, corrCoeffs=self.corrCoeffs)
            h_c = Nu_c * k_c / self.Dh_ch

            # ---- gas (shell-side) properties from the manifold at (h_removed [, Yc]) ----
            if combustion_on and fr:
                T_g, rho_g, mu_g, k_g, cp_g, X_H2O, X_CO2, omega_Yc = self._fpv.state(h_removed, Yc)
            elif combustion_on:
                T_g, rho_g, mu_g, k_g, cp_g, X_H2O, X_CO2 = self._gas_at(h_removed)
            else:
                T_g = T_g_gox
                p_gox = max(self.hotgasProp.p0, 1e4)
                rho_g = self._thermo.density("Oxygen", T_g, p_gox)
                mu_g = self._thermo.viscosity("Oxygen", T_g, p_gox)
                k_g = self._thermo.conductivity("Oxygen", T_g, p_gox)
                cp_g = self._thermo.cp("Oxygen", T_g, p_gox)
                omega_Yc = 0.0
            U_g = mdot_g / (rho_g * self.Ap_cc) if mdot_g > 0 else 0.0
            Re_g = rho_g * U_g * self.Dh_cc / mu_g if U_g > 0 else 1.0
            Pr_g = cp_g * mu_g / k_g
            Re_sh = rho_g * U_g * (self.Dh_ch + 2 * s_w) / mu_g if U_g > 0 else 1.0

            # ---- combined in-node fixed point on T_wg -----------------------------
            # h_g depends on T_wg two ways: the Kays-Crawford (T_bulk/T_wall)^n convective
            # correction AND the radiation coefficient. The steady solver breaks this with a
            # one-node lag; here the wall is a field (no "previous node"), so iterate in-node.
            # Radiation uses the 1-D tabulated emissivity (frozen composition) — no backend
            # calls — and fluxes_at_Tbar takes the precomputed h_g_rad and does one closed-form
            # solve. Both loops from the earlier draft are merged into this single one.
            node = OneDimensionalSteadyConduction_ShellnHelicalTube(
                h_g=1.0, h_c=h_c, T_c=T_c, T_g=T_g, s_w=s_w, Dh_ch=self.Dh_ch,
                f_kw_at_T=self.func_conductivity_HX,
                T_wg_0=Tbar_vec[i], T_wc_0=Tbar_vec[i], T_c_check_0=T_c, dx=dx,
                rad_enabled=False)  # radiation supplied externally as h_g_rad
            hot_flowing = bool(U_g > 0)
            rad_on = combustion_on and hot_flowing and (self._rad_Teval is not None)
            T_wg_est = float(Tbar_vec[i])
            for _fp in range(6):
                if hot_flowing:
                    Nu_g, h_g = dispatch_nu_shell(
                        self.combustorProp.Nusselt_shell, Re_sh=Re_sh, Re_g=Re_g, Pr_g=Pr_g,
                        k_g=k_g, U_g=U_g, rho_g=rho_g, mu_g=mu_g, coil_pitch=self.coil_pitch,
                        Dh_cc=self.Dh_cc, Dh_ch=self.Dh_ch, D_coil=self.D_coil,
                        thickness_wall=s_w, T_bulk=float(T_g), T_wall=T_wg_est,
                        Nusselt_correction=self.combustorProp.Nusselt_correction,
                        error_factor=1.0, corrCoeffs=self.corrCoeffs)
                else:
                    h_g = 1.0  # inert / no-flow floor
                h_g_rad = self._h_g_rad(T_g, T_wg_est, h_removed) if rad_on else 0.0
                node.h_g = h_g
                r = node.fluxes_at_Tbar(Tbar_vec[i], h_g_rad=h_g_rad)
                if abs(r["T_wg"] - T_wg_est) < 1e-4:
                    break
                T_wg_est = r["T_wg"]

            dq_hot[i] = r["dq_hot__dx"]; dq_cold[i] = r["dq_cold__dx"]
            T_wg_a[i] = r["T_wg"]; T_wc_a[i] = r["T_wc"]
            T_g_a[i] = T_g; T_c_a[i] = T_c; h_g_a[i] = h_g; h_c_a[i] = h_c

            # ---- advance fluids one node ----
            if not prescribed_Tc:
                # He gains the cold-side flux (co-flow: T_c is marched forward here;
                # counter-flow's own backward march lives in _relax_counter_flow)
                dT_c = (dq_cold[i] / self.N_ch) * dx / (bc["mdot_c"] / self.N_ch * cp_c)
                T_c += dT_c
                # Guard against explicit-march blowup at tiny mdot_c (local NTU =
                # h_c*P*dx/(mdot_c*cp_c) -> large as mdot_c->0).
                T_c = float(np.clip(T_c, 100.0, 4000.0))
                # simple isothermal-wall-consistent pressure drop (kept 1st-order, as steady)
                dp_c = -f_c * rho_c * U_c ** 2 / (2 * self.Dh_ch) * dx
                p_c += dp_c
                rho_c = self._thermo.density(cool, T_c, max(p_c, 1e4))
                U_c = bc["mdot_c"] / (rho_c * self.A_ch * self.N_ch)

            # gas loses the hot-side flux: advance along the manifold by accumulating the
            # specific enthalpy removed (no Cantera call — T_g etc. are read from the
            # table at the top of the next node). Per-unit-mass, so valid for any mdot_g.
            # Finite-rate: also transport the progress variable dYc/dx = omega_Yc / U_g.
            if combustion_on and mdot_g > 0:
                h_removed += dq_hot[i] * dx / mdot_g
                if fr and U_g > 0:
                    Yc += omega_Yc * dx / U_g
            elif mdot_g > 0:
                T_g_gox -= dq_hot[i] * dx / (mdot_g * cp_g)
                T_g_gox = float(np.clip(T_g_gox, 95.0, 1200.0))

        return dict(dq_hot__dx=dq_hot, dq_cold__dx=dq_cold, T_wg=T_wg_a, T_wc=T_wc_a,
                    T_g=T_g_a, T_c=T_c_a, h_g=h_g_a, h_c=h_c_a,
                    T_c_out=T_c, T_g_out=T_g)

    # ------------------------------------------------------------------
    def _wall_rhs(self, t, Tbar_vec):
        return self._wall_rate(t, Tbar_vec)[0]

    def _wall_rate(self, t, Tbar_vec):
        """One fluid pass -> (dTbar/dt, lambda) per node. `lambda` is the local
        stability eigenvalue -d(RHS_i)/d(Tbar_i) ~ (sum of film conductances)/
        (thermal mass) — used by the linearly-implicit fixed-step integrator to
        stay unconditionally stable at large dt. Overestimating lambda is safe
        (over-damped, still stable); the He film (h_c ~ 1e5) dominates it, so the
        small radiation contribution to the hot side is negligibly omitted here.
        """
        bc = self._bc_at(t)
        res = self.fluid_pass(Tbar_vec, bc)
        cp_w = np.array([self.func_cp_HX(T - 273.15) for T in Tbar_vec])
        denom = self.rho_w * cp_w * self.A_wall
        dTbar_dt = (res["dq_hot__dx"] - res["dq_cold__dx"]) / denom
        P_h = np.pi * (self.Dh_ch + 2 * self.combustorProp.thickness_coil_wall)
        P_c = np.pi * self.Dh_ch
        lam = (np.asarray(res["h_g"]) * P_h + np.asarray(res["h_c"]) * P_c) / denom
        return dTbar_dt, lam

    # ------------------------------------------------------------------
    def solve_transient(self, verbose=True):
        tp = self.transientProp
        if verbose:
            self.time_scale_audit()
        Tbar0 = np.full(self.N, tp.T_wall_initial)
        t_eval = np.linspace(0.0, tp.t_end, tp.n_save)
        if verbose:
            print(f"Integrating wall ODE: {self.N} nodes, t_end={tp.t_end}s, "
                  f"method={tp.solver_method}, max_step={tp.max_step}s ...")
        if tp.solver_method.lower() in ("fixed_step", "euler", "explicit_euler"):
            # Bounded-cost production path (shared with the shell-and-tube transient):
            # one fluid pass per time step, so counter-flow's per-RHS relaxation cost
            # is paid once per step instead of the hundreds of times BDF probes it.
            sol = self._solve_fixed_step(Tbar0, t_eval)
        else:
            ivp_options = {}
            if tp.use_sparse_jacobian and tp.solver_method in ("BDF", "Radau"):
                ivp_options["jac_sparsity"] = np.eye(self.N, dtype=bool)
            sol = solve_ivp(self._wall_rhs, (0.0, tp.t_end), Tbar0, method=tp.solver_method,
                            t_eval=t_eval, max_step=tp.max_step, rtol=1e-5, atol=1e-2,
                            **ivp_options)
            if not sol.success:
                print(f"  WARNING: solve_ivp did not fully succeed: {sol.message}")
        self.sol = sol
        self._build_time_series(sol)
        if verbose:
            self._print_transient_summary()
        return sol

    def solve_transient_core(self, verbose=True):
        """Run helical transient wall + conserved helium mass/internal energy.

        Hot gas remains quasi-steady through the existing helical `fluid_pass`
        physics. The coolant temperature profile is reconstructed from conserved
        `(m, U, V)` and passed into the wall/hot-side flux evaluator; coolant
        mass and internal energy are then advanced by the finite-volume kernel.
        """

        tp = self.transientProp
        if verbose:
            self.time_scale_audit()
            print(
                "Integrating helical transient core: "
                f"{self.N} nodes, t_end={tp.t_end}s, max_step={tp.max_step}s, "
                f"chemistry={self._chem_mode}, flow={self.combustorProp.flow_config}"
            )

        t_base = fixed_time_grid(
            t_end=float(tp.t_end),
            max_step=float(tp.max_step),
            t_eval=np.linspace(0.0, float(tp.t_end), int(tp.n_save)),
            schedules=(
                tp.schedule_mass_flow_c,
                tp.schedule_mass_flow_g,
                getattr(tp, "schedule_mass_flow_lox", None),
                getattr(tp, "schedule_mass_flow_diesel", None),
                tp.schedule_T_c_in,
                tp.schedule_p_c_in,
                getattr(tp, "schedule_p_c_out", None),
                tp.schedule_T_lox_in,
                tp.schedule_ignition_state,
                tp.schedule_OF,
            ),
        )
        momentum_model = str(getattr(tp, "coolant_momentum_model", "quasi_steady"))
        if momentum_model not in ("quasi_steady", "low_mach"):
            raise ValueError("coolant_momentum_model must be 'quasi_steady' or 'low_mach'")
        flow = self._helical_flow_direction()

        volume = np.full(self.N, self.A_ch * self.N_ch * self.dx, dtype=float)
        T_c0 = np.full(
            self.N,
            float(self.coolantProp.T_in if tp.T_coolant_initial is None else tp.T_coolant_initial),
        )
        p_in0 = float(_interp_schedule(tp.schedule_p_c_in, 0.0, self.coolantProp.p_in))
        mdot_ref = max(
            self._schedule_max_abs(tp.schedule_mass_flow_c, self.coolantProp.mass_flow_c),
            getattr(tp, "transient_coolant_mdot_floor", 1e-9),
        )
        dp_nominal = self._helical_nominal_coolant_dp(T_c0[0], p_in0, mdot_ref)
        p_outlet0 = self._helical_outlet_pressure_at(0.0, max(p_in0 - dp_nominal, 1.0e3))
        if momentum_model == "low_mach":
            p_initial = self._helical_boundary_pressure_profile(
                inlet_pressure=p_in0,
                outlet_pressure=p_outlet0,
                flow_direction=flow,
            )
        else:
            p_initial = self._helical_initial_pressure_profile(p_in0, dp_nominal)
        mass, internal_energy = initial_mass_energy_from_TP(
            T_c0,
            p_initial,
            volume,
            self.coolantProp.coolant,
        )

        t = np.asarray(t_base, dtype=float)
        n_time = len(t)
        progress = (
            TransientProgressPrinter.from_config(tp, total_steps=max(n_time - 1, 0))
            if verbose else
            None
        )
        outlet_pressure = max(float(p_outlet0), 1.0e3)
        resistance = self._helical_face_resistance(dp_nominal, mass / volume, mdot_ref)
        inertance = self._helical_face_inertance()
        global_mdot_cap = max(2.0 * mdot_ref, 1e-12)

        T_wall = np.zeros((n_time, self.N))
        T_coolant = np.zeros((n_time, self.N))
        coolant_mass = np.zeros((n_time, self.N))
        coolant_U = np.zeros((n_time, self.N))
        coolant_pressure = np.zeros((n_time, self.N))
        coolant_density = np.zeros((n_time, self.N))
        coolant_h = np.zeros((n_time, self.N))
        face_mdot = np.zeros((n_time, self.N + 1))
        heat_wall_to_coolant_W = np.zeros((n_time, self.N))
        energy_residual = np.zeros(n_time)
        mass_residual = np.zeros(n_time)
        hot_heat_added = np.zeros(n_time)
        adv_in = np.zeros(n_time)
        adv_out = np.zeros(n_time)

        extra_fields = {
            key: np.full((n_time, self.N), np.nan)
            for key in ("T_wg", "T_wc", "T_g", "dq_hot__dx", "dq_cold__dx", "h_g", "h_c")
        }

        T_wall[0] = np.full(self.N, float(tp.T_wall_initial))
        coolant_mass[0] = mass
        coolant_U[0] = internal_energy
        state0 = coolprop_state_from_mass_energy(mass, internal_energy, volume, self.coolantProp.coolant)
        T_coolant[0] = state0.temperature
        coolant_pressure[0] = state0.pressure
        coolant_density[0] = state0.density
        coolant_h[0] = state0.specific_enthalpy_J_kg

        last_step = None
        for j in range(n_time - 1):
            tj = float(t[j])
            dt = float(t[j + 1] - t[j])
            bc = self._bc_at(tj)
            state = coolprop_state_from_mass_energy(
                coolant_mass[j],
                coolant_U[j],
                volume,
                self.coolantProp.coolant,
            )
            p_inlet = max(float(bc["p_c_in"]), 1.0e3)
            outlet_pressure = max(
                self._helical_outlet_pressure_at(tj, outlet_pressure),
                1.0e3,
            )
            mdot_cmd = max(float(bc["mdot_c"]), 0.0)
            mdot_floor = getattr(tp, "transient_coolant_mdot_floor", 1e-9)
            if momentum_model == "low_mach":
                p_transport = self._helical_boundary_pressure_profile(
                    inlet_pressure=p_inlet,
                    outlet_pressure=outlet_pressure,
                    flow_direction=flow,
                )
                faces = self._helical_low_mach_faces(
                    face_old=face_mdot[j],
                    temperature=state.temperature,
                    pressure_profile=p_transport,
                    density=state.density,
                    dt=dt,
                    inlet_pressure=p_inlet,
                    outlet_pressure=outlet_pressure,
                    flow_direction=flow,
                    mdot_reference=mdot_ref,
                    mdot_floor=mdot_floor,
                )
                mdot_cap = max(mdot_cmd, mdot_floor)
                faces = np.clip(faces, -min(global_mdot_cap, mdot_cap), min(global_mdot_cap, mdot_cap))
                faces = _limit_face_mdot_for_inventory(
                    coolant_mass[j],
                    faces,
                    dt,
                    internal_energy_J=coolant_U[j],
                    specific_enthalpy_J_kg=state.specific_enthalpy_J_kg,
                )
            else:
                p_transport = state.pressure
                faces = self._helical_quasi_steady_faces(
                    state.pressure,
                    state.density,
                    resistance,
                    mdot_inlet=mdot_cmd,
                    outlet_pressure=outlet_pressure,
                    flow_direction=flow,
                    mdot_floor=mdot_floor,
                )
                faces = np.clip(faces, -global_mdot_cap, global_mdot_cap)
            if momentum_model != "low_mach" and mdot_cmd <= mdot_floor and j > 0:
                outlet_face = -1 if flow == 1 else 0
                if abs(faces[outlet_face]) <= mdot_floor:
                    faces[outlet_face] = 0.5 * face_mdot[j, outlet_face]
            if momentum_model != "low_mach" and mdot_cmd <= mdot_floor:
                faces = _limit_face_mdot_for_inventory(
                    coolant_mass[j],
                    faces,
                    dt,
                    internal_energy_J=coolant_U[j],
                    specific_enthalpy_J_kg=state.specific_enthalpy_J_kg,
                )
            mdot_effective = max(float(np.mean(np.abs(faces))), mdot_floor)
            bc_eval = dict(bc)
            bc_eval["mdot_c"] = mdot_effective

            flux = self._march_fluids(
                None,
                float(np.mean(p_transport)),
                T_wall[j],
                bc_eval,
                T_c_profile=state.temperature,
                p_c_profile=p_transport,
            )
            for key in extra_fields:
                extra_fields[key][j] = flux[key]

            cp_w = np.array([self.func_cp_HX(Ti - 273.15) for Ti in T_wall[j]])
            wall_capacity = self.rho_w * cp_w * self.A_wall * self.dx
            hot_heat_W = np.asarray(flux["dq_hot__dx"], dtype=float) * self.dx
            delta_wc = np.asarray(flux["T_wc"], dtype=float) - state.temperature
            cold_heat_W = np.asarray(flux["dq_cold__dx"], dtype=float) * self.dx
            conductance = np.divide(
                cold_heat_W,
                delta_wc,
                out=np.zeros_like(cold_heat_W),
                where=np.abs(delta_wc) > 1.0e-9,
            )
            conductance = np.maximum(conductance, 0.0)
            h_in = self._thermo.enthalpy(
                self.coolantProp.coolant, float(bc["T_c_in"]), max(float(bc["p_c_in"]), 1.0e3)
            )
            step = semi_implicit_wall_compressible_coolant_step(
                T_wall[j],
                wall_capacity,
                coolant_mass[j],
                coolant_U[j],
                state.temperature,
                state.specific_enthalpy_J_kg,
                faces,
                hot_heat_W,
                conductance,
                dt,
                inlet_enthalpy_J_kg=h_in,
                outlet_backflow_enthalpy_J_kg=h_in,
                mass_floor=1.0e-12,
            )
            m_candidate, U_candidate = enforce_density_bounds(
                step.coolant.mass_new,
                step.coolant.internal_energy_new_J,
                volume,
            )
            if momentum_model == "low_mach":
                provisional_U = enforce_internal_energy_floor(
                    m_candidate,
                    U_candidate,
                    volume,
                    self.coolantProp.coolant,
                    clip=False,
                )
                provisional_state = coolprop_state_from_mass_energy(
                    m_candidate,
                    provisional_U,
                    volume,
                    self.coolantProp.coolant,
                )
                p_projected = self._helical_boundary_pressure_profile(
                    inlet_pressure=p_inlet,
                    outlet_pressure=outlet_pressure,
                    flow_direction=flow,
                )
                m_new, U_new = _coolant_mass_energy_from_TP_profile(
                    provisional_state.temperature,
                    p_projected,
                    volume,
                    self.coolantProp.coolant,
                )
            else:
                m_new = m_candidate
                U_new = enforce_internal_energy_floor(
                    m_new,
                    U_candidate,
                    volume,
                    self.coolantProp.coolant,
                    clip=False,
                )
            new_state = coolprop_state_from_mass_energy(
                m_new,
                U_new,
                volume,
                self.coolantProp.coolant,
            )

            T_wall[j + 1] = step.T_wall_new
            coolant_mass[j + 1] = m_new
            coolant_U[j + 1] = U_new
            T_coolant[j + 1] = new_state.temperature
            coolant_pressure[j + 1] = new_state.pressure
            coolant_density[j + 1] = new_state.density
            coolant_h[j + 1] = new_state.specific_enthalpy_J_kg
            face_mdot[j + 1] = faces
            heat_wall_to_coolant_W[j + 1] = step.heat_wall_to_coolant_W
            energy_residual[j + 1] = step.total_energy_residual_J
            mass_residual[j + 1] = step.coolant.mass_residual_kg
            hot_heat_added[j + 1] = step.hot_heat_added_J
            adv_in[j + 1] = step.coolant.advective_energy_in_J
            adv_out[j + 1] = step.coolant.advective_energy_out_J
            last_step = step
            if progress is not None:
                outlet_index = self.N - 1 if flow == 1 else 0
                progress.update(
                    step=j + 1,
                    time_s=float(t[j + 1]),
                    T_wall=T_wall[j + 1],
                    T_coolant_outlet=float(new_state.temperature[outlet_index]),
                    p_coolant_outlet=float(new_state.pressure[outlet_index]),
                    T_gas_outlet=flux.get("T_g_out"),
                )

        # Fill final diagnostic fields from the final state for dashboard consistency.
        final_bc = self._bc_at(float(t[-1]))
        final_state = coolprop_state_from_mass_energy(coolant_mass[-1], coolant_U[-1], volume, self.coolantProp.coolant)
        final_flux = self._march_fluids(
            None,
            float(np.mean(final_state.pressure)),
            T_wall[-1],
            {**final_bc, "mdot_c": max(float(np.mean(np.abs(face_mdot[-1]))), getattr(tp, "transient_coolant_mdot_floor", 1e-9))},
            T_c_profile=final_state.temperature,
            p_c_profile=(
                final_state.pressure
                if momentum_model != "low_mach" else
                self._helical_boundary_pressure_profile(
                    inlet_pressure=max(float(final_bc["p_c_in"]), 1.0e3),
                    outlet_pressure=max(
                        self._helical_outlet_pressure_at(float(t[-1]), outlet_pressure),
                        1.0e3,
                    ),
                    flow_direction=flow,
                )
            ),
        )
        for key in extra_fields:
            extra_fields[key][-1] = final_flux[key]

        outlet_index = self.N - 1 if flow == 1 else 0
        fields = {
            "Tbar": T_wall,
            "T_c": T_coolant,
            "coolant_mass_kg": coolant_mass,
            "coolant_internal_energy_J": coolant_U,
            "p_c": coolant_pressure,
            "rho_c_state": coolant_density,
            "h_c_state": coolant_h,
            "face_mdot_c": face_mdot,
            "heat_wall_to_coolant_W": heat_wall_to_coolant_W,
            **extra_fields,
        }
        scalars = {
            "T_c_out": T_coolant[:, outlet_index],
            "T_g_out": extra_fields["T_g"][:, -1],
            "Q_hot_kW": np.sum(extra_fields["dq_hot__dx"], axis=1) * self.dx / 1.0e3,
            "Q_cold_kW": np.sum(heat_wall_to_coolant_W, axis=1) / 1.0e3,
            "dT_wall_max": np.max(extra_fields["T_wg"] - extra_fields["T_wc"], axis=1),
            "T_wall_max": np.max(T_wall, axis=1),
            "T_wall_min": np.min(T_wall, axis=1),
            "T_c_min": np.min(T_coolant, axis=1),
            "T_c_max": np.max(T_coolant, axis=1),
            "mdot_c": np.array([self._bc_at(float(ti))["mdot_c"] for ti in t]),
            "mdot_c_effective": np.mean(np.abs(face_mdot), axis=1),
            "mdot_c_inlet_face": face_mdot[:, 0] if flow == 1 else -face_mdot[:, -1],
            "mdot_c_outlet_face": face_mdot[:, -1] if flow == 1 else -face_mdot[:, 0],
            "mdot_g": np.array([self._bc_at(float(ti))["mdot_g"] for ti in t]),
            "energy_residual_J": energy_residual,
            "coolant_mass_residual_kg": mass_residual,
            "hot_heat_added_J": hot_heat_added,
            "advective_energy_in_J": adv_in,
            "advective_energy_out_J": adv_out,
            "He_residence_s": np.sum(coolant_mass, axis=1)
            / np.maximum(np.abs(np.array([self._bc_at(float(ti))["mdot_c"] for ti in t])), 1.0e-9),
            "tau_wall_min_s": np.array([self._bc_change_timescale(float(ti)) for ti in t]),
            "He_outlet_reliable": np.ones(n_time),
        }
        self.core_result = SimpleNamespace(
            integration=SimpleNamespace(
                t=t,
                T_wall=T_wall,
                T_coolant=T_coolant,
                coolant_mass_kg=coolant_mass,
                coolant_internal_energy_J=coolant_U,
                coolant_pressure_Pa=coolant_pressure,
                face_mdot_kg_s=face_mdot,
                last_step=last_step,
            )
        )
        self.sol = SimpleNamespace(
            t=t,
            y=T_wall.T,
            success=True,
            message="helical transient_core completed",
            nfev=max(len(t) - 1, 0),
        )
        self.time_series = dict(t=t, x=np.arange(self.N) * self.dx, fields=fields, scalars=scalars)
        if verbose:
            self._print_transient_summary()
        return self.core_result

    def _solve_fixed_step(self, Tbar0, t_eval):
        """Bounded-cost LINEARLY-IMPLICIT (semi-implicit Euler) wall integrator —
        one fluid pass per time step, but unconditionally stable at large dt.

        Plain forward Euler is UNSTABLE here: the He-cooled wall time constant is
        ~0.04-0.09 s (h_c ~ 1e5 W/m2K), far below a practical dt=0.25 s, so an
        explicit step blows up (verified: coolant exits hotter than the gas inlet,
        ~30x energy imbalance). The stiff part is the LOCAL wall<->film coupling,
        which is diagonal, so we treat it implicitly per node:

            Tbar_i^{n+1} = Tbar_i^n + dt * R_i(Tbar^n) / (1 + dt * lambda_i)

        The damping factor dt*lambda/(1+dt*lambda) in (0,1) for all dt>0 => stable
        for any step, one fluid pass per step (lambda comes from the same pass).
        Schedule breakpoints are inserted so ignition/flow jumps aren't stepped
        across. (The shell-and-tube transient's _solve_fixed_step should get the
        same treatment — its wall is also He-film-dominated.)
        """
        grid = self._fixed_time_grid(t_eval)
        y = np.zeros((self.N, len(grid)))
        y[:, 0] = Tbar0
        nfev = 0
        for j in range(len(grid) - 1):
            t = float(grid[j])
            dt = float(grid[j + 1] - grid[j])
            R, lam = self._wall_rate(t, y[:, j])
            nfev += 1
            y[:, j + 1] = np.clip(y[:, j] + dt * R / (1.0 + dt * lam), 80.0, 4000.0)
        return SimpleNamespace(t=grid, y=y, success=True,
                               message="fixed_step (linearly-implicit) completed", nfev=nfev)

    def _fixed_time_grid(self, t_eval):
        tp = self.transientProp
        base = np.arange(0.0, tp.t_end + 0.5 * tp.max_step, tp.max_step)
        points = [base, np.asarray(t_eval, dtype=float), np.array([0.0, tp.t_end])]
        for name in ("schedule_mass_flow_c", "schedule_mass_flow_g", "schedule_mass_flow_lox",
                     "schedule_mass_flow_diesel", "schedule_T_c_in", "schedule_p_c_in",
                     "schedule_T_lox_in", "schedule_ignition_state", "schedule_OF"):
            schedule = getattr(tp, name, None)
            if schedule:
                points.append(np.array([float(t) for t, _ in schedule]))
        return np.unique(np.clip(np.concatenate(points), 0.0, tp.t_end))

    # ------------------------------------------------------------------
    def _build_time_series(self, sol):
        """Post-process: one fluid_pass per stored time to log full per-node fields."""
        times = sol.t
        n_t = len(times)
        fields = {k: np.zeros((n_t, self.N)) for k in
                  ("Tbar", "T_wg", "T_wc", "T_g", "T_c", "dq_hot__dx", "dq_cold__dx", "h_g", "h_c")}
        scalars = {k: np.zeros(n_t) for k in
                   ("T_c_out", "T_g_out", "Q_hot_kW", "Q_cold_kW", "dT_wall_max",
                    "mdot_c", "mdot_g", "He_residence_s", "tau_wall_min_s", "He_outlet_reliable")}
        He_inv = self._He_inventory()
        for j, t in enumerate(times):
            Tbar = sol.y[:, j]
            bc = self._bc_at(t)
            res = self.fluid_pass(Tbar, bc)
            fields["Tbar"][j] = Tbar
            for k in ("T_wg", "T_wc", "T_g", "T_c", "dq_hot__dx", "dq_cold__dx", "h_g", "h_c"):
                fields[k][j] = res[k]
            scalars["T_c_out"][j] = res["T_c_out"]
            scalars["T_g_out"][j] = res["T_g_out"]
            scalars["Q_hot_kW"][j] = np.sum(res["dq_hot__dx"]) * self.dx / 1e3
            scalars["Q_cold_kW"][j] = np.sum(res["dq_cold__dx"]) * self.dx / 1e3
            scalars["dT_wall_max"][j] = np.max(res["T_wg"] - res["T_wc"])
            scalars["mdot_c"][j] = bc["mdot_c"]
            scalars["mdot_g"][j] = bc["mdot_g"]
            res_t = He_inv / max(bc["mdot_c"], 1e-9)
            scalars["He_residence_s"][j] = res_t
            # He-outlet reliability (section 4.6): the quasi-steady-fluid assumption for
            # He fails when its residence time is not small compared to how fast the
            # DRIVING CONDITIONS change — NOT compared to tau_wall. Under steady flow
            # (post-ramp) nothing changes, so the outlet is reliable no matter how the
            # residence compares to tau_wall. Use the BC rate-of-change timescale
            # tau_BC = |mdot_c / (dmdot_c/dt)| (== inf when the schedule is flat).
            tau_BC = self._bc_change_timescale(t)
            scalars["tau_wall_min_s"][j] = tau_BC   # (kept key name; now stores tau_BC)
            unreliable = (res_t > tau_BC) and self.transientProp.flag_He_outlet_when_residence_gt_tau
            scalars["He_outlet_reliable"][j] = 0.0 if unreliable else 1.0
        self.time_series = dict(t=times, x=np.arange(self.N) * self.dx, fields=fields, scalars=scalars)

    # ------------------------------------------------------------------
    def _bc_change_timescale(self, t, eps=0.05):
        """Fastest relative rate-of-change timescale min|X/(dX/dt)| across the mass-flow
        schedules at time t. Returns +inf when all driving schedules are locally flat
        (post-ramp steady operation). Used for the He-outlet reliability flag (4.6)."""
        tau = np.inf
        for sched, base in ((self.transientProp.schedule_mass_flow_c, self.coolantProp.mass_flow_c),
                            (self.transientProp.schedule_mass_flow_g, self.hotgasProp.mass_flow_g)):
            if sched is None:
                continue
            x0 = _interp_schedule(sched, t - eps, base)
            x1 = _interp_schedule(sched, t + eps, base)
            dxdt = (x1 - x0) / (2 * eps)
            xmid = 0.5 * (x0 + x1)
            if abs(dxdt) > 1e-12:
                tau = min(tau, abs(xmid / dxdt))
        return tau

    # ------------------------------------------------------------------
    def _helical_flow_direction(self):
        if self.combustorProp.flow_config == "co":
            return 1
        if self.combustorProp.flow_config == "counter":
            return -1
        raise ValueError("flow_config must be 'co' or 'counter'")

    def _helical_nominal_coolant_dp(self, T: float, p: float, mdot_total: float) -> float:
        rho = self._thermo.density(self.coolantProp.coolant, float(T), float(p))
        mu = self._thermo.viscosity(self.coolantProp.coolant, float(T), float(p))
        velocity = abs(float(mdot_total)) / max(rho * self.A_ch * self.N_ch, 1.0e-30)
        Re = max(rho * velocity * self.Dh_ch / mu, 1.0)
        f = dispatch_friction_coil(
            self.combustorProp.friction_coil,
            Re=Re,
            Dh=self.Dh_ch,
            Rc=self.D_coil / 2,
            roughness=self.combustorProp.channel_roughness,
            x=10e10,
            error_factor=1.0,
            corrCoeffs=self.corrCoeffs,
        )
        return max(float(f * rho * velocity**2 / (2.0 * self.Dh_ch) * self.L_ch_max), 1.0)

    def _helical_initial_pressure_profile(self, inlet_pressure: float, pressure_drop: float) -> np.ndarray:
        frac = np.linspace(0.0, 1.0, self.N)
        if self._helical_flow_direction() == 1:
            return float(inlet_pressure) - float(pressure_drop) * frac
        return float(inlet_pressure) - float(pressure_drop) * frac[::-1]

    def _helical_outlet_pressure_at(self, t: float, default: float) -> float:
        schedule = getattr(self.transientProp, "schedule_p_c_out", None)
        fixed = getattr(self.transientProp, "transient_coolant_outlet_pressure", None)
        if schedule is not None:
            return float(_interp_schedule(schedule, t, default))
        if fixed is not None:
            return float(fixed)
        return float(default)

    def _helical_boundary_pressure_profile(
        self,
        *,
        inlet_pressure: float,
        outlet_pressure: float,
        flow_direction: int,
    ) -> np.ndarray:
        """Return cell-centered pressure field between hydraulic boundaries."""

        frac = (np.arange(self.N, dtype=float) + 0.5) / max(self.N, 1)
        if int(flow_direction) == 1:
            left = float(inlet_pressure)
            right = float(outlet_pressure)
        else:
            left = float(outlet_pressure)
            right = float(inlet_pressure)
        return left + (right - left) * frac

    def _helical_face_inertance(self) -> np.ndarray:
        """Return low-Mach face inertance coefficients for the coil pipe."""

        area = max(float(self.A_ch * self.N_ch), 1.0e-30)
        dx = float(self.dx)
        inertance = np.empty(self.N + 1, dtype=float)
        inertance[0] = 0.5 * dx / area
        inertance[1:-1] = dx / area
        inertance[-1] = 0.5 * dx / area
        return inertance

    def _helical_low_mach_faces(
        self,
        *,
        face_old: np.ndarray,
        temperature: np.ndarray,
        pressure_profile: np.ndarray,
        density: np.ndarray,
        dt: float,
        inlet_pressure: float,
        outlet_pressure: float,
        flow_direction: int,
        mdot_reference: float,
        mdot_floor: float,
    ) -> np.ndarray:
        """Advance a pressure-driven helical through-flow momentum state.

        This is a low-Mach pipe closure for the coupled validation use case:
        pressure is assumed to equilibrate over the coil much faster than the
        wall/coolant thermal transient, so the pipe uses one implicit
        inertance/friction state for the through-flow. Coolant mass and energy
        remain distributed finite-volume states.
        """

        old = np.asarray(face_old, dtype=float)
        p = np.asarray(pressure_profile, dtype=float)
        rho = np.asarray(density, dtype=float)
        if old.shape != (self.N + 1,) or p.shape != (self.N,) or rho.shape != (self.N,):
            raise ValueError("helical low-Mach arrays have inconsistent shapes")
        if dt <= 0.0:
            return old.copy()

        flow = 1 if int(flow_direction) >= 0 else -1
        old_through = max(float(flow) * float(np.nanmean(old)), 0.0)
        drive = max(float(inlet_pressure) - float(outlet_pressure), 0.0)
        if drive <= 0.0:
            return np.zeros(self.N + 1, dtype=float)

        resistance = self._helical_lumped_resistance_over_density(
            temperature=temperature,
            pressure_profile=pressure_profile,
            density=rho,
            mdot=max(old_through, min(abs(float(mdot_reference)), 0.2), float(mdot_floor)),
        )
        inertance = max(float(self.L_ch_max) / max(float(self.A_ch * self.N_ch), 1.0e-30), 1.0e-30)
        mdot_new = _implicit_quadratic_momentum_update(
            old_mdot=old_through,
            pressure_drive=drive,
            inertance=inertance,
            resistance_over_density=resistance,
            dt=dt,
        )
        return np.full(self.N + 1, flow * max(float(mdot_new), 0.0), dtype=float)

    def _helical_lumped_resistance_over_density(
        self,
        *,
        temperature: np.ndarray,
        pressure_profile: np.ndarray,
        density: np.ndarray,
        mdot: float,
    ) -> float:
        """Return whole-coil Darcy coefficient `DeltaP / mdot^2`."""

        T_mean = float(np.nanmean(np.asarray(temperature, dtype=float)))
        p_mean = max(float(np.nanmean(np.asarray(pressure_profile, dtype=float))), 1.0e3)
        rho_mean = max(float(np.nanmean(np.asarray(density, dtype=float))), 1.0e-9)
        area = max(float(self.A_ch * self.N_ch), 1.0e-30)
        mu = self._thermo.viscosity(self.coolantProp.coolant, T_mean, p_mean)
        Re = max(abs(float(mdot)) * self.Dh_ch / max(mu * area, 1.0e-30), 1.0)
        f = dispatch_friction_coil(
            self.combustorProp.friction_coil,
            Re=Re,
            Dh=self.Dh_ch,
            Rc=self.D_coil / 2,
            roughness=self.combustorProp.channel_roughness,
            x=10e10,
            error_factor=1.0,
            corrCoeffs=self.corrCoeffs,
        )
        return float(f * self.L_ch_max / (2.0 * self.Dh_ch * rho_mean * area * area))

    def _helical_face_resistance_over_density(
        self,
        *,
        temperature: np.ndarray,
        pressure_profile: np.ndarray,
        density_face: np.ndarray,
        mdot_reference: float,
        mdot_floor: float,
    ) -> np.ndarray:
        """Return Darcy pipe coefficients `DeltaP / mdot^2` per face."""

        T = np.asarray(temperature, dtype=float)
        p = np.asarray(pressure_profile, dtype=float)
        area = max(float(self.A_ch * self.N_ch), 1.0e-30)
        dx = float(self.dx)
        length = np.empty(self.N + 1, dtype=float)
        length[0] = 0.5 * dx
        length[1:-1] = dx
        length[-1] = 0.5 * dx
        T_face = np.empty(self.N + 1, dtype=float)
        p_face = np.empty(self.N + 1, dtype=float)
        T_face[0] = T[0]
        T_face[1:-1] = 0.5 * (T[:-1] + T[1:])
        T_face[-1] = T[-1]
        p_face[0] = p[0]
        p_face[1:-1] = 0.5 * (p[:-1] + p[1:])
        p_face[-1] = p[-1]

        mdot = max(abs(float(mdot_reference)), float(mdot_floor))
        out = np.empty(self.N + 1, dtype=float)
        for i in range(self.N + 1):
            mu = self._thermo.viscosity(
                self.coolantProp.coolant, float(T_face[i]), max(float(p_face[i]), 1.0e3)
            )
            Re = max(mdot * self.Dh_ch / max(mu * area, 1.0e-30), 1.0)
            f = dispatch_friction_coil(
                self.combustorProp.friction_coil,
                Re=Re,
                Dh=self.Dh_ch,
                Rc=self.D_coil / 2,
                roughness=self.combustorProp.channel_roughness,
                x=10e10,
                error_factor=1.0,
                corrCoeffs=self.corrCoeffs,
            )
            rho_i = max(float(density_face[i]), 1.0e-9)
            out[i] = f * float(length[i]) / (2.0 * self.Dh_ch * rho_i * area * area)
        return out

    @staticmethod
    def _cap_helical_inlet_face(
        faces: np.ndarray,
        *,
        mdot_cap: float,
        flow_direction: int,
    ) -> np.ndarray:
        """Limit only the upstream valve face; internal/outlet residual flow may continue."""

        capped = np.asarray(faces, dtype=float).copy()
        cap = max(float(mdot_cap), 0.0)
        if int(flow_direction) == 1:
            capped[0] = min(max(capped[0], 0.0), cap)
        else:
            capped[-1] = -min(max(-capped[-1], 0.0), cap)
        return capped

    def _helical_face_resistance(self, pressure_drop: float, density: np.ndarray, mdot_reference: float) -> np.ndarray:
        rho_mean = max(float(np.nanmean(density)), 1.0e-9)
        mdot_ref = max(abs(float(mdot_reference)), 1.0e-9)
        dp_face = max(float(pressure_drop), 1.0) / max(self.N, 1)
        return np.full(self.N, rho_mean * dp_face / mdot_ref**2)

    def _helical_quasi_steady_faces(
        self,
        pressure: np.ndarray,
        density: np.ndarray,
        resistance: np.ndarray,
        *,
        mdot_inlet: float,
        outlet_pressure: float,
        flow_direction: int,
        mdot_floor: float,
    ) -> np.ndarray:
        p = np.asarray(pressure, dtype=float)
        rho = np.asarray(density, dtype=float)
        face = np.zeros(self.N + 1, dtype=float)
        for i in range(self.N - 1):
            rho_face = 0.5 * (rho[i] + rho[i + 1])
            face[i + 1] = self._orifice_mdot(p[i] - p[i + 1], rho_face, resistance[i])

        mdot_cmd = max(float(mdot_inlet), 0.0)
        if mdot_cmd > mdot_floor:
            return np.full(self.N + 1, mdot_cmd if flow_direction == 1 else -mdot_cmd, dtype=float)

        if flow_direction == 1:
            face[-1] = max(
                self._orifice_mdot(p[-1] - float(outlet_pressure), rho[-1], resistance[-1]),
                0.0,
            )
        else:
            face[0] = min(
                self._orifice_mdot(float(outlet_pressure) - p[0], rho[0], resistance[-1]),
                0.0,
            )
        return face

    @staticmethod
    def _orifice_mdot(dp: float, rho: float, resistance: float) -> float:
        if dp == 0.0:
            return 0.0
        return float(np.sign(dp) * np.sqrt(max(float(rho) * abs(float(dp)) / float(resistance), 0.0)))

    @staticmethod
    def _refine_time_grid_max_step(t_grid, max_dt: float) -> np.ndarray:
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

    @staticmethod
    def _schedule_max_abs(schedule, default: float) -> float:
        values = [abs(float(default))]
        if schedule:
            for row in schedule:
                if row is not None and len(row) >= 2:
                    values.append(abs(float(row[1])))
        return max(values)

    # ------------------------------------------------------------------
    def _He_inventory(self):
        """He mass held in the coil void [kg] at inlet conditions (for residence time)."""
        rho = self._thermo.density(self.coolantProp.coolant, self.coolantProp.T_in, self.coolantProp.p_in)
        V = self.A_ch * self.N_ch * self.L_ch_max
        return rho * V

    # ------------------------------------------------------------------
    def time_scale_audit(self):
        print("=" * 60)
        print("TRANSIENT TIME-SCALE AUDIT")
        print("=" * 60)
        s_w = self.combustorProp.thickness_coil_wall
        k = self.func_conductivity_HX(500.0)
        cp = self.func_cp_HX(500.0)
        alpha = k / (self.rho_w * cp)
        tau_diff = s_w ** 2 / (np.pi ** 2 * alpha)
        print(f"  wall: delta={s_w*1e3:.2f} mm, alpha={alpha:.2e} m2/s")
        print(f"  tau_diff,delta (profile shape) = {tau_diff*1e3:.1f} ms  "
              f"(quasi-static profile OK if << ramp)")
        He_inv = self._He_inventory()
        print(f"  He coil inventory = {He_inv*1e3:.1f} g")
        # tau_wall across the mdot_c schedule extremes
        tp = self.transientProp
        for label, mdot in [("min", 1e-3), ("full", self.coolantProp.mass_flow_c)]:
            res_t = He_inv / max(mdot, 1e-9)
            print(f"  He residence @ mdot_c={mdot*1e3:.1f} g/s = {res_t:.2f} s")
        print(f"  (warn: during early ramp He residence may exceed tau_wall -> "
              f"He-outlet T flagged unreliable, section 4.6)")
        print("=" * 60)

    # ------------------------------------------------------------------
    def _print_transient_summary(self):
        ts = self.time_series
        sc = ts["scalars"]
        print()
        print("=" * 60)
        print("TRANSIENT RESULTS")
        print("=" * 60)
        print(f"  Config: {self.combustorProp.HX_config} ({self.combustorProp.flow_config}-flow) | "
              f"{self.N} nodes | t_end={self.transientProp.t_end}s")
        i_end = -1
        print(f"  Final wall T range: {ts['fields']['Tbar'][i_end].min():.1f} -> "
              f"{ts['fields']['Tbar'][i_end].max():.1f} K")
        print(f"  Final He outlet: {sc['T_c_out'][i_end]:.1f} K | "
              f"Final gas outlet: {sc['T_g_out'][i_end]:.1f} K")
        print(f"  Final Q_hot: {sc['Q_hot_kW'][i_end]:.1f} kW | "
              f"Q_cold: {sc['Q_cold_kW'][i_end]:.1f} kW")
        print(f"  Peak wall dT (any node,time): {sc['dT_wall_max'].max():.1f} K "
              f"at t={ts['t'][np.argmax(sc['dT_wall_max'])]:.2f}s")
        n_flag = int(np.sum(sc['He_outlet_reliable'] == 0.0))
        if n_flag:
            print(f"  NOTE: He-outlet T flagged unreliable at {n_flag}/{len(ts['t'])} "
                  f"early-ramp snapshots (residence > tau_wall, section 4.6)")
        # settle check: compare final to steady probe when it was requested.
        if np.any(np.isfinite(self.Tbar_steady)):
            dT_settle = np.max(np.abs(ts['fields']['Tbar'][i_end] - self.Tbar_steady))
            print(f"  Wall vs steady solver at t_end: max|dTbar| = {dT_settle:.1f} K "
                  f"(-> {'settled' if dT_settle < 5 else 'still transient / mismatch'})")
        else:
            print("  Wall vs steady solver at t_end: skipped (no steady reference probe)")
        print("=" * 60)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    from .input_data import (coolantProp, hotgasProp, combustorProp, numericalProp,
                             system_requirements, transientProp)

    solver = transient_solver(coolantProp=coolantProp(), hotgasProp=hotgasProp(),
                              combustorProp=combustorProp(), numericalProp=numericalProp(),
                              system_requirements=system_requirements(),
                              transientProp=transientProp())
    solver.solve_transient()
    # dashboard (WP3c)
    try:
        from .model_data_process.data_plotting_transient import TransientDashboard
        out = TransientDashboard(solver.time_series, meta=dict(
            config=solver.combustorProp.HX_config, material=solver.combustorProp.material_HX,
            Tbar_steady=solver.Tbar_steady.tolist())).to_html()
        print(f"Dynamic dashboard written: {out}")
    except Exception as e:
        print(f"(dashboard skipped: {e})")
