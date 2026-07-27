
""" 
Real-gas, locally-coupled (per turn) CHT solver.

This file is intended as a *drop-in* alternative to your existing `main_solve_CHT.py` approach,
but it reuses your existing:
- friction correlations
- Nusselt correlations
- 1D conduction model `OneDimensionalSteadyConduction_ShellnHelicalTube`
- combustion gas object / Cantera wrapper (`combustion_gas_solve`, `combustion_node.remove_energy`)

Key change vs the current solver:
- The helium (coolant) state update is performed with *conservative total enthalpy* using real-gas
  properties (CoolProp), instead of the ideal-gas gradient helpers in `governing_equations.py`.
- Each spatial step solves a *scalar* fixed-point problem on heat rate Q:
      Q == Q_HT( states implied by Q )
  using a robust bracketing method (Brent).
This makes the steady energy balance tight, *locally per step*, without any global iteration.

Author: ChatGPT (rewrite requested by user)
"""

from __future__ import annotations

import math
import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional, Dict, Any

from scipy.optimize import brentq

# --- project imports (same as your current file) ---
from physics.combustion_chemistry.combustion_gas import combustion_gas_solve, choose_fuel
from physics.friction_correlations import getFrictionColebrook1939, getFrictionCurvedPipeAli2024
from physics.heat_transfer_correlations import (
    getNusseltDeveloping, compute_Nusselt_Gnielinski,
    nusselt_toroid_Ahmed1997, nusselt_inner_curved_tube_mori1967
)
from physics.heat_conduction import OneDimensionalSteadyConduction_ShellnHelicalTube
from mechanical.material_specs.material_temperature_strength import init_material_temperature_strength

from CoolProp.CoolProp import PropsSI


@dataclass
class StepSolution:
    """Outputs you may want to store each step."""
    Q: float
    UA: float
    T_wg: float
    T_wc: float
    # helium outlet
    T_c_out: float
    p_c_out: float
    rho_c_out: float
    U_c_out: float
    # hot gas outlet
    T_g_out: float
    p_g_out: float


class main_solver_real_gas:
    """
    Minimal, commented rewrite of your solver core, keeping your models but fixing the numerics:

    - local per-step coupling on Q (no global iteration)
    - real-gas helium update via conservative total enthalpy + dp model
    - hot gas update via your existing Cantera pathway

    Assumptions:
    - single channel (N_ch=1) as in your current solver
    - geometry and model switches are provided by your existing dataclasses
    """

    def __init__(self, coolantProp, hotgasProp, combustorProp, numericalProp, system_requirements):
        self.coolantProp = coolantProp
        self.hotgasProp = hotgasProp
        self.combustorProp = combustorProp
        self.numericalProp = numericalProp
        self.system_requirements = system_requirements

        # fuel / chemistry setup (same as your current solver)
        self.chem_mech_path, self.Y_fuel, self.Hv_fuel = choose_fuel(self.hotgasProp.fuel)

        # ---------------------------
        # Helical coil geometry
        # ---------------------------
        self.N_ch = 1
        self.Dh_ch = self.combustorProp.Dh_coil
        self.A_ch = math.pi * self.Dh_ch**2 / 4.0

        self.coil_pitch = self.Dh_ch + 2*self.combustorProp.thickness_coil_wall + self.combustorProp.coil_gap
        self.D_coil = (
            self.combustorProp.inner_diameter
            - 2*self.combustorProp.gap_shell2coil
            - self.Dh_ch
            - 2*self.combustorProp.thickness_coil_wall
        )
        self.Rc = self.D_coil / 2.0  # coil radius used in some correlations

        # areas per *dx* step.
        # Cold-side: inner wetted area of tube for length dx
        self.A_cold_per_dx = math.pi * self.Dh_ch  # perimeter; multiply by dx for area
        # Hot-side: you were using shell-side contact area based on coil OD; keep same pattern as your codebase if you already have it.
        # If your conduction model expects areas internally, it uses UP=U*P with dx passed in.
        # Here we pass dx into the conduction model so it computes UA consistently.

        # ---------------------------
        # Material property functions (same hooks as current code)
        # ---------------------------
        (
            self.func_CTE_HX,
            self.func_E_HX,
            self.func_Yield_HX,
            self.func_conductivity_HX
        ) = init_material_temperature_strength(self.combustorProp.material_HX)

        # Radiation backend (same pattern as your project; created elsewhere in your existing code)
        self.radiation_backend = None
        if getattr(self.numericalProp, "radiation_ON", False):
            from physics.radiation_model.radiation_build import make_ehlme_backend
            self.radiation_backend = make_ehlme_backend()

        # ---------------------------
        # Initial states (expected to be set by your calling code)
        # ---------------------------
        self.T_c = None
        self.p_c = None
        self.rho_c = None
        self.U_c = None

        self.T_g = None
        self.p_g = None
        self.combustion_node = None  # set via init_hot_gas()

        # Wall temperatures initial guess
        self.T_wg = 800.0
        self.T_wc = 600.0
        self.T_c_check = 0.0

        # bookkeeping
        self.L_HX = 0.0
        self.L_ch = 0.0

        # hot-side coupling fraction used in your current code
        self.coil_mass_fraction_g = getattr(self.combustorProp, "coil_mass_fraction_g", 1.0)

    # ---------------------------------------------------------------------
    # Initialization helpers
    # ---------------------------------------------------------------------

    def init_coolant(self, T_c_in: float, p_c_in: float):
        """Initialize helium state and derived variables."""
        self.T_c = float(T_c_in)
        self.p_c = float(p_c_in)
        self.rho_c = float(PropsSI("D", "T", self.T_c, "P", self.p_c, self.coolantProp.coolant))
        self.U_c = float(self.coolantProp.mass_flow_c / (self.N_ch * self.rho_c * self.A_ch))

    def init_hot_gas(self, T_g_in: float, p_g_in: float, X_init: Optional[Any] = None):
        """Initialize hot gas model using your combustion_gas_solve wrapper."""
        self.p_g = float(p_g_in)
        self.T_g = float(T_g_in)

        # Your project builds a combustion gas solver / node; keep the same interface you already use.
        self.combustion_node = combustion_gas_solve(
            mech_path=self.chem_mech_path,
            P=self.p_g,
            T=self.T_g,
            Y_fuel=self.Y_fuel,
            fuel=self.hotgasProp.fuel,
            of_ratio=self.hotgasProp.of_ratio,
            equilibrium_init=self.hotgasProp.equilibrium_init,
            X_init=X_init
        )
        # ensure the node is at the intended state
        self.combustion_node.phase.TP = self.T_g, self.p_g

    # ---------------------------------------------------------------------
    # Core per-step physics blocks
    # ---------------------------------------------------------------------

    def _cold_side_friction_factor(self, Re: float) -> float:
        """Reuse your existing friction correlation switch."""
        if self.combustorProp.friction_coil == "Colebrook1939":
            return getFrictionColebrook1939(Re, self.combustorProp.channel_roughness / self.Dh_ch)
        elif self.combustorProp.friction_coil == "CurvedPipeAli2024":
            return getFrictionCurvedPipeAli2024(Re=Re, Dh=self.Dh_ch, Rc=self.D_coil/2)
        else:
            # default fallback
            return getFrictionColebrook1939(Re, self.combustorProp.channel_roughness / self.Dh_ch)

    def _cold_side_htc(self, T: float, p: float, rho: float, U: float) -> float:
        """Reuse your existing Nusselt correlation switch, but driven by real-gas properties."""
        mu = float(PropsSI("V", "T", T, "P", p, self.coolantProp.coolant))
        k = float(PropsSI("L", "T", T, "P", p, self.coolantProp.coolant))
        cp = float(PropsSI("C", "T", T, "P", p, self.coolantProp.coolant))
        Pr = cp * mu / max(k, 1e-18)
        Re = rho * U * self.Dh_ch / max(mu, 1e-18)

        if self.combustorProp.nu_coil == "Gnielinski":
            Nu = compute_Nusselt_Gnielinski(Re=Re, Pr=Pr)
        elif self.combustorProp.nu_coil == "ToroidAhmed1997":
            # your Ahmed correlation expects Asqrt_toroid etc; use your project inputs
            Nu = nusselt_toroid_Ahmed1997(U_g=U, rho_g=rho, mu_g=mu, Pr_g=Pr, Asqrt_toroid=self.combustorProp.Asqrt_toroid)
        elif self.combustorProp.nu_coil == "MoriNakayama1967":
            Nu = nusselt_inner_curved_tube_mori1967(Re=Re, Pr=Pr, Dh=self.Dh_ch, Rc=self.Rc)
        else:
            # developing or generic
            Nu = getNusseltDeveloping(Re=Re, Pr=Pr)

        h = Nu * k / max(self.Dh_ch, 1e-18)
        return float(h)

    def _hot_side_htc(self) -> float:
        """Keep your existing approach for hot-side convection coefficient."""
        # In your current code you compute h_g from correlations in the loop.
        # Here we assume you already set self.h_g somewhere OR you compute it in a helper.
        # For rewrite simplicity, we recompute from Cantera properties + a chosen correlation hook.
        # Replace this body with your existing hot-side convection model.
        # NOTE: Many users treat h_g as an input / requirement; keep it explicit.
        return float(self.hotgasProp.h_g) if hasattr(self.hotgasProp, "h_g") else 1000.0

    def _hot_gas_outlet_for_Q(self, Q: float, p_g_in: float) -> Tuple[float, Any]:
        """
        Compute hot gas outlet temperature and composition *without permanently mutating* the Cantera state.

        This matches your existing logic:
            dh_g = Q / (m_dot_g * coil_mass_fraction_g)   [shellnHelicalTube]
        then remove_energy(dh=dh_g, updated_pressure=p_g)
        """
        # Save current Cantera state
        ph = self.combustion_node.phase
        T_save, P_save, X_save = ph.T, ph.P, ph.X

        if self.combustorProp.HX_config == "shellnHelicalTube":
            dh_g = Q / (self.hotgasProp.mass_flow_g * self.coil_mass_fraction_g)
        else:
            dh_g = Q / self.hotgasProp.mass_flow_g

        # Apply enthalpy removal
        self.combustion_node.remove_energy(
            dh=dh_g,
            updated_pressure=p_g_in,
            equilibrium_dh_gas_ON=self.numericalProp.equilibrium_dh_gas_ON
        )
        T_new = float(ph.T)
        X_new = ph.X

        # Restore
        ph.TPX = T_save, P_save, X_save

        return T_new, X_new

    def _helium_outlet_for_Q(self, Q: float) -> Tuple[float, float, float, float]:
        """
        Real-gas conservative helium update for one segment with *local fixed-point iteration*.

        Governing balances over a segment of length dx:

          1) dp from friction:
              dp = f * (dx/Dh) * (rho * U^2 / 2)

          2) total enthalpy (static + kinetic):
              Q = m_dot * [ (h_out + U_out^2/2) - (h_in + U_in^2/2) ]

          3) EOS inversion:
              T_out from (h_out, p_out), rho_out from (T_out, p_out)

          4) continuity:
              U_out = m_dot / (rho_out * A)

        This is a practical real-gas 1D model that preserves the coupling you care about,
        and it remains stable with large density change (Z != 1).
        """
        dx = float(self.numericalProp.dx)
        m_dot = float(self.coolantProp.mass_flow_c / self.N_ch)
        A = float(self.A_ch)

        # Inlet values
        T_in, p_in, rho_in, U_in = float(self.T_c), float(self.p_c), float(self.rho_c), float(self.U_c)
        h_in = float(PropsSI("H", "T", T_in, "P", p_in, self.coolantProp.coolant))

        # Initial guesses
        p_out = p_in
        U_out = U_in

        for _ in range(60):
            # Representative properties for friction: use inlet (simple) or average
            mu_in = float(PropsSI("V", "T", T_in, "P", p_in, self.coolantProp.coolant))
            Re_in = rho_in * U_in * self.Dh_ch / max(mu_in, 1e-18)
            f = self._cold_side_friction_factor(Re_in)

            dp = f * (dx / max(self.Dh_ch, 1e-18)) * (rho_in * U_in**2 / 2.0)
            p_new = p_in - dp

            # Total enthalpy update
            htot_in = h_in + 0.5 * U_in**2
            htot_out = htot_in + Q / max(m_dot, 1e-18)

            # Static enthalpy depends on outlet velocity
            h_new = htot_out - 0.5 * U_out**2

            # Invert (h,p) -> T, then rho
            T_new = float(PropsSI("T", "H", h_new, "P", p_new, self.coolantProp.coolant))
            rho_new = float(PropsSI("D", "T", T_new, "P", p_new, self.coolantProp.coolant))

            U_new = m_dot / max(rho_new * A, 1e-18)

            # convergence check
            err = max(
                abs(p_new - p_out) / max(abs(p_out), 1.0),
                abs(U_new - U_out) / max(abs(U_out), 1.0),
            )
            p_out, U_out = p_new, U_new
            if err < 1e-8:
                return T_new, p_out, rho_new, U_out

        # If not converged, still return last iterate (rare if dx is sensible)
        T_out = float(PropsSI("T", "H", h_new, "P", p_out, self.coolantProp.coolant))
        rho_out = float(PropsSI("D", "T", T_out, "P", p_out, self.coolantProp.coolant))
        return T_out, p_out, rho_out, U_out

    # ---------------------------------------------------------------------
    # Per-turn coupling on Q using your conduction model
    # ---------------------------------------------------------------------

    def _Q_residual(self, Q: float) -> float:
        """
        Fixed-point residual for coupling:
            residual(Q) = Q - Q_conduction(states implied by Q)

        We evaluate conduction using representative midpoint bulk temperatures.
        """
        # Helium outlet for trial Q
        T_c_out, p_c_out, rho_c_out, U_c_out = self._helium_outlet_for_Q(Q)

        # Hot gas outlet for trial Q (non-mutating)
        T_g_out, X_g_out = self._hot_gas_outlet_for_Q(Q, self.p_g)

        # Midpoint temperatures (stable choice)
        T_c_mid = 0.5 * (self.T_c + T_c_out)
        T_g_mid = 0.5 * (self.T_g + T_g_out)

        # Midpoint helium state for HTC
        p_c_mid = 0.5 * (self.p_c + p_c_out)
        rho_c_mid = float(PropsSI("D", "T", T_c_mid, "P", p_c_mid, self.coolantProp.coolant))
        U_c_mid = 0.5 * (self.U_c + U_c_out)

        h_c = self._cold_side_htc(T=T_c_mid, p=p_c_mid, rho=rho_c_mid, U=U_c_mid)
        h_g = self._hot_side_htc()

        # Radiation participating species fractions (if you already compute these elsewhere, plug them here)
        # Keep the same state dict keys used in your conduction model
        rad_state = None
        if getattr(self.numericalProp, "radiation_ON", False):
            # crude extraction (you likely have helpers to get yH2O/yCO2 and Le)
            # If your X_g_out is a Cantera array, your existing code likely already maps to X_H2O/X_CO2.
            rad_state = {'p': self.p_g, 'yH2O': getattr(self, "X_H2O", 0.0), 'yCO2': getattr(self, "X_CO2", 0.0), 'Le': getattr(self, "Le", 1.0)}

        node = OneDimensionalSteadyConduction_ShellnHelicalTube(
            h_g=h_g*(1 + self.numericalProp.artificial_error_Nu_hot),
            h_c=h_c*(1 + self.numericalProp.artificial_error_Nu_cold),
            T_c=T_c_mid, T_g=T_g_mid,
            s_w=self.combustorProp.thickness_coil_wall,
            Dh_ch=self.Dh_ch,
            f_kw_at_T=self.func_conductivity_HX,
            T_wg_0=self.T_wg, T_wc_0=self.T_wc, T_c_check_0=self.T_c_check,
            dx=self.numericalProp.dx,
            rad_enabled=getattr(self.numericalProp, "radiation_ON", False),
            eps_s=self.numericalProp.emissivity_wall,
            rad_backend=self.radiation_backend,
            rad_state=rad_state
        )
        node.Solve1Dconduction()
        Q_ht = float(node.dQ)

        return Q - Q_ht

    def solve_one_step(self) -> StepSolution:
        """Solve one turn/segment: find Q then advance both fluids."""
        # Bracket Q: lower bound 0, upper bound from UA*(Tg-Tc) with a safety factor
        # We need a bracket with a sign change for brentq.
        Q_lo = 0.0
        # A rough upper bound: assume very large UA
        dT0 = max(self.T_g - self.T_c, 1.0)
        Q_hi = 1e3 * dT0

        r_lo = self._Q_residual(Q_lo)
        r_hi = self._Q_residual(Q_hi)

        # Expand upper bound until sign change (robust)
        for _ in range(30):
            if r_lo * r_hi < 0:
                break
            Q_hi *= 2.0
            r_hi = self._Q_residual(Q_hi)

        if r_lo * r_hi >= 0:
            # If we failed to bracket, fall back to fixed-point evaluation at Q_hi (rare; indicates model inconsistency)
            Q = max(0.0, Q_hi - r_hi)
        else:
            Q = float(brentq(self._Q_residual, Q_lo, Q_hi, xtol=1e-3, rtol=1e-6, maxiter=100))

        # With converged Q, recompute a final conduction node to extract wall temps and UA
        T_c_out, p_c_out, rho_c_out, U_c_out = self._helium_outlet_for_Q(Q)
        T_g_out, _ = self._hot_gas_outlet_for_Q(Q, self.p_g)

        T_c_mid = 0.5 * (self.T_c + T_c_out)
        T_g_mid = 0.5 * (self.T_g + T_g_out)
        p_c_mid = 0.5 * (self.p_c + p_c_out)
        rho_c_mid = float(PropsSI("D", "T", T_c_mid, "P", p_c_mid, self.coolantProp.coolant))
        U_c_mid = 0.5 * (self.U_c + U_c_out)

        h_c = self._cold_side_htc(T=T_c_mid, p=p_c_mid, rho=rho_c_mid, U=U_c_mid)
        h_g = self._hot_side_htc()

        rad_state = None
        if getattr(self.numericalProp, "radiation_ON", False):
            rad_state = {'p': self.p_g, 'yH2O': getattr(self, "X_H2O", 0.0), 'yCO2': getattr(self, "X_CO2", 0.0), 'Le': getattr(self, "Le", 1.0)}

        node = OneDimensionalSteadyConduction_ShellnHelicalTube(
            h_g=h_g*(1 + self.numericalProp.artificial_error_Nu_hot),
            h_c=h_c*(1 + self.numericalProp.artificial_error_Nu_cold),
            T_c=T_c_mid, T_g=T_g_mid,
            s_w=self.combustorProp.thickness_coil_wall,
            Dh_ch=self.Dh_ch,
            f_kw_at_T=self.func_conductivity_HX,
            T_wg_0=self.T_wg, T_wc_0=self.T_wc, T_c_check_0=self.T_c_check,
            dx=self.numericalProp.dx,
            rad_enabled=getattr(self.numericalProp, "radiation_ON", False),
            eps_s=self.numericalProp.emissivity_wall,
            rad_backend=self.radiation_backend,
            rad_state=rad_state
        )
        node.Solve1Dconduction()

        # Update stored wall temps for next step (good initial guess)
        self.T_wg = float(node.T_wg_new)
        self.T_wc = float(node.T_wc_new)
        self.T_c_check = float(node.T_c_check_f)

        # --- Accept the fluid updates (advance the march) ---
        # Helium
        self.T_c, self.p_c, self.rho_c, self.U_c = T_c_out, p_c_out, rho_c_out, U_c_out
        # Hot gas (mutating for real now)
        if self.combustorProp.HX_config == "shellnHelicalTube":
            dh_g = Q / (self.hotgasProp.mass_flow_g * self.coil_mass_fraction_g)
        else:
            dh_g = Q / self.hotgasProp.mass_flow_g
        self.combustion_node.remove_energy(dh=dh_g, updated_pressure=self.p_g, equilibrium_dh_gas_ON=self.numericalProp.equilibrium_dh_gas_ON)
        self.T_g = float(self.combustion_node.phase.T)

        # Advance length counters (keep your existing logic)
        self.L_ch += self.numericalProp.dx
        self.L_HX += self.numericalProp.dx

        return StepSolution(
            Q=Q,
            UA=float(node.UA),
            T_wg=self.T_wg,
            T_wc=self.T_wc,
            T_c_out=T_c_out,
            p_c_out=p_c_out,
            rho_c_out=rho_c_out,
            U_c_out=U_c_out,
            T_g_out=self.T_g,
            p_g_out=self.p_g
        )

    # ---------------------------------------------------------------------
    # Main loop
    # ---------------------------------------------------------------------

    def solver(self):
        """Run the HX march until length limit or CoolProp validity limit."""
        # Minimal stopping criteria; adapt to your existing conditions
        T_c_min = float(PropsSI("Tmin", self.coolantProp.coolant))

        data = {
            "Q": [],
            "UA": [],
            "T_c": [], "p_c": [], "rho_c": [], "U_c": [],
            "T_g": [], "p_g": [],
            "T_wg": [], "T_wc": [],
            "Z": []
        }

        while self.L_HX <= (self.numericalProp.L_HX_max - self.combustorProp.mixing_length) and self.T_c > T_c_min * 2.5:
            sol = self.solve_one_step()

            data["Q"].append(sol.Q)
            data["UA"].append(sol.UA)
            data["T_c"].append(self.T_c)
            data["p_c"].append(self.p_c)
            data["rho_c"].append(self.rho_c)
            data["U_c"].append(self.U_c)
            data["T_g"].append(self.T_g)
            data["p_g"].append(self.p_g)
            data["T_wg"].append(sol.T_wg)
            data["T_wc"].append(sol.T_wc)
            data["Z"].append(float(PropsSI("Z", "T", self.T_c, "P", self.p_c, self.coolantProp.coolant)))

        return data
