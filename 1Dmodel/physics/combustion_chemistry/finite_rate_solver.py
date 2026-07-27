""" 
@ author : Raphaël Aubry
"""
#%%
import cantera as ctr
import numpy as np
from pathlib import Path

ctr.suppress_thermo_warnings()

#%%

base_path = Path(__file__).parent
JETA_chem_path = base_path / "A2highT.yaml"

import cantera as ct
import numpy as np

class FiniteRateGasTracker:
    """
    UA-free finite-rate chemistry tracker for a flowing control volume.
    - Energy removal is applied explicitly by lowering mixture enthalpy (Δh) each substep.
    - Residence time per cell: self.tau_conv = dx / U (U from mdot, rho, area).
    - Reactor volume is consistent with geometry: V = A * dx.
    """

    def __init__(self, fuel, oxidizer, OF=None, OF_st=None,
                 T_init=300.0, p_init=1e5, mass_flow=0.1, Dh=0.05, dx=0.01,
                 raise_on_negative_T=True, mech_file=JETA_chem_path):
        
        self.gas = ct.Solution(mech_file)

        self.gas.set_equivalence_ratio(phi=OF_st/OF, fuel=fuel, oxidizer=oxidizer)  # use your oxidizer string if not pure O2
        self.gas.TP = float(T_init), float(p_init)

        # --- geometry & reactor ---
        self.Dh = float(Dh)
        self.dx = float(dx)
        self.A  = 0.25 * np.pi * self.Dh**2
        self.V  = self.A * self.dx
        self.reactor = ct.IdealGasReactor(self.gas, energy='on', volume=self.V)
        self.net = ct.ReactorNet([self.reactor])

        # --- flow & controls ---
        self.mdot = float(mass_flow)
        self.raise_on_negative_T = bool(raise_on_negative_T)

        # initialize attributes
        self._assign_props()

    def _reaction_time_scales(self):
        """Species-based proxy: tau_k ≈ (Y_k * rho)/|wdot_k| (mass basis)."""
        rho = self.gas.density
        Y   = np.clip(self.gas.Y, 1e-16, None)
        wdot = self.gas.net_production_rates             # kmol/m^3/s
        MW    = self.gas.molecular_weights               # kg/kmol
        wdot_m = wdot * MW                               # kg/m^3/s
        tau = np.where(np.abs(wdot_m) > 0.0, Y * rho / np.abs(wdot_m), np.inf)
        finite = np.isfinite(tau)
        tau_min = np.min(tau[finite]) if np.any(finite) else np.inf
        tau_med = np.median(tau[finite]) if np.any(finite) else np.inf
        return tau_min, tau_med

    def _assign_props(self):
        """Assign current thermophysical properties to attributes (single place)."""
        self.T   = self.gas.T
        self.p   = self.gas.P
        self.rho = self.gas.density
        self.mu  = self.gas.viscosity
        self.k   = self.gas.thermal_conductivity
        self.cp  = self.gas.cp_mass
        self.cv = self.gas.cv_mass
        self.gamma = self.cp / self.cv
        self.W = self.gas.mean_molecular_weight
        self.Y   = self.gas.Y

    def step_with_energy_removal(self, Q_cell_W):
        """
        Remove Q_cell_W [W] over the cell residence time by decrementing enthalpy.
        No UA, no wall; all cooling is explicit via Δh each substep.
        """
        # residence time
        self.U = self.mdot / max(self.rho * self.A, 1e-12)
        self.tau_conv = self.dx / max(self.U, 1e-12)

        # target energy in this cell
        E_target = max(float(Q_cell_W), 0.0) * self.tau_conv

        # substep integration for stability
        N = max(10, int(np.ceil(self.tau_conv / max(self.tau_conv/20.0, 1e-6))))
        dt = self.tau_conv / N
        E_removed = 0.0

        for _ in range(N):
            # advance finite-rate chemistry (no cooling during this tiny dt)
            self.net.advance(self.net.time + dt)

            # apply explicit energy removal for this substep
            dE_sub = min(Q_cell_W * dt, E_target - E_removed)
            if dE_sub > 0.0:
                rho_now = self.gas.density
                h_now   = self.gas.enthalpy_mass
                dh      = dE_sub / max(rho_now * self.V, 1e-18)  # J/kg
                self.gas.HP = h_now - dh, self.gas.P             # enforce Δh at constant p
                self.reactor.syncState()                          # keep reactor consistent
                E_removed += dE_sub

            # safety: negative T is physically invalid
            if self.raise_on_negative_T and self.gas.T < 0.0:
                raise ValueError(f"Negative gas temperature encountered: {self.gas.T:.6g} K")

            if E_removed >= E_target:
                break

        # diagnostics & attributes
        self.nu   = self.mu / max(self.rho, 1e-12)
        self.alpha= self.k / max(self.rho * self.cp, 1e-12)
        
        self.tau_visc = self.Dh**2 / max(self.nu, 1e-16)
        self.tau_th   = self.Dh**2 / max(self.alpha, 1e-16)
        self.tau_chem_min, self.tau_chem_med = self._reaction_time_scales()
        self.Da_min = self.tau_conv / self.tau_chem_min if np.isfinite(self.tau_chem_min) and self.tau_chem_min > 0 else 0.0
        self.Da_med = self.tau_conv / self.tau_chem_med if np.isfinite(self.tau_chem_med) and self.tau_chem_med > 0 else 0.0

        # final props
        self._assign_props()
