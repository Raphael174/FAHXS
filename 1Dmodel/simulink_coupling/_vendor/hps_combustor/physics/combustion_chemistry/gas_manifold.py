"""
Shared 1-D equilibrium/frozen gas-cooling manifold — the C0 slice of the FGM
(DESIGN_PLAN_shellntube_transient.md section 4b/C0), factored out so both the
helical transient solver and the shell-and-tube transient solver can build the
"no Cantera in the march" gas-property table without duplicating the logic.

For finite-rate chemistry, reuse `fpv_manifold.build_fpv_manifold`/`FPVManifold`
directly (already generic over any Cantera gas object) — not duplicated here.

@ author : Raphaël Aubry
"""
import numpy as np


def build_equilibrium_manifold(gas, T_inlet, p_inlet, Y_inlet, mode="equilibrium",
                               n_h=250, T_floor=340.0):
    """
    Tabulate gas state vs specific enthalpy removed [J/kg] at fixed inlet
    composition/O/F/pressure. `mode`: "equilibrium" (re-equilibrate at each
    level) or "frozen" (composition held at Y_inlet throughout).

    Mutates `gas` locally; restores it to (T_inlet, p_inlet, Y_inlet) on return.
    """
    gas.TPY = T_inlet, p_inlet, Y_inlet
    h0 = float(gas.enthalpy_mass)
    eq = (mode == "equilibrium")

    hr, T, rho, mu, k, cp, xh2o, xco2 = ([] for _ in range(8))
    iH2O = gas.species_index("H2O"); iCO2 = gas.species_index("CO2")

    def _rec(h_removed):
        hr.append(h_removed); T.append(float(gas.T)); rho.append(float(gas.density))
        mu.append(float(gas.viscosity)); k.append(float(gas.thermal_conductivity))
        cp.append(float(gas.cp)); xh2o.append(float(gas.X[iH2O])); xco2.append(float(gas.X[iCO2]))

    _rec(0.0)
    dh = 2200.0 * (T_inlet - T_floor) / n_h
    cur = h0
    for _ in range(3 * n_h):
        cur -= dh
        gas.HPY = cur, p_inlet, Y_inlet if not eq else gas.Y
        if eq:
            gas.equilibrate('HP')
        _rec(h0 - cur)
        if gas.T < T_floor:
            break

    gas.TPY = T_inlet, p_inlet, Y_inlet  # restore inlet state

    return dict(h_grid=np.array(hr), T=np.array(T), rho=np.array(rho), mu=np.array(mu),
               k=np.array(k), cp=np.array(cp), xH2O=np.array(xh2o), xCO2=np.array(xco2))


class EquilibriumManifold:
    """Runtime interpolator over a built equilibrium/frozen manifold."""

    def __init__(self, m):
        self.m = m
        self.h = m["h_grid"]

    def at(self, h_removed):
        """Return (T, rho, mu, k, cp, xH2O, xCO2) at the given enthalpy removed."""
        m = self.m
        return (float(np.interp(h_removed, self.h, m["T"])),
                float(np.interp(h_removed, self.h, m["rho"])),
                float(np.interp(h_removed, self.h, m["mu"])),
                float(np.interp(h_removed, self.h, m["k"])),
                float(np.interp(h_removed, self.h, m["cp"])),
                float(np.interp(h_removed, self.h, m["xH2O"])),
                float(np.interp(h_removed, self.h, m["xCO2"])))
